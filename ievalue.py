import sys
import shlex
import subprocess
from itertools import chain
from typing import Optional
from ievalue_db import IevalueDB, DatabaseData, HitData, DbIdx, HitIdx


def get_delta_sizes(dbs: list[str]) -> list[DatabaseData]:
    db_sizes = []
    for db in dbs:
        db_residue, _ = read_database_parts(db)
        db_sizes.append((db, db_residue))
    return db_sizes


def read_database_parts(db_name: str) -> tuple[int, list[str]]:
    try:
        result = subprocess.run(
            ["blastdbcmd", "-db", db_name, "-info", "-exact_length"], stdout=subprocess.PIPE, check=True
        )
    except subprocess.CalledProcessError:
        raise ValueError("TODO write error")

    db_info = result.stdout.decode("utf-8").split("\n")

    capture_db_size = False
    capture_hits = False
    db_residue = -1
    parts = []
    for line in db_info:
        if line != "":
            if line.startswith("Database:"):
                capture_db_size = True
            elif capture_db_size:
                data = line.strip().split(" ")
                db_residue = int(data[2].replace(",", ""))
                capture_db_size = False

            if line == "Volumes:":
                capture_hits = True
            elif capture_hits:
                parts.append(line.strip().split("/")[-1])

    if db_residue == -1:
        raise ValueError("TODO write error")

    return db_residue, parts


def get_delta_db(
    blast_database: str, blast_type: str, out_file_name: str, prev_dbs: list[Optional[str]]
) -> tuple[str, list[str]]:
    _, curr_dbs = read_database_parts(blast_database)
    delta_parts = list(set(curr_dbs).difference(set(prev_dbs)))

    if len(delta_parts) == 0:
        print("No new data")
        return "", []

    if len(delta_parts) > 1:
        delta_dbs = " ".join(delta_parts)
        delta_dbs = f'"{delta_dbs}"'
    else:
        delta_dbs = delta_parts[0]

    # TODO try/except
    command = [
        "blastdb_aliastool",
        "-dbtype",
        blast_type,
        "-dblist",
        delta_dbs,
        "-out",
        f"{out_file_name}-delta",
        "-title",
        f"{out_file_name}-delta",
    ]
    print(" ".join(command))
    subprocess.run(command, check=True)

    return delta_dbs, delta_parts


def do_delta_blast(blast_kwargs: dict[str, str], delta_dbs: str) -> None:
    delta_kwargs = blast_kwargs.copy()

    blast_program = delta_kwargs.pop("blast_program")
    delta_kwargs.pop("blast_type")

    delta_kwargs["db"] = delta_dbs
    delta_kwargs["out"] = f"{delta_kwargs['out']}-delta"
    command = list(chain.from_iterable([[f"-{kw}", arg] for kw, arg in delta_kwargs.items()]))
    command.insert(0, blast_program)
    print(" ".join(command))
    subprocess.run(command, check=True)


def deconstruct_call(program_call: str) -> dict[str, str]:
    # TODO add verification that needed kwargs are in the call
    blast_kwargs = {}
    blast_program, *kwargs = [kwarg.strip() for kwarg in program_call.split("-")]

    blast_kwargs["blast_program"] = blast_program
    if blast_program == "blastn":
        blast_kwargs["blast_type"] = "nuc"
    elif blast_program == "blastp":
        blast_kwargs["blast_type"] = "prot"
    else:
        raise ValueError(f"Unknown program call: {blast_program}")

    for kwarg in kwargs:
        # TODO unsure if shlex.split is needed here, but it won't hurt for now
        kw, arg = shlex.split(kwarg)
        blast_kwargs[kw] = arg

    return blast_kwargs


def get_updated_evalues(hits: list[HitData], total_residues: int, partial_residues: int) -> list[HitData]:
    updates = []
    for hit in hits:
        scaling_factor = (1.0 * total_residues) / (1.0 * partial_residues)
        evalue = hit[HitIdx.EVALUE] * scaling_factor
        updated_hit = (hit[HitIdx.QUERY], hit[HitIdx.HIT], evalue, hit[HitIdx.THE_REST])
        updates.append(updated_hit)

    return updates


def parse_delta_db(out_file_name: str) -> list[HitData]:
    delta_file_name = out_file_name + "-delta"
    hits = []
    with open(delta_file_name, "r") as f:
        for line in f:
            result = line.strip().split("\t")
            hits.append(
                (
                    result[0],
                    result[1],
                    float(f"{float(result[10]):.5}"),
                    "\t".join(result[2:10] + [result[11]]),
                )
            )
    return hits


def write_updated_output(query_file_name: str, out_file_name: str, db) -> None:
    out_file_name = out_file_name.split(".")[0] + ".m8"
    with open(query_file_name, "r") as query_f, open(out_file_name, "w") as out_f:
        for line in query_f:
            if line.startswith(">"):
                query = line.strip()[1:]
                # values = query, hit, evalue
                for hit in db.get_query(query):
                    # TODO DON'T HARD CODE THIS, GET IT FROM KWARGS
                    evalue = hit[HitIdx.EVALUE]
                    if evalue <= 10:
                        the_rest = hit[HitIdx.THE_REST].split("\t")
                        all_values = [query, hit[HitIdx.HIT]] + the_rest[:-1] + [f"{evalue:.3}"] + [the_rest[-1]]
                        print("\t".join(all_values), file=out_f)


def make_updated_evalue_dict(delta_updated_evalue_hits: list[HitData]) -> dict[str, float]:
    updated_evalue_lookup = {}
    for hit in delta_updated_evalue_hits:
        updated_evalue_lookup[f"{hit[HitIdx.QUERY]}{hit[HitIdx.HIT]}{hit[HitIdx.THE_REST]}"] = hit[HitIdx.EVALUE]
    return updated_evalue_lookup


def main(program_call: str) -> None:
    db = IevalueDB()
    blast_kwargs = deconstruct_call(program_call)

    # TODO only search sequences that have not been searched on the exact previous databases
    prev_dbs, prev_residues = [], 0
    if prev_data := db.get_db_info():
        prev_dbs = [db_data[DbIdx.DATABASE] for db_data in prev_data]  # type: ignore
        prev_residues = int(sum([db_data[DbIdx.RESIDUE] for db_data in prev_data]))  # type: ignore # int for type checker

    # dbs that have not been searched on yet
    delta_dbs, delta_parts = get_delta_db(
        blast_kwargs["db"], blast_kwargs["blast_type"], blast_kwargs["out"], prev_dbs  # type: ignore
    )
    if not delta_dbs:
        return

    delta_data = get_delta_sizes(delta_parts)

    # add the size of the residues for each database
    db.add_database_record(delta_data)

    # update evalues from the old hits
    delta_residue = int(sum([delta_db[DbIdx.RESIDUE] for delta_db in delta_data]))  # int for type checker
    total_residues = prev_residues + delta_residue

    # TODO better explain the math here

    # update the evalues
    if prev_residues:
        db.update_old_evalues(total_residues, prev_residues)

    # perform blast on just the new dbs
    do_delta_blast(blast_kwargs, delta_dbs)

    # get the new results and add them to the database
    delta_hits = parse_delta_db(blast_kwargs["out"])

    delta_updated_evalue_hits = get_updated_evalues(delta_hits, total_residues, delta_residue)
    db.insert_hits(delta_updated_evalue_hits)

    write_updated_output(blast_kwargs["query"], blast_kwargs["out"], db)


if __name__ == "__main__":
    print(sys.argv[1])
    main(sys.argv[1])
