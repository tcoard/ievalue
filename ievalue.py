import sys
import shlex
import subprocess
from itertools import chain
from typing import Optional, cast
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
                # parts.append(line.strip().split("/")[-1])
                parts.append(line.strip())

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


def run_blast(blast_kwargs: dict[str, str], delta_dbs: str, is_delta: bool) -> None:
    delta_kwargs = blast_kwargs.copy()

    blast_program = delta_kwargs.pop("blast_program")
    delta_kwargs.pop("blast_type")

    delta_kwargs["db"] = delta_dbs
    if is_delta:
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
        evalue = cast(float, hit[HitIdx.EVALUE]) * scaling_factor  # float for type checker
        updated_hit = (hit[HitIdx.QUERY], hit[HitIdx.HIT], evalue, hit[HitIdx.THE_REST])
        updates.append(updated_hit)

    return updates


def parse_delta_db(out_file_name: str, is_delta: bool) -> list[HitData]:
    if is_delta:
        out_file_name = out_file_name + "-delta"
    hits = []
    with open(out_file_name, "r") as f:
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


def main(program_call: str, path: str) -> None:
    if path and path[-1] != "/":
        path += "/"
    db = IevalueDB(path)
    blast_kwargs = deconstruct_call(program_call)
    blast_kwargs["out"] = path + blast_kwargs["out"]

    # TODO only search sequences that have not been searched on the exact previous databases
    prev_dbs, prev_residues = [], 0
    if prev_data := db.get_db_info():
        prev_data = cast(list[DatabaseData], prev_data) # for the type checker
        prev_dbs = [db_data[DbIdx.DATABASE] for db_data in prev_data]
        prev_residues = int(sum([db_data[DbIdx.RESIDUE] for db_data in prev_data]))

    # dbs that have not been searched on yet
    delta_dbs, delta_parts = get_delta_db(
        blast_kwargs["db"], blast_kwargs["blast_type"], blast_kwargs["out"], prev_dbs  # type: ignore
    )
    if not delta_dbs:
        return

    delta_data = get_delta_sizes(delta_parts)

    # add the size of the residues for each database
    db.add_database_record(delta_data)

    if prev_residues:
        # update evalues from the old hits
        delta_residue = cast(int, sum([delta_db[DbIdx.RESIDUE] for delta_db in delta_data]))  # int for type checker
        total_residues = prev_residues + delta_residue
        db.update_old_evalues(total_residues, prev_residues)

        # perform blast on just the new dbs
        run_blast(blast_kwargs, delta_dbs, is_delta=True)

        # get the new results and add them to the database
        delta_hits = parse_delta_db(blast_kwargs["out"], is_delta=True)
        updated_delta_hits = get_updated_evalues(delta_hits, total_residues, delta_residue)
        db.insert_hits(updated_delta_hits)

        write_updated_output(blast_kwargs["query"], blast_kwargs["out"], db)

    else:
        run_blast(blast_kwargs, delta_dbs, is_delta=False)
        delta_hits = parse_delta_db(blast_kwargs["out"], is_delta=False)
        db.insert_hits(delta_hits)


if __name__ == "__main__":
    if len(sys.argv) >= 2:
        _path = ""
        if len(sys.argv) == 3:
            _path = sys.argv[2]
        elif len(sys.argv) > 3:
            print("Error: Too many input arguments")
        main(sys.argv[1], _path)
    else:
        print("Error: Need blast kwargs")
