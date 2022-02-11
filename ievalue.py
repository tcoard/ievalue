import sys
import shlex
import subprocess
from typing import Optional, cast, Any, Union
from ievalue_db import IevalueDB, DatabaseData, HitData, DbIdx, HitIdx

# TODO:
# * change ValueError to something better

# TODO:
# Verify numbers and add more


def get_needed_defaults(program_args: list[str]) -> float:
    program = program_args[0]
    needed_defaults: dict[str, Any] = {}
    defaults: dict[str, dict[str, Any]] = {
        "blast": {"evalue_cutoff": 10},
        "mmseqs": {"evalue_cutoff": 1.000e-03},
        "diamond": {"evalue_cutoff": 0.001},
    }
    for needed_arg in defaults[program]:
        val, _ = normalize_args(program_args, needed_arg)
        if val:
            needed_defaults[needed_arg] = float(val)
        else:
            needed_defaults[needed_arg] = defaults[program][needed_arg]

    return needed_defaults["evalue_cutoff"]


def normalize_args(program_args: list[str], arg: str) -> tuple[Optional[str], Optional[int]]:
    # this will be a normalization layer from application specific keywords to normalized keyword
    program = program_args[0]
    needed = ["db", "out", "query"]
    normalized_keywords: dict[str, dict[str, Union[str, list[str], int]]] = {
        "blastp": {"db": "-db", "out": "-out", "query": "-query", "evalue_cutoff": "-evalue"},
        "diamond": {"db": "-d", "out": "-o", "query": "-q", "evalue_cutoff": "-e"},
        "mmseqs": {"db": 3, "out": 4, "query": 2, "evalue_cutoff": ["-e", "--evalue"]},
    }
    idx = normalized_keywords[program][arg]
    if isinstance(idx, str):
        try:
            idx = program_args.index(idx) + 1
        except ValueError:
            if arg not in needed:
                return None, None
            raise ValueError(f"Could not find {arg}")

    elif isinstance(idx, list):
        for variation in idx:
            try:
                idx = program_args.index(variation) + 1
                break
            except KeyError:
                pass
        else:
            # since this is not for any needed keywords, we can skip this for now
            return None, None

    value = program_args[idx]
    return value, idx


def get_delta_sizes(dbs: list[str], search_type: str) -> list[DatabaseData]:
    db_sizes: list[DatabaseData] = []
    for db in dbs:
        db_residue, _ = read_database_parts(db, search_type)
        db_sizes.append((db, db_residue))
    return db_sizes


def read_blast_database_parts(db_name: str, search_type: str) -> tuple[int, list[str]]:
    db_residue = 0
    parts = []

    command = ["blastdbcmd", "-db", db_name, "-info", "-exact_length"]
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError:
        raise ValueError(f"Blast Database Information Query Failed: {command}")

    db_info = result.stdout.decode("utf-8").split("\n")
    capture_db_size = False
    capture_hits = False
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
    return db_residue, parts


def read_mmseqs_database_parts(db_name: str, search_type: str) -> tuple[int, list[str]]:
    db_residue = 0
    parts = []
    with open(f"{db_name}.source", "r") as f:
        for line in f:
            parts.append(line.strip().split("\t")[1])

    with open(db_name, "r") as f:
        for line in f:
            if not line.startswith(">"):
                db_residue += len(line.strip())
    return db_residue, parts


def read_diamond_database_parts(db_name: str, search_type: str) -> tuple[int, list[str]]:
    db_residue = 0
    parts = []
    parts.append(f"{db_name}.fa")
    command = ["diamond", "dbinfo", "--db", f"{db_name}.dmnd"]
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError:
        raise ValueError(f"Diamond Database Information Query Failed: {command}")

    db_info = result.stdout.decode("utf-8").split("\n")
    for line in db_info:
        line = line.strip()
        if line.startswith("Letters"):
            db_residue = int(line.split(" ")[-1])
    return db_residue, parts


def read_database_parts(db_name: str, search_type: str) -> tuple[int, list[str]]:
    db_residue = 0
    parts: list[str] = []
    if search_type == "blastp":
        db_residue, parts = read_blast_database_parts(db_name, search_type)
    elif search_type == "mmseqs":
        db_residue, parts = read_mmseqs_database_parts(db_name, search_type)
    elif search_type == "diamond":
        db_residue, parts = read_diamond_database_parts(db_name, search_type)

    if db_residue == 0:
        raise ValueError("Error finding database information")

    return db_residue, parts


def get_blast_delta_db(
    blast_database: str, search_type: str, out_file_name: str, prev_dbs: list[Optional[str]]
) -> tuple[str, list[str]]:
    _, curr_dbs = read_database_parts(blast_database, search_type)
    delta_parts = list(set(curr_dbs).difference(set(prev_dbs)))

    if len(delta_parts) == 0:
        print("No new data")
        return "", []

    if len(delta_parts) > 1:
        delta_dbs = " ".join(delta_parts)
        delta_dbs = f'"{delta_dbs}"'
    else:
        delta_dbs = delta_parts[0]

    command = [
        "blastdb_aliastool",
        "-dbtype",
        search_type,
        "-dblist",
        delta_dbs,
        "-out",
        f"{out_file_name}-delta",
        "-title",
        f"{out_file_name}-delta",
    ]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError:
        raise ValueError(f"Failed to read blast database: {command}")

    return delta_dbs, delta_parts


def run_program(program_args, delta_dbs: str, is_delta: bool) -> None:
    delta_args = program_args.copy()

    _, db_idx = normalize_args(delta_args, "db")
    delta_args[db_idx] = delta_dbs
    if is_delta:
        _, out_idx = normalize_args(delta_args, "out")
        delta_args[out_idx] += "-delta"

    try:
        subprocess.run(delta_args, check=True)
    except subprocess.CalledProcessError:
        raise ValueError(f"Running Program Failed: {delta_args}")


def deconstruct_call(program_call: str, path: str) -> tuple[str, str, str, str, float, list[str]]:
    # TODO add verification that needed kwargs are in the call
    program_args = shlex.split(program_call)
    program = program_args[0]
    query, _ = normalize_args(program_args, "query")
    db, _ = normalize_args(program_args, "db")
    out, out_idx = normalize_args(program_args, "out")
    out = path + out
    program_args[out_idx] = out
    evalue = get_needed_defaults(program_args)

    return program, query, db, out, evalue, program_args


def calculate_updated_evalues(hits: list[HitData], total_residues: int, partial_residues: int) -> list[HitData]:
    updates: list[HitData] = []
    for hit in hits:
        scaling_factor = (1.0 * total_residues) / (1.0 * partial_residues)
        evalue = cast(float, hit[HitIdx.EVALUE]) * scaling_factor  # float for type checker
        updated_hit = (hit[HitIdx.QUERY], hit[HitIdx.HIT], evalue, hit[HitIdx.THE_REST])
        updates.append(updated_hit)  # type: ignore

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


def write_updated_output(query_file_name: str, out_file_name: str, evalue_cutoff: float, db_client) -> None:
    out_file_name = out_file_name.split(".")[0] + ".m8"
    with open(query_file_name, "r") as query_f, open(out_file_name, "w") as out_f:
        to_print = ""
        for line in query_f:
            if line.startswith(">"):
                query = line.strip()[1:]
                # values = query, hit, evalue
                for hit in db_client.get_query(query):
                    evalue = hit[HitIdx.EVALUE]
                    if evalue <= evalue_cutoff:
                        the_rest = hit[HitIdx.THE_REST].split("\t")
                        all_values = [query, hit[HitIdx.HIT]] + the_rest[:-1] + [f"{evalue:.3}"] + [the_rest[-1]]
                        to_print += "\t".join(all_values) + "\n"
        out_f.write(to_print)


def main(program_call: str, path: str) -> None:
    if path and path[-1] != "/":
        path += "/"
    db_client = IevalueDB(path)
    program, query, db, out, evalue_cutoff, program_args = deconstruct_call(program_call, path)


    prev_dbs, prev_residues = [], 0
    if prev_data := db_client.get_db_info():
        prev_data = cast(list[DatabaseData], prev_data)  # for the type checker
        prev_dbs = cast(list[Optional[str]], [db_data[DbIdx.DATABASE] for db_data in prev_data])
        prev_residues = int(sum([db_data[DbIdx.RESIDUE] for db_data in prev_data]))

    delta_dbs = ""
    delta_parts: list[str] = []
    if program == "blastp":
        # dbs that have not been searched on yet
        delta_dbs, delta_parts = get_blast_delta_db(db, program, out, prev_dbs)
        if not delta_dbs:
            return
    else:
        delta_dbs = db
        delta_parts = [
            db,
        ]

    delta_data = get_delta_sizes(delta_parts, program)

    # add the size of the residues for each database
    db_client.add_database_record(delta_data)

    if prev_residues:
        # update evalues from the old hits
        delta_residue = cast(int, sum([delta_db[DbIdx.RESIDUE] for delta_db in delta_data]))  # int for type checker
        total_residues = prev_residues + delta_residue
        db_client.update_old_evalues(total_residues, prev_residues)

        # perform blast on just the new dbs
        run_program(program_args, delta_dbs, is_delta=True)

        # get the new results and add them to the database
        delta_hits = parse_delta_db(out, is_delta=True)
        updated_delta_hits = calculate_updated_evalues(delta_hits, total_residues, delta_residue)
        db_client.insert_hits(updated_delta_hits)

        write_updated_output(query, out, evalue_cutoff, db_client)

    else:
        run_program(program_args, delta_dbs, is_delta=False)
        delta_hits = parse_delta_db(out, is_delta=False)
        db_client.insert_hits(delta_hits)


if __name__ == "__main__":
    if len(sys.argv) >= 2:
        _path = ""
        if len(sys.argv) == 3:
            _path = sys.argv[2]
        elif len(sys.argv) > 3:
            print("Error: Too many input arguments")
        main(sys.argv[1], _path)
    else:
        raise ValueError("Error: Need program kwargs")
