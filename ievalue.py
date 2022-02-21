import shlex
import subprocess
import argparse
from collections import defaultdict
from typing import Optional, cast, Any, Union
from ievalue_db import IevalueDB, DatabaseData, HitData, DbIdx, HitIdx

# TODO:
# * change ValueError to something better


def get_needed_defaults(program_args: list[str]) -> tuple[float, int]:
    program = program_args[0]
    needed_defaults: dict[str, Any] = {}
    defaults: dict[str, dict[str, Any]] = {
        "blast": {"evalue_cutoff": 10.0, "max_target_seqs": 500},
        "mmseqs": {"evalue_cutoff": 1.000e-03, "max_target_seqs": 300},
        "diamond": {"evalue_cutoff": 0.001, "max_target_seqs": 25},
    }
    for needed_arg in defaults[program]:
        val, _ = normalize_args(program_args, needed_arg)
        if val:
            needed_defaults[needed_arg] = type(defaults[program][needed_arg])(val)
        else:
            needed_defaults[needed_arg] = defaults[program][needed_arg]

    return needed_defaults["evalue_cutoff"], needed_defaults["max_target_seqs"]


def normalize_args(program_args: list[str], arg: str) -> tuple[Optional[str], Optional[int]]:
    # this will be a normalization layer from application specific keywords to normalized keyword
    program = program_args[0]
    needed = ["db", "out", "query"]
    normalized_keywords: dict[str, dict[str, Union[str, list[str], int]]] = {
        "blastp": {
            "db": "-db",
            "out": "-out",
            "query": "-query",
            "evalue_cutoff": "-evalue",
            "max_target_seqs": "-max_target_seqs",
        },
        "diamond": {
            "db": ["--db", "-d"],
            "out": ["--out", "-o"],
            "query": ["--query", "-q"],
            "evalue_cutoff": ["--evalue", "-e"],
            "max_target_seqs": ["max-target-seqs", "-k"],
        },
        "mmseqs": {
            "db": 3,
            "out": 4,
            "query": 2,
            "evalue_cutoff": ["--evalue", "-e"],
            "max_target_seqs": "--max-seqs",
        },
    }
    kwarg = normalized_keywords[program][arg]
    idx = None
    value = None
    if isinstance(kwarg, str):
        try:
            idx = program_args.index(kwarg) + 1
            value = program_args[idx]
        except ValueError:
            if arg in needed:
                raise ValueError(f"Could not find needed input: {arg}")
    elif isinstance(kwarg, list):
        for variation in kwarg:
            try:
                idx = program_args.index(variation) + 1
                value = program_args[idx]
                break
            except ValueError:
                pass
        else:
            if arg in needed:
                raise ValueError(f"Could not find needed input: {arg}")
    elif isinstance(kwarg, int):
        idx = kwarg
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


def deconstruct_call(program_call: str, path: str) -> tuple[str, str, str, str, float, int, list[str]]:
    program_args = shlex.split(program_call)
    program = program_args[0]
    query, _ = cast(tuple[str, int], normalize_args(program_args, "query"))
    db, _ = cast(tuple[str, int], normalize_args(program_args, "db"))
    out, out_idx = cast(tuple[str, int], normalize_args(program_args, "out"))
    out = path + out
    program_args[out_idx] = out
    evalue, max_seqs = get_needed_defaults(program_args)

    return program, query, db, out, evalue, max_seqs, program_args


def calculate_updated_evalues(
    hits: list[HitData], total_residues: int, partial_residues: int, evalue_cutoff: float
) -> list[HitData]:
    updates = defaultdict(list)
    scaling_factor = (1.0 * total_residues) / (1.0 * partial_residues)
    print('aaa', scaling_factor)
    for hit in hits:
        evalue = cast(float, hit[HitIdx.EVALUE]) * scaling_factor  # float for type checker
        if evalue <= evalue_cutoff:
            updated_hit = (hit[HitIdx.QUERY], hit[HitIdx.HIT], evalue, hit[HitIdx.THE_REST])
            updates[str(hit[HitIdx.QUERY])].append(cast(tuple[str, str, float, str], updated_hit))

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


def write_updated_output(
    query_file_name: str, out_file_name: str, evalue_cutoff: float, max_seqs: int, db_client, path, delta_data, updated_delta_hits, prev_data
) -> None:
    out_file_name = out_file_name.split(".")[0] + ".m8"
    hits_to_keep = []
    with open(query_file_name, "r") as query_f, open(out_file_name, "w") as out_f:
        to_print = ""
        for line in query_f:
            if line.startswith(">"):
                query = line.strip()[1:]
                # values = query, hit, evalue
                # if query == "WP_013226005.1":
                #     print(db_client.get_query(query, max_seqs))
                hits = updated_delta_hits[query]
                hits.extend(db_client.get_query(query))
                hits = sorted(hits, key=lambda x: x[2])
                for hit in hits[:max_seqs]:
                    evalue = hit[HitIdx.EVALUE]
                    if evalue <= evalue_cutoff:
                        hits_to_keep.append(hit)
                        the_rest = hit[HitIdx.THE_REST].split("\t")
                        all_values = [query, hit[HitIdx.HIT]] + the_rest[:-1] + [f"{evalue:.3}"] + [the_rest[-1]]
                        to_print += "\t".join(all_values) + "\n"
        out_f.write(to_print)

    db_client.close()
    db_client.del_old()
    new_db_client = IevalueDB(path)
    new_db_client.insert_hits(hits_to_keep)
    new_db_client.add_database_record(prev_data)
    new_db_client.add_database_record(delta_data)
    new_db_client.close()


def main(program_call: str, path: str, clean_db: bool) -> None:
    if path and path[-1] != "/":
        path += "/"
    db_client = IevalueDB(path)
    program, query, db, out, evalue_cutoff, max_seqs, program_args = deconstruct_call(program_call, path)

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
    print(db_client.get_db_info(), delta_data)

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
        updated_delta_hits = calculate_updated_evalues(delta_hits, total_residues, delta_residue, evalue_cutoff)
        # db_client.insert_hits(updated_delta_hits)
        if clean_db:
            db_client.clean_db(evalue_cutoff, max_seqs)
        write_updated_output(query, out, evalue_cutoff, max_seqs, db_client, path, delta_data, updated_delta_hits, prev_data)

    else:
        run_program(program_args, delta_dbs, is_delta=False)
        delta_hits = parse_delta_db(out, is_delta=False)
        db_client.insert_hits(delta_hits)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("program_call", type=str)  # TODO, required=True)
    parser.add_argument("--path", type=str, default="")
    parser.add_argument("--clean_db", action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    main(args.program_call, args.path, args.clean_db)
