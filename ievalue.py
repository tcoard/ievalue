import shlex
import os
import subprocess
import argparse
import configparser
from collections import defaultdict
from typing import Optional, cast, Any
from ievalue_db import (
    IevalueDB,
    DatabaseData,
    HitData,
    DATABASE_IDX,
    RESIDUE_IDX,
    QUERY_IDX,
    HIT_IDX,
    EVALUE_IDX,
    THE_REST_IDX,
)

# TODO:
# * change ValueError to something better


def get_needed_defaults(program_args: list[str], default_values, normalized_keywords) -> tuple[float, int]:
    needed_defaults: dict[str, Any] = {}
    # defaults: dict[str, dict[str, Any]] = {
    #     "blastp": {"evalue_cutoff": 10.0, "max_target_seqs": 500},
    #     "mmseqs": {"evalue_cutoff": 1.000e-03, "max_target_seqs": 300},
    #     "diamond": {"evalue_cutoff": 0.001, "max_target_seqs": 25},
    # }
    for needed_arg in default_values:
        val, _ = normalize_args(program_args, needed_arg, normalized_keywords)
        if val:
            needed_val = val
        else:
            needed_val = default_values[needed_arg]

        if needed_val.isnumeric():
            if "." in default_values[needed_arg]:
                needed_val = float(needed_val)  # type: ignore
            else:
                needed_val = int(needed_val)  # type: ignore

        needed_defaults[needed_arg] = needed_val

    return needed_defaults["evalue_cutoff"], needed_defaults["max_target_seqs"]


def normalize_args(program_args: list[str], arg: str, normalized_keywords) -> tuple[Optional[str], Optional[int]]:
    # this will be a normalization layer from application specific keywords to normalized keyword
    needed = ["db", "out", "query"]
    kwarg = normalized_keywords[arg]
    if kwarg.isnumeric():
        kwarg = int(kwarg)
    else:
        # if there is only one element, this will make a list of one element
        kwarg = kwarg.split(",")

    idx = None
    value = None
    if isinstance(kwarg, list):
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


def get_delta_sizes(dbs: list[str], search_type: str, database_parts_function) -> list[DatabaseData]:
    db_sizes: list[DatabaseData] = []
    for db in dbs:
        db_residue, _ = read_database_parts(db, search_type, database_parts_function)
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


def read_database_parts(db_name: str, search_type: str, database_parts_function) -> tuple[int, list[str]]:
    db_residue = 0
    parts: list[str] = []
    db_residue, parts = database_parts_function(db_name, search_type)

    if db_residue == 0:
        raise ValueError("Error finding database information")

    return db_residue, parts


def get_blast_delta_db(
    blast_database: str, search_type: str, out_file_name: str, prev_dbs: list[Optional[str]], database_parts_function
) -> tuple[str, list[str]]:
    _, curr_dbs = read_database_parts(blast_database, search_type, database_parts_function)
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
        "prot",  # hard coding until we support nucleotides
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


def run_program(program_args, delta_dbs: str, is_delta: bool, normalized_keywords) -> None:
    delta_args = program_args.copy()

    _, db_idx = normalize_args(delta_args, "db", normalized_keywords)
    delta_args[db_idx] = delta_dbs
    if is_delta:
        _, out_idx = normalize_args(delta_args, "out", normalized_keywords)
        delta_args[out_idx] += "-delta"

    try:
        subprocess.run(delta_args, check=True)
    except subprocess.CalledProcessError:
        raise ValueError(f"Running Program Failed: {delta_args}")


def deconstruct_call(
    program_call: str, path: str, default_values, normalized_keywords
) -> tuple[str, str, str, str, float, int, list[str]]:
    program_args = shlex.split(program_call)
    program = program_args[0]
    query, _ = cast(tuple[str, int], normalize_args(program_args, "query", normalized_keywords))
    db, _ = cast(tuple[str, int], normalize_args(program_args, "db", normalized_keywords))
    out, out_idx = cast(tuple[str, int], normalize_args(program_args, "out", normalized_keywords))
    out = path + out
    program_args[out_idx] = out
    evalue, max_seqs = get_needed_defaults(program_args, default_values, normalized_keywords)

    return program, query, db, out, evalue, max_seqs, program_args


def calculate_updated_evalues(
    hits: list[HitData], total_residues: int, partial_residues: int, evalue_cutoff: float
) -> dict[str, list[HitData]]:
    updates = defaultdict(list)
    scaling_factor = (1.0 * total_residues) / (1.0 * partial_residues)
    for hit in hits:
        evalue = hit[EVALUE_IDX] * scaling_factor
        if evalue <= evalue_cutoff:
            updated_hit = (hit[QUERY_IDX], hit[HIT_IDX], evalue, hit[THE_REST_IDX])
            updates[hit[QUERY_IDX]].append(updated_hit)

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
    query_file_name: str,
    out_file_name: str,
    evalue_cutoff: float,
    max_seqs: int,
    db_client: IevalueDB,
    updated_delta_hits: dict[str, list[HitData]],
) -> list[HitData]:
    out_file_name = out_file_name.split(".")[0] + ".m8"
    hits_to_keep = []
    with open(query_file_name, "r") as query_f, open(out_file_name, "w") as out_f:
        to_print = ""
        for line in query_f:
            if line.startswith(">"):
                query = line.strip()[1:]

                # combine the new and old hits and keep the top <max_seqs>
                # that have an evalue_cutoff <= <evalue_cutoff>
                hits = updated_delta_hits[query]
                previous_hits = db_client.get_query(query)
                hits.extend(previous_hits)
                hits = sorted(hits, key=lambda x: x[EVALUE_IDX])
                for hit in hits[:max_seqs]:
                    evalue = hit[EVALUE_IDX]
                    if evalue <= evalue_cutoff:
                        hits_to_keep.append(hit)
                        # format for the text file
                        the_rest = hit[THE_REST_IDX].split("\t")
                        all_values = [query, hit[HIT_IDX]] + the_rest[:-1] + [f"{evalue:.3}"] + [the_rest[-1]]
                        to_print += "\t".join(all_values) + "\n"
        out_f.write(to_print)

    return hits_to_keep


def main(program_call: str, path: str) -> None:
    config = configparser.ConfigParser()
    # move this later
    program = program_call.split()[0]
    config.read(f"{os.path.abspath(os.path.dirname(__file__))}/configs/{program}.ini")
    print(f"{os.path.abspath(os.path.dirname(__file__))}/configs/{program}.ini")
    print(config.sections())
    default_values = config["DEFAULT_VALUES"]
    normalized_keywords = config["ARG_NAMES_OR_LOC"]
    database_parts_function = config["ARG_NAMES_OR_LOC"]
    database_parts_function = globals()[config["OTHER"]["database_parts_function"]]
    # output_file_type = config["OTHER"]["m8"]

    if path and path[-1] != "/":
        path += "/"
    db_client = IevalueDB(path)
    program, query, db, out, evalue_cutoff, max_seqs, program_args = deconstruct_call(
        program_call, path, default_values, normalized_keywords
    )

    prev_dbs, prev_residues = [], 0
    if prev_data := db_client.get_db_info():
        prev_data = cast(list[DatabaseData], prev_data)  # for the type checker
        prev_dbs = cast(list[Optional[str]], [db_data[DATABASE_IDX] for db_data in prev_data])
        prev_residues = int(sum([db_data[RESIDUE_IDX] for db_data in prev_data]))

    delta_dbs = ""
    delta_parts: list[str] = []
    if program == "blastp":
        # dbs that have not been searched on yet
        delta_dbs, delta_parts = get_blast_delta_db(db, program, out, prev_dbs, database_parts_function)
        if not delta_dbs:
            return
    else:
        delta_dbs = db
        delta_parts = [
            db,
        ]

    delta_data = get_delta_sizes(delta_parts, program, database_parts_function)

    if prev_residues:
        # update evalues from the old hits
        delta_residue = sum([delta_db[RESIDUE_IDX] for delta_db in delta_data])
        total_residues = prev_residues + delta_residue
        db_client.update_old_evalues(total_residues, prev_residues)

        # perform blast on just the new dbs
        run_program(program_args, delta_dbs, is_delta=True, normalized_keywords=normalized_keywords)

        # get the new results and add them to the database
        delta_hits = parse_delta_db(out, is_delta=True)
        updated_delta_hits = calculate_updated_evalues(delta_hits, total_residues, delta_residue, evalue_cutoff)
        hits_to_keep = write_updated_output(query, out, evalue_cutoff, max_seqs, db_client, updated_delta_hits)

        print("Output is ready, cleaning up database")

        db_client.close()
        db_client.del_old()
        new_db_client = IevalueDB(path)
        new_db_client.insert_hits(hits_to_keep)
        new_db_client.add_database_record(prev_data + delta_data)  # type: ignore

    else:
        run_program(program_args, delta_dbs, is_delta=False, normalized_keywords=normalized_keywords)
        delta_hits = parse_delta_db(out, is_delta=False)
        db_client.insert_hits(delta_hits)
        db_client.add_database_record(delta_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("program_call", type=str)
    parser.add_argument("--path", type=str, default="")
    args = parser.parse_args()
    main(args.program_call, args.path)
