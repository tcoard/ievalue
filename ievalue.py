import sys
import shlex
import subprocess
from itertools import chain
from ievalue_db import IevalueDB
from Bio.Blast import NCBIXML

class NoNewData(Exception):
    pass

def read_database_parts(db_name: str) -> list[str]:
    # TODO this code is from iblast, look into optimizing it
    try:
        result = subprocess.run(
            ["blastdbcmd", "-db", db_name, "-info", "-exact_length"], stdout=subprocess.PIPE, check=True
        )
    except subprocess.CalledProcessError:
        # TODO
        return []

    db_info = result.stdout.decode("utf-8").split("\n")

    start_capture = False
    parts = []
    for line in db_info:
        if line != "":
            if line == "Volumes:":
                start_capture = True
            elif start_capture:
                parts.append(line.strip().split("/")[-1])
    return parts


def get_delta_db(blast_database: str, blast_type: str, out_file_name: str, prev_dbs: list[str]) -> tuple[str, list[dict[str,str]]]:

    curr_dbs = read_database_parts(blast_database)
    delta_parts = list(set(curr_dbs).difference(set(prev_dbs)))

    if len(delta_parts) == 0:
        raise NoNewData("There are no new databases to query")

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

    delta_parts = [{"database_name": i for i in delta_parts}]
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


def update_old_evalues(db) -> None:
    updates = []
    for query, hit, evalue in db.get_all_hits():
        evalue = str(float(evalue) - 1)  # placeholder for actual evalue update
        updates.append({"query": query, "hit": hit, "evalue": evalue})
    db.update_evalues(updates)


def parse_delta_db(out_file_name: str, db) -> None:
    out_file_name += "-delta"
    hits = []
    with open(out_file_name, "r") as f:
        blast_record_iter = NCBIXML.parse(f)
        for record in blast_record_iter:
            for hit in record.descriptions:
                hits.append({"query": record.query, "hit": hit.title.split(" ")[-1], "evalue": hit.e})
    db.insert_hits(hits)


def write_updated_output(out_file_name: str, db) -> None:
    out_file_name = out_file_name.split('.')[0] + ".tsv"
    with open(out_file_name, 'w') as f:
         # values = query, hit, evalue
        for values in db.get_all_hits():
            print('\t'.join(values), file=f)

def main(program_call: str) -> None:
    db = IevalueDB()

    blast_kwargs = deconstruct_call(program_call)
    # TODO only search sequences that have not been searched on the exact previous databases
    prev_dbs = db.get_prev_dbs()

    # dbs that have not been searched on yet
    try:
        delta_dbs, delta_parts = get_delta_db(blast_kwargs["db"], blast_kwargs["blast_type"], blast_kwargs["out"], prev_dbs)
    except NoNewData as e:
        print(e)
        return
    db.add_database_record(delta_parts)

    # update evalues from the old hits
    update_old_evalues(db)

    # perform blast on just the new dbs
    do_delta_blast(blast_kwargs, delta_dbs)

    # get the new results and add them to the database
    parse_delta_db(blast_kwargs["out"], db)

    write_updated_output(blast_kwargs["out"], db)


if __name__ == "__main__":
    print(sys.argv[1])
    main(sys.argv[1])
