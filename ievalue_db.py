import sqlite3
import os
from enum import IntEnum
from typing import Union


# TODO save previous database, so that if things mess up part way through, we can revert to that
# maybe even automatically revert

DatabaseData = Union[tuple[str, int], tuple[str, int, str]]

HitData = tuple[str, str, float, str]


# These are so that we can save time by returning the default list of tuples from sqlite,
# but won't accidentally call the wrong index.
class DbIdx(IntEnum):
    DATABASE = 0
    RESIDUE = 1
    FASTAS = 1


class HitIdx(IntEnum):
    QUERY = 0
    HIT = 1
    EVALUE = 2
    THE_REST = 3
    HIT_PK = 4


class IevalueDB:
    def __init__(self, db_path: str) -> None:
        database = db_path + "ievalue_metadata.db"
        self.database = database
        self._con = sqlite3.connect(database)
        self._con.execute("PRAGMA synchronous = OFF")
        self._con.execute("PRAGMA journal_mode = OFF")
        self._cur = self._con.cursor()
        # If the database does not exist, initialize it
        if os.path.getsize(database) == 0:
            print(f"Creating initial database: {database}")
            self._create_db()

    def _create_db(self) -> None:
        self._cur.execute(
            """CREATE TABLE databases
                       (database text NOT NULL, residue int NOT NULL)"""
        )
        self._cur.execute(
            """CREATE TABLE hits
                       (query text NOT NULL, hit text NOT NULL, evalue float NOT NULL, the_rest TEXT NOT NULL)"""  # WITHOUT ROWID"""
        )
        self._cur.execute(
            """CREATE INDEX hits_idx on hits (query)"""
        )
        self._con.commit()
        # (query text NOT NULL, hit text NOT NULL, evalue float NOT NULL, the_rest TEXT PRIMAY KEY)"""# WITHOUT ROWID"""
        # ;(query text NOT NULL, hit text NOT NULL, evalue float NOT NULL, the_rest text NOT NULL, CONSTRAINT hit_pk PRIMARY KEY (query, hit, the_rest)) WITHOUT ROWID"""

    def get_all_hits(self) -> list[HitData]:
        return self._cur.execute("SELECT * from hits").fetchall()

    def get_query(self, query: str) -> list[HitData]:
        return self._cur.execute("SELECT * FROM hits WHERE query = ?", (query,)).fetchall()

    def insert_hits(self, hits: list[HitData]) -> None:
        self._cur.executemany("INSERT INTO hits VALUES (?, ?, ?, ?)", hits)
        self._con.commit()

    def update_old_evalues(self, total_residues, prev_residues) -> None:
        # we shouldn't have to worry about collision
        # because samples shouldn't be run more than once
        # but this should be kept in mind
        scaling_factor = 1.0 * total_residues / prev_residues
        self._cur.execute("UPDATE hits SET evalue = evalue * ?", (scaling_factor,))
        self._con.commit()

    def get_db_info(self) -> Union[list[None], list[DatabaseData]]:
        return self._cur.execute("SELECT * FROM databases").fetchall()

    def add_database_record(self, new_databases: list[DatabaseData]) -> None:
        self._cur.executemany("INSERT INTO databases VALUES (?, ?)", new_databases)
        self._con.commit()

    def clean_db(self, evalue_cutoff, max_seqs):
        # self._cur.execute("SELECT * FROM hits where (query, hit) in (SELECT query, hit FROM hits ORDER by evalue LIMIT 10)").fetchall()

        clean_evalue_command = "DELETE FROM hits where evalue > ?"
        self._cur.execute(clean_evalue_command, (evalue_cutoff,))

        clean_low_hits_command = """
            DELETE
            FROM hits
            where (query, hit) in (
                SELECT rs.query, rs.hit
                FROM (
                    SELECT query, hit,
                    Rank() over (
                        Partition BY query
                        ORDER BY evalue ASC
                    )
                    AS Rank FROM hits
                )
                rs WHERE Rank > ?
            )
        """
        self._cur.execute(clean_low_hits_command, (max_seqs,))
        self._con.commit()

        # breakpoint()
        #data = self._cur.execute("SELECT * FROM hits WHERE hit IN (SELECT * FROM hits ORDER by evalue) LIMIT 10")

        # self._cur.executemany("DELETE FROM Table hits WHERE query NOT IN (SELECT TOP 10 ID FROM Table)")
        # self._con.commit()

    def del_old(self):
        os.remove(self.database)

    def close(self) -> None:
        self._con.close()

    def __del__(self) -> None:
        self._con.close()
