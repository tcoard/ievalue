import sqlite3
import os.path
from enum import IntEnum
from typing import Union, Literal

DATABASE = "ievalue_metadata.db"

DatabaseData = tuple[str, int]
HitData = tuple[str, str, float, str]

# These are so that we can save time by returning the default list of tuples from sqlite,
# but won't accidentally call the wrong index.
class DbIdx(IntEnum):
    DATABASE = 0
    RESIDUE = 1


class HitIdx(IntEnum):
    QUERY = 0
    HIT = 1
    EVALUE = 2
    THE_REST = 3


class IevalueDB:
    def __init__(self) -> None:
        self._con = sqlite3.connect(DATABASE)
        self._cur = self._con.cursor()
        # If the database does not exist, initialize it
        if os.path.getsize(DATABASE) == 0:
            print(f"Creating initial database: {DATABASE}")
            self._create_db()

    def _create_db(self) -> None:
        self._cur.execute(
            """CREATE TABLE databases
                       (database text, residue int)"""
        )
        self._cur.execute(
            """CREATE TABLE hits
                       (query text, hit text, evalue float, the_rest text)"""
        )
        self._con.commit()

    def get_all_hits(self) -> list[HitData]:
        return self._cur.execute("SELECT * from hits").fetchall()

    def get_query(self, query: str) -> list[HitData]:
        return self._cur.execute("SELECT * FROM hits WHERE query = ? ORDER BY evalue", (query,)).fetchall()

    def insert_hits(self, hits: list[HitData]) -> None:
        self._cur.executemany("INSERT INTO hits VALUES (?, ?, ?, ?)", hits)
        self._con.commit()

    def update_evalues(self, updates: list[HitData]) -> None:
        # we shouldn't have to worry about collision
        # because samples shouldn't be run more than once
        # but this should be kept in mind
        self._cur.executemany(
            """UPDATE hits
            SET evalue = ?
            WHERE query = ? and hit = ? and the_rest = ?;
            """,
            (
                [
                    [update[HitIdx.EVALUE], update[HitIdx.QUERY], update[HitIdx.HIT], update[HitIdx.THE_REST]]
                    for update in updates
                ]
            ),
        )
        self._con.commit()

    def get_db_info(self) -> Union[list[None], list[DatabaseData]]:
        return self._cur.execute("SELECT * FROM databases").fetchall()

    def add_database_record(self, new_databases: list[DatabaseData]) -> None:
        self._cur.executemany("INSERT INTO databases VALUES (?, ?)", new_databases)
        self._con.commit()

    def __del__(self) -> None:
        self._con.close()
