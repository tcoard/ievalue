import sqlite3
import os.path
from typing import Optional, Union
from ievalue_types import DatabaseData, HitData

DATABASE = "ievalue_metadata.db"


def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        # we could do something here like if the first value is x, set this as a particular data type
        d[col[0]] = row[idx]
    return d


class IevalueDB:
    def __init__(self) -> None:
        self._con = sqlite3.connect(DATABASE)
        self._con.row_factory = dict_factory
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
        return self._cur.execute(
            "SELECT hit, evalue, the_rest FROM hits WHERE query = :query ORDER BY evalue", {"query": query}
        ).fetchall()

    def insert_hits(self, hits: list[HitData]) -> None:
        self._cur.executemany("INSERT INTO hits VALUES (:query, :hit, :evalue, :the_rest)", hits)
        self._con.commit()

    def update_evalues(self, updates: list[HitData]) -> None:
        # we shouldn't have to worry about collision
        # because samples shouldn't be run more than once
        # but this should be kept in mind
        self._cur.executemany(
            """UPDATE hits
            SET evalue = :evalue
            WHERE query = :query and hit = :hit and the_rest = :the_rest;
            """,
            updates,
        )
        self._con.commit()

    def get_db_info(self) -> list[Optional[DatabaseData]]:
        return self._cur.execute("SELECT * FROM databases").fetchall()

    def add_database_record(self, new_databases: list[DatabaseData]) -> None:
        self._cur.executemany("INSERT INTO databases VALUES (:database_name, :residue)", new_databases)
        self._con.commit()

    def __del__(self) -> None:
        self._con.close()
