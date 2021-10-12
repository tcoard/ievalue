import sqlite3
import os.path

DATABASE = "ievalue_metadata.db"


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
                       (database text)"""
        )
        self._cur.execute(
            """CREATE TABLE hits
                       (query text, hit text, evalue text)"""
        )
        self._con.commit()

    def get_all_hits(self) -> list:
        return self._cur.execute("SELECT * from hits").fetchall()

    def insert_hits(self, hits: list[dict[str, str]]) -> None:
        # we shouldn't have to worry about collision
        # because samples shouldn't be run more than once
        # but this should be kept in mind
        self._cur.executemany("INSERT INTO hits VALUES (:query, :hit, :evalue)", hits)
        self._con.commit()

    def update_evalues(self, updates: list[dict[str, str]]) -> None:
        # we shouldn't have to worry about collision
        # because samples shouldn't be run more than once
        # but this should be kept in mind
        self._cur.executemany(
            """UPDATE hits
            SET evalue = :evalue
            WHERE query = :query and hit = :hit;
            """,
            updates,
        )
        self._con.commit()

    def get_prev_dbs(self) -> list[str]:
        results = self._cur.execute("SELECT * FROM databases").fetchall()
        dbs = [db_tup[0] for db_tup in results]
        return dbs

    def add_database_record(self, new_databases: list[dict[str, str]]) -> None:
        self._cur.executemany("INSERT INTO databases VALUES (:database_name)", new_databases)
        self._con.commit()

    def __del__(self) -> None:
        self._con.close()
