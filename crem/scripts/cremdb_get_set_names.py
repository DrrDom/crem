#!/usr/bin/env python3
import argparse
import re
import sqlite3


def get_set_names(db_path):
    conn = sqlite3.connect(db_path)

    tables = [
        row[0] for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'radius%'"
        )
    ]
    tables.sort(key=lambda t: int(re.search(r'\d+', t).group()))

    reserved = {'env_id', 'core_smi_id'}
    result = {}

    for table in tables:
        cols = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
        result[table] = sorted(cols - reserved)

    conn.close()
    return result


def main():
    parser = argparse.ArgumentParser(description="List set names from each radius table in a CReM database")
    parser.add_argument("-i", "--input", required=True, help="Path to the SQLite database")
    args = parser.parse_args()

    for table, names in get_set_names(args.input).items():
        print(f"{table}: {names}")


if __name__ == "__main__":
    main()
