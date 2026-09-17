#!/usr/bin/env python3
"""
Print the schema version, fragment sets and property columns of CReM databases.

Several databases may be given at once; each is reported as a block headed by the path
as it was typed. Only PRAGMAs and sqlite_master are read, so the command is instant on
databases of any size, and every file is opened read-only.

Usage:
    cremdb_info -i fragments.db
    cremdb_info -i *.db --json

This command supersedes cremdb_get_set_names, which reports set names only and will be
removed in a future release.
"""

import argparse
import json
import sqlite3
import sys

from crem.db import get_db_info
from crem.scripts.cremdb_create import DB_SCHEMA_VERSION

# Exceptions that mean "this input is unusable" rather than "the tool is broken": with
# several databases on the command line, one bad file must not hide the others.
_INPUT_ERRORS = (OSError, ValueError, sqlite3.Error)


def get_set_names(db_path):
    """Fragment set names per radius table: ``{'radius1': ['chembl'], ...}``.

    Kept as the shape returned by the deprecated cremdb_get_set_names, which imports it
    from here. v0 databases have no sets, so their lists are empty.
    """
    return dict(get_db_info(db_path)['radius_tables'])


def version_label(version):
    """Human-readable schema version, e.g. ``2 (current)``."""
    if version == 0:
        return "0 (legacy) - no fragment sets"
    if version == 1:
        return "1 (deprecated, support may be dropped in a future release)"
    if version == DB_SCHEMA_VERSION:
        return f"{version} (current)"
    return f"{version} (unknown - newer than this CReM, which writes v{DB_SCHEMA_VERSION})"


def format_info(info):
    """Render one database's info as an indented block, without the path heading."""
    rows = [("schema version", version_label(info['version']))]

    if info['has_sets']:
        for table, set_names in info['radius_tables'].items():
            rows.append((table, ", ".join(set_names) if set_names else "(none)"))
    else:
        # v0 radius tables hold no set columns, so listing them by name is all there is
        # to say; one row keeps the block short when a database has many radii.
        rows.append(("radius tables", ", ".join(info['radius_tables'])))

    # Only tables that actually carry property columns are listed; a freshly built
    # database has none anywhere, and a line of "(none)" per table would be noise.
    with_props = [(table, cols) for table, cols in info['properties'].items() if cols]
    rows.append(("properties",
                 "\n".join(f"{table}: {', '.join(cols)}" for table, cols in with_props)
                 if with_props else "(none)"))

    if info['is_shard']:
        rows.append(("note", "unmerged stride shard (cremdb_create --parallel-shards)"))

    width = max(len(label) for label, _ in rows)
    lines = []
    for label, value in rows:
        head, *rest = value.split("\n")
        lines.append(f"  {label.ljust(width)} : {head}")
        # continuation lines (one per property table) line up under the first value
        lines.extend(f"  {' ' * width}   {line}" for line in rest)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Print the schema version, fragment sets and property columns of "
                    "one or more CReM fragment databases.")
    parser.add_argument("-i", "--input", metavar="FILENAME", required=True, nargs="+",
                        help="one or more SQLite DBs with CReM fragments")
    parser.add_argument("--json", action="store_true",
                        help="print the same information as a JSON list")
    args = parser.parse_args()

    blocks, records, failed = [], [], False

    for db_path in args.input:
        try:
            info = get_db_info(db_path)
        except _INPUT_ERRORS as e:
            failed = True
            records.append({'path': db_path, 'error': str(e)})
            blocks.append(f"{db_path}\n  error: {e}")
        else:
            records.append(info)
            blocks.append(f"{db_path}\n{format_info(info)}")

    if args.json:
        print(json.dumps(records, indent=2))
    else:
        print("\n\n".join(blocks))

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
