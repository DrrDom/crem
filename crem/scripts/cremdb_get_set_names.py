#!/usr/bin/env python3
"""
Deprecated. Use cremdb_info, which reports the schema version and property columns as
well, and accepts several databases at once.

This script is a thin wrapper kept so that existing pipelines keep working: its stdout
format is unchanged, and get_set_names remains importable from here. It will be removed
in a future CReM release.
"""

import argparse
import sys

# re-exported for code that imports get_set_names from this module
from crem.scripts.cremdb_info import get_set_names

_DEPRECATION_MSG = (
    "cremdb_get_set_names is deprecated and will be removed in a future CReM release. "
    "Use cremdb_info instead - it reports the database version and property columns too, "
    "and accepts several databases at once.\n"
)


def main():
    parser = argparse.ArgumentParser(
        description="Deprecated, use cremdb_info. List set names from each radius table "
                    "in a CReM database")
    parser.add_argument("-i", "--input", required=True, help="Path to the SQLite database")
    args = parser.parse_args()

    # stderr rather than warnings.warn: DeprecationWarning is suppressed by Python's
    # default filters outside __main__, and console entry points are wrapper modules, so
    # the warning would be invisible to most users. Writing to stderr also leaves stdout
    # byte-identical for pipelines that parse it.
    sys.stderr.write(_DEPRECATION_MSG)

    for table, names in get_set_names(args.input).items():
        print(f"{table}: {names}")


if __name__ == "__main__":
    main()
