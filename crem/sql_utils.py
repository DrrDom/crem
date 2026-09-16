"""SQL helpers shared by the build, merge, convert and query paths.

Deliberately free of imports so that the scripts which keep their module import
RDKit-free (cremdb_merge, cremdb_add_prop) can use it at the top level.
"""


def quote_ident(name):
    """Return `name` as a quoted SQLite identifier.

    Set names and fragment property names become columns of the radius{N} / frags
    tables, so they reach SQL as identifiers rather than as bound values, and an
    unquoted one that happens to be a keyword ('all', 'order', 'index', ...) is a
    syntax error. SQLite keeps adding keywords, so every interpolated identifier is
    quoted instead of being checked against a keyword list.

    Quoting is purely syntactic - "all" and all name the same column - so databases
    built before this was introduced are read and extended unchanged.
    """
    return '"' + str(name).replace('"', '""') + '"'
