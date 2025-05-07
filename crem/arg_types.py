import os
from multiprocessing import cpu_count


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def filepath_type(x):
    if x:
        return os.path.abspath(x)
    else:
        return x


def str_lower_type(x):
    if x:
        return x.lower()
    else:
        return x


def similarity_value_type(x):
    return max(0, min(1, float(x)))
