import os
import sys
import datetime
from multiprocessing import Pool, current_process


def fib(n):
    if n == 0: return 0
    elif n == 1: return 1
    else: return fib(n-1)+fib(n-2)


def calc(n):
    res = fib(n)
    return str(datetime.datetime.now()), str(os.getpid()), str(current_process()), n, res


def get_data(n=100):
    for _ in range(n):
        yield 35


ncpu = int(sys.argv[1])
fname = sys.argv[2]

p = Pool(ncpu)

with open(fname, 'wt') as f:
    f.write(str(os.getpid()) + '\n')
    for res in p.imap(calc, get_data()):
        f.write("\t".join(map(str, res)) + '\n')
        f.flush()
