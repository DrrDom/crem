# Pavel Polishchuk, 2017

from subprocess import Popen, PIPE
import os
from random import shuffle


def divide(n, iterable):
    # https://github.com/erikrose/more-itertools/blob/dfde6c75eb520605ac6a2d374bcdb72e7ce2d894/more_itertools/more.py
    """Divide the elements from *iterable* into *n* parts, maintaining
    order.
        >>> group_1, group_2 = divide(2, [1, 2, 3, 4, 5, 6])
        >>> list(group_1)
        [1, 2, 3]
        >>> list(group_2)
        [4, 5, 6]
    If the length of *iterable* is not evenly divisible by *n*, then the
    length of the returned iterables will not be identical:
        >>> children = divide(3, [1, 2, 3, 4, 5, 6, 7])
        >>> [list(c) for c in children]
        [[1, 2, 3], [4, 5], [6, 7]]
    If the length of the iterable is smaller than n, then the last returned
    iterables will be empty:
        >>> children = divide(5, [1, 2, 3])
        >>> [list(c) for c in children]
        [[1], [2], [3], [], []]
    This function will exhaust the iterable before returning and may require
    significant storage. If order is not important, see :func:`distribute`,
    which does not first pull the iterable into memory.
    """
    if n < 1:
        raise ValueError('n must be at least 1')

    seq = tuple(iterable)
    q, r = divmod(len(seq), n)

    ret = []
    for i in range(n):
        start = (i * q) + (i if i < r else r)
        stop = ((i + 1) * q) + (i + 1 if i + 1 < r else r)
        ret.append(iter(seq[start:stop]))

    return ret


if __name__ == '__main__':

    nnodes = 2

    py_dir = '/home/pavlop/python/crem/pbs'
    # base_dir = '/home/pavlop/imtm/crem/generation'
    # dirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    # dirs = ['/home/pavlop/imtm/crem/generation/orgelm_sa3_benzene_r1_freq0']

    # dirs = []

    base_dir = '/home/pavlop/imtm/crem/drugbank/pains'
    dirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

    # for d in os.listdir(base_dir):
    #     if os.path.isdir(os.path.join(base_dir, d)):
    #         for dd in os.listdir(os.path.join(base_dir, d)):
    #             if os.path.isdir(os.path.join(base_dir, d, dd)):
    #                 dirs.append(os.path.join(base_dir, d, dd))

    for d in dirs:

        print(d)

        job_dir = os.path.join(d, 'jobs')
        if not os.path.exists(job_dir):
            os.mkdir(job_dir)

        files = [f for f in os.listdir(d) if f.endswith('.txt')]
        shuffle(files)
        files = divide(nnodes, files)

        for n in range(nnodes):

            pbs_name = os.path.join(job_dir, 'sim_%i.pbs' % n)

            script = """
            #!/bin/bash
            #PBS -l select=1:ncpus=32
            #PBS -k eo
            source activate rdkit-1709
            cd %s
            """ % (d, )

            for f in files[n]:
                script += 'python3 %s/similarity2.py %s > %s\n' % \
                          (py_dir, f, f.replace('.txt', '.diversity'))

            with open(pbs_name, "wt") as f:
                f.write(script)
            p = Popen(['qsub', pbs_name], stdout=PIPE, encoding='utf8')
            outs, errs = p.communicate(timeout=30)
            id = outs.strip()
            print(id, pbs_name)
