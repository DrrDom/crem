# Pavel Polishchuk, 2017

import random
from subprocess import Popen, PIPE
import time
import os
import datetime


def chunks(l, n):
    """ Yield n successive chunks from l.
        https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
    """
    if n > len(l):
        n = len(l)
    newn = int(1.0 * len(l) / n + 0.5)
    for i in range(0, n-1):
        yield l[i*newn:i*newn+newn]
    yield l[n*newn-newn:]


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


def not_finished(job_ids):
    p = Popen(['qstat'], stdout=PIPE, encoding='utf8')
    outs, errs = p.communicate(timeout=30)
    lines = outs.split('\n')
    if len(lines) > 2:
        for line in lines[2:]:
            job = [s for s in line.strip().split(' ') if s]  # Job_id Name User Time_Use S Queue
            if job and job[0] in job_ids and job[4] != 'E':
                return True
    return False


def get_jobs(job_ids, unfinished=True):
    unfinished_jobs = []
    p = Popen(['qstat'], stdout=PIPE, encoding='utf8')
    outs, errs = p.communicate(timeout=30)
    lines = outs.split('\n')
    if len(lines) > 2:
        for line in lines[2:]:
            job = [s for s in line.strip().split(' ') if s]  # Job_id Name User Time_Use S Queue
            if job and job[0] in job_ids and job[4] != 'E':
                unfinished_jobs.append(job[0])
    if unfinished:
        return unfinished_jobs
    else:
        return list(set(job_ids).difference(unfinished_jobs))


if __name__ == '__main__':

    nmols = 10
    niter = 150
    # freq = 10

    db_fname = '/home/pavlop/imtm/crem/db/r3_newmolcontext_1709_freq.db'

    output_dir = '/bigdisk1/pavlop/crem/benzene_r1_freq10'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    job_dir = os.path.join(output_dir, 'jobs')
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    stdout_dir = os.path.join(output_dir, 'out')
    if not os.path.exists(stdout_dir):
        os.mkdir(stdout_dir)

    stderr_dir = os.path.join(output_dir, 'err')
    if not os.path.exists(stderr_dir):
        os.mkdir(stderr_dir)

    input_mols = [('c1ccccc1', '0')]

    for i in range(niter):

        new_mols = []

        job_ids = {}   # job_id: output_fname

        for j, (smi, mol_id) in enumerate(reversed(input_mols)):
            pbs_name = os.path.join(job_dir, 'i%s_j%s.pbs' % (str(i).zfill(3), str(j).zfill(2)))
            output_fname = os.path.join(output_dir, 'out_i%s_j%s.txt' % (str(i).zfill(3), str(j).zfill(2)))
            script = """
            #!/bin/bash
            #PBS -l select=1:ncpus=16
            #PDS -d %s
            source activate rdkit-1709
            python3 /home/pavlop/python/crem/pbs/iterative_test4_mutate.py -i "%s" --id %s --iteration %i --job_id %i -d %s -c %i -o %s
            """ % (output_dir, smi, mol_id, i, j, db_fname, 16, output_fname)
            with open(pbs_name, "wt") as f:
                f.write(script)
            p = Popen(['qsub', pbs_name], stdout=PIPE, encoding='utf8')
            outs, errs = p.communicate(timeout=30)
            id = outs.strip()
            job_ids[id] = output_fname
            # print("job id: %s was submitted" % id)

        time.sleep(10)

        uniq_smi = dict()  # smi: (id/name, mw)

        while True:

            unfinished_jobs = get_jobs(job_ids, unfinished=True)
            finished_jobs = list(set(job_ids).difference(unfinished_jobs))
            # read files with enumerated mols
            for k in finished_jobs:
                with open(job_ids[k]) as f:
                    f.readline()  # header
                    for line in f:
                        items = line.split('\t')
                        mw = float(items[4])
                        if mw <= 500:
                            uniq_smi[items[0]] = (items[1], mw)
                del job_ids[k]

            if not unfinished_jobs:
                break

            time.sleep(min(i * 2, 60))

        if len(uniq_smi) > nmols:
            mw = sorted((v[1], (k, v[0])) for k, v in uniq_smi.items())  # v[0] - name, v[1] - mw
            input_mols = []
            for c in divide(nmols, mw):
                input_mols.append(random.choice(list(c))[1])
        else:
            input_mols = [(k, v[0]) for k, v in uniq_smi.items()]

        print("[%s] Iteration: %i, mols: %i" % (str(datetime.datetime.now()), i + 1, len(uniq_smi)))