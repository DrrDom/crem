#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 19-08-2018
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2018
# license         : 
#==============================================================================

import argparse
import os
from subprocess import Popen, PIPE
import time
import re


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

    parser = argparse.ArgumentParser(description='Generate files for fragment database.')
    parser.add_argument('-o', '--out', metavar='OUTPUT_DIR', required=True,
                        help='output dir to store files.')
    parser.add_argument('-k', '--keep_ids', required=True,
                        help='path to the file with mol ids to keep at SB generation.')
    parser.add_argument('-r', '--radius', required=False, nargs='*', default=[1, 2, 3],
                        help='list of radius values. Default: [1,2,3].')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "out": output_dir = os.path.abspath(v)
        if o == "keep_ids": fname = os.path.abspath(v)
        if o == "radius": radius = list(map(int, v))

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    job_dir = os.path.join(output_dir, 'jobs')
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    job_ids = []

    for r in radius:
        pbs_name = os.path.join(job_dir, '%s_r%i.pbs' % (os.path.basename(output_dir), r))
        script = """
        #!/usr/bin/env bash
        #PBS -l select=1:ncpus=32
        #PBS -k oe
        
        RADIUS=%i
        
        cd %s
        source ~/anaconda3/bin/activate rdkit-1709
        
        python3 ~/python/crem/frag_to_env_mp.py -i ~/imtm/crem/frags.txt -o r${RADIUS}.txt -k %s -r ${RADIUS} -c 32 -v
        sort r${RADIUS}.txt | uniq -c > r${RADIUS}_c.txt
        """ % (r, output_dir, fname)
        with open(pbs_name, "wt") as f:
            f.write(script)
        p = Popen(['qsub', pbs_name], stdout=PIPE, encoding='utf8')
        outs, errs = p.communicate(timeout=30)
        id = outs.strip()
        job_ids.append(id)
        print(id)

    time.sleep(60)

    while get_jobs(job_ids, unfinished=True):
        time.sleep(60)

    pbs_name = os.path.join(job_dir, 'db_gen_%s.pbs' % (os.path.basename(output_dir)))
    script = """
    #!/usr/bin/env bash
    #PBS -l select=1:ncpus=32
    #PBS -k oe

    cd %s
    source ~/anaconda3/bin/activate rdkit-1709

    """ % output_dir

    for fname in sorted(os.listdir(output_dir)):
        if re.search('^r[1-5]_c\.txt$', fname):
            script += 'python3 ~/python/crem/import_env_to_db.py -i %s -o replacements.db -r %s -c -n 32\n' % (fname, fname[1])
    with open(pbs_name, "wt") as f:
        f.write(script)
    p = Popen(['qsub', pbs_name], stdout=PIPE, encoding='utf8')
    outs, errs = p.communicate(timeout=30)
    id = outs.strip()
    print(id)
