#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 20-08-2018
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2018
# license         : 
#==============================================================================

import argparse
import os
import re
from subprocess import Popen


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create fragment DB from sorted text files with context and fragments.')
    parser.add_argument('-d', '--dir', metavar='DIR', required=True,
                        help='dir where input files are stored.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "dir": output_dir = os.path.abspath(v)

    job_dir = os.path.join(output_dir, 'jobs')
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    pbs_name = os.path.join(job_dir, 'db_gen_%s.pbs' % (os.path.basename(output_dir)))
    script = """
    #!/usr/bin/env bash
    #PBS -l select=1:ncpus=32
    #PBS -k oe

    cd %s

    """ % output_dir

    for fname in sorted(os.listdir(output_dir)):
        if re.search('^r[1-5]_c\.txt$', fname):
            script += 'python3 ~/python/crem/import_env_to_db.py -i %s -o replacements.db -t radius%s -c\n' % (fname, fname[1])

    with open(pbs_name, "wt") as f:
        f.write(script)
    p = Popen(['qsub', pbs_name], encoding='utf8')


