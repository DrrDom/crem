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
from subprocess import Popen


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate files for fragment database.')
    parser.add_argument('-o', '--out', metavar='OUTPUT_DIR', required=True,
                        help='output dir to store files.')
    parser.add_argument('-k', '--keep_ids', required=True,
                        help='path to the file with mol ids to keep at SB generation.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "out": output_dir = os.path.abspath(v)
        if o == "keep_ids": fname = os.path.abspath(v)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    job_dir = os.path.join(output_dir, 'jobs')
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    for r in [1, 2, 3]:
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
        p = Popen(['qsub', pbs_name], encoding='utf8')

