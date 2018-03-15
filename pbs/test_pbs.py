from subprocess import Popen, PIPE, TimeoutExpired, call
import time


def not_finished(job_ids):
    p = Popen(['qstat'], stdout=PIPE, encoding='utf8')
    outs, errs = p.communicate(timeout=30)
    lines = outs.split('\n')
    if len(lines) > 2:
        for line in lines[2:]:
            job = [s for s in line.strip().split(' ') if s]  # Job_id Name User Time_Use S Queue
            if job[0] in job_ids:
                return True
    return False


job_ids = []

smiles = ['C1=CC=CC=C1', 'C1COCC1']

for i, smi in enumerate(smiles):
    job_name = "pyjob_%i" % i
    pbs_name = "job_%i.pbs" % i
    script = """
    #!/bin/bash
    #PBS -l select=1:ncpus=3:mem=1gb,walltime=1:00:00
    #PBS -o=/home/pavlop/python/crem/%s.o
    #PBS -e=/home/pavlop/python/crem/%s.e
    #PBS -N %s
    cd /home/pavlop/python/crem
    source activate rdkit-1709
    python3 test_job_pbs.py job_output_%i.txt %s
    ### python3 -c "f = open('job_output_%i.txt', 'wt'); f.write('%s'); f.close()"
    """ % (job_name, job_name, job_name, i, smi, i, smi)
    with open(pbs_name, "wt") as f:
        f.write(script)
    p = Popen(['qsub', pbs_name], stdout=PIPE, encoding='utf8')
    outs, errs = p.communicate(timeout=30)
    job_ids.append(outs.strip())
    print("job id: %s was submitted" % job_ids[-1])

time.sleep(10)
while not_finished(job_ids):
    time.sleep(30)

print("All done")




# p = Popen(['qstat'], stdout=PIPE, encoding='utf8')
#
# try:
#     outs, errs = p.communicate(timeout=30)
#
#     lines = outs.split('\n')
#
#     if len(lines) > 2:
#         for line in lines[2:]:
#             job = [s for s in line.strip().split(' ') if s]  # Job_id Name User Time_Use S Queue
#             if job:
#                 print('job name: %s, user: %s, status: %s' % (job[1], job[2], job[4]))
#
# except TimeoutExpired:
#
#     p.kill()
#     outs, errs = p.communicate()
#
#
