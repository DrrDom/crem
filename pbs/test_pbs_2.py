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

ncpus = 30
nnodes = 1
njobs = 1

for i in range(njobs):
    job_name = "pyjob_%i" % i
    pbs_name = "job_%i.pbs" % i
    script = """
    #!/bin/bash
    #PBS -l select=1:ncpus=30:mem=1gb
    cd /home/pavlop/python/crem
    python3 test_job_pbs_2.py %i job_output_%i.txt
    """ % (ncpus, i)
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

