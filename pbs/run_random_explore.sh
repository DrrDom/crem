#!/usr/bin/env bash

cd ~/imtm/crem/comb_smarts/

#screen -d -m -S crem_o_r1 bash -c 'python ~/python/crem/pbs/random_explore.py -o ~/imtm/crem/gen/orgelm_r1_f0 -d ~/imtm/crem/db/orgelm/replacements.db -r 1'
#screen -d -m -S crem_o_r2 bash -c 'python ~/python/crem/pbs/random_explore.py -o ~/imtm/crem/gen/orgelm_r2_f0 -d ~/imtm/crem/db/orgelm/replacements.db -r 2'
#screen -d -m -S crem_o_r3 bash -c 'python ~/python/crem/pbs/random_explore.py -o ~/imtm/crem/gen/orgelm_r3_f0 -d ~/imtm/crem/db/orgelm/replacements.db -r 3'
screen -d -m -S crem_pains_r1 bash -c 'python ~/python/crem/pbs/random_explore.py -o gen/pains_r1_f0 -d db/pains/replacements.db -r 1'
screen -d -m -S crem_pains_r2 bash -c 'python ~/python/crem/pbs/random_explore.py -o gen/pains_r2_f0 -d db/pains/replacements.db -r 2'
screen -d -m -S crem_pains_r3 bash -c 'python ~/python/crem/pbs/random_explore.py -o gen/pains_r3_f0 -d db/pains/replacements.db -r 3'
screen -d -m -S crem_pains_r4 bash -c 'python ~/python/crem/pbs/random_explore.py -o gen/pains_r4_f0 -d db/pains/replacements.db -r 4'
screen -d -m -S crem_pains_r5 bash -c 'python ~/python/crem/pbs/random_explore.py -o gen/pains_r5_f0 -d db/pains/replacements.db -r 5'
