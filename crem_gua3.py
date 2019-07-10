#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 26-06-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

# from __future__ import print_function

import argparse
import heapq
import json
import os
import random
from time import time
from typing import List, Optional

import joblib
import numpy as np

import sys
sys.path.append('/home/pavlop/python/guacamol')
# sys.path.append('/home/pavel/python/guacamol')
db_fname = '/home/pavlop/imtm/crem/update/db/orgelm/replacements.db'
# db_fname = '/home/pavel/QSAR/crem/update/db/lilly_pains_sc2/replacements.db'
smi_fname = '/home/pavlop/python/guacamol/guacamol/data/guacamol_v1_all.smi'
# smi_fname = '/home/pavel/python/guacamol/guacamol/data/guacamol_v1_all.smi'
total_ncpu = 31
# total_ncpu = 2


from guacamol.assess_goal_directed_generation import assess_goal_directed_generation
from guacamol.goal_directed_generator import GoalDirectedGenerator
from guacamol.scoring_function import ScoringFunction
from guacamol.utils.chemistry import canonicalize
from guacamol.utils.helpers import setup_default_logger
from joblib import delayed
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from crem import mutate_mol2, mutate_mol

from pprint import pprint
from multiprocessing import Pool, cpu_count

def make_mating_pool(population_mol: List[Mol], population_scores, offspring_size: int):
    """
    Given a population of RDKit Mol and their scores, sample a list of the same size
    with replacement using the population_scores as weights
    Args:
        population_mol: list of RDKit Mol
        population_scores: list of un-normalised scores given by ScoringFunction
        offspring_size: number of molecules to return
    Returns: a list of RDKit Mol (probably not unique)
    """
    # scores -> probs
    sum_scores = sum(population_scores)
    population_probs = [p / sum_scores for p in population_scores]
    mating_pool = np.random.choice(population_mol, p=population_probs, size=offspring_size, replace=True)
    return mating_pool


def select_new_pool(ids, scores, num):
    # sum_scores = sum(scores)
    # probs = [p / sum_scores for p in scores]
    probs = [i / sum(range(len(ids))) for i in reversed(range(len(ids)))]
    return np.random.choice(ids, p=probs, size=num, replace=False)


def score_mol(mol, score_fn):
    return score_fn(Chem.MolToSmiles(mol))


class CREM_Generator(GoalDirectedGenerator):

    def __init__(self, smi_file, selection_size, radius, replacements, 
                 max_size, min_size, max_inc, min_inc,
                 generations, ncpu, random_start, output_dir):
        self.pool = joblib.Parallel(n_jobs=ncpu)
        self.smiles = self.load_smiles_from_file(smi_file)
        self.selection_size = selection_size
        self.radius = radius
        self.min_size = min_size
        self.max_size = max_size
        self.min_inc = min_inc
        self.max_inc = max_inc
        self.replacements = replacements
        self.replacements_baseline = replacements
        self.generations = generations
        self.random_start = random_start
        self.patience = 2
        self.patience2 = 10
        self.patience_restart = 24
        self.task = 0
        self.output_dir = output_dir

    def load_smiles_from_file(self, smi_file):
        with open(smi_file) as f:
            return list(line.strip() for line in f)
            # return self.pool(delayed(canonicalize)(s.strip()) for s in f)

    def top_k(self, smiles, scoring_function, k):
        joblist = (delayed(scoring_function.score)(s) for s in smiles)
        scores = self.pool(joblist)
        scored_smiles = list(zip(scores, smiles))
        scored_smiles = sorted(scored_smiles, key=lambda x: x[0], reverse=True)
        return [smile for score, smile in scored_smiles][:k]

    def generate(self, mols):
        res = self.pool(delayed(mutate_mol2)(mol, db_name=db_fname,
                                             radius=self.radius, min_size=self.min_size,
                                             max_size=self.max_size,
                                             min_rel_size=0, max_rel_size=1,
                                             min_inc=self.min_inc, max_inc=self.max_inc,
                                             max_replacements=self.replacements,
                                             replace_cycles=False,
                                             protected_ids=None, min_freq=0,
                                             return_rxn=False, return_rxn_freq=False,
                                             ncores=1) for mol in mols)
        res = set(m[0] for sublist in res for m in sublist)
        return list(res)

    def get_params(self, score):
        # get min_inc, max_inc, max_replacements
        if score > 0.8:
            return -4, 4, self.replacements_baseline + 1000
        if score > 0.6:
            return -5, 5, self.replacements_baseline + 200
        if score > 0.4:
            return -6, 6, self.replacements_baseline + 100
        if score > 0.3:
            return -7, 7, self.replacements_baseline
        if score > 0.2:
            return -8, 8, self.replacements_baseline
        return -9, 9, self.replacements_baseline

    def generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int,
                                     starting_population: Optional[List[str]] = None) -> List[str]:

        self.task += 1
        restart = starting_population is None

        if number_molecules > self.selection_size:
            self.selection_size = number_molecules
            print(f'Benchmark requested more molecules than expected: new population is {number_molecules}')

        if starting_population is None:
            print('selecting initial population...')
            if self.random_start:
                starting_population = np.random.choice(self.smiles, self.selection_size)
            else:
                starting_population = self.top_k(self.smiles, scoring_function, self.selection_size)

        # select initial population
        population_smiles = heapq.nlargest(self.selection_size, starting_population, key=scoring_function.score)
        population_mols = [Chem.MolFromSmiles(s) for s in population_smiles]
        population_scores = self.pool(delayed(score_mol)(m, scoring_function.score) for m in population_mols)

        # evolution: go go go!!
        t0 = time()

        patience = 0
        patience2 = 0
        patience2_score = 0
        patience_restart = 0

        used_smiles = set(population_smiles)   # smiles already used for mutation
        best_smiles = set()
        ref_score = max(population_scores)

        self.min_inc, self.max_inc, self.replacements = self.get_params(max(population_scores))

        for generation in range(self.generations):

            new_smiles = set(self.generate(population_mols))
            new_mols = [Chem.MolFromSmiles(s) for s in new_smiles]
            new_scores = self.pool(delayed(score_mol)(m, scoring_function.score) for m in new_mols)

            best_smiles = list(best_smiles) + list(zip(new_scores, new_smiles))
            best_smiles = set(best_smiles)
            best_smiles = sorted(best_smiles, key=lambda x: x[0], reverse=True)[:number_molecules]

            # update parameters
            if max(new_scores) <= ref_score and max(new_scores) < 1:
                if patience >= self.patience:
                    self.min_inc -= 1
                    self.max_inc += 1
                    patience = 0
                    used_smiles = set()
                    ref_score = max(new_scores)
                else:
                    patience += 1
            else:
                patience = 0
                ref_score = max(new_scores)
                self.min_inc, self.max_inc, self.replacements = self.get_params(best_smiles[0][0])

            if sum(score for score, smi in best_smiles) <= patience2_score:
                patience_restart += 1
                if patience2 >= self.patience2:
                    self.min_inc -= 15
                    self.max_inc += 15
                    self.replacements = 10000
                    patience2 = 0
                    patience = -5
                else:
                    patience2 += 1
            else:
                patience2_score = sum(score for score, smi in best_smiles)
                patience2 = 0
                patience_restart = 0

            # select next population
            new_tuples = []
            for score, mol, smi in zip(new_scores, new_mols, new_smiles):
                if smi not in population_smiles:
                    new_tuples.append((score, mol, smi))
            new_tuples += list(zip(population_scores, population_mols, population_smiles))
            new_tuples = sorted(new_tuples, key=lambda x: x[0], reverse=True)

            ids = []
            for i, (score, mol, smi) in enumerate(new_tuples):
                if smi not in used_smiles:
                    ids.append(i)
                    if len(ids) == self.selection_size:
                        break
                    
            #ids = select_new_pool(ids, [new_tuples[i][0] for i in ids], self.selection_size)

            population_scores = [new_tuples[i][0] for i in ids]
            population_mols = [new_tuples[i][1] for i in ids]
            population_smiles = [new_tuples[i][2] for i in ids]

            used_smiles.update(population_smiles)

            # stats
            gen_time = time() - t0
            t0 = time()

            print(f'{generation: >5} | '
                  f'best avg: {round(sum(score for score, smi in best_smiles) / len(best_smiles), 3)} | '
                  f'max: {np.max(population_scores):.3f} | '
                  f'avg: {np.mean(population_scores):.3f} | '
                  f'min: {np.min(population_scores):.3f} | '
                  f'std: {np.std(population_scores):.3f} | '
                  f'sum: {np.sum(population_scores):.3f} | '
                  f'min_inc: {self.min_inc} | '
                  f'max_inc: {self.max_inc} | '
                  f'replacements: {self.replacements} | '
                  f'patience: {patience} | '
                  f'patience2: {patience2} | '
                  f'{gen_time:.2f} sec')

            sys.stdout.flush()

            # everything is perfect
            if all(score == 1 for score, smi in best_smiles):
                break

            # optimization get stuck - no improvement for a long time
            if patience_restart >= self.patience_restart:
                if restart:
                    # restart optimization from another random seeds
                    starting_population = np.random.choice(self.smiles, self.selection_size)
                    population_smiles = heapq.nlargest(self.selection_size, starting_population, key=scoring_function.score)
                    population_mols = [Chem.MolFromSmiles(s) for s in population_smiles]
                    population_scores = self.pool(delayed(score_mol)(m, scoring_function.score) for m in population_mols)
                    ref_score = max(population_scores)
                    patience_restart = 0
                else:
                    # give up
                    break

        # finally
        with open(os.path.join(self.output_dir, f'{self.task}.smi'), 'wt') as f:
            for score, smi in best_smiles:
                f.write(f'{smi}\t{round(score, 3)}\n')
        return [smi for score, smi in best_smiles]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles_file', default=smi_fname)
    parser.add_argument('--selection_size', type=int, default=5)
    parser.add_argument('--radius', type=int, default=3)
    parser.add_argument('--replacements', type=int, default=100)
    parser.add_argument('--min_size', type=int, default=0)
    parser.add_argument('--max_size', type=int, default=10)
    parser.add_argument('--min_inc', type=int, default=-7)
    parser.add_argument('--max_inc', type=int, default=7)
    parser.add_argument('--generations', type=int, default=1000)
    parser.add_argument('--ncpu', type=int, default=total_ncpu)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output_dir', type=str, default=None)
    parser.add_argument('--suite', default='v2')

    args = parser.parse_args()

    np.random.seed(args.seed)

    setup_default_logger()

    if args.output_dir is None:
        args.output_dir = os.path.dirname(os.path.realpath(__file__))

    # save command line args
    with open(os.path.join(args.output_dir, 'goal_directed_params.json'), 'w') as jf:
        json.dump(vars(args), jf, sort_keys=True, indent=4)

    optimiser = CREM_Generator(smi_file=args.smiles_file,
                               selection_size=args.selection_size,
                               radius=args.radius,
                               min_size=args.min_size,
                               max_size=args.max_size,
                               min_inc=args.min_inc,
                               max_inc=args.max_inc,
                               replacements=args.replacements,
                               generations=args.generations,
                               ncpu=args.ncpu,
                               random_start=True,
                               output_dir=args.output_dir)

    json_file_path = os.path.join(args.output_dir, 'goal_directed_results.json')
    assess_goal_directed_generation(optimiser, json_output_file=json_file_path, benchmark_version=args.suite)


if __name__ == "__main__":
    main()
