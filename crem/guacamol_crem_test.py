#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 26-06-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

import argparse
import json
import os
import sys
from time import time
from typing import List, Optional

import joblib
import numpy as np
import pandas as pd

from guacamol.assess_goal_directed_generation import assess_goal_directed_generation
from guacamol.goal_directed_generator import GoalDirectedGenerator
from guacamol.scoring_function import ScoringFunction
from guacamol.utils.helpers import setup_default_logger
from joblib import delayed
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from .crem import mutate_mol2


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


def score_mol(mol, score_fn):
    return score_fn(Chem.MolToSmiles(mol))


class CREM_Generator(GoalDirectedGenerator):

    def __init__(self, smi_file, selection_size, db_fname, radius,
                 replacements, max_size, min_size, max_inc, min_inc,
                 generations, ncpu, random_start, output_dir):
        self.pool = joblib.Parallel(n_jobs=ncpu)
        self.smiles = self.load_smiles_from_file(smi_file)
        self.N = selection_size
        self.db_fname = db_fname
        self.radius = radius
        self.min_size = min_size
        self.max_size = max_size
        self.min_inc = min_inc
        self.max_inc = max_inc
        self.replacements = replacements
        self.replacements_baseline = replacements
        self.generations = generations
        self.random_start = random_start
        self.patience1 = 3
        self.patience2 = 10
        self.patience3 = 33
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

    def generate(self, smiles):
        mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles]
        res = self.pool(delayed(mutate_mol2)(mol, db_name=self.db_fname,
                                             radius=self.radius, min_size=self.min_size,
                                             max_size=self.max_size,
                                             min_rel_size=0, max_rel_size=1,
                                             min_inc=self.min_inc, max_inc=self.max_inc,
                                             max_replacements=self.replacements,
                                             replace_cycles=False,
                                             protected_ids=None, min_freq=0,
                                             return_rxn=False, return_rxn_freq=False,
                                             ncores=1) for mol in mols)
        res = set(m for sublist in res for m in sublist)
        return list(res)

    def set_params(self, score):
        # get min_inc, max_inc, max_replacements
        self.replacements = self.replacements_baseline
        if score > 0.8:
            self.min_inc = -4
            self.max_inc = 4
        elif score > 0.7:
            self.min_inc = -5
            self.max_inc = 5
        elif score > 0.6:
            self.min_inc = -6
            self.max_inc = 6
        elif score > 0.5:
            self.min_inc = -7
            self.max_inc = 7
        elif score > 0.4:
            self.min_inc = -8
            self.max_inc = 8
        elif score > 0.3:
            self.min_inc = -9
            self.max_inc = 9
        else:
            self.min_inc = -10
            self.max_inc = 10

    def get_scores(self, scoring_function, smiles):
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        return self.pool(delayed(score_mol)(m, scoring_function.score) for m in mols)

    def generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int,
                                     starting_population: Optional[List[str]] = None) -> List[str]:

        self.task += 1

        if number_molecules > self.N:
            self.N = number_molecules
            print(f'Benchmark requested more molecules than expected: new population is {number_molecules}')

        # select initial population
        if starting_population is None:
            print('selecting initial population...')
            if self.random_start:
                population = pd.DataFrame(np.random.choice(self.smiles, self.N), columns=['smi'])
            else:
                population = pd.DataFrame(self.top_k(self.smiles, scoring_function, self.N), columns=['smi'])
        else:
            population = pd.DataFrame(starting_population, columns=['smi'])
        population['score'] = self.get_scores(scoring_function, population['smi'])

        # evolution: go go go!!
        t0 = time()
        time_start = t0

        patience1 = 0
        patience2 = 0
        patience3 = 0

        best = population.copy().drop_duplicates(subset='smi')
        ref_score = np.mean(best['score'].iloc[:number_molecules])
        self.set_params(max(best['score']))
        used_smiles = set(population['smi'])   # smiles already used for mutation

        for generation in range(self.generations):

            if ref_score == 1:
                break

            population = pd.DataFrame(list(set(self.generate(population['smi']))), columns=['smi'])
            population['score'] = self.get_scores(scoring_function, population['smi'])
            population.sort_values(by='score', ascending=False, inplace=True)
            population.drop_duplicates(subset='smi', inplace=True)

            best = best.append(population).\
                drop_duplicates(subset='smi').\
                sort_values(by='score', ascending=False).\
                head(self.N)
            cur_score = np.mean(best['score'].iloc[:number_molecules])

            if cur_score > ref_score:
                ref_score = cur_score
                population = population.head(self.N)
                self.set_params(max(population['score']))
                used_smiles.update(population['smi'])
                patience1 = 0
                patience2 = 0
                patience3 = 0
            else:
                patience1 += 1
                patience2 += 1
                patience3 += 1
                if patience3 >= self.patience3:
                    if starting_population is None and self.random_start:
                        patience1 = 0
                        patience2 = 0
                        patience3 = 0
                        population = pd.DataFrame(np.random.choice(self.smiles, self.N), columns=['smi']).drop_duplicates(subset='smi')
                        population['score'] = self.get_scores(scoring_function, population['smi'])
                        population.sort_values(by='score', ascending=False, inplace=True)
                        self.set_params(max(population['score']))
                        used_smiles = set(population['smi'])
                    else:
                        break
                else:
                    population = population.head(self.N)
                    used_smiles.update(population['smi'])
                    if patience2 >= self.patience2:
                        patience1 = 0
                        patience2 = 0
                        self.min_inc -= 10
                        self.max_inc += 10
                        self.replacements += 500
                        used_smiles = set(population['smi'])
                    elif patience1 >= self.patience1:
                        patience1 = 0
                        self.min_inc -= 1
                        self.max_inc += 1
                        self.replacements += 100
                        used_smiles = set(population['smi'])

            # stats
            gen_time = time() - t0
            t0 = time()
            print(f'{generation: >5} | '
                  f'best avg: {np.round(np.mean(best["score"].iloc[:number_molecules]), 3)} | '
                  f'max: {np.max(population["score"]):.3f} | '
                  f'avg: {np.mean(population["score"]):.3f} | '
                  f'min: {np.min(population["score"]):.3f} | '
                  f'std: {np.std(population["score"]):.3f} | '
                  f'sum: {np.sum(population["score"]):.3f} | '
                  f'min_inc: {self.min_inc} | '
                  f'max_inc: {self.max_inc} | '
                  f'repl: {self.replacements} | '
                  f'p1: {patience1} | '
                  f'p2: {patience2} | '
                  f'p3: {patience3} | '
                  f'{gen_time:.2f} sec')
            sys.stdout.flush()

            if t0 - time_start > 18000:   # execution time > 5hr
                break

        # finally
        best.round({'score': 3}).to_csv(os.path.join(self.output_dir, f'{self.task}.smi'), sep="\t", header=False, index=False)
        return best['smi'][:number_molecules]


def entry_point():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles_file', type=str)
    parser.add_argument('--db_fname', type=str)
    parser.add_argument('--selection_size', type=int, default=10)
    parser.add_argument('--radius', type=int, default=3)
    parser.add_argument('--replacements', type=int, default=1000)
    parser.add_argument('--min_size', type=int, default=0)
    parser.add_argument('--max_size', type=int, default=10)
    parser.add_argument('--min_inc', type=int, default=-7)
    parser.add_argument('--max_inc', type=int, default=7)
    parser.add_argument('--generations', type=int, default=1000)
    parser.add_argument('--ncpu', type=int, default=1)
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
                               db_fname=args.db_fname,
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
    entry_point()
