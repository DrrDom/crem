# Pavel Polishchuk, 2017

import sys
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMMPA
from mol_context import get_canon_context_core


cycle_pattern = re.compile("(?<!:)[1-9]+")


def smiles_to_smarts(smi):

    mol = Chem.MolFromSmiles(smi)

    if mol is None:
        sys.stderr.write("Can't generate mol for: %s\n" % smi)
        return None

    # change the isotope to 42
    for atom in mol.GetAtoms():
        atom.SetIsotope(42)

    # print out the smiles - all the atom attributes will be fully specified
    smarts = Chem.MolToSmiles(mol, isomericSmiles=True)
    # remove the 42 isotope labels
    smarts = re.sub(r'\[42', "[", smarts)
    # now have a fully specified SMARTS - simples!

    return smarts


def __fragment_mol(mol, radius=3, return_ids=True, keep_stereo=False, protected_ids=None):
    """
    INPUT:
        mol - Mol
        radius - integer, number of bonds to cut context

    OUTPUT:
        list of tuples (env_smi, core_smi, tuple of core atom ids)
        ('C[*:1].C[*:2]', 'CC(C(=O)O)c1ccc(CC([*:1])[*:2])c(Br)c1', (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
        ('Cc(c)c(cc)[*:1]', '[H][*:1]', (25,))

    If input mol has explicit hydrogens the output will contain also fragments where core = [H][*:1].
    Smiles of fragments with heavy atoms will contain only heavy atoms
    """

    def get_atom_prop(molecule, prop="Index"):
        res = []
        for a in molecule.GetAtoms():
            try:
                res.append(a.GetIntProp(prop))
            except KeyError:
                continue
        return tuple(sorted(res))

    # if chosen
    if protected_ids and not return_ids:
        return_ids = True

    output = []

    # set original atom idx to keep them in fragmented mol
    if return_ids:
        for atom in mol.GetAtoms():
            atom.SetIntProp("Index", atom.GetIdx())

    # heavy atoms
    frags = rdMMPA.FragmentMol(mol, pattern="[!#1]!@!=!#[!#1]", maxCuts=4, resultsAsMols=True, maxCutBonds=30)
    # hydrogen atoms
    frags += rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)

    for i, (core, chains) in enumerate(frags):
        if core is None:  # single cut
            components = list(Chem.GetMolFrags(chains, asMols=True))
            ids_0 = get_atom_prop(components[0]) if return_ids else tuple()
            ids_1 = get_atom_prop(components[1]) if return_ids else tuple()
            if Chem.MolToSmiles(components[0]) != '[H][*:1]':  # context cannot be H
                env, frag = get_canon_context_core(components[0], components[1], radius, keep_stereo)
                output.append((env, frag, ids_1))
            if Chem.MolToSmiles(components[1]) != '[H][*:1]':  # context cannot be H
                env, frag = get_canon_context_core(components[1], components[0], radius, keep_stereo)
                output.append((env, frag, ids_0))
        else:   # multiple cuts
            # there are no checks for H needed because H can be present only in single cuts
            env, frag = get_canon_context_core(chains, core, radius, keep_stereo)
            output.append((env, frag, get_atom_prop(core)))

    return output  # list of tuples (env smiles, core smiles, list of atom ids)


def __frag_replace(mol, frag_sma, replace_sma, frag_ids=None):
    """
    INPUT
        mol:         mol,
        frag_sma:    SMARTS of a fragment,
        replace_sma: SMARTS of a replacement (from DB),
        frag_ids:    atom ids of a fragment if you need to make exact replacement, if None all possible matches will be
                     replaced
    OUTPUT
        list of mols with replaced fragment.
        Each output mol has a new field named 'transformation' with information about reaction SMARTS applied
    """

    frag_sma = frag_sma.replace('*', '!#1')    # to avoid map H in mol with explicit H (lead to wrong replacement)
    rxn_sma = "%s>>%s" % (frag_sma, replace_sma)
    rxn = AllChem.ReactionFromSmarts(rxn_sma)

    if frag_ids:
        ids = set(frag_ids)
        # extend atom ids on neighbour atoms
        for i in frag_ids:
            a = mol.GetAtomWithIdx(i)
            ids.update(na.GetIdx() for na in a.GetNeighbors())
        # protect untouched atoms
        for a in mol.GetAtoms():
            if a.GetAtomicNum() > 1 and a.GetIdx() not in ids:
                a.SetProp('_protected', '1')
            else:
                a.ClearProp('_protected')

    ps = rxn.RunReactants([mol])

    products = set()
    for y in ps:
        for p in y:
            e = Chem.SanitizeMol(p, catchErrors=True)
            if e:
                sys.stderr.write("Molecule %s caused sanitization error %i" % (Chem.MolToSmiles(p, isomericSmiles=True), e))
                sys.stderr.flush()
            else:
                smi = Chem.MolToSmiles(Chem.RemoveHs(p), isomericSmiles=True)
                if smi not in products:
                    products.add(smi)
                    yield smi, rxn_sma


def __get_replacements(db_cur, env, min_atoms, max_atoms, radius):
    if radius == 3:
        db_cur.execute("""SELECT core_smi, core_sma
                          FROM radius3
                          WHERE env IN (SELECT env FROM radius3 WHERE env = ?)
                                AND 
                                freq >= 10 AND
                                core_num_atoms BETWEEN ? AND ?""", (env, min_atoms, max_atoms))
    elif radius == 2:
        db_cur.execute("""SELECT core_smi, core_sma
                          FROM radius2
                          WHERE env IN (SELECT env FROM radius2 WHERE env = ?)
                                AND 
                                freq >= 10 AND
                                core_num_atoms BETWEEN ? AND ?""", (env, min_atoms, max_atoms))
    return db_cur.fetchall()


def mutate_mol(mol, db_cur, radius=3, min_size=1, max_size=10, min_rel_size=0, max_rel_size=1, min_inc=-2, max_inc=2,
               replace_cycles=False, protected_ids=None):
    """
    INPUT
        mol:      mol
        db_cur:   cursor of SQLite3 DB with fragment replacements. DB should contain tables named "radius2", "radius3".
        radius:   integer. How far molecular context should be considered while looking for interchangeable fragments.
        min_size, max_size: integer. Size of a fragment to replace.
        min_rel_size, max_rel_size: float. Relative size of a fragment to the whole mol
                                    (in terms of a number of heavy atoms)
        min_inc, max_inc:   integer. Relative minimum and maximum size of new fragments which will replace
                            the existed one. -2 and 2 mean that the existed fragment with N atoms will be replaced
                            with fragments from a DB having from N-2 to N+2 atoms.
        replace_cycles:     looking for replacement of a fragment containing cycles irrespectively of the fragment size

    OUTPUT
        list of unique mols

    Supply mol with explicit Hs if H replacement is desired
    """

    def get_core_sma(db_cur, env, core_smi, radius):
        if radius == 3:
            db_cur.execute("SELECT core_sma FROM radius3 WHERE env = ? AND core_smi = ?", (env, core_smi))
        elif radius == 2:
            db_cur.execute("SELECT core_sma FROM radius2 WHERE env = ? AND core_smi = ?", (env, core_smi))
        res = db_cur.fetchone()
        return res[0] if res else None

    products = set()

    f = __fragment_mol(mol, radius, protected_ids=protected_ids)

    mol_hac = mol.GetNumHeavyAtoms()

    for env, core, ids in f:

        num_heavy_atoms = Chem.MolFromSmiles(core).GetNumHeavyAtoms()
        hac_ratio = num_heavy_atoms / mol_hac

        if (min_size <= num_heavy_atoms <= max_size and min_rel_size <= hac_ratio <= max_rel_size) \
                or (replace_cycles and cycle_pattern.search(core)):

            # frag_sma = get_core_sma(db_cur, env, core, radius)
            frag_sma = smiles_to_smarts(core)

            # if frag_sma:

            min_atoms = num_heavy_atoms + min_inc
            max_atoms = num_heavy_atoms + max_inc

            rep = __get_replacements(db_cur, env, min_atoms, max_atoms, radius)
            for core_smi, core_sma in rep:
                if core_smi != core:
                    for smi, rxn in __frag_replace(mol, frag_sma, core_sma, ids):
                        # smi = Chem.MolToSmiles(new_mol, isomericSmiles=True)
                        if smi not in products:
                            products.add(smi)
                            yield smi, rxn
                            # products[smi] = new_mol

    # return list(products.values())


from pprint import pprint
# mm = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1Br)C(C)C(=O)O')
# mm = Chem.AddHs(mm)
# f = __fragment_mol(mm, 3)
# pprint(f)


# frag_sma = '[H][*:1]'
# frag_ids = (32, )
# replace_sma = '[CH3][*:1]'
#
# r = __frag_replace(mm, frag_sma, replace_sma, frag_ids)
#
#
# for item in r:
#     print(item[0])

#
# print(Chem.MolToSmiles(mm, isomericSmiles=True))
# for a in mm.GetAtoms():
#     print(a.GetIdx(), a.GetSymbol())
#
# for item in r:
#     print(Chem.MolToSmiles(item, isomericSmiles=True))
#     for a in item.GetAtoms():
#         print(a.GetIdx(), a.GetSymbol())
#     item = Chem.RemoveHs(item)
#     print(Chem.MolToSmiles(item, isomericSmiles=True))
#     for a in item.GetAtoms():
#         print(a.GetIdx(), a.GetSymbol())


# mm = Chem.MolFromSmiles('CC1=CC=C(CC(=O)NC2CCCC2)C=C1Br')
# mm = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1Br)C(C(O)=O)c1ccc[nH]1')
# # mm = Chem.MolFromSmiles('c1ccccc1c1ccccc1')
# print(Chem.MolToSmiles(mm))
# pprint([(a.GetSymbol(), a.GetNumExplicitHs(), a.GetNumImplicitHs()) for a in mm.GetAtoms()])
# mm = Chem.AddHs(mm)
# print(Chem.MolToSmiles(mm))
# pprint([(a.GetSymbol(), a.GetNumExplicitHs(), a.GetNumImplicitHs()) for a in mm.GetAtoms()])
# mm = Chem.RemoveHs(mm)
# print(Chem.MolToSmiles(mm))
# pprint([(a.GetSymbol(), a.GetNumExplicitHs(), a.GetNumImplicitHs()) for a in mm.GetAtoms()])

# f = __fragment_mol(mm, 2)
#
# pprint(f)

# for a in mm.GetAtoms():
#     print(a.GetIdx(), a.GetSymbol(), a.GetIsAromatic())

# r = __frag_replace(mm, '[CH2]1[CH2][CH2][CH]([*:1])[CH2]1', '[CH3][C]1([CH3])[S][CH]2[CH]([*:1])[C](=[O])[N]2[CH]1[C]([NH2])=[O]', (9,10,11,12,13))
# print(r)
# for item in r:
#     print(Chem.MolToSmiles(item))
#     print(item.GetProp('transformation'))
