# CReM - chemically reasonable mutations

**CReM** is an open-source Python framework to generate chemical structures using a fragment-based approach.

The main idea behind is similar to matched molecular pairs considering context that fragments in the identical context are interchangeable. Therefore, one can create a database of interchangeable fragments and use it for generation of chemically valid structures.

**Features:**  
1) Generation of a custom fragment database  
2) Three modes of structure generation: MUTATE, GROW, LINK  
3) Context radius to consider for replacement  
4) Fragment size to replace and the size of a replacing fragment  
5) Protection of atoms from modification (e.g. scaffold protection)  
6) Replacements with fragments occurred in a fragment database with certain minimal frequency  
7) Make randomly chosen replacements up to the specified number  

**Limitations and known issues**
1) New ring systems cannot be constructed from fragments, thus representativeness of ring systems in generated structures depends on a used fragment database. We are working on that issue.
2) Very large molecules will not be processed by CReM. If a molecule has more than 30 non-ring single bonds it will not be MUTATED. If a molecule has more than 100 hydrogen atoms it will not be processed by GROW and LINK.
3) Canonicalisation of contexts depends on RDKit SMILES representation. Thus, changing in RDKit SMILES representation may affect fragment databases and make impossible to use a database prepared with previous RDKit version from code running under later RDKit versions.  

## Documentation

https://crem.readthedocs.io/en/latest/

## Web app

To play with a tool online.  
https://crem.imtm.cz/

## Installation

Several command line utilities will be installed to create fragment databases and `crem` module will become available in Python imports to generate structures.

From pypi package
```text
pip install crem
```

Manually from repository
```text
git clone https://github.com/DrrDom/crem
cd crem
python3 setup.py sdist bdist_wheel
pip install dist/crem-0.1-py3-none-any.whl
```

Uninstall
```text
pip uninstall crem
```

## Dependencies

`crem` requires `rdkit>=2017.09`. To run the guacamol test `guacamol` should be installed.

## Generation of a fragment database

This step is required if you want to generate a custom fragment database. You can download precompiled databases obtained by fragmentation of the whole ChEMBL by links provided on this page - http://www.qsar4u.com/pages/crem.php.  

For convenience there is the bash script crem_create_frag_db.sh which includes all steps below. It takes three positional arguments: input file with SMILES, output directory where intermediate files and a final database will be stored and number of CPUs to use (this is optional, default value is 1).
```text
crem_create_frag_db.sh input.smi fragdb_dir 32
```
 
Fragmentation of input structures:
```text
fragmentation -i input.smi -o frags.txt -c 32 -v
```

Convert fragments to standardized representation of a core and a context of a given radius:
```text
frag_to_env -i frags.txt -o r3.txt -r 3 -c 32 -v
```

Remove duplicated lines in the output file and count frequency of occurrence of fragemnt-context pairs. These (`sort` and `uniq`) are `bash` utilities but since Win10 is Linux-friendly that should not be a big issue for Win users to execute them
```text
sort r3.txt | uniq -c > r3_c.txt
```

Create DB and import the file to a database table
```text
env_to_db -i r3_c.txt -o fragments.db -r 3 -c -v
```

Last three steps should be executed for each radius. All tables can be stored in the same database.

## Structure generation

Import necessary functions from the main module
```python
from crem.crem import mutate_mol, grow_mol, link_mols
from rdkit import Chem
```

Create a molecute and **mutate** it. Only one heavy atom will be substituted. Default radius is 3.
```python
m = Chem.MolFromSmiles('c1cc(OC)ccc1C')  # methoxytoluene
mols = list(mutate_mol(m, db_name='replacements.db', max_size=1))
```
output example
```text
['CCc1ccc(C)cc1',
 'CC#Cc1ccc(C)cc1',
 'C=C(C)c1ccc(C)cc1',
 'CCCc1ccc(C)cc1',
 'CC=Cc1ccc(C)cc1',
 'CCCCc1ccc(C)cc1',
 'CCCOc1ccc(C)cc1',
 'CNCCc1ccc(C)cc1',
 'COCCc1ccc(C)cc1',
 ...
 'Cc1ccc(C(C)(C)C)cc1']
```


Add hydrogens to the molecule to **mutate hydrogens** as well
```python
mols = list(mutate_mol(Chem.AddHs(m), db_name='replacements.db', max_size=1))
```
output
```text
['CCc1ccc(C)cc1',
 'CC#Cc1ccc(C)cc1',
 'C=C(C)c1ccc(C)cc1',
 'CCCc1ccc(C)cc1',
 'Cc1ccc(C(C)C)cc1',
 'CC=Cc1ccc(C)cc1',
 ...
 'COc1ccc(C)cc1C',
 'C=Cc1cc(C)ccc1OC',
 'COc1ccc(C)cc1Cl',
 'COc1ccc(C)cc1CCl']
```

**Grow** molecule. Only hydrogens will be replaced. Hydrogens should not be added explicitly.
```python
mols = list(grow_mol(m, db_name='replacements.db'))
```
output
```text
['COc1ccc(C)c(Br)c1',
 'COc1ccc(C)c(C)c1',
 'COc1ccc(C)c(Cl)c1',
 'COc1ccc(C)c(OC)c1',
 'COc1ccc(C)c(N)c1',
 ...
 'COc1ccc(CCN)cc1']
```

Create the second molecule and **link** it to toluene
```python
m2 = Chem.MolFromSmiles('NCC(=O)O')  # glycine
mols = list(link_mols(m, m2, db_name='replacements.db'))
```
output
```text
['Cc1ccc(OCC(=O)NCC(=O)O)cc1',
 'Cc1ccc(OCCOC(=O)CN)cc1',
 'COc1ccc(CC(=N)NCC(=O)O)cc1',
 'COc1ccc(CC(=O)NCC(=O)O)cc1',
 'COc1ccc(CC(=S)NCC(=O)O)cc1',
 'COc1ccc(CCOC(=O)CN)cc1']
```

You can vary the size of a linker and specify the distance between two attachment points in a linking fragment. There are many other arguments available in these functions, look at their **docstrings** for details.

##### Additional filters to control fragments chosen for replacing

An example of a filtering function which will keep only fragments containing a specific atom to be chosen for replacing.

```python
from collections import defaultdict
from functools import partial
from rdkit import Chem

def filter_function(row_ids, cur, radius, atom_number):

    """
    The first three arguments should be always the same as shown in the example. These parameters will be passed to a function from a main function, e.g. from mutate_mol. All other arguments are user-defined. The function should return the list of row ids of fragments which will be used for replacing. 

    :param row_id: a list of row ids from CReM database of those fragments which satisfy other selection criteria
    :param cur: cursor of CReM database
    :param radius: radius of a context 
    :param atom_number: an atomic number, fragments with this number will be discarded
    :return list of remaining row ids
    """

    # this part may be kept intact, it collects from DB SMILES of fragments with given row ids
    # since fragments may occur multiple times (due to different contexts) the results are collected in a dict
    if not row_ids:
        return []
    batch_size = 32000  # SQLite has a limit on a number of passed values to a query
    row_ids = list(row_ids)
    smis = defaultdict(list)  # {smi_1: [rowid_1, rowid_5, ...], ...}
    for start in range(0, len(row_ids), batch_size):
        batch = row_ids[start:start + batch_size]
        sql = f"SELECT rowid, core_smi FROM radius{radius} WHERE rowid IN ({','.join('?' * len(batch))})"
        for i, smi in cur.execute(sql, batch).fetchall():
            smis[smi].append(i)

    output_row_ids = []
    for smi, ids in smis.items():
        for a in Chem.MolFromSmiles(smi).GetAtoms():
            if a.GetAtomicNum() == atom_number:
                output_row_ids.extend(ids)
    return output_row_ids

# only F-containing fragments will be chosen for replacing
mol = Chem.MolFromSmiles('c1ccccc1C')
mols = mutate_mol(mol, db_name='replacements.db', filter_func=partial(filter_function, atom_number=9), max_size=1, max_inc=3)
```
output
```text
['Fc1ccccc1', 
 'FC(F)(F)c1ccccc1',
 'FC(F)Oc1ccccc1',
 'FC(F)Sc1ccccc1']
```

##### Custom sampling of randomly chosen fragments

If not all possible derivatives are necessary to generate a user may limit the number of returned compounds by `max_replacements` argument which enable uniform sampling of a desired number of molecules. To enable custom sampling and bias selection to more desired chemotypes there is an argument `sample_func` which takes a function implementing the selection. This function should take four necessary arguments: row_ids (list or set of row_ids from the fragment database), cursor of that fragment database, radius (int) and the number of returned items (int). An example of such a function is implemented in `utils` module and showed below. Other biases can be introduced via this option. Please not, using of complex sampling functions may slow down structure generation.
```python
def sample_csp3(row_ids, cur, radius, n):
    """
    Performs random selection of fragments proportionally to a squared fraction of sp3 carbon atoms.
    :param row_ids: the list of row ids of fragments to consider
    :param cur: cursor to the fragment database
    :param radius: context radius
    :param n: the number of fragments to select
    :return: the list of row ids of selected fragments
    """
    d = defaultdict(list)
    for rowid, core_smi, _, _ in _get_replacements(cur, radius, row_ids):
        d[core_smi].append(rowid)
    smis = list(d.keys())
    values = [rdMolDescriptors.CalcFractionCSP3(Chem.MolFromSmiles(smi)) ** 2 for smi in smis]
    values = [v + 1e-8 for v in values]
    values = np.array(values) / sum(values)
    selected_smiles = np.random.choice(smis, n, replace=False, p=values).tolist()
    ids = []
    for smi in selected_smiles:
        ids.extend(d[smi])
    ids = random.sample(ids, n)
    return ids
```
Example of `sample_func` application. F atom is replaced and 10 derivatives is returned biased by the fraction of sp3 carbons and not.
```python
m = Chem.MolFromSmiles('c1ccccc1F')

res = list(mutate_mol(m, 'replacements_sa2_f5.db',
                      radius=3, min_inc=0, max_inc=10, max_replacements=10,
                      replace_ids=[6]))
values = sorted(round(rdMolDescriptors.CalcFractionCSP3(Chem.MolFromSmiles(smi)), 4) for smi in res)
print(res)
print(values)

res = list(mutate_mol(m, 'replacements_sa2_f5.db',
                      radius=3, min_inc=0, max_inc=10, max_replacements=10,
                      replace_ids=[6], sample_func=sample_csp3))
values = sorted(round(rdMolDescriptors.CalcFractionCSP3(Chem.MolFromSmiles(smi)), 4) for smi in res)
print(res)
print(values)
```
output
```text
# uniform sampling
['c1ccc(-c2ccc(-c3csnn3)cc2)cc1', 'c1ccc(COc2cccnc2)cc1', 'CCCc1ccc(OCc2ccccc2)cc1', 'CC(C)(O)C(=O)NCCc1ccccc1', 'CN(C)C(=O)CNC(=O)OCc1ccccc1', 'COc1ccc(-c2ccccc2)cc1C(N)=O', 'Fc1ccccc1-c1ccccc1', 'Nc1cccnc1Sc1ccccc1', 'O=C(Nc1ccc(F)c(F)c1)c1ccccc1', 'O=C(COC(=O)c1ccco1)c1ccccc1']
[0.0, 0.0, 0.0, 0.0, 0.0714, 0.0769, 0.0833, 0.25, 0.3333, 0.4167]

# sampling biased by the fraction of sp3 carbons
['c1ccc(CNCc2ccncc2)cc1', 'Cc1cccc(CSc2ccccc2)n1', 'CC(=Cc1ccccc1)CN1CCN(C)CC1', 'CC(C)N(CCOc1ccccc1)C(C)C', 'CCN(CCC#N)C(=O)Nc1ccccc1', 'CSCc1ccccc1', 'O=C(Cc1ccccc1)NCCc1ccoc1', 'O=C(CCc1ccccc1)NC1CCCCC1', 'O=C(CN1CCCC1)NCCc1ccccc1', 'O=C(CSc1ccccc1)NC1CC1']
[0.1538, 0.1538, 0.2143, 0.25, 0.3333, 0.3636, 0.4667, 0.5, 0.5333, 0.5714]
```
In the latter case there are molecules having the greater fraction of sp3-carbon atoms.


##### Iterative enumeration

For convenience there is a function `enumerate_compounds` in `utils` module (added in version 0.2.6). It performs iterative growing (scaffold decoration) or mutation (analog enumeration) of a supplied molecule. More details are in docstring of the function.

Example. Enumerate derivatives of 1-chloro-3-methylbenzene at positions 2 and 4 of the ring and at the methyl group at the same time. In this case one should choose `scaffold` mode, 3 iterations, specify atom ids (0-based indices) where fragments can be attached and set `protect_added_frag=True` to restrict enumeration only to selected positions.

```python
from crem.utils import enumerate_compounds

mol = Chem.MolFromMolBlock("""
  Mrv1922 05242309182D          

  8  8  0  0  0  0            999 V2000
   -3.2813    1.3161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9957    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9957    0.0786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2813   -0.3339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5668    0.0786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5668    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2813   -1.1589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8523    1.3161    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
  4  7  1  0  0  0  0
  6  8  1  0  0  0  0
M  END
""")

mols = enumerate_compounds(mol, 'replacements_sa2.db', mode='scaffold', n_iterations=3,
                           radius=3, max_replacements=2, replace_ids=[2,4,6], protect_added_frag=True, 
                           return_smi=True)
``` 
output
```
['COc1c(C)cccc1Cl', 
'Cc1cc(Cl)ccc1Cl', 
'COc1ccc(Cl)c(OC)c1C', 
'COc1c(Cl)cccc1CF', 
'Cc1c(Cl)ccc(Cl)c1C', 
'CSCc1cc(Cl)ccc1Cl', 
'COc1ccc(Cl)c(OC)c1CC#N', 
'COCc1c(OC)ccc(Cl)c1OC', 
'COc1c(Cl)cccc1C(F)F', 
'COc1c(Cl)ccc(CO)c1CF', 
'Cc1c(Cl)ccc(Cl)c1CC#N', 
'Cc1c(Cl)ccc(Cl)c1CCN', 
'CSCc1c(Cl)ccc(Cl)c1C', 
'CSCc1c(Cl)ccc(Cl)c1Cl']
```

##### Multiprocessing
All functions have an argument `ncores` and can make mupltile replacement in one molecule in parallel. If you want to process several molecules in parallel you have to write your own code. However, the described functions are generators and cannot be used with `multiprocessing` module. Therefore, three complementary functions `mutate_mol2`, `grow_mol2` and `link_mols2` were created. They return the list with results and can be pickled and used with `multiprocessing.Pool` or other tools.

Example:
```python
from multiprocessing import Pool
from functools import partial
from crem.crem import mutate_mol2
from rdkit import Chem

p = Pool(2)
input_smi = ['c1ccccc1N', 'NCC(=O)OC', 'NCCCO']
input_mols = [Chem.MolFromSmiles(s) for s in input_smi]

res = list(p.imap(partial(mutate_mol2, db_name='replacements.db', max_size=1), input_mols))
```

`res` would be a list of lists with SMILES of generated molecules

## Bechmarks

##### Guacamol

|task|SMILES LSTM*|SMILES GA*|Graph GA*|Graph MCTS*|CReM
|---|:---:|:---:|:---:|:---:|:---:|
|Celecoxib rediscovery|**1.000**|0.732|**1.000**|0.355|**1.000**
|Troglitazone rediscovery|**1.000**|0.515|**1.000**|0.311|**1.000**
|Thiothixene rediscovery|**1.000**|0.598|**1.000**|0.311|**1.000**
|Aripiprazole similarity|**1.000**|0.834|**1.000**|0.380|**1.000**
|Albuterol similarity|**1.000**|0.907|**1.000**|0.749|**1.000**
|Mestranol similarity|**1.000**|0.79|**1.000**|0.402|**1.000**
|C11H24|**0.993**|0.829|0.971|0.410|0.966
|C9H10N2O2PF2Cl|0.879|0.889|**0.982**|0.631|0.940
|Median molecules 1|**0.438**|0.334|0.406|0.225|0.371
|Median molecules 2|0.422|0.38|0.432|0.170|**0.434**
|Osimertinib MPO|0.907|0.886|0.953|0.784|**0.995**
|Fexofenadine MPO|0.959|0.931|0.998|0.695|**1.000**
|Ranolazine MPO|0.855|0.881|0.92|0.616|**0.969**
|Perindopril MPO|0.808|0.661|0.792|0.385|**0.815**
|Amlodipine MPO|0.894|0.722|0.894|0.533|**0.902**
|Sitagliptin MPO|0.545|0.689|**0.891**|0.458|0.763
|Zaleplon MPO|0.669|0.413|0.754|0.488|**0.770**
|Valsartan SMARTS|0.978|0.552|0.990|0.04|**0.994**
|Deco Hop|0.996|0.970|**1.000**|0.590|**1.000**
|Scaffold Hop|0.998|0.885|**1.000**|0.478|**1.000**
|total score|17.341|14.398|17.983|9.011|17.919

## License
BSD-3

## Citation
CReM: chemically reasonable mutations framework for structure generation    
Pavel Polishchuk  
*Journal of Cheminformatics* **2020**, 12, (1), 28  
https://doi.org/10.1186/s13321-020-00431-w

Control of Synthetic Feasibility of Compounds Generated with CReM  
Pavel Polishchuk  
*Journal of Chemical Information and Modeling* **2020**, 60, 6074-6080  
https://dx.doi.org/10.1021/acs.jcim.0c00792
