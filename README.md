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
mols = list(grow_mol(m, db_name='replacements_sc2.db'))
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
