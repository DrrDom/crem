import sys
from rdkit import Chem

if __name__ == '__main__':

    fname = sys.argv[1]
    smi = sys.argv[2]
    mol = Chem.MolFromSmiles(smi)

    print(smi + ' OK')

    with open(fname, 'wt') as f:
        f.write(Chem.MolToSmiles(mol) + '\n')
        f.write(smi + "\n")
        f.write(str(sys.argv))
        f.write('\n\n')
        f.write(sys.version)
        f.write('\n\n')
        if smi == 'C1COCC1':
            raise Exception('my test error')
        f.write("Success")

