import os
import datetime
import pandas as pd
import re
from collections import defaultdict


base_dir = '/home/pavlop/imtm/crem/generation'
dirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]

# dirs = ['/home/pavel/python/crem/test/data_collect']

for dirname in dirs:

    print(datetime.datetime.now(), dirname)

    files = [f.rsplit('.', 1)[0] for f in os.listdir(dirname) if os.path.isfile(os.path.join(dirname, f)) and re.match('out_i[0-9]+_j[0-9]+.sa', f)]

    files_dict = defaultdict(list)
    for f in files:
        files_dict[int(re.search('(?<=i)[0-9]+', f).group(0))].append(f)

    for i, (k, v) in enumerate(sorted(files_dict.items())):

        df = []

        for f in v:
            a = pd.read_table(os.path.join(dirname, f + '.txt'), index_col=1)
            a = a.drop(['transformation'], axis=1)
            b = pd.read_table(os.path.join(dirname, f + '.sa'), index_col=1)
            b = b.drop(['SMILES'], axis=1)
            df.append(pd.concat([a, b], axis=1))

        if len(df) > 1:
            df = df[0].append(df[1:])
        else:
            df = df[0]
        df = df.drop_duplicates(['Iteration', 'SMILES'])

        if i == 0:
            df.to_csv(os.path.join(dirname, 'res.out'), sep="\t", index_label='Name')
        else:
            df.to_csv(os.path.join(dirname, 'res.out'), sep="\t", header=False, mode='a')
