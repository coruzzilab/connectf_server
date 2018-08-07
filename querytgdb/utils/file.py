from collections import defaultdict
from contextlib import closing
from string import whitespace
from typing import Dict, TextIO, Tuple

import pandas as pd


def get_gene_lists(f: TextIO) -> Tuple[pd.DataFrame, Dict]:
    gene_to_name = defaultdict(set)
    name_to_gene = defaultdict(set)

    with closing(f) as gene_file:
        list_name = 'default'
        for line in gene_file:
            line = line.strip()
            if line.startswith('>'):
                list_name = line.strip('>' + whitespace)
            else:
                gene_to_name[line].add(list_name)
                name_to_gene[list_name].add(line)

    df = pd.DataFrame(
        ((key, ', '.join(val), len(val)) for key, val in gene_to_name.items()),
        columns=['TARGET', 'User List', 'User List Count']
    )

    df = df.set_index('TARGET')

    return df, name_to_gene
