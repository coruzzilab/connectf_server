from collections import OrderedDict
from contextlib import closing
from string import whitespace
from typing import TextIO, Tuple

import pandas as pd


def get_gene_lists(f: TextIO) -> Tuple[pd.DataFrame, OrderedDict]:
    gene_to_name = OrderedDict()
    name_to_gene = OrderedDict()

    with closing(f) as gene_file:
        list_name = 'default'
        for line in gene_file:
            line = line.strip()
            if line.startswith('>'):
                list_name = line.strip('>' + whitespace)
            else:
                gene_to_name.setdefault(line, set()).add(list_name)
                name_to_gene.setdefault(list_name, set()).add(line)

    df = pd.DataFrame(
        ((key, ', '.join(val), len(val)) for key, val in gene_to_name.items()),
        columns=['TARGET', 'User List', 'User List Count']
    )

    df = df.set_index('TARGET')

    return df, name_to_gene
