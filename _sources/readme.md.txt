[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)
[![PyPI-Server](https://img.shields.io/pypi/v/MultiAssayExperiment.svg)](https://pypi.org/project/MultiAssayExperiment/)
![Unit tests](https://github.com/BiocPy/MultiAssayExperiment/actions/workflows/pypi-test.yml/badge.svg)

# MultiAssayExperiment

Container class to represent and manage multi-omics genomic experiments. Follows Bioconductor's [MAE R/Package](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html).

## Install

Package is published to [PyPI](https://pypi.org/project/multiassayexperiment/)

```shell
pip install multiassayexperiment
```

## Usage

First create mock sample data

```python
from random import random

import numpy as np
from biocframe import BiocFrame
from genomicranges import GenomicRanges
from iranges import IRanges

nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
gr = GenomicRanges(
    seqnames=[
            "chr1",
            "chr2",
            "chr2",
            "chr2",
            "chr1",
            "chr1",
            "chr3",
            "chr3",
            "chr3",
            "chr3",
        ] * 20,
    ranges=IRanges(range(100, 300), range(110, 310)),
    strand = ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"] * 20,
    mcols=BiocFrame({
        "score": range(0, 200),
        "GC": [random() for _ in range(10)] * 20,
    })
)

col_data_sce = BiocFrame({"treatment": ["ChIP", "Input"] * 3},
    row_names=[f"sce_{i}" for i in range(6)],
)

col_data_se = BiocFrame({"treatment": ["ChIP", "Input"] * 3},
    row_names=[f"se_{i}" for i in range(6)],
)

sample_map = BiocFrame({
    "assay": ["sce", "se"] * 6,
    "primary": ["sample1", "sample2"] * 6,
    "colname": ["sce_0", "se_0", "sce_1", "se_1", "sce_2", "se_2", "sce_3", "se_3", "sce_4", "se_4", "sce_5", "se_5"]
})

sample_data = BiocFrame({"samples": ["sample1", "sample2"]}, row_names= ["sample1", "sample2"])
```

Now we can create an instance of an MAE -

```python
from multiassayexperiment import MultiAssayExperiment
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

tsce = SingleCellExperiment(
    assays={"counts": counts}, row_data=gr.to_pandas(), column_data=col_data_sce
)

tse2 = SummarizedExperiment(
    assays={"counts": counts.copy()},
    row_data=gr.to_pandas().copy(),
    column_data=col_data_se.copy(),
)

mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    column_data=sample_data,
    sample_map=sample_map,
    metadata={"could be": "anything"},
)
```

    ## output
    class: MultiAssayExperiment containing 2 experiments
    [0] sce: SingleCellExperiment with 200 rows and 6 columns
    [1] se: SummarizedExperiment with 200 rows and 6 columns
    column_data columns(1): ['samples']
    sample_map columns(3): ['assay', 'primary', 'colname']
    metadata(1): could be

For more use cases, checkout the [documentation](https://biocpy.github.io/MultiAssayExperiment/).

<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
