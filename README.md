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
import pandas as pd
import numpy as np
from genomicranges import GenomicRanges

nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
gr = GenomicRanges(
    {
        "seqnames": [
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
        ]
        * 20,
        "starts": range(100, 300),
        "ends": range(110, 310),
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"] * 20,
        "score": range(0, 200),
        "GC": [random() for _ in range(10)] * 20,
    }
)

col_data_sce = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["sce"] * 6,
)

col_data_se = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["se"] * 6,
)

sample_map = pd.DataFrame(
    {
        "assay": ["sce", "se"] * 6,
        "primary": ["sample1", "sample2"] * 6,
        "colname": ["sce", "se"] * 6,
    }
)

sample_data = pd.DataFrame({"samples": ["sample1", "sample2"]})
```

Now we can create an instance of an MAE -

```python
from multiassayexperiment import MultiAssayExperiment
from singlecellexperiment import SingleCellExperiment
from summarizedExperiment import SummarizedExperiment

tsce = SingleCellExperiment(
    assays={"counts": counts}, row_data=df_gr, col_data=col_data_sce
)

tse2 = SummarizedExperiment(
    assays={"counts": counts.copy()},
    row_data=df_gr.copy(),
    col_data=col_data_se.copy(),
)

mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    col_data=sample_data,
    sample_map=sample_map,
    metadata={"could be": "anything"},
)
```

For more use cases, checkout the [documentation](https://biocpy.github.io/MultiAssayExperiment/).

<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
