# Tutorial

`MultiAssayExperiment` is a container class to represent multiple experiments represented as `SingleCellExperiment`, `SummarizedExperiment` or `RangeSummarizedExperiment`. For more detailed description checkout the [Bioc MultiAssayExperiment R package](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html))

## Mock data 

we first create a mock dataset of 200 rows and 6 columns, also adding a cell annotations, sample mapping and sample data across experiments.

```python
import pandas as pd
import numpy as np
from genomicranges import GenomicRanges

nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
df_gr = pd.DataFrame(
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

gr = GenomicRanges.fromPandas(df_gr)

colData_sce = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["sce"] * 6,
)

colData_se = pd.DataFrame(
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

### `MultiAssayExperiment`


Lets first create various experiment classes

```python
from singlecellexperiment import SingleCellExperiment
from summarizedExperiment import SummarizedExperiment

tsce = SingleCellExperiment(
    assays={"counts": counts}, rowData=df_gr, colData=colData_sce
)

tse2 = SummarizedExperiment(
    assays={"counts": counts.copy()},
    rowData=df_gr.copy(),
    colData=colData_se.copy(),
)
```

Now that we have all the pieces together, we can now create an MAE,


```python
from multiassayexperiment import MultiAssayExperiment

mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    colData=sample_data,
    sampleMap=sample_map,
    metadata={"could be": "anything"},
)
```

### Accessors

Multiple methods are available to access various slots of a `MultiAssayExperiment` object

```python
tse.assays()
tse.colData()
tse.SampleMap()
tse.experiments()
tse.metadata()
```

### Access specific sets

You can also specify a specific experiment to access

```python
tse.experiment("sce")
```

## Export and import AnnData objects

coming soon...