# Tutorial

Container class to represent multiple experiments and assays performed over a set of samples. For more detailed description checkout the [MultiAssayExperiment Bioc/R package](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html))

# Construct an `MultiAssayExperiment`

An MAE contains three main entities

- Primary information (`coldata`): Bio-specimen information on which experiments were run. represented as a Pandas `DataFrame`.
- Experiments (`experiments`): genomic data from each experiment. represented as `SingleCellExperiment`, `SummarizedExperiment`, `RangeSummarizedExperiment`.
- Sample Mapping (`sampleMap`): Mapping from biospecimens to samples/columns in each experiment. in the context of single cell, these are cells.

Lets create these objects

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

Then, create various experiment classes,

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

maeObj = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    colData=sample_data,
    sampleMap=sample_map,
    metadata={"could be": "anything"},
)
```

To make your life easier, we also provide methods to naively create sample mapping from experiments.

**_This is not a recommended approach, but if you don't have sample mapping, then it doesn't matter._**

```python
maeObj = mae.makeMAE(experiments={"sce": tsce, "se": tse2})
```

## Import `MuData` and `AnnData` as `MultiAssayExperiment`

If you have a dataset stored as `MuData`, these can be easily converted to an MAE using the `fromMuData` method.

Lets first construct AnnData objects and then an MAE

```python
import multiassayexperiment as mae
import numpy as np
from anndata import AnnData

np.random.seed(1)

n, d, k = 1000, 100, 10

z = np.random.normal(loc=np.arange(k), scale=np.arange(k) * 2, size=(n, k))
w = np.random.normal(size=(d, k))
y = np.dot(z, w.T)

adata = AnnData(y)
adata.obs_names = [f"obs_{i+1}" for i in range(n)]
adata.var_names = [f"var_{j+1}" for j in range(d)]

d2 = 50
w2 = np.random.normal(size=(d2, k))
y2 = np.dot(z, w2.T)

adata2 = AnnData(y2)
adata2.obs_names = [f"obs_{i+1}" for i in range(n)]
adata2.var_names = [f"var2_{j+1}" for j in range(d2)]
```

we can now construct a `MuData` object and convert that to an MAE

```python
mdata = MuData({"rna": adata, "spatial": adata2})

maeObj = mae.fromMuData(mudata=mdata)
```

Methods are also available to convert an `AnnData` object to `MAE`.

```python
maeObj = mae.readH5AD("tests/data/adata.h5ad")
```

# Accessors

Multiple methods are available to access various slots of a `MultiAssayExperiment` object

```python
maeObj.assays
maeObj.colData
maeObj.sampleMap
maeObj.experiments
maeObj.metadata
```

## Access experiments

if you want to access a specific experiment

```python
# access a specific experiment
maeObj.experiment(experiment_name)
```

This does not include the sample data stored in the MAE. If you want to include this information

***Note: This creates a copy of the experiment object.***

```python
expt_with_sampleData = maeObj.experiment(experiment_name, withSampleData=True)
```

# Slice a `MultiAssayExperiment`

`MultiAssayExperiment` allows subsetting by `rows`, `columns`, and `experiments`. Samples are automatically sliced during this operation.

The structure for slicing,

```
maeObj[rows, columns, experiments]
```

- rows, columns: accepts either a slice, list of indices or a dictionary to specify slices per experiment.
- experiments: accepts a list of experiment names to subset to.

## Slice by row and column slices

```python
maeObj[1:5, 0:4]
```

## Slice by rows, columns, experiments

```python
maeObj[1:5, 0:4, ["spatial"]]
```

## Specify slices per experiment

You can specify slices by experiment, rest of the experiments are not sliced.

```python
maeObj[{"rna": slice(0,10)}, {"spatial": slice(0,5)}, ["spatial"]]
```

Checkout other methods that perform similar operations - `subsetByRows`, `subsetByColumns` & `subsetByExperiments`.

# Helper methods

## completedCases

This method returns a boolean vector that specifies which biospecimens have data across all experiments.

```python
maeObj.completedCases()
```

## replicated

replicated identifies biospecimens that have multiple observations per experiment.

```python
maeObj.replicated()
```
