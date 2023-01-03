from typing import Union, MutableMapping, Optional, Tuple, Sequence

from singlecellexperiment.SingleCellExperiment import SingleCellExperiment
from summarizedexperiment.SummarizedExperiment import SummarizedExperiment
from summarizedexperiment.RangeSummarizedExperiment import RangeSummarizedExperiment
from biocframe import BiocFrame

import pandas as pd
from collections import OrderedDict

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class MultiAssayExperiment:
    def __init__(
        self,
        experiments: MutableMapping[
            str,
            Union[
                SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment,
            ],
        ],
        colData: pd.DataFrame,
        sampleMap: pd.DataFrame,
        metadata: Optional[MutableMapping] = None,
    ) -> None:
        """Class for managing multi-modal and multi-sample genomic experiments

        Args:
            experiments (MutableMapping[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]): 
                dictionary of experiments
            colData (pd.DataFrame]): sample data. 
            sampleMap (pd.DataFrame): Mappings between sample data across experiments. 
                Must contain columns `assay`, `primary` and `colname`.
                For more info, checkout [MAE docs](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html).
            metadata (MutableMapping, optional): study level metadata. Defaults to None.
        """
        self._experiments = experiments
        self._coldata = colData
        self._sampleMap = sampleMap
        self._metadata = metadata

        self._validate()

    def _validate(self):
        """Internal method to validate the object

        Raises:
            ValueError: when attributes don't match expectations
        """
        if not isinstance(self._experiments, dict):
            raise TypeError("experiments must be an instance of dict")

        if not (
            isinstance(self._sampleMap, BiocFrame)
            or isinstance(self._sampleMap, pd.DataFrame)
        ):
            raise TypeError(
                "sampleMap must be either a pandas dataframe or a biocframe object"
            )

        if not (
            isinstance(self._coldata, BiocFrame)
            or isinstance(self._coldata, pd.DataFrame)
        ):
            raise TypeError(
                "colData must be either a pandas dataframe or a biocframe object"
            )

        if not set(["assay", "primary", "colname"]).issubset(
            set(list(self._sampleMap.columns))
        ):
            raise ValueError(
                f"Sample data does not contain required columns: `assay`, `primary` and `colname`"
            )

        # check if unique samples is same as in sample data
        smaps_list = list(self._sampleMap["primary"])
        smap_uniq_length = len(set(smaps_list))

        if self._coldata.shape[0] != smap_uniq_length:
            raise ValueError(
                f"SampleMap and SampleData do not match: provided {smap_uniq_length}, needs to be {self._coldata.shape[0]}"
            )

        # check if coldata has index
        if self._coldata.index is None:
            raise ValueError(
                "SampleData must contain an index with all sample names (primary column) from SampleMap"
            )

        missing = set(smaps_list).difference(set(self._coldata.index.tolist()))
        if len(missing) > 0:
            raise ValueError(
                f"SampleData contains missing samples from SampleMap: {missing}"
            )

        # check if all assay names are in experiments
        smap_unique_assaynames = set(self._sampleMap["assay"])

        if not smap_unique_assaynames.issubset(set(list(self._experiments.keys()))):
            raise ValueError(
                f"Not all assays {smap_unique_assaynames} in sampleMap map to experiments: {list(self._experiments.keys())}"
            )

        # check if colnames exist
        agroups = self._sampleMap.groupby(["assay"])
        for group, rows in agroups:
            if group not in self._experiments:
                raise ValueError(f"Experiment {group} does not exist")

            gcolData = self._experiments[group].colData

            if not set(rows["colname"].unique().tolist()).issubset(
                set(gcolData.index.tolist())
            ):
                raise ValueError(
                    f"Assay {group} does not contain all columns in sampleMap"
                )

    @property
    def experiments(
        self,
    ) -> MutableMapping[
        str,
        Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment,],
    ]:
        """Get experiments.

        Returns:
            MutableMapping[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]: all experiments
        """

        return self._experiments

    @experiments.setter
    def experiments(
        self,
        expts: MutableMapping[
            str,
            Union[
                SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment
            ],
        ],
    ):
        """Set new experiments.

        Args:
            expts (MutableMapping[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]): new experiments dictionary to set.
        """
        if not isinstance(expts, dict):
            raise TypeError("expts must be a dictionary like object")

        self._experiments = expts
        self._validate()

    def experiment(
        self, name: str
    ) -> Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]:
        """Get experiment by name

        Args:
            name (str): experiment name

        Raises:
            ValueError: if experiment name does not exist

        Returns:
            Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]: experiment
        """
        if name not in self._experiments:
            raise ValueError(f"Experiment {name} does not exist")

        return self._experiments[name]

    @property
    def sampleMap(self) -> pd.DataFrame:
        """Get sample map between experiments and sample metadata

        Returns:
            pd.DataFrame: sample map dataframe
        """
        return self._sampleMap

    @sampleMap.setter
    def sampleMap(self, sampleMap: pd.DataFrame):
        """Set new sample map
        """
        if not isinstance(sampleMap, pd.DataFrame):
            raise TypeError("sample mapping must be a pandas dataframe")

        self._sampleMap = sampleMap
        self._validate()

    @property
    def colData(self) -> pd.DataFrame:
        """Get sample metadata

        Returns:
            pd.DataFrame: sample metadata
        """
        return self._coldata

    @colData.setter
    def colData(self, colData: pd.DataFrame):
        """Set new sample metadata
        """
        if not isinstance(colData, pd.DataFrame):
            raise TypeError("sample metadata must be a pandas dataframe")

        self._coldata = colData
        self._validate()

    @property
    def assays(
        self,
    ) -> MutableMapping[
        str,
        Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment,],
    ]:
        """Get experiments

        Returns:
            MutableMapping[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]: all experiments
        """
        return self.experiments

    @property
    def metadata(self) -> Optional[MutableMapping]:
        """Get metadata.

        Returns:
            Optional[MutableMapping]: metadata if available
        """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: MutableMapping):
        """Set metadata.

        Args:
            metadata (MutableMapping): new metadata tobject
        """
        self._metadata = metadata

    def _subsetExpt(
        self, subset: Union[str, Sequence[str]]
    ) -> MutableMapping[
        str,
        Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment,],
    ]:
        """Internal method to subset experiments.

        Args:
            subset (Sequence[str]): list of experiments to keep.

        Returns:
            MutableMapping[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]: a dictionary with subset experiments
        """
        if isinstance(subset, str):
            subset = [subset]

        if not all([isinstance(x, str) for x in subset]):
            raise ValueError("all experiment slices must be strings")

        new_expt = OrderedDict()

        for texpt in subset:
            if texpt not in self._experiments:
                raise ValueError(
                    f"experiment {texpt} does not exist. should be {list(self._experiments.keys())}"
                )
            new_expt[texpt] = self._experiments[texpt]

        return new_expt

    def _slice(
        self,
        args: Tuple[
            Union[Sequence[int], slice],
            Optional[Union[Sequence[int], slice]],
            Optional[Sequence[str]],
        ],
    ) -> Tuple[
        MutableMapping[
            str,
            Union[
                SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment,
            ],
        ],
        pd.DataFrame,
        pd.DataFrame,
    ]:
        """Internal method to slice `MAE` by index.

        Args:
            args (Tuple[Union[Sequence[int], slice], Optional[Union[Sequence[int], slice]], Optional[Sequence[str]]]): indices to slice. tuple can
                contains slices along dimensions (rows, columns, experiments).

        Raises:
            ValueError: Too many slices

        Returns:
             Tuple[MutableMapping[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]], pd.DataFrame, pd.DataFrame]: 
                sliced row, cols and assays.
        """

        if len(args) == 0:
            raise ValueError("Arguments must contain atleast one slice")

        print(args)

        rowIndices = args[0]
        colIndices = None
        exptIndices = None

        print(len(args))

        if len(args) > 1:
            colIndices = args[1]
        
        if len(args) > 2:
            print("in exptindices")
            exptIndices = args[2]

            if exptIndices is not None:
                if isinstance(exptIndices, str):
                    exptIndices = [exptIndices]

                if not all([isinstance(x, str) for x in exptIndices]):
                    raise ValueError("all assay slices must be strings")
        
        if len(args) > 3:
            raise ValueError("contains too many slices")

        subset_expts = self._experiments.copy()

        print(rowIndices, colIndices, exptIndices)

        if exptIndices is not None:
            subset_expts = self._subsetExpt(exptIndices)

        if rowIndices is not None:
            for expname, expt in subset_expts.items():
                subset_expts[expname] = expt[rowIndices, :]

        if colIndices is not None:
            for expname, expt in subset_expts.items():
                subset_expts[expname] = expt[:, colIndices]

        subset_colnames = []
        for expname, expt in subset_expts.items():
            subset_colnames.extend(expt.colnames)

        # filter sampleMap
        subset_sampleMap = self._sampleMap[
            (self._sampleMap["assay"].isin(list(subset_expts.keys())))
            & (self._sampleMap["colname"].isin(subset_colnames))
        ]

        # filter coldata
        subset_coldata = self._coldata[
            self._coldata.index.isin(subset_sampleMap["primary"].unique().tolist())
        ]

        return (subset_expts, subset_sampleMap, subset_coldata)

    def subsetByExperiments(
        self, subset: Union[str, Sequence[str]]
    ) -> "MultiAssayExperiment":
        """Subset by experiment(s).

        Args:
            subset (Union[str, Sequence[str]]): experiment or experiments list to subset

        Returns:
            MultiAssayExperiment: a new MAE with the subset.
        """
        expt, smap, sdata = self._slice(args=(None, None, subset))
        return MultiAssayExperiment(expt, sdata, smap, self._metadata)

    def subsetByRow(
        self, subset: Tuple[Union[Sequence[int], slice]]
    ) -> "MultiAssayExperiment":
        """Subset by rows.

        Args:
            subset (Tuple[Union[Sequence[int], slice]]): row indices or slice to subset.

        Returns:
            MultiAssayExperiment: a new MAE with the subset.
        """
        expt, smap, sdata = self._slice(args=(subset, None, None))
        return MultiAssayExperiment(expt, sdata, smap, self._metadata)

    def subsetByColumn(
        self, subset: Tuple[Union[Sequence[int], slice]]
    ) -> "MultiAssayExperiment":
        """Subset by column.

        Args:
            subset (Tuple[Union[Sequence[int], slice]]): column indices or slice to subset.

        Returns:
            MultiAssayExperiment: a new MAE with the subset.
        """
        expt, smap, sdata = self._slice(args=(None, subset, None))
        return MultiAssayExperiment(expt, sdata, smap, self._metadata)

    def __getitem__(
        self,
        args: Tuple[
            Union[Sequence[int], slice],
            Optional[Union[Sequence[int], slice]],
            Optional[Sequence[str]],
        ],
    ) -> "MultiAssayExperiment":
        """Subset a `MultiAssayExperiment`. supports a tuple specifying slices along (rows, columns and experiments).

        Args:
            args (Tuple[Union[Sequence[int], slice], Optional[Union[Sequence[int], slice]], Optional[str]]): indices to slice. tuple can
                contains slices along dimensions (row, column, experiments)

        Raises:
            ValueError: Too many slices

        Returns:
            MultiAssayExperiment: new sliced `MultiAssayExperiment` object
        """
        expt, smap, sdata = self._slice(args=args)
        return MultiAssayExperiment(expt, sdata, smap, self._metadata)
