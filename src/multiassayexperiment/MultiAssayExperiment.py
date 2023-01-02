from typing import Union, MutableMapping, Optional

from singlecellexperiment.SingleCellExperiment import SingleCellExperiment
from summarizedexperiment.SummarizedExperiment import SummarizedExperiment
from summarizedexperiment.RangeSummarizedExperiment import RangeSummarizedExperiment
from biocframe import BiocFrame

import pandas as pd

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
        smap_uniq_length = len(set(list(self._sampleMap["primary"])))

        if self._coldata.shape[0] != smap_uniq_length:
            raise ValueError(
                f"SampleMap and SampleData do not match: provided {smap_uniq_length}, needs to be {self._coldata.shape[0]}"
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
