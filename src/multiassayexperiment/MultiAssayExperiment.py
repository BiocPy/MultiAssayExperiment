from typing import Any, Dict, Union
from singlecellexperiment.SingleCellExperiment import SingleCellExperiment
from summarizedexperiment.SummarizedExperiment import SummarizedExperiment
from summarizedexperiment.RangeSummarizedExperiment import (
    RangeSummarizedExperiment,
)
import pandas as pd
import logging

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class MultiAssayExperiment:
    def __init__(
        self,
        experiments: Dict[
            str,
            Union[
                SingleCellExperiment,
                SummarizedExperiment,
                RangeSummarizedExperiment,
            ],
        ],
        colData: pd.DataFrame,
        sampleMap: pd.DataFrame,
        metadata: Any = None,
    ) -> None:
        """Class for managing multiple experiments

        Args:
            experiments (Dict[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]): dictionary of experiments
            colData (pd.DataFrame, optional): sample data. Defaults to None.
            sampleMap (pd.DataFrame, optional): Mappings between sample data across experiments. Defaults to None.
            metadata (Any, optional): study level metadata. Defaults to None.
        """

        if not set(["assay", "primary", "colname"]).issubset(
            set(sampleMap.columns.tolist())
        ):
            logging.error(
                f"Sample Data does not contain columns: `assay`, `primary` and `colname`"
            )
            raise Exception(
                f"Sample Data does not contain columns: `assay`, `primary` and `colname`"
            )

        # check if unique samples is same as in sample data
        smap_uniq_length = len(sampleMap["primary"].unique())

        if colData.shape[0] != smap_uniq_length:
            logging.error(
                f"SampleMap {smap_uniq_length} and SampleData {colData.shape[0]} contain incorrect samples"
            )
            raise Exception(
                f"SampleMap {smap_uniq_length} and SampleData {colData.shape[0]} contain incorrect samples"
            )

        # check if all assay names are in experiments
        smap_unique_assaynames = sampleMap["assay"].unique().tolist()

        if not set(smap_unique_assaynames).issubset(
            set(list(experiments.keys()))
        ):
            logging.error(
                f"Not all assays {smap_unique_assaynames} in SampleMap map to experiments: {list(experiments.keys())}"
            )
            raise Exception(
                f"Not all assays {smap_unique_assaynames} in SampleMap map to experiments: {list(experiments.keys())}"
            )

        # check if colnames exist
        agroups = sampleMap.groupby(["assay"])

        for group, rows in agroups:
            if group not in experiments:
                logging.error(f"Experiment {group} does not exist")
                raise Exception(f"Experiment {group} does not exist")

            gcolData = experiments[group].colData()

            if not set(rows["colname"].unique().tolist()).issubset(
                set(gcolData.index.tolist())
            ):
                logging.error(
                    f"Assay {group} does not contain all columns in sampleMap"
                )
                raise Exception(
                    f"Assay {group} does not contain all columns in sampleMap"
                )

        self._experiments = (experiments,)
        self._coldata = (colData,)
        self._sampleMap = (sampleMap,)
        self._metadata = metadata

    def experiments(
        self,
    ) -> Dict[
        str,
        Union[
            SingleCellExperiment,
            SummarizedExperiment,
            RangeSummarizedExperiment,
        ],
    ]:
        """Accessor for experiments in the object

        Returns:
            Dict[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]: all experiments
        """

        return self._experiments

    def experiment(
        self, name: str
    ) -> Union[
        SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment
    ]:
        """Access an experiment by name

        Args:
            name (str): experiment name

        Raises:
            Exception: if experiment name does not exist

        Returns:
            Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]: experiment
        """
        if name not in self._experiments:
            logging.error(f"Experiment {name} does not exist")
            raise Exception(f"Experiment {name} does not exist")

        return self._experiments[name]

    def sampleMap(self) -> pd.DataFrame:
        """Access sample map between experimentsand sample metadata

        Returns:
            pd.DataFrame: sample map dataframe
        """
        return self._sampleMap

    def colData(self) -> pd.DataFrame:
        """Access sample metadata

        Returns:
            pd.DataFrame: sample metadata
        """
        return self._coldata

    def assays(
        self,
    ) -> Dict[
        str,
        Union[
            SingleCellExperiment,
            SummarizedExperiment,
            RangeSummarizedExperiment,
        ],
    ]:
        """Accessor for experiments in the object

        Returns:
            Dict[str, Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]: all experiments
        """

        return self.experiments()

    def metadata(self) -> Any:
        """Accessor for metadata

        Returns:
            Any: metadata slot
        """
        return self._metadata

    @staticmethod
    def fromMuData(data) -> "MultiAssayExperiment":
        """Transform MuData object to MAE representation

        Args:
            data (MuData): MuData object

        Returns:
            MultiAssayExperiment: MAE representation
        """
        raise NotImplementedError

    @staticmethod
    def fromAnnData(data) -> "MultiAssayExperiment":
        """Transform AnnData object to MAE representation

        Args:
            data (AnnData): MuData object

        Returns:
            MultiAssayExperiment: MAE representation
        """
        raise NotImplementedError
