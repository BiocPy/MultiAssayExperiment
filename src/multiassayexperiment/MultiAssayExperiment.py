from collections import OrderedDict
from copy import deepcopy
from typing import MutableMapping, Optional, Sequence

import pandas as pd
from mudata import MuData
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment
from summarizedexperiment._type_checks import is_bioc_or_pandas_frame, is_list_of_type
from summarizedexperiment.BaseSE import BaseSE

from .types import SlicerArgTypes, SlicerResult, SlicerTypes, StrOrListStr

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class MultiAssayExperiment:
    def __init__(
        self,
        experiments: MutableMapping[str, BaseSE],
        colData: pd.DataFrame,
        sampleMap: pd.DataFrame,
        metadata: Optional[MutableMapping] = None,
    ) -> None:
        """Class for managing multi-modal and multi-sample genomic experiments

        Args:
            experiments (MutableMapping[str, BaseSE]):
                dictionary of experiments
            colData (pd.DataFrame]): sample data.
            sampleMap (pd.DataFrame): Mappings between sample data across experiments.
                Must contain columns `assay`, `primary` and `colname`.
                For more info, checkout [MAE docs](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html).
            metadata (MutableMapping, optional): study level metadata. Defaults to None.
        """
        self._validate_experiments(experiments)
        self._validate_sampleMap(
            sampleMap=sampleMap, colData=colData, experiments=experiments
        )
        self._sampleMap = sampleMap
        self._colData = colData
        self._experiments = experiments

        self._metadata = metadata

    def _validate_experiments(self, experiments: MutableMapping[str, BaseSE]):
        """Internal method to validate experiments.

        Raises:
            TypeError: if experiments is not a dict.
        """
        if not isinstance(experiments, dict):
            raise TypeError("experiments must be an instance of dict")

    def _validate_colData(self, colData: pd.DataFrame):
        """Internal method to validate coldata.

        Args:
            colData (pd.DataFrame): column data.

        Raises:
            TypeError: if object is neither a dataframe nor biocframe.
        """
        if not is_bioc_or_pandas_frame(colData):
            raise TypeError(
                "colData must be either a pandas dataframe or a biocframe object"
            )

    def _validate_sampleMap_with_colData(
        self, sampleMap: pd.DataFrame, colData: pd.DataFrame
    ):
        """Internal method to validate sample mapping and coldata.

        Args:
            sampleMap (pd.DataFrame): sample mapping.
            colData (pd.DataFrame): column data.

        Raises:
            ValueError: if any of the checks fail.
        """
        # check if unique samples is same as in sample data
        smapsList = list(sampleMap["primary"])
        smapUniqLength = len(set(smapsList))

        if colData.shape[0] != smapUniqLength:
            raise ValueError(
                f"SampleMap and SampleData do not match: provided {smapUniqLength}, needs to be {colData.shape[0]}"
            )

        # check if coldata has index
        if colData.index is None:
            raise ValueError(
                "SampleData must contain an index with all sample names (primary column) from SampleMap"
            )

        missing = set(smapsList).difference(set(colData.index.tolist()))
        if len(missing) > 0:
            raise ValueError(
                f"SampleData contains missing samples from SampleMap: {missing}"
            )

    def _validate_sampleMap_with_Expts(
        self, sampleMap: pd.DataFrame, experiments: MutableMapping[str, BaseSE]
    ):
        """Internal method to validate sample map and experiments

        Args:
            sampleMap (pd.DataFrame): sample mapping.
            experiments (MutableMapping[str, BaseSE]): experiments.

        Raises:
            ValueError: if any of the checks fail.
        """
        # check if all assay names are in experiments
        smapUniqueAssaynames = set(sampleMap["assay"].unique())
        UniqueExperimentname = set(list(experiments.keys()))

        if not UniqueExperimentname.issubset(smapUniqueAssaynames):
            raise ValueError(
                f"Not all primary assays {smapUniqueAssaynames} in `sampleMap` map to experiments: {list(experiments.keys())}"
            )

        # check if colnames exist
        agroups = sampleMap.groupby(["assay"])
        for group, rows in agroups:
            if group not in experiments:
                raise ValueError(f"Experiment {group} does not exist")

            gcolData = experiments[group].colData

            if not set(rows["colname"].unique().tolist()).issubset(
                set(gcolData.index.tolist())
            ):
                raise ValueError(
                    f"Assay {group} does not contain all columns in sampleMap"
                )

    def _validate_sampleMap(
        self,
        sampleMap: pd.DataFrame,
        colData: pd.DataFrame,
        experiments: MutableMapping[str, BaseSE],
    ):
        """Validate sample mapping.

        Args:
            sampleMap (pd.DataFrame): sample mapping.
            colData (pd.DataFrame): column data.
            experiments (MutableMapping[str, BaseSE]): experiments.

        Raises:
            TypeError, ValueError: any of the checks fail.
        """
        if not is_bioc_or_pandas_frame(sampleMap):
            raise TypeError(
                "sampleMap must be either a pandas dataframe or a biocframe object"
            )

        if not set(["assay", "primary", "colname"]).issubset(
            set(list(sampleMap.columns))
        ):
            raise ValueError(
                "Sample data does not contain required columns: `assay`, `primary` and `colname`"
            )

        self._validate_sampleMap_with_colData(sampleMap, colData)
        self._validate_sampleMap_with_Expts(sampleMap, experiments)

    def _validate(self):
        """Internal method to validate the object

        Raises:
            ValueError: when attributes don't match expectations
        """
        self._validate_experiments(self._experiments)
        self._validate_colData(self._colData)
        self._validate_sampleMap(self._sampleMap, self._colData, self._experiments)

    @property
    def experiments(
        self,
    ) -> MutableMapping[str, BaseSE]:
        """Get experiments.

        Returns:
            MutableMapping[str, BaseSE]: all experiments
        """

        return self._experiments

    @experiments.setter
    def experiments(
        self,
        experiments: MutableMapping[str, BaseSE],
    ):
        """Set new experiments.

        Args:
            experiments (MutableMapping[str, BaseSE]): new experiments to set.
        """

        self._validate_experiments(experiments)
        self._validate_sampleMap_with_Expts(self._sampleMap, experiments)
        self._experiments = experiments

    def experiment(self, name: str, withSampleData: bool = False) -> BaseSE:
        """Get experiment by name.

        if withSampleData is True, a copy of the experiment object is returned.

        Args:
            name (str): experiment name.
            with_sampleData (bool, optional): include sample data in returned object?
                Defaults to False.

        Raises:
            ValueError: if experiment name does not exist.

        Returns:
            BaseSE: experiment.
        """
        if name not in self._experiments:
            raise ValueError(f"Experiment {name} does not exist")

        expt = self._experiments[name]

        if withSampleData is True:
            expt = deepcopy(self._experiments[name])

            subset_map = self.sampleMap[self.sampleMap["assay"] == name]
            subset_map = subset_map.set_index("colname")

            expt_colData = expt.colData
            new_colData = pd.concat([subset_map, expt_colData], axis=1)
            expt.colData = new_colData

        return expt

    @property
    def sampleMap(self) -> pd.DataFrame:
        """Get sample map between experiments and sample metadata.

        Returns:
            pd.DataFrame: sample map.
        """
        return self._sampleMap

    @sampleMap.setter
    def sampleMap(self, sampleMap: pd.DataFrame):
        """Set new sample mapping.

        Args:
            sampleMap (pd.DataFrame): new sample map.
        """
        self._validate_sampleMap(sampleMap, self._colData, self._experiments)
        self._sampleMap = sampleMap

    @property
    def colData(self) -> pd.DataFrame:
        """Get sample metadata.

        Returns:
            pd.DataFrame: sample metadata.
        """
        return self._colData

    @colData.setter
    def colData(self, colData: pd.DataFrame):
        """Set sample metadata.

        Args:
            colData (pd.DataFrame): new metadata.
        """
        self._validate_colData(colData)
        self._validate_sampleMap_with_colData(self._sampleMap, colData)
        self._colData = colData

    @property
    def assays(
        self,
    ) -> MutableMapping[str, BaseSE,]:
        """Get experiments.

        Alias to the experiment property.

        Returns:
            MutableMapping[str, BaseSE]: all experiments.
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

    def _subsetExpt(self, subset: StrOrListStr) -> MutableMapping[str, BaseSE]:
        """Internal method to subset experiments.

        Args:
            subset (Sequence[str]): list of experiments to keep.

        Returns:
            MutableMapping[str, BaseSE]: a dictionary with subset experiments.
        """
        if isinstance(subset, str):
            subset = [subset]

        if not is_list_of_type(subset, str):
            raise ValueError("all provided experiment names must be strings")

        newExpt = OrderedDict()

        for texpt in subset:
            if texpt not in self.experiments:
                raise ValueError(
                    f"experiment {texpt} does not exist. should be {list(self.experiments.keys())}"
                )
            newExpt[texpt] = self.experiments[texpt]

        return newExpt

    def _slice(
        self,
        args: SlicerArgTypes,
    ) -> SlicerResult:
        """Internal method to slice `MAE` by index.

        Args:
            args (SlicerArgTypes): indices to slice. tuple can
                contains slices along dimensions (rows, columns, experiments).

        Raises:
            ValueError: Too many or too few slices.

        Returns:
            SlicerResult: sliced row, cols and assays.
        """

        if len(args) == 0:
            raise ValueError("Arguments must contain atleast one slice")

        rowIndices = args[0]
        colIndices = None
        exptIndices = None

        if len(args) > 1:
            colIndices = args[1]

        if len(args) > 2:
            exptIndices = args[2]

            if exptIndices is not None:
                if isinstance(exptIndices, str):
                    exptIndices = [exptIndices]

                if not all([isinstance(x, str) for x in exptIndices]):
                    raise ValueError("all assay slices must be strings")

        if len(args) > 3:
            raise ValueError("contains too many slices")

        subsetExpts = self.experiments.copy()

        if exptIndices is not None:
            subsetExpts = self._subsetExpt(exptIndices)

        if rowIndices is not None:
            if isinstance(rowIndices, dict):
                incorrect = set(list(rowIndices.keys())).difference(
                    list(subsetExpts.keys())
                )
                if len(incorrect) > 0:
                    raise ValueError(f"Incorrect experiment name provided: {incorrect}")

                for expname, expt in subsetExpts.items():
                    if expname in rowIndices:
                        subsetExpts[expname] = expt[rowIndices[expname], :]
                    else:
                        subsetExpts[expname] = expt

            elif isinstance(rowIndices, slice) or all(
                isinstance(x, int) for x in rowIndices
            ):
                for expname, expt in subsetExpts.items():
                    subsetExpts[expname] = expt[rowIndices, :]
            else:
                raise TypeError(
                    "slice for rows is not an expected type. It should be either a dict"
                )

        if colIndices is not None:
            if isinstance(colIndices, dict):
                incorrect = set(list(colIndices.keys())).difference(
                    list(subsetExpts.keys())
                )
                if len(incorrect) > 0:
                    raise ValueError(f"Incorrect experiment name provided: {incorrect}")

                for expname, expt in subsetExpts.items():
                    if expname in colIndices:
                        subsetExpts[expname] = expt[:, colIndices[expname]]
                    else:
                        subsetExpts[expname] = expt

            elif isinstance(colIndices, slice) or all(
                isinstance(x, int) for x in colIndices
            ):
                for expname, expt in subsetExpts.items():
                    subsetExpts[expname] = expt[:, colIndices]
            else:
                raise TypeError(
                    "slice for columns is not an expected type. It should be either a dict"
                )

        # filter sampleMap
        subsetColnames = []
        subsetSampleMap = pd.DataFrame()
        for expname, expt in subsetExpts.items():
            subsetColnames.extend(expt.colnames)
            subsetSampleMap = pd.concat(
                [
                    subsetSampleMap,
                    self.sampleMap[
                        (self.sampleMap["assay"] == expname)
                        & (self.sampleMap["colname"].isin(expt.colnames))
                    ],
                ]
            )

        # filter coldata
        subsetColdata = self.colData[
            self.colData.index.isin(subsetSampleMap["primary"].unique().tolist())
        ]

        return SlicerResult(subsetExpts, subsetSampleMap, subsetColdata)

    def subsetByExperiments(self, subset: StrOrListStr) -> "MultiAssayExperiment":
        """Subset by experiment(s).

        Args:
            subset (StrOrListStr): experiment or experiments list to subset.

        Returns:
            MultiAssayExperiment: a new `MultiAssayExperiment` with the subset.
        """
        sresult = self._slice(args=(None, None, subset))
        return MultiAssayExperiment(
            sresult.experiments, sresult.colData, sresult.sampleMap, self.metadata
        )

    def subsetByRow(self, subset: SlicerTypes) -> "MultiAssayExperiment":
        """Subset by rows.

        Args:
            subset (SlicerTypes): column indices or slice to subset.

        Returns:
            MultiAssayExperiment: a new `MultiAssayExperiment` with the subset.
        """
        sresult = self._slice(args=(subset, None, None))
        return MultiAssayExperiment(
            sresult.experiments, sresult.colData, sresult.sampleMap, self.metadata
        )

    def subsetByColumn(self, subset: SlicerTypes) -> "MultiAssayExperiment":
        """Subset by column.

        Args:
            subset (SlicerTypes): column indices or slice to subset.

        Returns:
            MultiAssayExperiment: a new `MultiAssayExperiment` with the subset.
        """
        sresult = self._slice(args=(None, subset, None))
        return MultiAssayExperiment(
            sresult.experiments, sresult.colData, sresult.sampleMap, self.metadata
        )

    def __getitem__(self, args: SlicerArgTypes) -> "MultiAssayExperiment":
        """Subset a `MultiAssayExperiment`.

        supports a tuple specifying slices along (rows, columns and experiments).

        Args:
            args (SlicerArgTypes): indices to slice. tuple can
                contains slices along dimensions (row, column, experiments)

        Raises:
            ValueError: Too many or too few slices.

        Returns:
            MultiAssayExperiment: new sliced `MultiAssayExperiment` object.
        """
        sresult = self._slice(args=args)
        return MultiAssayExperiment(
            sresult.experiments, sresult.colData, sresult.sampleMap, self.metadata
        )

    def __str__(self) -> str:
        pattern = (
            f"Class MultiAssayExperiment with {len(self.experiments.keys())} experiments and {len(self.colData)} samples \n"
            f"  experiments: "
        )

        for expname, expt in self.experiments.items():
            pattern = f"{pattern} \n" f"    {expname}: {str(expt)}"
        return pattern

    def completeCases(self) -> Sequence[bool]:
        """Identify samples that have data across all experiments.

        Returns:
            Sequence[bool]: a list, True if sample is present in all experiments.
        """
        vec = []
        for x in self.colData.index.tolist():
            subset = self.sampleMap[self.sampleMap["primary"] == x]

            vec.append(len(subset["assay"].unique()) == len(self.experiments.keys()))

        return vec

    def replicated(self) -> MutableMapping[str, MutableMapping[str, Sequence[bool]]]:
        """Identify samples with replicates within each experiment.

        Returns:
            MutableMapping[str, MutableMapping[str, Sequence[bool]]]: return true for replicates.
        """
        replicates = {}
        allSamples = self.colData.index.tolist()
        for expname, expt in self.experiments.items():
            if expname not in replicates:
                replicates[expname] = {}

                for s in allSamples:
                    replicates[expname][s] = []

            colnames = expt.colnames
            smap = self.sampleMap[self.sampleMap["assay"] == expname]

            for x in colnames:
                colmap = smap[smap["colname"] == x]
                for s in allSamples:
                    replicates[expname][s].append(s in colmap["primary"])

        return replicates

    def addExperiment(
        self,
        name: str,
        experiment: BaseSE,
        sampleMap: pd.DataFrame,
        colData: Optional[pd.DataFrame] = None,
    ):
        """Add an new experiment to MAE.
            Note: you have to provide information about new samples and a sample map.

        Args:
            name (str): Name of the experiment
            experiment (Union[SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment, ]): The experiment to add
            sampleMap (pd.DataFrame): sample map to append to the MAE
            colData (pd.DataFrame, optional): Sample data to append to the MAE. Defaults to None.
        """

        if name in self.experiments:
            raise ValueError(
                f"an experiment with {name} already exists, provide a different name"
            )

        self._validate_colData(colData)

        new_experiments = self._experiments.copy()
        new_experiments[name] = experiment

        new_colData = colData
        if new_colData is not None:
            new_colData = pd.concat([self.colData, colData], axis=0)

        new_sampleMap = pd.concat([self.sampleMap, sampleMap], axis=0)

        self._validate_experiments(new_experiments)
        self._validate_sampleMap(
            sampleMap=new_sampleMap, colData=new_colData, experiments=new_experiments
        )

        self._experiments = new_experiments
        self._sampleMap = new_sampleMap
        self._colData = new_colData

    def toMuData(self) -> MuData:
        """Transform `SingleCellExperiment` object to `MuData`.

        Returns:
            MuData: MuData representation
        """

        exptsList = OrderedDict()

        for expname, expt in self.experiments.items():
            if isinstance(expt, SingleCellExperiment):
                obj, adatas = expt.toAnnData(alts=True)

                exptsList[expname] = obj

                if adatas is not None:
                    for aname, aexpt in adatas.items():
                        exptsList[f"{expname}_{aname}"] = aexpt
            elif isinstance(expt, SummarizedExperiment):
                exptsList[expname] = expt.toAnnData()
            else:
                print(f"Experiment: {expname} is not supported!")

        return MuData(exptsList)
