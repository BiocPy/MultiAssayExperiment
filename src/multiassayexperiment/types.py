from collections import namedtuple
from typing import MutableMapping, Optional, Sequence, Tuple, Union

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

StrOrListStr = Union[str, Sequence[str]]
SlicerTypes = Union[
    MutableMapping[str, Union[Sequence[int], slice]],
    Union[Sequence[int], slice],
]

SlicerArgTypes = Tuple[
    Optional[SlicerTypes],
    Optional[SlicerTypes],
    Optional[Sequence[str]],
]

SlicerResult = namedtuple("SlicerResult", ["experiments", "sampleMap", "colData"])
