from collections import namedtuple
from typing import Dict, List, Optional, Tuple, Union

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

StrOrListStr = Union[str, List[str]]
SlicerTypes = Union[
    Dict[str, Union[List[int], slice]],
    List[int],
    slice,
]

SlicerArgTypes = Tuple[
    Optional[SlicerTypes],
    Optional[SlicerTypes],
    Optional[List[str]],
]

SlicerResult = namedtuple("SlicerResult", ["experiments", "sample_map", "col_data"])
