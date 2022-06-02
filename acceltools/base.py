from typing import Union

from accel.base.box import Box
from accel.base.mols import Mols


class ToolBox:
    def __init__(self, value: Union[Box, Mols]):
        if isinstance(value, Box):
            self.mols: Mols = value.mols
        elif isinstance(value, Mols):
            self.mols: Mols = value
        else:
            self.mols: Mols = Mols()
