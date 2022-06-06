from typing import Iterable, Union

from accel import Mol
from accel.base.box import Box
from accel.base.mols import Mols


class ToolBox:
    def __init__(self, box: Union[Box, Mols, Iterable[Mol], Mol]):
        if isinstance(box, Box):
            self._mols: Mols = box.mols
        elif isinstance(box, Mols):
            self._mols: Mols = box
        elif isinstance(box, Mol):
            self._mols: Mols = Mols().append(box)
        else:
            self._mols: Mols = Box(box).mols

    @property
    def mols(self):
        return self._mols.has_state(True)

    @property
    def allmols(self) -> Mols:
        return self._mols.has_state(None)
