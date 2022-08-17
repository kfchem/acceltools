from typing import Iterable, Union

from accel.base.box import Box
from accel.base.systems import System, Systems


class ToolBox:
    def __init__(self, contents: Union[Box, Systems, Iterable[System], System]):
        if isinstance(contents, Box):
            self.contents: Systems = contents.get()
        elif isinstance(contents, Systems):
            self.contents: Systems = contents
        elif isinstance(contents, System):
            self.contents: Systems = Systems()
            self.contents.append(contents)
        else:
            self.contents: Systems = Box(contents).get()

    def get(self) -> Systems:
        return self.contents.has_state(True)
