import csv
from copy import deepcopy
from pathlib import Path
from typing import Iterable, Union

import matplotlib.pyplot as plt
import numpy as np
from accel.base.box import Box
from accel.base.mols import Mol, Mols
from accel.base.tools import change_dir, make_dir
from accel.util.log import logger

from acceltools.base import ToolBox


def get_average(mols: Mols, averaged_mols: Mols = None, ecd_key="ecd") -> Mols:
    if averaged_mols is None:
        averaged_mols = Box().bind(mols).get_average()
    for ave in averaged_mols:
        confs = mols.has_state().has_label(ave.name)
        for idx in range(len(confs)):
            for state_number, state_dict in confs[idx].data["ecd"].items():
                w_state = deepcopy(state_dict)
                w_state["R_velocity"] *= confs[idx].data["distribution"]
                w_state["R_length"] *= confs[idx].data["distribution"]
                ave.data["ecd"][f"{idx}_{confs[idx].name}_{state_number}"] = w_state
    return averaged_mols


class EcdBox(ToolBox):
    def __init__(self, box: Union[Box, Mols, Iterable[Mol], Mol]):
        self.expt: list[tuple[float]] = []
        self.ecd_key: str = "ecd"
        self.curve_key: str = "ecd_curve"
        self.nm_ev: float = 1239.84193
        self.const: float = 22.96
        self.calc_start: float = 100
        self.calc_stop: float = 800
        self.calc_step: float = 0.1
        super().__init__(box)

    def load_expt(self, filepath: Union[Path, str], x_column: int = 1, y_column: int = 2, start_row: int = 2):
        with Path(filepath).open() as f:
            ls = [_l for _l in csv.reader(f)]
        self.expt = []
        for _l in ls[(start_row - 1) :]:
            self.expt.append((float(_l[x_column - 1]), float(_l[y_column - 1])))
        return self

    def get_average(self) -> Mols:
        return get_average(self.mols, ecd_key=self.ecd_key)

    def calc_curve(self, half_width: float = 0.19, shift: float = 0.0, scale: float = 1.0, key: str = "R_velocity"):
        x_values = np.arange(self.calc_start, self.calc_stop, self.calc_step)
        for _c in self.mols:
            y_values = np.zeros(len(x_values))
            for _d in _c.data[self.ecd_key]:
                y_values += (
                    _d["energy"]
                    * _d[key]
                    * np.exp(-1 * np.square(((self.nm_ev / x_values) - _d["energy"]) / half_width))
                )
            y_values /= self.const * half_width * np.sqrt(np.pi)
            y_values *= scale
            _c.data[self.curve_key] = [(x + shift, y) for x, y in zip(x_values, y_values)]
            logger.info(f"ECD of {_c.name} was calculated using half-width {half_width}")
        return self

    def plot(
        self,
        directory: Path = None,
        start: float = 200,
        stop: float = 400,
        max_strength: float = None,
        with_expt: bool = True,
        with_ent: bool = False,
        transparent: bool = False,
    ):
        if max_strength is None:
            abs_max = 0.0
            for _c in self.mols:
                x_vals = [abs(xy[1]) for xy in _c.data[self.curve_key] if start <= xy[0] and xy[0] <= stop]
                abs_max = max([abs_max] + x_vals)
            max_strength = abs_max

        make_dir(directory)
        for _c in self.mols:
            if max_strength == 0:
                max_strength = max([abs(xy[1]) for xy in _c.data[self.curve_key] if start <= xy[0] and xy[0] <= stop])
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlim(start, stop)
            ax.set_ylim(-max_strength, max_strength)
            ax.plot(
                [xy[0] for xy in _c.data[self.curve_key]],
                [0.0 for _ in _c.data[self.curve_key]],
                color="black",
                linewidth=0.5,
                linestyle="-",
            )
            ax.plot(
                [xy[0] for xy in _c.data[self.curve_key]],
                [xy[1] for xy in _c.data[self.curve_key]],
                color="black",
                linewidth=1.5,
                linestyle="-",
                label=_c.name,
            )
            if with_ent:
                ax.plot(
                    [xy[0] for xy in _c.data[self.curve_key]],
                    [-1 * xy[1] for xy in _c.data[self.curve_key]],
                    color="black",
                    linewidth=1.5,
                    linestyle="--",
                    label=_c.name,
                )
            if with_expt:
                ax.plot(
                    [xy[0] for xy in self.expt],
                    [-1 * xy[1] for xy in self.expt],
                    color="black",
                    linewidth=1.5,
                    linestyle="-.",
                    label=_c.name,
                )
            _p = change_dir(_c.path, directory).with_suffix(".png")
            plt.savefig(_p, transparent=transparent, dpi=600)
            logger.info(f"ECD spectra of {_c.name} plotted")
            plt.close()
        return self
