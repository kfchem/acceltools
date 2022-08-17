import csv
from copy import deepcopy
from pathlib import Path
from typing import Iterable, Union

import matplotlib.pyplot as plt
import numpy as np
from accel.base.box import Box
from accel.base.systems import System, Systems
from accel.util.log import logger

from acceltools.base import ToolBox


def get_average(
    systems: Systems, averaged_systems: Systems = None, ecd_key="ecd", weighted_keys=["R_velocity", "R_length", "f"]
) -> Systems:
    if averaged_systems is None:
        averaged_systems = Box().bind(systems).get_average()
    for ave in averaged_systems:
        ecd_dict = {}
        confs = systems.has_state().has_label(ave.name)
        for idx in range(len(confs)):
            for state_number, state_dict in confs[idx].data[ecd_key].items():
                w_state = deepcopy(state_dict)
                for _key in weighted_keys:
                    if _key not in w_state:
                        logger.error(f"weighted key {_key} not found in data {ecd_key}")
                        continue
                    w_state[_key] *= confs[idx].distribution
                ecd_dict[f"{idx}_{confs[idx].name}_{state_number}"] = w_state
        ave.data["ecd"] = ecd_dict
    return averaged_systems


class EcdBox(ToolBox):
    def __init__(self, contents: Union[Box, Systems, Iterable[System], System]):
        self.expt: list[tuple[float]] = []
        self.expt_uv: list[tuple[float]] = []
        self.ecd_key: str = "ecd"
        self.curve_key: str = "ecd_curve"
        self.curve_key_uv: str = "uv_curve"
        self.nm_ev: float = 1239.84193
        self.const: float = 22.97
        self.const_uv: float = 28700.0
        self.calc_start: float = 100
        self.calc_stop: float = 800
        self.calc_step: float = 0.1
        super().__init__(contents)

    def load_expt(self, filepath: Union[Path, str], x_column: int = 1, y_column: int = 2, start_row: int = 2):
        with Path(filepath).open() as f:
            ls = [_l for _l in csv.reader(f)]
        self.expt = []
        for _l in ls[(start_row - 1) :]:
            self.expt.append((float(_l[x_column - 1]), float(_l[y_column - 1])))
        return self

    def get_average(self) -> Systems:
        return get_average(self.get(), ecd_key=self.ecd_key)

    def calc_curve(self, half_width: float = 0.19, shift: float = 0.0, scale: float = 1.0, key: str = "R_velocity"):
        x_values = np.arange(self.calc_start, self.calc_stop, self.calc_step)
        for _c in self.get():
            y_values = np.zeros(len(x_values))
            for state_num in _c.data[self.ecd_key]:
                _d: dict[str, Union[float, str]] = _c.data[self.ecd_key][state_num]
                y_values += (
                    _d["energy"]
                    * _d[key]
                    * np.exp(-1 * np.square(((self.nm_ev / x_values) - _d["energy"]) / half_width))
                )
            y_values /= self.const * half_width * np.sqrt(np.pi)
            y_values *= scale
            _c.data[self.curve_key] = [(x + shift, y) for x, y in zip(x_values, y_values)]
            logger.info(f"ECD of {_c.name} calculated: half-width {half_width}, shift {shift} and scale {scale}")
        return self

    def calc_curve_uv(self, half_width: float = 0.19, shift: float = 0.0, scale: float = 1.0, key: str = "f"):
        x_values = np.arange(self.calc_start, self.calc_stop, self.calc_step)
        for _c in self.get():
            y_values = np.zeros(len(x_values))
            for state_num in _c.data[self.ecd_key]:
                _d: dict[str, Union[float, str]] = _c.data[self.ecd_key][state_num]
                y_values += (
                    _d["energy"]
                    * _d[key]
                    * np.exp(-1 * np.square(((self.nm_ev / x_values) - _d["energy"]) / half_width))
                )
            y_values = y_values * self.const / (half_width * np.sqrt(np.pi))
            y_values *= scale
            _c.data[self.curve_key_uv] = [(x + shift, y) for x, y in zip(x_values, y_values)]
            logger.info(f"UV of {_c.name} calculated: half-width {half_width}, shift {shift} and scale {scale}")
        return self

    def write_img(
        self,
        directory: Path,
        start: float = 200,
        stop: float = 400,
        max_strength: float = None,
        with_expt: bool = True,
        with_ent: bool = False,
        with_bar: bool = False,
        transparent: bool = False,
    ):
        if max_strength is None:
            abs_max = 0.0
            for _c in self.get():
                x_vals = [abs(xy[1]) for xy in _c.data[self.curve_key] if start <= xy[0] and xy[0] <= stop]
                abs_max = max([abs_max] + x_vals)
            max_strength = abs_max * 1.1
        Path(directory).mkdir(exist_ok=True)
        for _c in self.get():
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
                    [xy[1] for xy in self.expt],
                    color="black",
                    linewidth=1.5,
                    linestyle="-.",
                    label=_c.name,
                )
            _p = Path(directory).joinpath(_c.name).with_suffix(".png")
            plt.savefig(_p, transparent=transparent, dpi=600)
            logger.info(f"ECD spectra of {_c.name} plotted")
            plt.close()
        return self

    def write_img_uv(
        self,
        directory: Path,
        start: float = 200,
        stop: float = 400,
        max_strength: float = None,
        with_expt: bool = True,
        with_bar: bool = False,
        transparent: bool = False,
    ):
        if max_strength is None:
            abs_max = 0.0
            for _c in self.get():
                x_vals = [abs(xy[1]) for xy in _c.data[self.curve_key_uv] if start <= xy[0] and xy[0] <= stop]
                abs_max = max([abs_max] + x_vals)
            max_strength = abs_max * 1.1
        Path(directory).mkdir(exist_ok=True)
        for _c in self.get():
            if max_strength == 0:
                max_strength = max(
                    [abs(xy[1]) for xy in _c.data[self.curve_key_uv] if start <= xy[0] and xy[0] <= stop]
                )
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlim(start, stop)
            ax.set_ylim(0, max_strength)
            ax.plot(
                [xy[0] for xy in _c.data[self.curve_key_uv]],
                [xy[1] for xy in _c.data[self.curve_key_uv]],
                color="black",
                linewidth=1.5,
                linestyle="-",
                label=_c.name,
            )
            if with_expt:
                ax.plot(
                    [xy[0] for xy in self.expt_uv],
                    [xy[1] for xy in self.expt_uv],
                    color="black",
                    linewidth=1.5,
                    linestyle="-.",
                    label=_c.name,
                )
            _p = Path(directory).joinpath(_c.name).with_suffix(".png")
            plt.savefig(_p, transparent=transparent, dpi=600)
            logger.info(f"UV spectra of {_c.name} plotted")
            plt.close()
        return self

    def write_csv(self, directory: Path):
        Path(directory).mkdir(exist_ok=True)
        for c in self.get():
            p = Path(directory).joinpath(c.name).with_suffix(".csv")
            with p.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["X (nm)", "Y (System-1cm-1)"])
                w.writerows(c.data[self.curve_key])
            logger.info(f"curve data was exported as {p.name}")
        return self

    def write_csv_uv(self, directory: Path):
        Path(directory).mkdir(exist_ok=True)
        for c in self.get():
            p = Path(directory).joinpath(c.name).with_suffix(".csv")
            with p.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["X (nm)", "Y (System-1cm-1)"])
                w.writerows(c.data[self.curve_key_uv])
            logger.info(f"curve data was exported as {p.name}")
        return self
