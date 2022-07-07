import copy
from pathlib import Path
from typing import Iterable, Union

import matplotlib.pyplot as plt
from accel.base.box import Box
from accel.base.mols import Mol, Mols
from accel.util.log import logger
from matplotlib.axes import Axes

from acceltools.base import ToolBox


def _marker_generator(maker: str = None):
    if maker is not None:
        while True:
            yield maker
    else:
        markers = ["o", ",", "v", "^", "<", ">", "8", "p", "*", "h", "H", "D", "d"]
        i = 0
        while True:
            i += 1
            if i > len(markers):
                i = 1
            yield markers[i - 1]


class PlotBox(ToolBox):
    def __init__(self, box: Union[Box, Mols, Iterable[Mol], Mol]):
        self.path: Path = None
        super().__init__(box)

    def diagram(
        self,
        filepath: Path,
        diagram_roles: list[str] = ["reactant", "ts", "product"],
        zero_role_index: int = 0,
        role_key: str = "diagram_role",
        connection_key: str = "diagram_connection",
        non_mimimum: bool = True,
    ):

        _ze = min(_c.energy for _c in self.mols.has_data(role_key, diagram_roles[zero_role_index]))

        fig = plt.figure(figsize=((2 * len(diagram_roles)) + 0.4, 4.8))
        ax: Axes = fig.add_subplot(1, 1, 1)
        ax.set_xlim(0, (len(diagram_roles) + 1) * 2)
        ax.set_xticks([])
        x_values = {_key: [(i * 2) + 1, (i * 2) + 2] for i, _key in enumerate(diagram_roles)}
        ploted_conf: list[Mol] = []
        non_mimimum_conf: list[list[Mol]] = []
        for _role in diagram_roles:
            for _confs in Box(self.mols.has_data(role_key, _role)).mols.labels.values():
                _confs_orderd = sorted(_confs, key=lambda t: t.energy)
                ploted_conf.append(_confs_orderd[0])
                non_mimimum_conf.append(_confs_orderd[1:])

        for _c, _non_min_cs in zip(ploted_conf, non_mimimum_conf):
            _key = _c.data[role_key]
            _energy = _c.energy - _ze
            _connect: list[list[str]] = copy.deepcopy(_c.data.get(connection_key))
            _role_idx = diagram_roles.index(_key)
            if _connect is not None and len(_connect) == 2:
                for _idx in range(len(_connect)):
                    if isinstance(_connect[_idx], str):
                        _connect[_idx] = [_connect[_idx]]
                for _tc in ploted_conf:
                    if _role_idx == 0:
                        pass
                    elif _tc.data[role_key] == diagram_roles[_role_idx - 1] and _tc.label in _connect[0]:
                        ax.plot(
                            [x_values[_key][0] - 1, x_values[_key][0]],
                            [_tc.energy - _ze, _energy],
                            color="black",
                            linewidth=0.1,
                            linestyle="--",
                        )
                    if _role_idx == len(diagram_roles) - 1:
                        pass
                    elif _tc.data[role_key] == diagram_roles[_role_idx + 1] and _tc.label in _connect[1]:
                        ax.plot(
                            [x_values[_key][1], x_values[_key][1] + 1],
                            [_energy, _tc.energy - _ze],
                            color="black",
                            linewidth=0.1,
                            linestyle="--",
                        )
            ax.plot(
                x_values[_key],
                [_energy, _energy],
                color="black",
                linewidth=1.5,
                linestyle="-",
                label=_c.path.stem,
            )
            ax.text(x_values[_key][1] + 0.1, _energy, _c.name, fontsize=7, va="center")
            ax.text(
                x_values[_key][0] - 0.1,
                _energy,
                f"{_energy:.1f}",
                fontsize=7,
                va="center",
                ha="right",
            )
            if non_mimimum:
                for _non_min in _non_min_cs:
                    _energy = _non_min.energy - _ze
                    ax.plot(
                        x_values[_key],
                        [_energy, _energy],
                        color="black",
                        linewidth=0.1,
                        linestyle="-",
                    )

        _png = Path(filepath).with_suffix(".png")
        plt.savefig(_png, transparent=False, dpi=300)
        plt.close()
        logger.info(f"{str(_png)} was ploted")
        self.path = _png
        return self

    def scatter(
        self,
        filepath: Path,
        x_key: str = None,
        y_key: str = None,
        size=10,
        marker=None,
        color=True,
        size_by_energy=False,
        color_by_energy=False,
    ):
        fig = plt.figure()
        ax: Axes = fig.add_subplot(1, 1, 1)
        ax.tick_params(direction="in")
        maker = _marker_generator(marker)
        _legend = []
        for label, _confs in self.mols.labels.items():
            _x = [_c.energy for _c in _confs]
            ax.set_xlabel("Energy")
            _y = [_c.energy for _c in _confs]
            ax.set_ylabel("Energy")
            if x_key is not None:
                _x = [_c.data.get(x_key) for _c in _confs]
                ax.set_xlabel(x_key)
            if y_key is not None:
                _y = [_c.data.get(y_key) for _c in _confs]
                ax.set_ylabel(y_key)
            if size_by_energy:
                energies = [_c.energy for _c in _confs]
                energy_min = min(energies)
                energy_max = max(energies) - energy_min
                energies = [1 - ((_e - energy_min) / energy_max) for _e in energies]
                energies = [(_e * 2) ** 6 for _e in energies]
                size_used = [float(size) * _e for _e in energies]
            else:
                size_used = size
            if color_by_energy:
                energies = [_c.energy for _c in _confs]
                if size_by_energy:
                    ax.scatter(_x, _y, s=size_used, marker=next(maker), c=energies, alpha=0.5, edgecolors="gray")
                else:
                    ax.scatter(_x, _y, s=size_used, marker=next(maker), c=energies, edgecolors="gray")
            elif color:
                if size_by_energy:
                    ax.scatter(_x, _y, s=size_used, marker=next(maker), alpha=0.5, edgecolors="gray")
                else:
                    ax.scatter(_x, _y, s=size_used, marker=next(maker))
            else:
                ax.scatter(_x, _y, s=size_used, marker=next(maker), c="white", edgecolors="black", linewidths=0.5)
            _legend.append(label)
        if len(_legend) != 1:
            ax.legend(_legend)
        _png = Path(filepath).with_suffix(".png")
        ax.set_title(_png.stem)
        plt.savefig(_png, transparent=False, dpi=300)
        plt.close()
        logger.info(f"{str(_png)} was ploted")
        self.path = _png
        return self

    def line(self):
        return self
