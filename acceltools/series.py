import math
from pathlib import Path
from statistics import mean
from typing import Iterable, Sequence, Set, Tuple, Union

from accel.base.atoms import Atoms
from accel.base.box import Box
from accel.base.systems import System, Systems
from accel.base.tools import change_dir, float_to_str
from accel.util.log import logger

from acceltools.base import ToolBox


def _get_angle_set(bonds: Set[Tuple[int]]) -> Set[Tuple[int]]:
    ang_set = set()
    for _b in bonds:
        for _al in [(_b[0], _b[1]), (_b[1], _b[0])]:
            for _pair_bond in [_pair_b for _pair_b in bonds if _al[0] in _pair_b and _b != _pair_b]:
                if _al[0] == _pair_bond[0]:
                    another_atom = _pair_bond[1]
                else:
                    another_atom = _pair_bond[0]
                if another_atom < _al[1]:
                    ang_set.add((another_atom, _al[0], _al[1]))
                elif another_atom > _al[1]:
                    ang_set.add((_al[1], _al[0], another_atom))
                else:
                    pass
    return ang_set


def _get_dihedral_set(bonds: Set[Tuple[int]]) -> Set[Tuple[int]]:
    dihed_set = set()
    for ang in _get_angle_set(bonds):
        for _al in [(ang[0], ang[1], ang[2]), (ang[2], ang[1], ang[0])]:
            for _pair_bond in [_pair_b for _pair_b in bonds if _al[0] in _pair_b and _al[1] not in _pair_b]:
                if _al[0] == _pair_bond[0]:
                    another_atom = _pair_bond[1]
                else:
                    another_atom = _pair_bond[0]
                if _al[0] < _al[1]:
                    dihed_set.add((another_atom, _al[0], _al[1], _al[2]))
                elif _al[0] > _al[1]:
                    dihed_set.add((_al[2], _al[1], _al[0], another_atom))
                else:
                    pass
    return dihed_set


def _remove_hydrogen(bonds: Set[Tuple[int]], atoms: Atoms) -> Set[Tuple[int]]:
    non_h_set = set()
    for _b in bonds:
        if "H" not in [atoms.get(_b[0]).symbol, atoms.get(_b[1]).symbol]:
            non_h_set.add(_b)
    return non_h_set


class SeriesBox(ToolBox):
    def __init__(self, contents: Union[Box, Systems, Iterable[System], System]):
        self.keys: list[str] = None
        super().__init__(contents)

    def modify_length(
        self,
        number_a: int,
        number_b: int,
        begin: float,
        end: float,
        step: float = 0.1,
        fix_a: bool = False,
        fix_b: bool = False,
        numbers_along_with_a: Sequence[int] = [],
        numbers_along_with_b: Sequence[int] = [],
    ):
        ret_list = []
        target = begin
        count = 0
        while target < end:
            mc = Box().bind(self.get()).duplicate()
            mc.modify_length(number_a, number_b, target, fix_a, fix_b, numbers_along_with_a, numbers_along_with_b)
            for c in mc.get():
                c.name = f"{c.name}_{count:03d}"
                ret_list.append(c)
            target += step
            count += 1
        return Systems().bind(ret_list)

    def calc_length(self, all=False, key: str = "L_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.get().labels.values():
            box = Box(confs)
            _atoms = confs.get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for bond in _bonds:
                tkey = f"{key}{_atoms.get(bond[0])}-{_atoms.get(bond[1])}"
                box.calc_length(bond[0], bond[1], tkey)
                key_list.append(tkey)
        self.keys = key_list
        return self

    def calc_angle(self, all=False, key: str = "A_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.get().labels.values():
            box = Box(confs)
            _atoms = confs.get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for angles in _get_angle_set(_bonds):
                tkey = f"{key}{_atoms.get(angles[0])}-{_atoms.get(angles[1])}-{_atoms.get(angles[2])}"
                box.calc_angle(angles[0], angles[1], angles[2], tkey)
                key_list.append(tkey)
        self.keys = key_list
        return self

    def calc_dihedral(self, all=False, key: str = "D_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.get().labels.values():
            box = Box(confs)
            _atoms = confs.get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for dhs in _get_dihedral_set(_bonds):
                tkey = f"{key}{_atoms.get(dhs[0])}-{_atoms.get(dhs[1])}-{_atoms.get(dhs[2])}-{_atoms.get(dhs[3])}"
                box.calc_dihedral(dhs[0], dhs[1], dhs[2], dhs[3], tkey)
                key_list.append(tkey)
        self.keys = key_list
        return self

    def calc_dihedral_xy(self, all=False, key: str = "D_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.get().labels.values():
            _atoms = confs.get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for dhs in _get_dihedral_set(_bonds):
                xkey = f"X{key}{_atoms.get(dhs[0])}-{_atoms.get(dhs[1])}-{_atoms.get(dhs[2])}-{_atoms.get(dhs[3])}"
                ykey = f"Y{key}{_atoms.get(dhs[0])}-{_atoms.get(dhs[1])}-{_atoms.get(dhs[2])}-{_atoms.get(dhs[3])}"
                for c in confs:
                    d_degree = c.atoms.get_dihedral(dhs[0], dhs[1], dhs[2], dhs[3])
                    c.data[xkey] = math.cos(math.radians(d_degree))
                    c.data[ykey] = math.sin(math.radians(d_degree))
                key_list.append(xkey)
                key_list.append(ykey)
        self.keys = key_list
        return self

    def write_trjxyz(self, output_dir: Path = None, order_key: str = None, centering: bool = True):
        for label, cs in self.get().labels.items():
            csls = cs.to_list()
            if order_key is not None:
                csls = [c for c in csls if c.data.get(order_key) is not None]
                csls = sorted(csls, key=lambda c: c.data[order_key])
            ls = []
            for c in csls:
                ls.append(f"{len(c.atoms)}\n")
                ls.append(f"{c.name}\n")
                if centering:
                    centering_vec = [mean([a.xyz[i] for a in c.atoms]) for i in range(3)]
                    try:
                        prec = max(max(len(float_to_str(a.xyz[i]).split(".")[1]) for a in c.atoms) for i in range(3))
                    except IndexError:
                        logger.error(f"could not resolve the precision in converting to xyz file of {c.name}")
                        prec = 10
                    centering_vec = [round(val, prec) for val in centering_vec]
                for a in c.atoms:
                    xyz = [a.x, a.y, a.z]
                    if centering:
                        xyz = [float_to_str(round(val - centering_vec[i], prec)) for i, val in enumerate(xyz)]
                    xyz = [float_to_str(val) for val in xyz]
                    ls.append(f"{a.symbol:<2} {xyz[0]:>15} {xyz[1]:>15} {xyz[2]:>15}\n")
            with change_dir(cs.get().path, output_dir, label).with_suffix(".xyz").open("w") as f:
                f.writelines(ls)
            logger.debug(f"{label}.xyz was exported")
        return self
