import math
from typing import Sequence, Set, Tuple

from accel import Box
from accel.base.atoms import Atoms
from accel.base.xyz import get_dihedral


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


class Series:
    def __init__(self, mulcos: Box) -> None:
        self.mulcos = mulcos

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
            mc = Box().bind(self.mulcos.pack()).get_duplicate()
            mc.modify_length(number_a, number_b, target, fix_a, fix_b, numbers_along_with_a, numbers_along_with_b)
            for _c in mc.pack():
                _c.name = f"{_c.name}_{count:03d}"
            ret_list.extend(mc.pack().to_list())
            target += step
            count += 1
        return Box(ret_list)

    def calc_length(self, all=False, key: str = "L_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.mulcos.pack().labels.values():
            label_mc = Box(confs)
            _atoms = label_mc.pack().get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for bond in _bonds:
                tkey = f"{key}{_atoms.get(bond[0])}-{_atoms.get(bond[1])}"
                label_mc.calc_length(bond[0], bond[1], tkey)
                key_list.append(tkey)
        return key_list

    def calc_angle(self, all=False, key: str = "A_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.mulcos.pack().labels.values():
            label_mc = Box(confs)
            _atoms = label_mc.pack().get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for angles in _get_angle_set(_bonds):
                tkey = f"{key}{_atoms.get(angles[0])}-{_atoms.get(angles[1])}-{_atoms.get(angles[2])}"
                label_mc.calc_angle(angles[0], angles[1], angles[2], tkey)
                key_list.append(tkey)
        return key_list

    def calc_dihedral(self, all=False, key: str = "D_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.mulcos.pack().labels.values():
            label_mc = Box(confs)
            _atoms = label_mc.pack().get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for dhs in _get_dihedral_set(_bonds):
                tkey = f"{key}{_atoms.get(dhs[0])}-{_atoms.get(dhs[1])}-{_atoms.get(dhs[2])}-{_atoms.get(dhs[3])}"
                label_mc.calc_dihedral(dhs[0], dhs[1], dhs[2], dhs[3], tkey)
                key_list.append(tkey)
        return key_list

    def calc_dihedral_xy(self, all=False, key: str = "D_", in_label=True, ignore_hydrogen=True):
        key_list = []
        for confs in self.mulcos.pack().labels.values():
            label_mc = Box(confs)
            _atoms = label_mc.pack().get().atoms
            _bonds = _atoms.bonds.keys()
            if ignore_hydrogen:
                _bonds = _remove_hydrogen(_bonds, _atoms)
            for dhs in _get_dihedral_set(_bonds):
                xkey = f"X{key}{_atoms.get(dhs[0])}-{_atoms.get(dhs[1])}-{_atoms.get(dhs[2])}-{_atoms.get(dhs[3])}"
                ykey = f"Y{key}{_atoms.get(dhs[0])}-{_atoms.get(dhs[1])}-{_atoms.get(dhs[2])}-{_atoms.get(dhs[3])}"
                for _c in label_mc.pack():
                    d_degree = get_dihedral(
                        _c.atoms.get(dhs[0]), _c.atoms.get(dhs[1]), _c.atoms.get(dhs[2]), _c.atoms.get(dhs[3])
                    )
                    _c.data[xkey] = math.cos(math.radians(d_degree))
                    _c.data[ykey] = math.sin(math.radians(d_degree))
                key_list.append(xkey)
                key_list.append(ykey)
        return key_list
