import csv
from collections.abc import MutableSequence
from copy import deepcopy
from pathlib import Path
from typing import Iterator, Union

import numpy as np
from accel.base.box import Box
from accel.base.mols import Mol, Mols
from accel.util.constants import Elements
from accel.util.datadict import Data
from accel.util.log import logger
from scipy import stats

from acceltools.base import ToolBox


class Peak:
    def __init__(self) -> None:
        self.val: float = None
        self.name: str = None
        self.numbers: list[int] = []
        self.dtopic: "Peak" = None
        self.shape: str = None
        self.is_sp3: bool = True
        self.jvalues: list[float] = []
        self.swapped_flag: bool = False
        self.root: "Peak" = None
        self.corrs: list["Peak"] = []

    @property
    def nuclei(self):
        return self._nuclei

    @nuclei.setter
    def nuclei(self, value):
        self._nuclei = Elements.canonicalize(value)

    def show(self):
        return self

    def duplicate(self):
        return deepcopy(self)


class Peaks(MutableSequence):
    def __init__(self, label: str = "") -> None:
        self._peaks: list[Peak] = []
        self.label: str = label

    def __str__(self):
        return f"{self.label}: {len(self._peaks)} peaks"

    def __getitem__(self, index):
        return self._peaks[index]

    def __setitem__(self, index, value):
        self._peaks[index] = value

    def __delitem__(self, index):
        del self._peaks[index]

    def __len__(self):
        return len(self._peaks)

    def __iter__(self) -> Iterator[Peak]:
        return super().__iter__()

    def insert(self, index: int, value):
        self._peaks.insert(index, value)

    def bind(self, peaks: list[Peak]):
        self._peaks = peaks
        return self

    def show(self):
        return self

    def copy(self) -> "Peaks":
        _new = Peaks()
        _new.label = self.label
        return _new

    def duplicate(self) -> "Peaks":
        return self.copy().bind([_pk.duplicate() for _pk in self._peaks])

    @property
    def vals(self) -> list[float]:
        return [_pk.val for _pk in self._peaks]

    @property
    def nuclei(self) -> dict[str, "Peaks"]:
        peaks_dict: dict[str, "Peaks"] = {}
        for _p in self:
            peaks_dict[_p.nuclei] = None
        for key in [peaks_dict.keys()]:
            peaks_dict[key] = self.copy().bind([_pk for _pk in self._peaks if _pk.nuclei == key])
        return peaks_dict

    def abs(self):
        for _p in self._peaks:
            _p.val = abs(_p.val)
        return self

    def __abs__(self):
        return self.duplicate().abs()

    def add(self, value: float):
        value = float(value)
        for _p in self._peaks:
            _p.val += value
        return self

    def __add__(self, value: float):
        return self.duplicate().add(value)

    def __radd__(self, value: float):
        return self.duplicate().add(value)

    def mul(self, value: float):
        value = float(value)
        for _p in self._peaks:
            _p.val *= value
        return self

    def __mul__(self, value: float):
        return self.duplicate().mul(value)

    def __rmul__(self, value: float):
        return self.duplicate().mul(value)

    def sub(self, value: float):
        value = float(value)
        for _p in self._peaks:
            _p.val -= value
        return self

    def __sub__(self, value: float):
        return self.duplicate().sub(value)

    def __rsub__(self, value: float):
        return self.duplicate().mul(-1).add(value)

    def div(self, value: float):
        value = float(value)
        for _p in self._peaks:
            _p.val = _p.val / value
        return self

    def __div__(self, value: float):
        return self.duplicate().div(value)

    def pow(self, value: float):
        value = float(value)
        for _p in self._peaks:
            _p.val = _p.val**value
        return self

    def __pow__(self, value: float):
        return self.duplicate().pow(value)

    def has_nuclei(self, symbol) -> "Peaks":
        symbol = Elements.canonicalize(symbol)
        return self.copy().bind([_pk for _pk in self._peaks if _pk.nuclei == symbol])

    def has_sp3(self, is_sp3=True) -> "Peaks":
        return self.copy().bind([_pk for _pk in self._peaks if _pk.is_sp3 == is_sp3])


def check_identity(peaks_a: Peaks, peaks_b: Peaks) -> bool:
    return True


def get_expt(experiment_csv_file: Union[str, Path], ver: int = 2) -> Peaks:
    _path = Path(experiment_csv_file)
    with _path.open() as f:
        _ls = [_l for _l in csv.reader(f)]
    _peaks = Peaks(_path.stem)
    if ver == 1:
        for _l in _ls:
            _pk = Peak(_l[0])
            _pk.name = _l[1]
            _pk.numbers = [int(i) for i in _l[2].split()]
            _pk.is_sp3 = bool(_l[5])
            _pk.nuclei = _l[6]
            _peaks.append(_pk)
        for _i, _l in enumerate(_ls):
            if _l[3] != "":
                _peaks[_i].dtopic = _peaks[int(_l[3]) - 1]
            if _l[4] != "":
                _peaks[_i].root = _peaks[int(_l[4]) - 1]
    elif ver == 2:
        index_list = []
        for _id, _l in enumerate(_ls):
            index_list.append(_l[0])
            _pk = Peak(_l[1])
            _pk.name = _l[2]
            _pk.numbers = [int(i) for i in _l[3].split()]
            _pk.is_sp3 = bool(_l[6])
            _pk.nuclei = _l[7]
            _peaks.append(_pk)
        for _id, _l in enumerate(_ls):
            if _l[4] != "":
                _peaks[_id].dtopic = _peaks(index_list.index(int(_l[4])))
            if _l[5] != "":
                _peaks[_id].root = _peaks(index_list.index(int(_l[5])))
    return _peaks


def get_tensor(mol: Mol, key: str = "isotropic") -> Peaks:
    _peaks = Peaks(mol.name)
    for _a in mol.atoms:
        _p = Peak()
        _p.val = _a.data[key]
        _p.numbers = [_a.number]
        _peaks.append(_p)
    return _peaks


def get_shifts(tensor: Peaks, shielding_constant: dict[str, float], sloping: bool = False) -> Peaks:
    calc_peaks = tensor.duplicate()
    for peak in calc_peaks:
        ref = shielding_constant.get(peak.nuclei)
        if ref is None:
            continue
        if sloping:
            peak.val = (ref - peak.val) / (1 - (ref / (10 ^ 6)))
        else:
            peak.val = ref - peak.val
    return calc_peaks


def get_shifts_by_csv(tensor: Peaks, reference_csv_file: Union[str, Path], sloping: bool = False) -> Peaks:
    refs: dict[str, float] = {}
    with Path(reference_csv_file).open() as f:
        for _l in csv.reader(f):
            refs[Elements.canonicalize(_l[0])] = float(_l[1]) + float(_l[2])
    logger.info(f"NMR: {Path(reference_csv_file).absolute()} was loaded")
    logger.info(f"reference: {refs}")
    return get_shifts(tensor=tensor, references=refs, sloping=sloping)


def get_scaled(peaks: Peaks, slope: float, intercept: float):
    # scaled = (input - intercept) / slope
    return (peaks.duplicate() - intercept) / slope


def get_assigned(peaks: Peaks, expt: Peaks) -> Peaks:
    assigned_peaks = expt.duplicate()
    assigned_peaks.label = peaks.label
    for peak in assigned_peaks:
        peak.val = np.mean([calc_peak.val for calc_peak in peaks if calc_peak.numbers[0] in peak.numbers])
    return assigned_peaks


def get_factors(peaks: Peaks, expt: Peaks):
    # expt = (peaks - intercept) / slope
    check_identity(expt, peaks)
    _slope, _intercept, _rval = stats.linregress(expt.vals, peaks.vals)
    _slope: float = _slope
    _intercept: float = _intercept
    _rval: float = _rval
    return _slope, _intercept, _rval


def get_swapped(assigned: Peaks, expt: Peaks) -> Peaks:
    swapped = assigned.duplicate()
    for e_peak in [_p for _p in expt if _p.dtopic is not None]:
        if e_peak.root in expt:
            continue
        if e_peak is not e_peak.dtopic.dtopic:
            logger.error(f"{e_peak.name} and {e_peak.dtopic.name} are not related")
            return None
        if e_peak.val > e_peak.dtopic.val:
            is_large = True
        elif e_peak.val == e_peak.dtopic.val:
            continue
        else:
            is_large = False
        a_peak = swapped[expt.index(e_peak)]
        a_peak_dtopic = swapped[expt.index(e_peak.dtopic)]
        if a_peak.swapped_flag is True:
            continue
        if is_large ^ (a_peak.val > a_peak_dtopic.val):
            a_peak.val, a_peak_dtopic.val = a_peak_dtopic.val, a_peak.val
            a_peak.swapped_flag = True
            a_peak_dtopic.swapped_flag = True
            logger.info(f"peak  {a_peak.name} and {a_peak.dtopic.name} are swaped")
    for e_peak in [_p for _p in expt if _p.root is not None]:
        if e_peak is not e_peak.dtopic.dtopic:
            logger.error(f"{e_peak.name} and {e_peak.dtopic.name} are not related")
            return None
        if e_peak.root.swapped_flag is False:
            continue
        a_peak = swapped[expt.index(e_peak)]
        a_peak_dtopic = swapped[expt.index(e_peak.dtopic)]
        if a_peak.swapped_flag is True:
            continue
        a_peak.val, a_peak_dtopic.val = a_peak_dtopic, a_peak.val
        a_peak.swapped_flag = True
        a_peak_dtopic.swapped_flag = True
        logger.info(f"peak  {a_peak.name} and {a_peak.dtopic.name} are swaped")
    return swapped


def export_peak(peaks: Peaks, filepath: Union[str, Path] = None) -> Path:
    keys = ["Shift", "Name", "Nuclei", "Number", "sp2"]
    vals_dicts = []
    for _p in peaks:
        dicts = {
            "Shift": _p.val,
            "Name": _p.name,
            "Nuclei": _p.nuclei,
            "Number": _p.numbers,
            "sp3": _p.is_sp3,
        }
        vals_dicts.append(dicts)
    if filepath is None:
        filepath = Path.cwd().joinpath(peaks.label)
    _p = Path(filepath).with_suffix(".csv")
    with _p.open("w", newline="") as f:
        _w = csv.DictWriter(f, keys)
        _w.writeheader()
        for _val in vals_dicts:
            _w.writerow(_val)
    return filepath


def get_t_probability(assigned: Peaks, expt: Peaks, mean: float = 0.0, stdev: float = 0.0, degree: float = 0.0):
    check_identity(assigned, expt)
    t_probs = 1.0
    for idx in range(len(expt)):
        _err = assigned[idx].val - expt[idx].val
        t_probs *= stats.t.sf(abs(_err - mean) / stdev, degree)
    logger.info(f"Survival function of Students t distribution was calculated: {t_probs}")
    return float(t_probs)


def get_n_probability(assigned: Peaks, expt: Peaks, mean: float = 0.0, stdev: float = 0.0):
    check_identity(assigned, expt)
    n_probs = 1.0
    for idx in range(len(expt)):
        _err = assigned[idx].val - expt[idx].val
        n_probs *= stats.norm.sf(abs(_err - mean) / stdev)
    logger.info(f"Survival function of normal distribution was calculated: {n_probs}")
    return float(n_probs)


class NmrBox(ToolBox):
    def __init__(self, value: Union[Box, Mols]):
        self.expt: Peaks = None
        self.ref: dict[str, float] = None
        self.tensor: list[Peaks] = []
        self._shift: list[Peaks] = []
        self._assigned: list[Peaks] = []
        self._shift_for_analysis: list[Peaks] = []
        self._assigned_for_analysis: list[Peaks] = []
        self.data = Data(self)
        self.analyzing = False
        super().__init__(value)

    @property
    def shift(self) -> list[Peaks]:
        if self.analyzing:
            return self._shift_for_analysis
        return self._shift

    @shift.setter
    def shift(self, value: list[Peaks]):
        if self.analyzing:
            self._shift_for_analysis = value
        else:
            self._shift = value

    @property
    def assigned(self) -> list[Peaks]:
        if self.analyzing:
            return self._assigned_for_analysis
        return self._assigned

    @assigned.setter
    def assigned(self, value: list[Peaks]):
        if self.analyzing:
            self._assigned_for_analysis = value
        else:
            self._assigned = value

    def init_analysis(self):
        self._shift_for_analysis = []
        self._assigned_for_analysis = []
        self.analyzing = True
        return self

    def start_analysis(self):
        self.analyzing = True
        return self

    def stop_analysis(self):
        self.analyzing = False
        return self

    def load_expt(self, experiment_csv_path: Union[str, Path], ver: int = 2):
        self.expt = get_expt(experiment_csv_path=experiment_csv_path, ver=ver)
        return self

    def load_tensor(self, tensor_key="isotropic"):
        for _c in self.mols:
            self.tensor.append(get_tensor(_c, key=tensor_key))
        return self

    def load_ref(self, reference: Union[str, Path, dict]):
        if isinstance(reference, (str, Path)):
            refs: dict[str, float] = {}
            with Path(reference).open() as f:
                for _l in csv.reader(f):
                    refs[Elements.canonicalize(_l[0])] = float(_l[1]) + float(_l[2])
            self.ref = refs
        elif isinstance(reference, dict):
            self.ref = reference
        return self

    def conv_to_shift(self, precise: bool = True):
        if self.ref is None:
            for tsr in self.tensor:
                self.shift.append(tsr * (-1))
        else:
            for tsr in self.tensor:
                self.shift.append(get_shifts(tsr, self.ref, precise))
        return self

    def conv_to_assigned(self):
        if len(self.shift) == 0:
            self.conv_to_shift()
        self.assigned: list[Peaks] = []
        for peaks in self.shift:
            self.assigned.append(get_assigned(peaks, self.expt))
        return self

    def swap_assigned(self):
        if len(self.assigned) == 0:
            self.conv_to_assigned()
        swapped: list[Peaks] = []
        for peaks in self.assigned:
            swapped.append(get_swapped(peaks, self.expt))
        self.assigned = swapped
        return self

    def check_assign(self):
        return self

    def scale_assigned(self, key="factor"):
        self.check_assign()
        for idx in range(len(self.assigned)):
            for nuc in self.expt.nuclei.keys():
                a_peaks = self.assigned[idx].has_nuclei(nuc)
                slope, intercept, rval = get_factors(a_peaks, self.expt.has_nuclei(nuc))
                a_peaks.sub(intercept).div(slope)
                self.data[f"{key}_slope_{self.assigned[idx].label}_{nuc}"] = slope
                self.data[f"{key}_intercept_{self.assigned[idx].label}_{nuc}"] = intercept
                self.data[f"{key}_rvalue_{self.assigned[idx].label}_{nuc}"] = rval
        return self

    def calc_mae(self, key="MAE"):
        self.check_assign()
        edict_for_label: dict[str, dict[str, float]] = {}
        for a_peaks in self.assigned:
            errors_for_nuc: dict[str, float] = {}
            for e_peak in self.expt:
                errors_for_nuc[e_peak.nuclei] += abs(e_peak.val - a_peaks[self.expt.index(e_peak)].val)
            for nuc in errors_for_nuc:
                errors_for_nuc[nuc] = errors_for_nuc[nuc] / len(self.expt.nuclei.get(nuc))
            edict_for_label[a_peaks.label] = errors_for_nuc
        for nuc in self.expt.nuclei.keys():
            self.data[key + "_" + nuc] = {_k: _v[nuc] for _k, _v in edict_for_label.items()}
        return self

    def analyze_mae(self, key="MAE"):
        self.init_analysis()
        self.conv_to_shift()
        self.conv_to_assigned()
        self.swap_assigned()
        self.calc_mae(key=key)
        return self.stop_analysis()

    def calc_rmse(self, key="RMSE"):
        self.check_assign()
        edict_for_label: dict[str, dict[str, float]] = {}
        for a_peaks in self.assigned:
            errors_for_nuc: dict[str, float] = {}
            for e_peak in self.expt:
                errors_for_nuc[e_peak.nuclei] += (e_peak.val - a_peaks[self.expt.index(e_peak)].val) ** 2
            for nuc in errors_for_nuc:
                errors_for_nuc[nuc] = (errors_for_nuc[nuc] / len(self.expt.nuclei.get(nuc))) ** 0.5
            edict_for_label[a_peaks.label] = errors_for_nuc
        for nuc in self.expt.nuclei.keys():
            self.data[key + "_" + nuc] = {_k: _v[nuc] for _k, _v in edict_for_label.items()}
        return self

    def analyze_rmse(self, key="RMSE"):
        self.init_analysis()
        self.conv_to_shift()
        self.conv_to_assigned()
        self.swap_assigned()
        self.calc_rmse(key=key)
        return self.stop_analysis()

    def calc_maxerr(self, key="MaxError"):
        self.check_assign()
        edict_for_label: dict[str, dict[str, float]] = {}
        for a_peaks in self.assigned:
            errors_for_nuc: dict[str, float] = {nuc: 0.0 for nuc in self.expt.nuclei.keys()}
            for e_peak in self.expt:
                _err = e_peak.val - a_peaks[self.expt.index(e_peak)].val
                errors_for_nuc[e_peak.nuclei] = max(_err, errors_for_nuc[e_peak.nuclei])
            edict_for_label[a_peaks.label] = errors_for_nuc
        for nuc in self.expt.nuclei.keys():
            self.data[key + "_" + nuc] = {_k: _v[nuc] for _k, _v in edict_for_label.items()}
        return self

    def analyze_maxerr(self, key="MaxError"):
        self.init_analysis()
        self.conv_to_shift()
        self.conv_to_assigned()
        self.swap_assigned()
        self.calc_maxerr(key=key)
        return self.stop_analysis()

    def analyze_cmae(self, key="CMAE"):
        self.init_analysis()
        self.conv_to_shift().conv_to_assigned().swap_assigned()
        for peaks in self.assigned:
            peaks.mul(-1)
        self.scale_assigned()
        self.calc_mae(key=key)
        return self.stop_analysis()

    def analyze_dp4(self, key="DP4", param: dict[str, tuple[float]] = {"C": (0.0, 0.0, 0.0), "H": (0.0, 0.0, 0.0)}):
        self.init_analysis()
        self.conv_to_shift().conv_to_assigned().swap_assigned()
        for peaks in self.assigned:
            peaks.mul(-1)
        self.scale_assigned()
        t_probs = {"C": [1.0 for _ in self.assigned], "H": [1.0 for _ in self.assigned]}
        for nuc in ("C", "H"):
            logger.info(f"DP4 parameter for {nuc}: {param[nuc]}")
            t_probs[nuc] = [
                get_t_probability(
                    _a.has_nuclei(nuc), self.expt.has_nuclei(nuc), param[nuc][0], param[nuc][1], param[nuc][2]
                )
                for _a in self.assigned
            ]
            self.data[f"{key}_{nuc}"] = [100.0 * _val / sum(t_probs[nuc]) for _val in t_probs[nuc]]
        t_probs["All"] = [t_probs["C"][idx] * t_probs["H"][idx] for idx in range(len(self.assigned))]
        self.data[f"{key}_All"] = [100.0 * _val / sum(t_probs["All"]) for _val in t_probs["All"]]
        return self.stop_analysis()

    def analyze_dp4p(
        self,
        key="DP4plus",
        scaled_param: dict[str, tuple[float]] = {"C": (0.0, 0.0, 0.0), "H": (0.0, 0.0, 0.0)},
        unscaled_sp3_param: dict[str, tuple[float]] = {"C": (0.0, 0.0, 0.0), "H": (0.0, 0.0, 0.0)},
        unscaled_sp2_param: dict[str, tuple[float]] = {"C": (0.0, 0.0, 0.0), "H": (0.0, 0.0, 0.0)},
        reference: dict[str, float] = {"C": 0.0, "H": 0.0},
        int_degree: bool = False,
    ):
        self.init_analysis().conv_to_shift(reference=reference).assign_shift().swap_assign()

        probs_us_sp3: dict[str, list[float]] = {}
        probs_us_sp2: dict[str, list[float]] = {}
        probs_scaled: dict[str, list[float]] = {}
        for nuc in ("C", "H"):
            logger.info(f"DP4+ parameter for unscaled {nuc}: {unscaled_sp3_param[nuc]}")
            probs_us_sp3[nuc] = [
                get_t_probability(
                    a_peaks.has_nuclei(nuc).has_sp3(True),
                    self.expt.has_nuclei(nuc).has_sp3(True),
                    unscaled_sp3_param[nuc][0],
                    unscaled_sp3_param[nuc][1],
                    unscaled_sp3_param[nuc][2],
                )
                for a_peaks in self.assigned
            ]
            logger.info(f"DP4+ parameter for unscaled sp2 {nuc}: {unscaled_sp2_param[nuc]}")
            probs_us_sp2[nuc] = [
                get_t_probability(
                    a_peaks.has_nuclei(nuc).has_sp3(False),
                    self.expt.has_nuclei(nuc).has_sp3(False),
                    unscaled_sp2_param[nuc][0],
                    unscaled_sp2_param[nuc][1],
                    unscaled_sp2_param[nuc][2],
                )
                for a_peaks in self.assigned
            ]

        self.scale_assigned()
        for nuc in ("C", "H"):
            logger.info(f"DP4+ parameter for scaled {nuc}: {scaled_param[nuc]}")
            probs_scaled[nuc] = [
                get_t_probability(
                    a_peaks.has_nuclei(nuc),
                    self.expt.has_nuclei(nuc),
                    scaled_param[nuc][0],
                    scaled_param[nuc][1],
                    scaled_param[nuc][2],
                )
                for a_peaks in self.assigned
            ]

        probs_unscaled: dict[str, list[float]] = {}
        probs_all: dict[str, list[float]] = {}

        for nuc in ("C", "H"):
            probs_unscaled[nuc] = [_sp3 * _sp2 for _sp3, _sp2 in zip(probs_us_sp3[nuc], probs_us_sp2[nuc])]
            probs_all[nuc] = [_us * _s for _us, _s in zip(probs_unscaled[nuc], probs_scaled[nuc])]

        probs_us_sp3["All"] = [_c * _h for _c, _h in zip(probs_us_sp3["C"], probs_us_sp3["H"])]
        probs_us_sp2["All"] = [_c * _h for _c, _h in zip(probs_us_sp2["C"], probs_us_sp2["H"])]
        probs_scaled["All"] = [_c * _h for _c, _h in zip(probs_scaled["C"], probs_scaled["H"])]
        probs_unscaled["All"] = [_c * _h for _c, _h in zip(probs_unscaled["C"], probs_unscaled["H"])]
        probs_all["All"] = [_c * _h for _c, _h in zip(probs_all["C"], probs_all["H"])]

        pct_unscaled: dict[str, list[float]] = {}
        pct_scaled: dict[str, list[float]] = {}
        pct_all: dict[str, list[float]] = {}
        for nuc in ("C", "H"):
            pct_unscaled[nuc] = [100.0 * val / sum(probs_unscaled[nuc]) for val in probs_unscaled[nuc]]
            pct_scaled[nuc] = [100.0 * val / sum(probs_scaled[nuc]) for val in probs_scaled[nuc]]
            pct_all[nuc] = [100.0 * val / sum(probs_all[nuc]) for val in probs_all[nuc]]
        pct_unscaled["All"] = [100.0 * val / sum(probs_unscaled[nuc]) for val in probs_unscaled[nuc]]
        pct_scaled["All"] = [100.0 * val / sum(probs_scaled[nuc]) for val in probs_scaled[nuc]]
        pct_all["All"] = [100.0 * val / sum(probs_all[nuc]) for val in probs_all[nuc]]

        logger.info(f"DP4+: unscaled: {pct_unscaled}")
        logger.info(f"DP4+: scaled: {pct_scaled}")
        logger.info(f"DP4+: all: {pct_all}")
        for nuc in ("C", "H", "All"):
            self.data[f"{key}_unscaled_{nuc}"] = pct_unscaled[nuc]
            self.data[f"{key}_scaled_{nuc}"] = pct_scaled[nuc]
            self.data[f"{key}_all_{nuc}"] = pct_all[nuc]
        return self.stop_analysis()

    def analyze_dice(
        self,
        key: str = "dice",
        scale: dict[str, tuple[float]] = {"C": (0.0, 0.0), "H": (0.0, 0.0), "N": (0.0, 0.0)},
        param: dict[str, tuple[float]] = {"C": (0.0, 0.0, 0.0), "H": (0.0, 0.0, 0.0), "N": (0.0, 0.0, None)},
    ):
        self.init_analysis().conv_to_shift().conv_to_assigned().swap_assigned()
        for peaks in self.assigned:
            peaks.mul(-1)

        for peaks in self.assigned:
            peak_num = len(peaks)
            for nuc, sl in scale.items():
                peak_num -= len(peaks.has_nuclei(nuc).sub(sl[1]).div(sl[0]))
            if peak_num != 0:
                logger.error("part of peaks not scaled: check scale")

        t_probs: dict[str, list[float]] = {"N": {1.0 for _ in self.assigned}}
        for nuc, par in param.items():
            if nuc == "N":
                t_probs[nuc] = [
                    get_n_probability(a_peaks.has_nuclei(nuc), self.expt.has_nuclei(nuc), par[0], par[1])
                    for a_peaks in self.assigned
                ]
            else:
                t_probs[nuc] = [
                    get_t_probability(a_peaks.has_nuclei(nuc), self.expt.has_nuclei(nuc), par[0], par[1], par[2])
                    for a_peaks in self.assigned
                ]
            self.data[f"{key}_{nuc}"] = [100.0 * _val / sum(t_probs[nuc]) for _val in t_probs[nuc]]
        t_probs["All"] = [
            t_probs["C"][idx] * t_probs["H"][idx] * t_probs["N"][idx] for idx in range(len(self.assigned))
        ]
        self.data[f"{key}_All"] = [100.0 * _val / sum(t_probs["All"]) for _val in t_probs["All"]]
        return self.stop_analysis()
