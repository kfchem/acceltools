import copy
import csv
from collections import defaultdict
from collections.abc import MutableSequence
from copy import deepcopy
from pathlib import Path
from statistics import mean
from typing import Dict, Iterator, List, Union

import numpy as np
from accel.base.box import Box
from accel.base.tools import make_dir
from accel.util.constants import Elements
from accel.util.datadict import Data
from accel.util.log import logger
from scipy import stats

from acceltools.parameters import Factor


def _r_str(val):
    if isinstance(val, dict):
        return {_k: _r_str(_v) for _k, _v in val.items()}
    elif isinstance(val, list):
        return [_r_str(_v) for _v in val]
    else:
        return str(val)


class CalcValues:
    def __init__(self, init_value: float):
        self.raw: float = init_value
        self.val: float = init_value
        self.tmp: float = init_value
        self.saved: Dict[str, float] = {}
        self.swap_flag: bool = False
        self.tmp_swap_flag: float = init_value


class Shift:
    def __init__(self, expt_value: float):
        self.expt = float(expt_value)
        self._nuclei: str = ""
        self.name = ""
        self.atom_numbers = []
        self.is_sp2 = False
        self.iso_shift: "Shift" = None
        self.root_shift: "Shift" = None
        self.calcs: Dict[str, CalcValues] = {}
        self.data = Data(self)
        self.cache = {}

    @property
    def nuclei(self):
        return self._nuclei

    @nuclei.setter
    def nuclei(self, value):
        self._nuclei = Elements.canonicalize(value)

    def duplicate(self, parent_obj=None):
        return deepcopy(self)


class Peaks(MutableSequence):
    def __init__(self, parent_obj=None):
        self._peaks: List[Shift] = []
        self._parent: Nmr = parent_obj
        super().__init__()

    def _new_peak(self, value) -> Shift:
        if isinstance(value, Shift):
            return value
        else:
            raise TypeError

    def __str__(self):
        return f"object including {len(self._peaks)} peaks"

    def __getitem__(self, index):
        return self._peaks[index]

    def __setitem__(self, index, value):
        _n = self._new_peak(value)
        self._peaks[index] = _n

    def __delitem__(self, index):
        del self._peaks[index]

    def __len__(self):
        return len(self._peaks)

    def __iter__(self) -> Iterator[Shift]:
        return super().__iter__()

    def get_keys(self) -> List[str]:
        keys = set()
        for _pk in self._peaks:
            keys = keys.union(set(_pk.calcs.keys()))
        return sorted(list(keys))

    def add(self, other: Union[Dict[str, float], float]) -> "Peaks":
        if isinstance(other, float):
            other = {_key: other for _key in self.get_keys()}
        for _key, _val in other.items():
            for _pk in self._peaks:
                _pk.calcs[_key].val = _pk.calcs[_key].val + float(_val)
        return self

    def sub(self, other: Union[Dict[str, float], float]) -> "Peaks":
        if isinstance(other, float):
            other = {_key: other for _key in self.get_keys()}
        for _key, _val in other.items():
            for _pk in self._peaks:
                _pk.calcs[_key].val = _pk.calcs[_key].val - float(_val)
        return self

    def mul(self, other: Union[Dict[str, float], float]) -> "Peaks":
        if isinstance(other, float):
            other = {_key: other for _key in self.get_keys()}
        for _key, _val in other.items():
            for _pk in self._peaks:
                _pk.calcs[_key].val = _pk.calcs[_key].val * float(_val)
        return self

    def div(self, other: Union[Dict[str, float], float]) -> "Peaks":
        if isinstance(other, float):
            other = {_key: other for _key in self.get_keys()}
        for _key, _val in other.items():
            for _pk in self._peaks:
                _pk.calcs[_key].val = _pk.calcs[_key].val / float(_val)
        return self

    def pow(self, other: Union[Dict[str, float], float]) -> "Peaks":
        if isinstance(other, float):
            other = {_key: other for _key in self.get_keys()}
        for _key, _val in other.items():
            for _pk in self._peaks:
                _pk.calcs[_key].val = _pk.calcs[_key].val ** float(_val)
        return self

    def scaling(self, factors: Dict[str, Factor]) -> "Peaks":
        self.mul({_k: _f.slope for _k, _f in factors.items()})
        self.add({_k: _f.intercept for _k, _f in factors.items()})
        return self

    def insert(self, index: int, value) -> Shift:
        _n = self._new_peak(value)
        self._peaks.insert(index, _n)
        return _n

    def duplicate(self) -> "Peaks":
        _n = Peaks(self._parent)
        _n._peaks = [_pk.duplicate() for _pk in self._peaks]
        return _n

    def has_nuclei(self, nuclei_symbol) -> "Peaks":
        nuclei_symbol = Elements.canonicalize(nuclei_symbol)
        _n = Peaks(self._parent)
        _n._peaks = [_pk for _pk in self._peaks if _pk.nuclei == nuclei_symbol]
        return _n

    def has_sp2(self, is_sp2=True) -> "Peaks":
        _n = Peaks(self._parent)
        _n._peaks = [_pk for _pk in self._peaks if _pk.is_sp2 == is_sp2]
        return _n

    def get_expt(self) -> List[float]:
        return [_pk.expt for _pk in self._peaks]

    def get_calcs(self) -> Dict[str, List[float]]:
        calcs_dict = {_k: list() for _k in self.get_keys()}
        for _pk in self._peaks:
            for _k in calcs_dict:
                calcs_dict[_k].append(_pk.calcs[_k].val)
        return calcs_dict

    def get_errors(self) -> Dict[str, List[float]]:
        calcs_dict = {_k: list() for _k in self.get_keys()}
        for _pk in self._peaks:
            for _k in calcs_dict:
                calcs_dict[_k].append(_pk.calcs[_k].val - _pk.expt)
        return calcs_dict

    def get_nucleis(self) -> Dict[str, "Peaks"]:
        keys = set()
        for _pk in self._peaks:
            keys.add(_pk.nuclei)
        ret_dict = {_k: self.has_nuclei(_k) for _k in sorted(list(keys))}
        return ret_dict

    def get_factors(self) -> Dict[str, Factor]:
        factors = {_k: Factor() for _k in self.get_keys()}
        for _key, _fac in factors.items():
            (
                _fac.slope,
                _fac.intercept,
                _fac.r_value,
                _fac.p_value,
                _fac.std_err,
            ) = stats.linregress(self.get_calcs()[_key], self.get_expt())

        abs_errs = {_k: 0.0 for _k in self.get_keys()}
        sqr_errs = {_k: 0.0 for _k in self.get_keys()}
        for _pk in self._peaks:
            for _key, _val in _pk.calcs.items():
                abs_errs[_key] += abs(_val.val - _pk.expt)
                sqr_errs[_key] += (_val.val - _pk.expt) ** 2
        for _key, _fac in factors.items():
            _fac.mae = abs_errs[_key] / len(self._peaks)
            _fac.rmse = (sqr_errs[_key] / len(self._peaks)) ** 0.5
        return factors

    def get_t_probability(self, param: Factor, int_degree=False):
        if len(self._peaks) == 0:
            return None
        t_probs = {_k: 1.0 for _k in self.get_keys()}
        for _pk in self._peaks:
            for _key, _cal in _pk.calcs.items():
                _err = _cal.val - _pk.expt
                if int_degree:
                    t_probs[_key] *= stats.t.sf(abs(_err - param.mean) / param.stdev, int(param.degree))
                else:
                    t_probs[_key] *= stats.t.sf(abs(_err - param.mean) / param.stdev, param.degree)
        logger.info(f"Survival function of Students t distribution was calculated: {t_probs}")
        return t_probs

    def get_n_probability(self, param: Factor):
        if len(self._peaks) == 0:
            return None
        n_probs = {_k: 1.0 for _k in self.get_keys()}
        for _pk in self._peaks:
            for _key, _cal in _pk.calcs.items():
                _err = _cal.val - _pk.expt
                n_probs[_key] *= stats.norm.sf(abs(_err - param.mean) / param.stdev)
        logger.info(f"Survival function of Students t distribution was calculated: {n_probs}")
        return n_probs


class NmrBase:
    def __init__(self, box: Box):
        self.peaks: Peaks = Peaks(self)
        self.box = box.get_duplicate()
        self.data = Data(self)
        self.ref: Dict[str, float] = {}

    def __str__(self):
        return f"NMR ({len(self.peaks)} peaks)"

    def save(self, keyword: str = None):
        if keyword is None:
            for _pk in self.peaks:
                for _key in _pk.calcs:
                    _pk.calcs[_key].tmp = _pk.calcs[_key].val
                    _pk.calcs[_key].tmp_swap_flag = _pk.calcs[_key].swap_flag
        else:
            for _pk in self.peaks:
                for _key in _pk.calcs:
                    _pk.calcs[_key].saved[keyword] = _pk.calcs[_key].val
        return self

    def restore(self):
        for _pk in self.peaks:
            for _key in _pk.calcs:
                _pk.calcs[_key].val = _pk.calcs[_key].tmp
        return self

    def reset(self):
        for _pk in self.peaks:
            for _key in _pk.calcs:
                _pk.calcs[_key].val = _pk.calcs[_key].raw
                _pk.calcs[_key].swap_flag = False
        return self

    def export_data(self, filepath):
        registered_dict = {}

        def _unpack_dict(parent_list: list, obj: list):
            if isinstance(obj, dict):
                for key, val in obj.items():
                    _unpack_dict(parent_list + [key], val)
            else:
                registered_dict[tuple(parent_list)] = obj

        _unpack_dict([], self.data._data)
        sorted_dict = defaultdict(dict)
        label_set = set()
        for key, val in registered_dict.items():
            sorted_dict[key[:-1]][key[-1]] = val
            label_set.add(key[-1])

        _p = Path(filepath).with_suffix(".csv")
        with _p.open("w", newline="") as f:
            _w = csv.DictWriter(f, [""] + sorted(list(label_set)))
            _w.writeheader()
            for _key, _val in sorted_dict.items():
                _val[""] = " | ".join(_key)
                _w.writerow(_val)
        logger.info(f"NMR: {_p.name} was created")
        return self

    def export_peaks(self, filepath):
        keys = ["Shift", "Name", "Nuclei", "Number", "sp2"]

        def keys_append(val: str):
            if val not in keys:
                keys.append(val)

        vals_dicts = []
        for peak in self.peaks:
            dicts = {
                "Shift": peak.expt,
                "Name": peak.name,
                "Nuclei": peak.nuclei,
                "Number": peak.atom_numbers,
                "sp2": peak.is_sp2,
            }
            for _k, _v in peak.calcs.items():
                dicts[_k + " | raw"] = _v.raw
                keys_append(_k + " | raw")
                dicts[_k + " | val"] = _v.val
                keys_append(_k + " | val")
                for _saved_key, _saved_val in _v.saved.items():
                    dicts[_k + " | " + _saved_key] = _saved_val
                    keys_append(_k + " | " + _saved_key)
            vals_dicts.append(dicts)
        _p = Path(filepath).with_suffix(".csv")
        with _p.open("w", newline="") as f:
            _w = csv.DictWriter(f, keys)
            _w.writeheader()
            for _val in vals_dicts:
                _w.writerow(_val)
        logger.info(f"NMR: {_p.name} was created")
        return self

    def load_expt(self, csv_path, ver=2):
        with Path(csv_path).open() as f:
            _ls = [_l for _l in csv.reader(f)]
        if ver == 1:
            for _l in _ls:
                _sh = Shift(_l[0])
                _sh.name = _l[1]
                _sh.atom_numbers = [int(i) for i in _l[2].split()]
                _sh.is_sp2 = bool(_l[5])
                _sh.nuclei = _l[6]
                self.peaks.append(_sh)
            for _i, _l in enumerate(_ls):
                if _l[3] != "":
                    self.peaks[_i].iso_shift = self.peaks[int(_l[3]) - 1]
                if _l[4] != "":
                    self.peaks[_i].root_shift = self.peaks[int(_l[4]) - 1]
        if ver == 2:
            index_list = []
            for _id, _l in enumerate(_ls):
                index_list.append(_l[0])
                _sh = Shift(_l[1])
                _sh.name = _l[2]
                _sh.atom_numbers = [int(i) for i in _l[3].split()]
                _sh.is_sp2 = bool(_l[6])
                _sh.nuclei = _l[7]
                self.peaks.append(_sh)
            for _id, _l in enumerate(_ls):
                if _l[4] != "":
                    _sh.iso_shift = self.peaks(index_list.index(int(_l[4])))
                if _l[5] != "":
                    _sh.root_shift = self.peaks(index_list.index(int(_l[5])))
        return self

    def load_ref(self, csv_path):
        with Path(csv_path).open() as f:
            for _l in csv.reader(f):
                self.ref[Elements.canonicalize(_l[0])] = float(_l[1]) + float(_l[2])
            logger.info(f"NMR: {Path(csv_path).absolute()} was loaded")
            logger.info(f"ref_data: {self.ref}")
        return self

    def apply_ref(self):
        for _pk in self.peaks:
            for _key in _pk.calcs:
                _pk.calcs[_key].val = self.ref[_pk.nuclei] - _pk.calcs[_key].val
        return self

    def assign(self, label=None):
        avg_mulcos = Box().bind(self.box.pack_average(keys_for_atoms=["isotropic"]))
        if label is None:
            _confs = avg_mulcos.pack()
        else:
            _confs = avg_mulcos.pack().labels[label]
        for _c in _confs:
            for _pk in self.peaks:
                if _pk.calcs.get(_c.name) is not None:
                    logger.error(f"confomer name {_c.name} is already registered")
                else:
                    _pk.calcs[_c.name] = CalcValues(
                        np.mean([_c.atoms.get(_num).data["isotropic"] for _num in _pk.atom_numbers])
                    )
        return self

    def swap_isos(self):
        for _pk in self.peaks:
            if not (_pk.iso_shift is not None and _pk.root_shift is None):
                continue
            if _pk != _pk.iso_shift.iso_shift:
                raise Exception
            if _pk.expt > _pk.iso_shift.expt:
                is_large = True
            elif (_pk.expt == _pk.iso_shift.expt) and (self.peaks.index(_pk) > self.peaks.index(_pk.iso_shift)):
                is_large = True
            else:
                is_large = False
            for _k, _val in _pk.calcs.items():
                if _val.swap_flag:
                    continue
                if is_large ^ (_val.val > _pk.iso_shift.calcs[_k].val):
                    _val.val, _pk.iso_shift.calcs[_k].val = (
                        _pk.iso_shift.calcs[_k].val,
                        _val.val,
                    )
                    _val.swap_flag = True
                    _pk.iso_shift.calcs[_k].swap_flag = True
                    logger.info(f"NMR: {_k}: calcd peak for {_pk.expt} and {_pk.iso_shift.expt} are swaped")

        for _pk in self.peaks:
            if not (_pk.iso_shift is not None and _pk.root_shift is not None):
                continue
            if _pk != _pk.iso_shift.iso_shift:
                raise Exception
            for _k, _val in _pk.calcs.items():
                if _val.swap_flag:
                    continue
                if _pk.root_shift.calcs[_k].swap_flag:
                    _val.val, _pk.iso_shift.calcs[_k].val = (
                        _pk.iso_shift.calcs[_k].val,
                        _val.val,
                    )
                    _val.swap_flag = True
                    _pk.iso_shift.calcs[_k].swap_flag = True
                    logger.info(
                        "NMR: {}: calcd peak for {} and {} are swaped according to root shift".format(
                            _k, _pk.expt, _pk.iso_shift.expt
                        )
                    )
        return self

    def swap_isos_invert(self):
        self.peaks.mul(-1.0)
        self.swap_isos()
        self.peaks.mul(-1.0)
        return self


class Nmr(NmrBase):
    def analyze_mae(self):
        self.save().reset()
        self.apply_ref().swap_isos()
        _data = dict()
        for _nuc, _pks in self.peaks.get_nucleis().items():
            _data[_nuc] = {_key: _fac.mae for _key, _fac in _pks.get_factors().items()}
            logger.info(f"MAE for {_nuc}: {_data[_nuc]}")
        self.data["MAE"] = _data
        return self.save("MAE").restore()

    def analyze_cmae(self):
        self.save().reset().swap_isos_invert()
        _data = dict()
        _data_factor = dict()
        for _nuc, _pks in self.peaks.get_nucleis().items():
            _data_factor[_nuc] = _pks.get_factors()
            _pks.scaling(_data_factor[_nuc])
            _data[_nuc] = {_key: _fac.mae for _key, _fac in _pks.get_factors().items()}
            logger.info(f"factors for {_nuc}: {_r_str(_data_factor[_nuc])}")
            logger.info(f"CMAE for {_nuc}: {_data[_nuc]}")
        self.data["CMAE"] = _data
        self.data["CMAE_FACTORS"] = _data_factor
        return self.save("CMAE").restore()

    def cal_averaged_distance(self, output_dir: Path, noe_path_dict: dict):
        for _mol_name, _p in noe_path_dict.items():
            with Path(_p).open() as f:
                _noes = {int(_l[0]): [float.fromhex(_v) for _v in _l[1:]] for _l in csv.reader(f)}
            with Path(_p).open() as f:
                _noes_index = [int(_l[0]) for _l in csv.reader(f)]
            shift_noes_list = []
            for _s in self.shifts:
                _s_noes = [_noes[cal_num] for cal_num in _s.calc_number if cal_num in _noes.keys()]
                if len(_s_noes) != 0:
                    _s_noes_ave = {_noes_index[i]: mean([_l[i] for _l in _s_noes]) for i in range(len(_s_noes[0]))}
                    shift_noes_list.append({"shift": _s, "noe": copy.deepcopy(_s_noes_ave)})
            already_swaped = []
            for _d in shift_noes_list:
                if _d["shift"].row_index in already_swaped:
                    continue
                if _mol_name in _d["shift"].swaped_data:
                    for _d_pair in shift_noes_list:
                        if _d["shift"].either_index == _d_pair["shift"].row_index:
                            _d["noe"], _d_pair["noe"] = _d_pair["noe"], _d["noe"]
                            logger.info(
                                "NOE data of {} and {} was swaped in {} 1st step".format(
                                    _d["shift"].label, _d_pair["shift"].label, _mol_name
                                )
                            )
                            already_swaped.append(_d["shift"].row_index)
                            already_swaped.append(_d_pair["shift"].row_index)
            _noes_T = {_id: [_d["noe"][_id] for _d in shift_noes_list] for _id in _noes_index}
            shift_noes_list_T = []
            for _s in self.shifts:
                _s_noes = [_noes_T[cal_num] for cal_num in _s.calc_number if cal_num in _noes_T.keys()]
                if len(_s_noes) != 0:
                    _s_noes_ave = [mean([_l[i] for _l in _s_noes]) for i in range(len(_s_noes[0]))]
                    shift_noes_list_T.append({"shift": _s, "noe": copy.deepcopy(_s_noes_ave)})
            already_swaped = []
            for _d in shift_noes_list_T:
                if _d["shift"].row_index in already_swaped:
                    continue
                if _mol_name in _d["shift"].swaped_data:
                    for _d_pair in shift_noes_list_T:
                        if _d["shift"].either_index == _d_pair["shift"].row_index:
                            _d["noe"], _d_pair["noe"] = _d_pair["noe"], _d["noe"]
                            logger.info(
                                "NOE data of {} and {} was swaped in {} 2nd step".format(
                                    _d["shift"].label, _d_pair["shift"].label, _mol_name
                                )
                            )
                            already_swaped.append(_d["shift"].row_index)
                            already_swaped.append(_d_pair["shift"].row_index)
            if output_dir is not None:
                make_dir(output_dir)
                _o = [["expt_shift", "swaped_calc_number", "label"] + [_d["shift"].label for _d in shift_noes_list]]
                for _d in shift_noes_list_T:
                    _dist_list = [_noe ** (-1 / 6) for _noe in _d["noe"]]
                    if _d["shift"].row_index in already_swaped:
                        _swaped_cal_num = self.shifts[_d["shift"].either_index - 1].calc_number
                    else:
                        _swaped_cal_num = _d["shift"].calc_number
                    _o.append([_d["shift"].expt_shift, _swaped_cal_num, _d["shift"].label] + _dist_list)

                _p = Path(output_dir).joinpath(_mol_name + "_NOE_distance").with_suffix(".csv")
                with _p.open("w", newline="") as f:
                    csv.writer(f).writerows(_o)
        return self

    def cal_averaged_jvalues(self, output_dir: Path, jvalue_path_dict: dict, scaling=1.0):
        logger.info(f"scaling factor: {scaling}")
        for _mol_name, _p in jvalue_path_dict.items():
            with Path(_p).open() as f:
                _jvalues = {int(_l[0]): [float(_v) for _v in _l[1:]] for _l in csv.reader(f)}
            with Path(_p).open() as f:
                _js_index = [int(_l[0]) for _l in csv.reader(f)]
            shift_js_list = []
            for _s in self.shifts:
                _s_js = [_jvalues[cal_num] for cal_num in _s.calc_number if cal_num in _jvalues.keys()]
                if len(_s_js) != 0:
                    _s_js_ave = {_js_index[i]: mean([_l[i] for _l in _s_js]) * scaling for i in range(len(_s_js[0]))}
                    shift_js_list.append({"shift": _s, "jvalue": copy.deepcopy(_s_js_ave)})
            already_swaped = []
            for _d in shift_js_list:
                if _d["shift"].row_index in already_swaped:
                    continue
                if _mol_name in _d["shift"].swaped_data:
                    for _d_pair in shift_js_list:
                        if _d["shift"].either_index == _d_pair["shift"].row_index:
                            _d["jvalue"], _d_pair["jvalue"] = (
                                _d_pair["jvalue"],
                                _d["jvalue"],
                            )
                            logger.debug(
                                "J vaule data of {} and {} was swaped in {} 1st step".format(
                                    _d["shift"].label, _d_pair["shift"].label, _mol_name
                                )
                            )
                            already_swaped.append(_d["shift"].row_index)
                            already_swaped.append(_d_pair["shift"].row_index)
            _js_T = {_id: [_d["jvalue"][_id] for _d in shift_js_list] for _id in _js_index}
            shift_js_list_T = []
            for _s in self.shifts:
                _s_js = [_js_T[cal_num] for cal_num in _s.calc_number if cal_num in _js_T.keys()]
                if len(_s_js) != 0:
                    _s_js_ave = [mean([_l[i] for _l in _s_js]) for i in range(len(_s_js[0]))]
                    shift_js_list_T.append({"shift": _s, "jvalue": copy.deepcopy(_s_js_ave)})
            already_swaped = []
            for _d in shift_js_list_T:
                if _d["shift"].row_index in already_swaped:
                    continue
                if _mol_name in _d["shift"].swaped_data:
                    for _d_pair in shift_js_list_T:
                        if _d["shift"].either_index == _d_pair["shift"].row_index:
                            _d["jvalue"], _d_pair["jvalue"] = (
                                _d_pair["jvalue"],
                                _d["jvalue"],
                            )
                            logger.debug(
                                "J value data of {} and {} was swaped in {} 2nd step".format(
                                    _d["shift"].label, _d_pair["shift"].label, _mol_name
                                )
                            )
                            already_swaped.append(_d["shift"].row_index)
                            already_swaped.append(_d_pair["shift"].row_index)
            if output_dir is not None:
                make_dir(output_dir)
                _o = [["expt_shift", "swaped_calc_number", "label"] + [_d["shift"].label for _d in shift_js_list]]
                for _d in shift_js_list_T:
                    if _d["shift"].row_index in already_swaped:
                        _swaped_cal_num = self.shifts[_d["shift"].either_index - 1].calc_number
                    else:
                        _swaped_cal_num = _d["shift"].calc_number
                    _o.append([_d["shift"].expt_shift, _swaped_cal_num, _d["shift"].label] + _d["jvalue"])

                _p = Path(output_dir).joinpath(_mol_name + "_J_values").with_suffix(".csv")
                with _p.open("w", newline="") as f:
                    csv.writer(f).writerows(_o)
        return self
