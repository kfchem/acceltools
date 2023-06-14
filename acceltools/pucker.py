from typing import Iterable, Union

import numpy as np
from accel.base.atoms import Atom
from accel.base.box import Box
from accel.base.systems import System, Systems
from accel.util.log import logger

from acceltools.base import ToolBox


class CPParams:
    def __init__(self, Q: float, theta: float, phi: float) -> None:
        self.Q = float(Q)
        self.theta = float(theta)
        self.phi = float(phi)


class PuckerParams:
    def __init__(self, cp: CPParams) -> None:
        self.cp: CPParams = cp


class CanonicalPuckers:
    # parameters from Tetrahderon 2001 57 477 and modified to match numbering.
    six: dict[str, PuckerParams] = {
        "4C1": PuckerParams(cp=CPParams(0.57, 180, 0)),
        "1C4": PuckerParams(cp=CPParams(0.57, 0, 0)),
        "25B": PuckerParams(cp=CPParams(0.76, 90, 240)),
        "B25": PuckerParams(cp=CPParams(0.76, 90, 60)),
        "36B": PuckerParams(cp=CPParams(0.76, 90, 120)),
        "B36": PuckerParams(cp=CPParams(0.76, 90, 300)),
        "14B": PuckerParams(cp=CPParams(0.76, 90, 0)),
        "B14": PuckerParams(cp=CPParams(0.76, 90, 180)),
        "2H3": PuckerParams(cp=CPParams(0.42, 129, 270)),
        "3H2": PuckerParams(cp=CPParams(0.42, 51, 90)),
        "3H4": PuckerParams(cp=CPParams(0.42, 51, 150)),
        "4H3": PuckerParams(cp=CPParams(0.42, 129, 330)),
        "4H5": PuckerParams(cp=CPParams(0.42, 129, 30)),
        "5H4": PuckerParams(cp=CPParams(0.42, 51, 210)),
        "5H6": PuckerParams(cp=CPParams(0.42, 51, 270)),
        "6H5": PuckerParams(cp=CPParams(0.42, 129, 90)),
        "6H1": PuckerParams(cp=CPParams(0.42, 129, 150)),
        "1H6": PuckerParams(cp=CPParams(0.42, 51, 330)),
        "1H2": PuckerParams(cp=CPParams(0.42, 51, 30)),
        "2H1": PuckerParams(cp=CPParams(0.42, 129, 210)),
        "5S1": PuckerParams(cp=CPParams(0.62, 88, 210)),
        "1S5": PuckerParams(cp=CPParams(0.62, 92, 30)),
        "6S2": PuckerParams(cp=CPParams(0.62, 92, 90)),
        "2S6": PuckerParams(cp=CPParams(0.62, 88, 270)),
        "1S3": PuckerParams(cp=CPParams(0.62, 88, 330)),
        "3S1": PuckerParams(cp=CPParams(0.62, 92, 150)),
        "2E": PuckerParams(cp=CPParams(0.45, 125, 240)),
        "E2": PuckerParams(cp=CPParams(0.45, 55, 60)),
        "3E": PuckerParams(cp=CPParams(0.45, 55, 120)),
        "E3": PuckerParams(cp=CPParams(0.45, 125, 300)),
        "4E": PuckerParams(cp=CPParams(0.45, 125, 360)),
        "E4": PuckerParams(cp=CPParams(0.45, 55, 180)),
        "5E": PuckerParams(cp=CPParams(0.45, 55, 240)),
        "E5": PuckerParams(cp=CPParams(0.45, 125, 60)),
        "6E": PuckerParams(cp=CPParams(0.45, 125, 120)),
        "E6": PuckerParams(cp=CPParams(0.45, 55, 300)),
        "1E": PuckerParams(cp=CPParams(0.45, 55, 360)),
        "E1": PuckerParams(cp=CPParams(0.45, 125, 180)),
    }
    five: dict[str, PuckerParams] = {}


CanonicalPuckers.six["1C4"].cp.phi


class PuckerBox(ToolBox):
    def __init__(self, contents: Union[Box, Systems, Iterable[System], System]):
        self.keys: list[str] = None
        super().__init__(contents)

    def calc_cp_six(
        self,
        number_1: int,
        number_2: int,
        number_3: int,
        number_4: int,
        number_5: int,
        number_6: int,
        key: str = "pucker_six_cp",
    ):
        # D. Cremer, J. A. Pople, J. Am. Chem. Soc. 1975, 97, 1354.
        numbers = [number_1, number_2, number_3, number_4, number_5, number_6]
        numbers = [int(n) for n in numbers]
        logger.debug(f"calclating puckering for {numbers}")
        for c in self.get():
            n_atoms: list[Atom] = [c.atoms.get(n) for n in numbers]
            N = len(n_atoms)
            n_xyzs = np.array([a.xyz for a in n_atoms])
            n_xyzs -= np.sum(n_xyzs, axis=0) / N
            Rp = np.sum(np.array([n_xyzs[j] * np.sin(2 * np.pi * j / N) for j in range(N)]), axis=0)
            Rpp = np.sum(np.array([n_xyzs[j] * np.cos(2 * np.pi * j / N) for j in range(N)]), axis=0)
            norm = np.cross(Rp, Rpp) / np.linalg.norm(np.cross(Rp, Rpp))
            zjs = np.array([np.dot(xyz, norm) for xyz in n_xyzs])
            m = 2
            qm_cos_phim = ((2 / N) ** 0.5) * sum([zjs[j] * np.cos(2 * np.pi * m * j / N) for j in range(N)])
            qm_sin_phim = -1 * ((2 / N) ** 0.5) * sum([zjs[j] * np.sin(2 * np.pi * m * j / N) for j in range(N)])
            qm = (qm_cos_phim**2 + qm_sin_phim**2) ** 0.5
            phim = np.arctan2(qm_sin_phim, qm_cos_phim)
            if phim < 0:
                phim += 2 * np.pi
            qeven = N ** (-1 / 2) * sum([((-1) ** j) * zjs[j] for j in range(N)])
            Q = (qm**2 + qeven**2) ** 0.5
            theta = np.arctan2(qm, qeven)
            if theta < 0:
                theta += 2 * np.pi
            c.data[f"{key}_Q"] = Q
            c.data[f"{key}_theta"] = np.rad2deg(theta)
            c.data[f"{key}_phi"] = np.rad2deg(phim)
            c.data[f"{key}_q2"] = np.rad2deg(qm)
            c.data[f"{key}_q3"] = np.rad2deg(qeven)
            logger.info(f"{c.name}: (Q, theta, phi) = ({Q:.2f}, {np.rad2deg(theta):.2f}, {np.rad2deg(phim):.2f})")
        return self

    def assign_canonical_pucker(
        self,
        key: str = "pucker_canonical",
        key_for_theta: str = "pucker_six_cp_theta",
        key_for_phi: str = "pucker_six_cp_phi",
    ):
        for c in self.get():
            dists = {}
            for pucker in CanonicalPuckers.six:
                theta1 = np.deg2rad(CanonicalPuckers.six[pucker].cp.phi)
                phi1 = np.deg2rad(90 - CanonicalPuckers.six[pucker].cp.theta)
                theta2 = np.deg2rad(c.data[key_for_phi])
                phi2 = np.deg2rad(90 - c.data[key_for_theta])
                dphi = abs(phi1 - phi2)
                dlam = abs(theta1 - theta2)
                temp1 = (np.sin(dphi / 2)) ** 2
                temp1 += np.cos(phi1) * np.cos(phi2) * ((np.sin(dlam / 2)) ** 2)
                temp1 = temp1**0.5
                dists[pucker] = (np.arcsin(temp1)) ** 2
            min_dist = dists["1C4"]
            min_dist_pucker: str = ""
            for pucker in dists:
                if dists[pucker] <= min_dist:
                    min_dist = dists[pucker]
                    min_dist_pucker = pucker
            c.data[key] = min_dist_pucker
        return self
