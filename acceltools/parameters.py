from typing import Union


class Factor:
    def __init__(self, factors: dict = {}):
        self.slope: Union[float, None] = None
        self.intercept: Union[float, None] = None
        self.r_value: Union[float, None] = None
        self.p_value: Union[float, None] = None
        self.std_err: Union[float, None] = None
        self.mae: Union[float, None] = None
        self.rmse: Union[float, None] = None
        self.max_dev: Union[float, None] = None
        self.mean: Union[float, None] = None
        self.stdev: Union[float, None] = None
        self.degree: Union[float, None] = None

        for _key in factors:
            self.__setattr__(_key, float(factors[_key]))

    def __str__(self):
        return "Factor " + str({_k: _v for _k, _v in self.__dict__.items() if _v is not None})


DP4_PAR = dict()
DP4_PAR["C"] = Factor({"mean": 0.0, "stdev": 2.306, "degree": 11.38})
DP4_PAR["H"] = Factor({"mean": 0.0, "stdev": 0.185, "degree": 14.18})

# parameters for mPW1PW91/6-31G+(d,p)-PCM//B3LYP/6-31G(d)
DP4P_PAR = {_key: dict() for _key in ("S", "U", "U2")}
DP4P_PAR["S"]["C"] = Factor({"mean": 0.0, "stdev": 1.557, "degree": 6.227})
DP4P_PAR["S"]["H"] = Factor({"mean": 0.0, "stdev": 0.104, "degree": 3.893})
DP4P_PAR["U2"]["C"] = Factor({"mean": -0.920, "stdev": 1.748, "degree": 5.364})
DP4P_PAR["U2"]["H"] = Factor({"mean": 0.347, "stdev": 0.118, "degree": 4.911})
DP4P_PAR["U"]["C"] = Factor({"mean": 2.909, "stdev": 1.600, "degree": 6.269})
DP4P_PAR["U"]["H"] = Factor({"mean": -0.018, "stdev": 0.112, "degree": 3.651})

DICE_SCALE = {"DMSO": {}, "CHCl3": {}, "MeOH": {}, "MeCN": {}}
DICE_SCALE["DMSO"]["C"] = Factor({"slope": -0.9770, "intercept": 189.14})
DICE_SCALE["DMSO"]["H"] = Factor({"slope": -0.9726, "intercept": 31.31})
DICE_SCALE["DMSO"]["N"] = Factor({"slope": -0.9776, "intercept": -126.77})
DICE_SCALE["CHCl3"]["C"] = Factor({"slope": -0.9758, "intercept": 188.78})
DICE_SCALE["CHCl3"]["H"] = Factor({"slope": -1.0089, "intercept": 31.50})
DICE_SCALE["MeOH"]["C"] = Factor({"slope": -0.9766, "intercept": 190.59})
DICE_SCALE["MeOH"]["H"] = Factor({"slope": -0.9923, "intercept": 31.41})
DICE_SCALE["MeCN"]["C"] = Factor({"slope": -0.9705, "intercept": 189.66})
DICE_SCALE["MeCN"]["H"] = Factor({"slope": -0.9872, "intercept": 31.36})

DICE_PAR = dict()
DICE_PAR["C"] = Factor({"mean": 0.0, "stdev": 2.038, "degree": 36.46})
DICE_PAR["H"] = Factor({"mean": 0.0, "stdev": 0.113, "degree": 3.523})
DICE_PAR["N"] = Factor({"mean": 0.0, "stdev": 4.78})
