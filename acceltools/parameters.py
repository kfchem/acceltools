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
