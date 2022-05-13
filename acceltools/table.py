from typing import List

import pandas as pd
from accel import Box
from accel.util.log import logger


class Table:
    def __init__(self, box: Box) -> None:
        self.box: Box = box

    def get_df(self, data_list: List[str] = []):
        df = pd.DataFrame()
        for _c in self.box.pack():
            ser_dict = {}
            for key in data_list:
                if key in [
                    "path",
                    "name",
                    "filetype",
                    "label",
                    "flag",
                    "cause",
                    "energy",
                    "atoms",
                    "data",
                    "cache",
                    "total_charge",
                    "multiplicity",
                ]:
                    ser_dict[key] = getattr(_c, key)
                else:
                    ser_dict[key] = _c.data.get(key)
            _ser = pd.Series(ser_dict, name=_c.name)
            df = df.append(_ser)
            logger.info(f"data of {_c.name} was added to dataframe")
        return df
