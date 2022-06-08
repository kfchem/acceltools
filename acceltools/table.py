from typing import List

import pandas as pd
from accel.util.log import logger

from acceltools.base import ToolBox


class TableBox(ToolBox):
    def get_df(self, data_list: List[str] = []):
        df = pd.DataFrame()
        for _c in self.mols:
            ser_dict = {}
            for key in data_list:
                if key in [
                    "path",
                    "name",
                    "filetype",
                    "label",
                    "flag",
                    "history",
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
            df = pd.concat([df, pd.DataFrame([_ser])])
            logger.info(f"data of {_c.name} was added to dataframe")
        return df
