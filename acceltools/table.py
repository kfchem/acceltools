from typing import List

import pandas as pd
from accel.util.log import logger

from acceltools.base import ToolBox


class TableBox(ToolBox):
    def get_df(self, data_list: List[str] = []):
        df = pd.DataFrame()
        for c in self.get():
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
                    ser_dict[key] = getattr(c, key)
                else:
                    ser_dict[key] = c.data.get(key)
            _ser = pd.Series(ser_dict, name=c.name)
            df = pd.concat([df, pd.DataFrame([_ser])])
            logger.info(f"data of {c.name} was added to dataframe")
        return df
