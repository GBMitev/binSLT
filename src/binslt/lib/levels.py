# %%
import pandas as pd
import numpy as np

from typing import Union
from binslt.utils.parse_files import StateLevels
import os

class Levels:
    
    def __init__(
        self,
        levels: pd.DataFrame,
        molecule_name: str,
        ) -> None:

        self.df = levels    
        self.molecule_name = molecule_name

    # def __repr__():
    #     representation = f"{molecule_name}: "
        
# %%
df = pd.read_csv("/home/gmitev/Documents/ExoMol/binSLT/prototyping/OH_levels.csv")
# %%

levels = Levels(df, "OH")
grouper = [i for i in df.columns if i not in ["L", "E", "NN"]]
levels.df.groupby(grouper, as_index = False).agg(count = ("L", "count"))["count"].max()



