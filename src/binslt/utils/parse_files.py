# %%
import pandas as pd
import numpy as np
import tarfile
import bz2
import io
import warnings

from binslt.utils.schema import Schema
from binslt.utils.grouping import Grouping

from typing import Optional
from tqdm import tqdm

class StateFile:
    """
    Class to store a states file loaded in from a tarball
    """
    def __init__(
            self, 
            molecule_name:str,
            schema:dict,
            box_length: float,
            tar_member: tarfile.TarInfo,
            tar_obj: tarfile.TarFile,
            ):
        
        self.molecule_name = molecule_name
        self.box_length = box_length
        self.tar_member = tar_member
        self.tar_obj = tar_obj        
        self.schema = schema

    def load(
            self,
            **kwargs
    ) -> pd.DataFrame:
        
        f = self.tar_obj.extractfile(self.tar_member)
        if f is None:
            raise FileNotFoundError(f"Cannot extract {self.tar_member.name}")
        
        names = list(self.schema.keys())
        
        with bz2.open(f) as bz:
            df = pd.read_csv(
                bz,
                sep = "\s+", 
                names = names, 
                dtype = self.schema,
                )
        
        return df

    def __repr__(self):
        representation = f"StateFile: {self.molecule_name} - <{self.box_length}>"
        return representation

class States:
    
    @staticmethod
    def _get_box_length(
        tar_info:tarfile.TarInfo, 
        molecule_name:str,
        ):

        name = tar_info.name

        prefix = f"{molecule_name}_"        
        extension = ".states.bz2"
        box_length = name.strip(prefix).strip(extension)

        try: 
            box_length = float(box_length)
            return box_length
        
        except:
            raise ValueError(f"{name} not in expected format: <molecule_name>_<box_length>.states.bz2, check your tarball.")

    def __init__(
            self, 
            tar_path: str,
            molecule_name:str,
            schema: Optional[dict] = None
            ):

        schema = schema if schema is not None else Schema.STATES.values()
        self.molecule_name = molecule_name

        self.tar_path = tar_path
        self.tar = tarfile.open(tar_path, "r")
        
        self.tar_members = self.tar.getmembers()
        
        self.states = {}
        for member in self.tar_members:
            
            box_length = self._get_box_length(
                member,
                molecule_name,
                )
            
            state_file = StateFile(
                molecule_name,
                schema,
                box_length,
                member,
                self.tar
            )
            
            box_length_key = str(float(box_length))

            if box_length_key in self.states.keys():
                raise ValueError(f"Box Length: {box_length_key} is duplicated in the dataset")
            
            else:
                self.states[box_length_key] = state_file
            
    def get(
            self, 
            box_length: float
            ):
        
        key = str(float(box_length))
        
        if key not in self.states.keys():
            raise ValueError(f"{box_length} not in States instance.")
        
        return self.states[key]

class StateLevels(States):

    def __init__(
            self, 
            tar_path,
            molecule_name,
            target_states:list[str],
            schema=None,
            ):
        
        if schema == None:
            Schema.STATES.values()

        super().__init__(tar_path, molecule_name, schema)

        self.target_states = target_states
        self.levels = None
        
    def load_levels(self):

        chunks = []
        for box_length, states in tqdm(self.states.items(), total=len(self.states)):
            states = states.load()
            states = states[states["Manifold"].isin(self.target_states)]
            states = states.assign(L = float(box_length))

            chunks.append(states)

        self.levels = pd.concat(
            chunks, 
            ignore_index=True
            )
        self.levels.sort_values(
            Grouping.STATES.values()+["L"], 
            inplace=True
            )

    def save_levels(self, path):
        path = path if path is not None else f"./{self.molecule_name}_levels.csv"
        self.levels.to_csv(path, index = False)