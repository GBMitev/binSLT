# %%
import pandas as pd
import numpy as np
import tarfile
import bz2
import io
import warnings

from binslt.utils.schema import Schema
from typing import Optional

class StateFile:
    """
    Class to store a states file loaded in from a tarball
    """
    def __init__(
            self, 
            box_length: float,
            tar_member: tarfile.TarInfo,
            tar_obj: tarfile.TarFile,
            ):
        
        self.box_length = box_length
        self.tar_member = tar_member
        self.tar_obj = tar_obj
        
        self.df = None
        self.name = None
        self._precision = len(str(self.box_length).split('.')[1]) 

    def load(
            self,
            **kwargs
    ) -> pd.DataFrame:
        
        f = self.tar_obj.extractfile(self.tar_member)
        if f is None:
            raise FileNotFoundError(f"Cannot extract {self.tar_member.name}")
        
        data = bz2.decompress(f.read())
        buffer = io.BytesIO(data)
        self.df = pd.read_csv(buffer, **kwargs)

    def _set_name(
            self, 
            name, 
            precision:Optional[int] = None,
            ):
    
        self.name = name
        if precision:
            self._precision = precision

    def __repr__(self):
        representation = f"StateFile: {self.name} - <{self.box_length:.{self._precision}f}>"
        return representation

class States:
    
    @staticmethod
    def _get_box_length(
        tar_info:tarfile.TarInfo, 
        molecule_name:str,
        prefix:str = None,
        extension:str = None,
        ):

        name = tar_info.name

        if not prefix:
            prefix = f"{molecule_name}_"
        
        if not extension:
            extension = ".states.bz2"

        box_length = name.strip(prefix).strip(extension)

        try: 
            box_length = float(box_length)
            return box_length
        
        except:
            raise ValueError(f"{name} not in expected format: <prefix><box_length><extension>")

    @staticmethod
    def _get_precision(box_lengths):
        precision = max([len(s.split('.')[1]) for s in box_lengths])
        return precision

    def __init__(
            self, 
            tar_path: str,
            molecule_name:str,
            prefix:Optional[str] = None,
            extension:Optional[str] = None,
            ):
    
        self.tar_path = tar_path
        self.tar = tarfile.open(tar_path, "r")
        
        self.tar_members = self.tar.getmembers()
        
        _states = {}
        for member in self.tar_members:
            if not member.isfile():
                warnings.warn(f"{member.name} is not a file, this will not be loaded")
                continue
            
            box_length = self._get_box_length(
                member,
                molecule_name,
                prefix,
                extension,
                )
            
            state_file = StateFile(
                box_length,
                member,
                self.tar
                )
            
            box_length_key = str(box_length)

            if box_length_key in _states.keys():
                raise ValueError(f"Box Length: {box_length_key} is duplicated in the dataset")
            
            else:
                _states[box_length_key] = state_file

        self._precision = self._get_precision(_states.keys())

        self.states = {}
        for box_length_key, state_file in _states.items():
            
            state_file._set_name(molecule_name, precision = self._precision)

            box_length = float(box_length_key)
            key = f"{box_length:.{self._precision}f}"
            
            self.states[key] = state_file

    def get(
            self, 
            box_length: float
            ):
        query = float(box_length)
        query_precision = len(str(query).split('.')[1])
        if query_precision > self._precision:
            raise ValueError(f"Box length: {query} does not exist.")
        
        box_length_key = f"{query:.{self._precision}f}"
        return self.states[box_length_key]

