from enum import Enum

class Grouping(Enum):

    STATES = (
        "Columns for grouping and ordering states files."
        ,
        [
            "Manifold",
            "tau"     ,
            "v"       ,
            "J"       ,
            "Lambda"  ,
            "Sigma"   ,
            "Omega"   ,
        ]
        )
    
    def describe(self):
        return self.value[0]
    
    def values(self):
        return self.value[1]