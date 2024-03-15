# %%
from .dependencies import *

def lorentzian(x: list,x0: float,Gamma: float):
    '''
    Returns Lorentzian distribution
    
    Inputs:
        x            = dependent variable values    : list  (float)
        x0           = mean                         : value (float)  
        Gamma        = HWHM                         : value (float)  
    
    Outputs:
        Lorentz      = Lorentzian Distribution      : list  (float)
    '''
    Numerator = (1/np.pi)*(0.5*Gamma)
    Denomenator = ((x-x0)**2+(0.5*Gamma)**2)
    
    Lorentz = Numerator/Denomenator 

    return Lorentz

def gaussian(x:list, mu:float, sigma:float):
    '''
    Returns Gaussian distribution
    
    Inputs:
        x            = dependent variable values    : list  (float)
        mu           = mean                         : value (float)  
        sigma        = standard deviation           : value (float)  
    
    Outputs:
        Gauss      = Gaussian Distribution          : list  (float)
    '''
    Gauss = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)
    return Gauss
# %%
