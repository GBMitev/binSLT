# %%
from .dependencies import *
from .histogram import get_centers
from .distributions import lorentzian

def fit_histogram(Count,Edges, formatted = True, guesses = None):
    '''
    Fits Lorentzian function to Histrogram Count and Centers and extracts halfwidth

    Inputs: 
        Count       = Count for bin         : np.1darray    (float)
        Edges       = Bin edges             : np.1darray    (float)
        formatted   = format of lifetime    : bool          
    Outputs:
        popt        = fitted parameters     : list          (float)
        lifetime    = lifetime of state     : value         (float or str)
    '''

    Centers = get_centers(Edges)
    
    if guesses == None:
        guesses = [0,1]
    else:
        guesses = guesses        

    popt = curve_fit(lorentzian,Centers,Count,p0=guesses)[0]
    FWHM = popt[1]

    gamma = FWHM*1.98630e-23
    hbar = 1.054571817e-34

    lifetime = hbar/(gamma)

    return [popt, Quantity(lifetime,"s")] if formatted == True else [popt, lifetime]
# %%
