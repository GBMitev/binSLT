# %%
from .dependencies import *

def diff(x, fx, sampling_frequency = 1):
    """
    Input: 
    x, fx: Input curve x and fx values - dtypes = np.array
    sampling frequency: number of data points between gradient pair - dtype = int
    
    Output:
    x: truncated x array accounting for lost terms - dtype = np.array
    f_prime: 1st differential f_prime(x) points - dtype = np.array
    """
    
    dx = x[sampling_frequency] - x[0] 
    f_prime = []
    
    for i in range(len(fx)-sampling_frequency):
        f_prime.append((fx[i+sampling_frequency]-fx[i])/dx)
        
    x = x[:(len(f_prime)-len(x))]
    return x, f_prime

def sav_gol(Lifetimes, window_length = 100, polyorder = 1, **kwargs):
    '''
    Applies Savitzky-Golay filter on Lifetime data:

    Inputs:
        Lifeitmes     = Lifetimes by bin number     : list  (float)     
        window_length = Savgol window length        : value (int) 
        polyorder     = Savgol polynomial order     : value (int)   
    Outputs:        
        Filtered      = Filtered data               : list  (float)
    '''
    Filtered = savgol_filter(Lifetimes, window_length=window_length, polyorder=polyorder, **kwargs)
    return [*Filtered]

def spline_smoothing(ActiveBins, Filtered, Bins, Degree = None):
    '''
    Nonsense
    '''

    Degree = len(Bins) if Degree is None else Degree

    SplineSmoothed = BSpline(*splrep(ActiveBins, Filtered, s=Degree))(Bins)
    return SplineSmoothed

def chi_squared(Obs, Calc, ddof=None, reduced = True):
    "ChiSquared test"
    Obs = np.array([*Obs])
    Calc = np.array([*Calc])

    Chi = sum(((Obs - Calc)**2)/Calc)

    Degrees_Of_Freedom =  len(Obs)-ddof-1 if ddof is not None else None

    p_value = 1-chi2.cdf(Chi, Degrees_Of_Freedom)
    
    #print("pvalue=",p_value)

    Chi_reduced = Chi/Degrees_Of_Freedom if Degrees_Of_Freedom is not None else Chi

    return (Chi_reduced, p_value) if reduced == True else (Chi, p_value)

def to_pickle(fname, cucumber):
    with open(fname, "wb") as handle:
        pickle.dump(cucumber, handle, protocol=pickle.HIGHEST_PROTOCOL)

def read_pickle(fname):
    with open(fname, "rb") as handle:
        cucumber = pickle.load(handle)
        return cucumber