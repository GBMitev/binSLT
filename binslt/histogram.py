# %%
from .dependencies import *

def get_histogram_data(E, bins, Centered = True, **kwargs):
    '''
    Returns Histogram plot data for distribution of Energy levels of a given state. 

    Inputs: 
        E           = Energies              : np.1darray    (float)    
        bins        = Number of Bins        : value         (float)   
        Centered    = Center on the mean    : bool (default = True)
    
    Outputs:
        Count       = Count for bin         : np.1darray    (float)
        Edges       = Bin edges             : np.1darray    (float)
        Mean        = Dataset mean          : value         (float)
    '''

    Mean = E.mean()    
    E = E-Mean if Centered==True else E
    density = kwargs.get("density", True)
    Count, Edges = np.histogram(E, bins=bins,density=density)

    return Count, Edges, Mean

def get_centers(Edges):
    '''
    Returns Centers of bins for given Edges

    Inputs: 
        Edges   = Bin edges from np.histogram   : np.1darray (float)
    Outputs:
        Centers = Bin centers from np.histogram : np.1darray (float)
    '''
    Centers = 0.5*(Edges[1:]+ Edges[:-1])
    return Centers 

def get_xrange(Edges):
    Centers = get_centers(Edges)
    return np.linspace(min(Centers), max(Centers), 1000)

def plot_histogram(count,edges, centered = True, **bar_kwargs):
    '''
    Plots Histogram from GetHistogram
    
    Inputs:
        Count       = Count for bin         : np.1darray    (float)
        Edges       = Bin edges             : np.1darray    (float)
        Centered    = Center on mean        : bool
    '''
    plt.figure(figsize = (9,7))
    plt.bar(edges[:-1], count,width=np.diff(edges),color="black", align="edge", **bar_kwargs)

    xlabel = r"$E-\langle E \rangle$ / cm $^{-1}$" if centered == True else r"$E$"
    plt.xlabel(xlabel, fontsize = 20)
    plt.ylabel("Count", fontsize = 20)

def measure_fwhm(Edges, Count, NumPoints=None):
    '''
    Measures the FWHM by interpolating Lorentzian distribution and finding intersection points.

    Inputs:
        Edges     = Bin edges of histogram pdf : np.1darray (float)
        Count     = Counts for histogram pdf   : np.1darray (float)
        NumPoints = Number of bins             : value      (int)
    Outputs: 
        Erange    = Interpolation x values     : np.1darray (float)
        Hrange    = Interpolation y values     : np.1darray (float)
        min_idx   = Index of frst intersection : value      (int)
        max_idx   = Index of scnd intersection : value      (int)
        FWHM_Interp_Measurement = MeasuredFWHM : value      (float)
    '''

    Centers = get_centers(Edges)
    
    NumPoints = 10000 if NumPoints is None else NumPoints

    Erange = np.linspace(min(Centers), max(Centers), NumPoints)

    H = interpolate.interp1d(Centers,Count)
    Hrange = H(Erange)

    InterpFWHM  = max(Hrange)/2 * np.ones(NumPoints)

    idx_Interp  = np.argwhere(np.diff(np.sign(InterpFWHM - Hrange))).flatten()

    if len(idx_Interp) > 2:
       raise ValueError("There are more than two intersection points, try lowering your bins, check your histogram.")
    
    min_idx = idx_Interp[0]
    max_idx = idx_Interp[1]
    FWHM_Interp_Measurement = abs(Erange[min_idx]-Erange[max_idx])

    return Erange, Hrange, min_idx, max_idx, FWHM_Interp_Measurement

def plot_gamma(Erange, Hrange, min_idx, max_idx):
    '''
    Plots FWHM of Lorentzian distribution as shaded area. 

    Inputs: 
        Erange    = Interpolation x values     : np.1darray (float)
        Hrange    = Interpolation y values     : np.1darray (float)
        min_idx   = Index of frst intersection : value      (int)
        max_idx   = Index of scnd intersection : value      (int)
    Oututs:
        None
    '''
    plt.figure(figsize=(9,7))

    plt.title("Measurement of FWHM for Lorentzian Distribution",fontsize = 20)

    plt.plot(Erange, Hrange, "k--", label = "Interpolated Distribution")
    #plt.plot(Erange[min_idx:max_idx+1], InterpFWHM[min_idx:max_idx+1], "b-", label = "FWHM")

    plt.xlabel(r"$E-\langle E \rangle$ / cm $^{-1}$",fontsize = 20)
    plt.ylabel("Count",fontsize = 20)
    #plt.xticks(fontsize = 20)
    #plt.yticks(fontsize = 20)

    #plt.legend(loc = "best", fontsize = 20)
    
    plt.axhline(min(Hrange), color="k" ,alpha = 0.5)
    plt.fill_between(Erange, Hrange, min(Hrange),
            where = (Erange >= Erange[min_idx]) & (Erange <= Erange[max_idx]),
            color = 'g')

def gamma_estimate(Edges, Count, **kwargs):
    '''
    Wrapper function to extract FWHM estimate and Plotting

    Inputs: 
        Edges     = Bin edges of histogram pdf : np.1darray (float)
        Count     = Counts for histogram pdf   : np.1darray (float)
    Outputs:
        FWHM_Interp_Measurement = MeasuredFWHM : value      (float)
    '''
    if "NumPoints" in kwargs.keys():
        NumPoints = kwargs["NumPoints"]
    else:
        NumPoints = None

    Erange, Hrange, min_idx, max_idx, FWHM_Interp_Measurement = measure_fwhm(Edges, Count, NumPoints)
    
    if "Plot" in kwargs.keys():
        Plot = kwargs["Plot"]
    else:
        Plot = False

    if Plot == True:
        plot_gamma(Erange, Hrange, min_idx, max_idx)
    
    return FWHM_Interp_Measurement

def add_histogram(Edges, Count, Centered = True):
    plt.bar(Edges[:-1], Count,width=np.diff(Edges),color="black", align="edge")