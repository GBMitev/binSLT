# %%
from .dependencies  import *
from .distributions import gaussian
from .wrangling     import filter, cutoff
from .histogram     import get_histogram_data, gamma_estimate, get_centers
from .fitting       import fit_histogram

def fwhm_to_lifetime(fwhm):
    """
    Converts energy profile width into predissociation lifetime

    Inputs:
        fwhm = energy profile width (cm-1) : value (flaot)
    Outputs:
        lifetime = energy level lifetime (s) : value (flaot)
    """
    gamma = fwhm*1.98630e-23
    hbar = 1.054571817e-34

    lifetime = hbar/(gamma)
    return lifetime

def lifetime_to_fwhm(lifetime):
    """
    Converts predissociation lifetime into energy profile width

    Inputs:
        lifetime = energy level lifetime (s) : value (float)
    Outputs:
        fwhm = energy profile width (cm-1) : value (float)
    """
    hbar = 1.054571817e-34
    gamma = hbar/lifetime

    fwhm  = gamma/1.98630e-23
    return fwhm

def lifetime(df,J,v,ef,nsigma,bins, moct="mean",LE=None):
    '''
    Returns lifetime for a given energy level with cutoff, nsigma

    Inputs:
        df          = Full data set             : pd.DataFrame  see wrangling.filter()
        J           = J quantum number          : value (float) see wrangling.filter()
        v           = v quantum number          : value (int)   see wrangling.filter()
        ef          = e/f quantum number        : value (str)   see wrangling.filter()
        nsigma      = Number of Std for cutoff                      : value (float) see wrangling.cutoff()
        bins        = Number of bins            : value (int)   see histogram.get_histogram_data()
        moct        = Measure of central tendancy - mean or median  : value (str)   see wrangling.cutoff()
        LE = L and E vals to use in place of wrangling.filter() results : 2D list (float) unpacks as L,E=LE
    Outputs:
        lifetime    = Lifetime (s)                 : quantiphy.Quantity (float)
        x0          = central energy offset (cm-1) : value (float)
        fwhm        = energy profile width  (cm-1) : value (float)   
        mean        = histogram center             : value (float)
    '''
    if LE is None:
        L,E = filter(df, J, v, ef)
        L,E = cutoff(L,E,NSigma=nsigma,moct=moct)
    else:
        L,E = LE
        L,E = cutoff(L,E,NSigma=nsigma,moct=moct)

    count, edges, mean = get_histogram_data(E, bins)

    fwhm_guess  = gamma_estimate(edges, count, Plot = False)
    x0_Guess    = 1

    guesses = [fwhm_guess,x0_Guess]
    
    popt, lifetime = fit_histogram(count, edges,guesses = guesses)
    
    x0     = popt[0]
    fwhm   = popt[1]

    return lifetime, x0, fwhm, mean

def lifetime_over_bins(df, J, v, ef, nsigma, bin_range, moct="mean", LE = None, progress_bar = False, bar_pos=0):
    '''
    Calculates the lifetime as a function of the number of bins in the histogram for a given cutoff NSigma. 

    As some fits fail, that bin number is removed from the ActiveBins object

    Inputs:
        df          = Full data set             : pd.DataFrame  see wrangling.filter()
        J           = J quantum number          : value (float) see wrangling.filter()
        v           = v quantum number          : value (int)   see wrangling.filter()
        ef          = e/f quantum number        : value (str)   see wrangling.filter()
        nsigma      = Number of Std for cutoff                      : value (float) see wrangling.cutoff()
        bin_range   = Range of bin number       : list  (float)   
        moct        = Measure of central tendancy - mean or median  : value (str)   see wrangling.cutoff()
        LE = L and E vals to use in place of wrangling.filter() results : 2D list (float) unpacks as L,E=LE see lifetime()
    Outputs:
        fit_info  = dataframe containing ActiveBins, Lifetimes, x0, FWHM, and Mean: pd.DataFrame
    '''
    Lifetimes  = [] 
    ActiveBins = []
    
    x0         = []
    FWHM       = []

    Mean       = []

    if progress_bar == True:
        bin_range = tqdm(bin_range, desc=f"NSigma = {nsigma}, J = {J}, v = {v}, e/f = {ef}", position = bar_pos)
    
    if LE is None:
        LE = filter(df, J, v, ef)

    for bin in bin_range:
        try:
            lifetime_output = lifetime(df,J,v,ef,nsigma,bin, moct=moct, LE=LE)

            Lifetimes .append(lifetime_output[0])
            ActiveBins.append(bin)

            x0        .append(lifetime_output[1])
            FWHM      .append(lifetime_output[2])
            Mean      .append(lifetime_output[3])
        except:
            pass
    
    fit_info = {"ActiveBins": ActiveBins, "Lifetimes":Lifetimes, "x0":x0, "FWHM":FWHM, "Mean":Mean}
    
    fit_info = pd.DataFrame(fit_info)
    return fit_info

def statistical_lifetime(lifetimes, bins, sqrt_bins = False, params=False, state_label = None):
    '''
    Uses the lifetime over the bins to produce a statistical value for the lifetime with an uncertainty by fitting a Gaussian to the spread of lifetimes. 
    Inputs:
        lifeitmes           = Lifetimes calculated in lifetime_over_bins()                   : list  (float)
        bin_range           = Range of bin number used in lifetime_over_bins()               : list  (float)   
        sqrt_bins           = Use sqrt(len(lifetimes)) as bin number, else len(lifetimes)/10 : value (bool)
        params              = Return paramters of gaussian and histogram properties          : value (bool)
        state_label         = If statistical_lifetime() fails, label to print with failure statement, if None, nothing is printed : value, list (any) 
    Outputs:
        if params == False:
            measured_lifetime    = Mean of Gaussian               : value (float)   
            uncertainty          = Standard deviation of Gaussian : value (float)   
            active_bins          = bin number for which a Lorentzian fit was possible : value (float)
        if params == True:
            measured_lifetime    = Mean of Gaussian                     : value (float)   
            uncertainty          = Standard deviation of Gaussian       : value (float)   
            popt                 = Optimized parameters for Gaussian    : list (float)
            count                = see histogram.get_histogram_data()   : list (float or int)
            edges                = see histogram.get_histogram_data()   : list (float)
            centers              = see histogram.get_centers()          : list (float)
    '''
    try:
        if len(lifetimes) < len(bins)/3:
            measured_lifetime = np.nan
            uncertainty      = np.nan
            return measured_lifetime, uncertainty, len(lifetimes)
        
        bin_number = int(np.floor(len(lifetimes)/10)) if sqrt_bins == False else int(np.floor(np.sqrt(len(lifetimes))))
        count, edges, _ = get_histogram_data(np.array(lifetimes),bin_number,Centered=False)

        centers = get_centers(edges)

        estimated_mean = np.mean(lifetimes)
        estimated_std  = np.std(lifetimes)

    
        popt = curve_fit(gaussian,centers, count, p0 = [estimated_mean, estimated_std])[0]

        measured_lifetime = popt[0]  
        uncertainty      = popt[1]  

        return (measured_lifetime, uncertainty, len(lifetimes)) if params ==False else (measured_lifetime, uncertainty, popt, count, edges, centers)

    except:
        if state_label is not None:
            print(f"Error unknown for state_label = {state_label}")
        return (np.nan, np.nan, len(lifetimes))