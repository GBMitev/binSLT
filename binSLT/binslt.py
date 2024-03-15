# %%
from .dependencies      import *
from .wrangling         import le_wrt_duo,filter, cutoff
from .lifetime          import lifetime_over_bins, statistical_lifetime
from .histogram         import get_histogram_data
from .fitting           import fit_histogram

class FitParams:
    '''
    class to hold fitting parameters for Lorentzian fits for a given sigma over bins for use in class FitInfos
    '''
    def __init__(self, lorentzian_fit_parameters, sigma): 
        '''
        stores lorentzian fit paramters for a given cutoff sigma

        Inputs:
            lorentzian_fit_parameters = x0 and FWHM for each fit over active_bins : list (float)
            sigma = cutoff sigma at which lorentzian_fit_paramters was produced   : value(float)
        '''
        self.df     = lorentzian_fit_parameters
        self.sigma  = sigma
class FitInfos:
    def __init__(self, list_of_fit_params):
        '''
        class to store fitting Lorentzian fits for varying sigmas over bins

        Inputs:
            lorentzian_parameters = class instance of LorentzianParameters : LorentzianParameters            
        '''
        dict = {f"{i.sigma}":i for i in list_of_fit_params}
        
        sigmas = [int(i.sigma) for i in list_of_fit_params]
        
        fit_infos = pd.DataFrame(data = zip(sigmas,list_of_fit_params), columns=["Sigma","FitParams"])
        
        self.fit_infos = fit_infos
    def get_fit_info(self, sigma):
        '''
        returns lorentzian fit parameters over active bins for a given cutoff sigma if it is contained in the FitInfos class instance atributes

        Inputs:
            sigma = cutoff sigma for which lorentzian parameters are requested : value (float or int, depending on sigma type)
        Outputs:
            fit_info = table containing lorentzian fit parameters for given cutoff sigma : pd.DataFrame 
        '''
        try:
            fit_info = self.fit_infos.loc[self.fit_infos["Sigma"] == sigma]["FitParams"].to_numpy()[0].df
            return fit_info
        except:
            raise ValueError(f"Could not return fit info for sigma = {sigma}, check if this sigma has been used.")

def slt_over_sigmas(df, sigma_range, bin_range, J,v,ef, Duo = None, lock = None,moct=None,try_without_Duo=False):
    '''
    Returns statistical lifetime over multiple cutoff sigmas using lifetime.lifetime_over_bins() for listed cutoff sigmas (argument: sigmas)

    Inputs:
    Positional arguments:
        df          = Full data set             : pd.DataFrame  see lifetime.lifetime_over_bins()
        sigma_range = range of cutoff sigmas over which to use lifetime.lifetime_over_bins() : list (float or int)
        bin_range   = Range of bin number       : list  (float) 
        J           = J quantum number          : value (float) see lifetime.lifetime_over_bins()
        v           = v quantum number          : value (int)   see lifetime.lifetime_over_bins()
        ef          = e/f quantum number        : value (str)   see lifetime.lifetime_over_bins()
    Key word arguments
        Duo = Duo *.states file used as reference see wrangling.le_wrt_duo() : pd.DataFrame 
        moct        = Measure of central tendancy - mean or median  : value (str) see lifetime.lifetime_over_bins()
        lock        = lock = only return geometries/energies where energies = duo_filter(*args, **kwargs) +/- lock : value (float) see wrangling.le_wrt_duo()
        try_without_duo = if wrangling.le_wrt_duo() fails, use wrangling.filter() : value (bool)
    '''
    slts        = []
    fit_infos   = []
    
    try:
        if Duo is not None:
            # print("Duo is not none")
            try:
                LE = le_wrt_duo(df, Duo,J,v,ef,lock=lock)
            except:
                if try_without_Duo==True:
                    LE = filter(df, J, v, ef)
                else:
                    average_slt = ufloat(np.nan, np.nan)
                    slts = [ufloat(np.nan,np.nan)]*len(sigma_range)
                    fit_infos = np.nan
                    output = [average_slt, *slts, fit_infos]
                    return output
        else:
            # print("Duo is None")
            LE = None

        # print("trying sigmas")
        for sigma in sigma_range:
            # print(f"trying sigma = {sigma}")
            fit_info  = lifetime_over_bins(df,J,v,ef,sigma,bin_range,moct=moct,LE=LE)        
            lifetimes = fit_info["Lifetimes"].to_numpy()
            fit_infos.append(FitParams(fit_info,sigma))

            # print("Getting statistical lifetime")
            stat_lifetime, uncertainty, _ = statistical_lifetime(lifetimes, bin_range)
            slts.append(ufloat(stat_lifetime, uncertainty))

        # print("unpacking data")
        data = [[val.n,val.s] for val in slts if np.isnan(val.n) == False and np.isnan(val.s)==False]


        if len(data) != 0:
            lifetimes, uncs = zip(*data)
            lifetimes, uncs = np.array(lifetimes), np.array(uncs)

            lifetime_ave    = np.mean(lifetimes)
            unc_ave         = 1/len(lifetimes)*np.sqrt(sum(uncs**2))

            average_slt = ufloat(lifetime_ave,unc_ave) 
        
        if len(data) == 0:
            average_slt = ufloat(np.nan, np.nan)

        fit_infos = FitInfos(fit_infos)
        output = [average_slt, *slts, fit_infos]
        return output
    except:
        print(f"Something Went Wrong: J = {J}, v = {v}, e/f = {ef}")

def slt_over_levels(df, sigma_range, bin_range, Duo=None, cores = None,lock=None,moct=None,try_without_Duo = False):
    '''
    returns statistical lifetimes over all levels in the full stabilization dataset, currently only supports doublet Sigma states

    Inputs:
    Positional arguments:
        df          = Full data set             : pd.DataFrame  
        sigma_range = range of cutoff sigmas over which to use lifetime.lifetime_over_bins() : list (float or int)
        bin_range   = Range of bin number       : list  (float) 
    Key word arguments:
        cores           = number of cores to use in parallel processing of statistical lifetimes : value (int)
        Duo             = Duo *.states file used as reference see wrangling.le_wrt_duo()         : pd.DataFrame 
        try_without_Duo = if wrangling.le_wrt_duo() fails, use wrangling.filter()                : value (bool)
        moct            = Measure of central tendancy - mean or median                           : value (str) see lifetime.lifetime_over_bins()
        lock            = lock = only return geometries/energies where energies = duo_filter(*args, **kwargs) +/- lock : value (float) see wrangling.le_wrt_duo()
    '''
    slt_df = df.groupby(["J","v","e/f","State","Lambda","Sigma"],as_index=False).agg({"E":"mean"})[["J","v","e/f","State","Lambda","Sigma"]]
    
    cores = cores if cores is not None else core.NB_PHYSICAL_CORES
    pandarallel.initialize(nb_workers = cores,progress_bar=True, verbose = 0)

    sigma_column_names = [str(i) for i in sigma_range]

    slt_df[["Lifetime",*sigma_column_names, "FitInfo"]] =  slt_df.parallel_apply(
        lambda x: slt_over_sigmas(df, sigma_range, bin_range
                                  ,x["J"],x["v"],x["e/f"]
                                  ,Duo=Duo,lock=lock,moct=moct,try_without_Duo=try_without_Duo)
        ,axis=1
        ,result_type = "expand")
    
    return slt_df

class BinSLT:
    def __init__(self, df, sigma_range, bin_range, cores = None, Duo = None, moct = None, lock = None, **kwargs):
        '''
        Establishes class instance of BinSLT and performs calculation of statistical lifetimes over all levels in full stabilization data set, 'df'.

        Inputs:
        Positional Arguments:
            df  = Full data set         : pd.DataFrame
            sigma_range = range of cutoff sigmas over which to use lifetime.lifetime_over_bins() : list (float or int)
            bin_range   = Range of bin number       : list  (float) 
        
        Keyword Arguments:
            cores = number of cores to use in parallel processing of statistical lifetimes
            Duo = Duo *.states file used as reference see wrangling.le_wrt_duo() : pd.DataFrame
            moct        = Measure of central tendancy - mean or median  : value (str) see lifetime.lifetime_over_bins()
            lock        = lock = only return geometries/energies where energies = duo_filter(*args, **kwargs) +/- lock : value (float) see wrangling.le_wrt_duo()
        
        **kwargs:
            try_without_duo = if wrangling.le_wrt_duo() fails, use wrangling.filter() : value (bool)
            blt_init = for use with read_blt() function
        '''
        self.df     = df
        if sigma_range is not None:
            self.sigmas_used = [str(i) for i in sigma_range]

        if Duo is not None:
            self.Duo = Duo

        try_without_Duo = kwargs.get("try_without_Duo",False)


        if "blt_init" not in kwargs.keys():
            self.blt = slt_over_levels(self.df, sigma_range, bin_range, cores = cores, Duo=Duo, lock=lock, moct=moct,try_without_Duo=try_without_Duo)
        elif kwargs["blt_init"] == True:
            pass

    def separate_blt(self, nullnt = None):
        '''
        returns DataFrame with the Lifetimes and Uncertainties is separate columns

        Inputs:
            nullnt = use self.get_nullnt() : value (bool)
        '''
        lt = self.blt if nullnt == None else self.get_nullnt()

        lt["Unc"] = lt.apply(lambda x:x["Lifetime"].s, axis = 1)
        lt["Lifetime"] = lt.apply(lambda x:x["Lifetime"].n, axis = 1)

        return lt[["J","v","e/f","State","Lambda","Sigma","Lifetime","Unc"]]

    def lifetime(self, J,v,ef):
        '''
        returns lifetime of a given level with quantum numbers J, v, ef
        '''
        df = self.blt[["J","v","e/f","Lifetime",*self.sigmas_used,"FitInfo"]]

        #Progressively filtering
        df = df[(df["J"]==J)&(df["v"]==v)&(df["e/f"]==ef)]
        tau = df["Lifetime"].to_numpy()
        return tau

    def get_fit_info(self,J, v, ef, sigma):
        '''
        returns lorentzian parameters over bins for a given sigma cutoff with quantum numbers J, v, ef
        '''
        df = self.blt[["J","v","e/f","FitInfo"]]
        df = df[
            (df["J"] == J)   & 
            (df["v"] == v)   &
            (df["e/f"] == ef)]
        if len(df) > 1:
            print("Warning")
        fit_info = df["FitInfo"].to_numpy()[0].get_fit_info(sigma)
        return fit_info
    
    def get_null(self):
        '''
        return quantum numbers for energy levels for which a statistical lifetime calculation failed
        '''
        LT = self.blt[["J","v","e/f","State","Lambda","Sigma","Lifetime"]]
        null = LT[np.isnan(unumpy.nominal_values(LT["Lifetime"]))==True].reset_index(drop = True)

        return null
    
    def get_nullnt(self): 
        '''
        return quantum numbers + lifetimes for energy levels for which a statistical lifetime calculation was successful
        '''  
        LT = self.blt[["J","v","e/f","State","Lambda","Sigma","Lifetime","FitInfo"]]
        nullnt = LT[np.isnan(unumpy.nominal_values(LT["Lifetime"]))==False].reset_index(drop = True)
        return nullnt
    
    def make_sandwich(self,fname):
        '''
        creates *.blt file at 'fname' to store class instance of BinSLT for future use.
        '''
        if fname[-4:] == ".blt":
            pass
        else:
            fname += ".blt"
            
        with open(fname, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)        
    
    def get_good_bins(self):
        print("will be developed at some point")

def read_blt(fname):
    '''
    reads *.blt file and establishes instance of BinSLT with atributes from given *.blt file.
    '''
    with open(fname, 'rb') as handle:
        Lvls = pickle.load(handle)
    Object = BinSLT(None, None, None, blt_init=True)
    Object.__dict__ = Lvls.__dict__
    return Object

def convergence_testing(df, Duo, J, v, ef, lock, load, bins,window,nsigma=None):
    '''
    Convergence testing for a given energy level with quantum numbers J, v, ef (currently only supports 2Sigma states)

    Inputs:
        df  = Full data set         : pd.DataFrame
        Duo = Duo *.states file used as reference see wrangling.le_wrt_duo() : pd.DataFrame
        J           = J quantum number          : value (float) see lifetime.lifetime_over_bins()
        v           = v quantum number          : value (int)   see lifetime.lifetime_over_bins()
        ef          = e/f quantum number        : value (str)   see lifetime.lifetime_over_bins()
        lock        = lock = only return geometries/energies where energies = duo_filter(*args, **kwargs) +/- lock : value (float) see wrangling.le_wrt_duo()
        load        = array of indecies over which to 'load' the energy profile : series (int)
        bins        = Number of bins            : value (int)   see histogram.get_histogram_data()
        window      = window for moving average and variance : value (int)
        nsigma      = cutoff sigma for use in wrangling.cutoff() : value (float)
    Outputs:
        percentage_loaded = how much of the energy profile has been 'loaded' as a percentage : np.1darray (float)
        lifetimes         = lifetime as a function of percentage_loaded : np.1darray (float)
        moving_average    = moving average lifetime as a function of percentage_loaded : np.1darray (float)
        moving_var        = moving variance of the lifetime as a function of percentage_loaded : np.1darray (float)
    '''

    try:
        nsigma = nsigma if nsigma is not None else np.inf
        L_init,E_init = cutoff(*le_wrt_duo(df, Duo, J, v, ef, lock),nsigma,moct="mean")
        lifetimes = []

        for i in load:
            try:
                L,E = L_init[:i],E_init[:i]
                count, edges, _ = get_histogram_data(E,bins)
                _, lt = fit_histogram(count, edges)
                lifetimes.append(lt*1e12)
            except:
                lifetimes.append(lifetimes[-1])
                
        lower = np.mean(lifetimes)-5*np.std(lifetimes)
        upper = np.mean(lifetimes)+5*np.std(lifetimes)

        load, lifetimes = np.transpose([[load[num], l] for num, l in enumerate(lifetimes) if lower<= l <=upper])
        
        moving_var = pd.Series(lifetimes).rolling(window).std()
        moving_average = pd.Series(lifetimes).rolling(window).mean()
        
        percentage_loaded = [100*i/max(load) for i in load]
        # print("done")
        return percentage_loaded,lifetimes, moving_average.to_numpy(),moving_var.to_numpy()
    except:
        return [np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan]