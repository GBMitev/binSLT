# %%
# Import needed packages
from binslt import binslt as blt
from binslt.dependencies import *
# initialize full dataset (df) and Duo stable energies (Duo)
df  = pd.read_csv("./full_data_set.csv")
Duo = pd.read_csv("./Duo_stable_energies.states",delim_whitespace = True, names = ["NN","E","g","J","tau","e/f","Manifold","v","Lambda","Sigma","Omega"])
# %%
# initialize sigma_range and bin_range
sigmas = [3,4,5]
bins   = np.arange(30,400,1)
# run BinSLT algorithm
lt = blt.BinSLT(
    df,                       # full dataset
    sigmas,                   # sigma_range
    bins,                     # bin_range
    Duo = Duo,                #stable Duo energies
    cores = 16,               # number of cores to use in calculation
    moct = "median",          # measure of central tendency for data wranling, see wrangling.cutoff()
    lock = 10,                # allowed energy deviation from Duo stable energy, see wrangling.le_wrt_duo()
    try_without_Duo = False,  # if le_wrt_duo() fails, try using wrangling.filter()
    blt_init = False          # initialise empty BinSLT class instance, used with binslt.read_blt()
    )

# save solution as *.blt file
lt.make_sandwich("calculated_predissociation_lifetimes") #.blt file extension added automatically if not present
# %%
# load a solution from a *.blt file
lt = blt.read_blt("calculated_predissociation_lifetimes.blt")

# see dataframe containing all information about lifetimes and fits
lt_blt = lt.blt

# see dataframe for which lifetimes have been computed
lt_nullnt = lt.get_nullnt()

# extract fit information for energy level with J = 0.5, v = 5, e/f = e, sigma_cutoff = 3sigma
fit_example = lt.get_fit_info(0.5,5,"e",3)
# %%