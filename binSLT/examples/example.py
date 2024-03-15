# %%
from binslt import binslt as blt
from binslt.dependencies import *
# %%
df  = pd.read_csv("path/to/full_data_set")
Duo = pd.read_csv("path/to/Duo_energies.states",delim_whitespace = True, names = ["NN","E","g","J","tau","e/f","Manifold","v","Lambda","Sigma","Omega"])
# # %%
sigmas = [3,4,5]
bins   = np.arange(30,400,1)
lt = blt.BinSLT(
  df,
  sigmas, 
  bins, 
  Duo = Duo, 
  cores = 16,
  moct = "median",
  lock = 10,
  try_without_Duo = False,
  blt_init = False
  )

lt.make_sandwich("calculated_predissociation_lifetimes.blt")
