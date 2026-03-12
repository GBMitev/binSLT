# %%
from .dependencies import *
#
def filter(df,J:float,v:int,ef:str, state=None):
    '''
    Returns geometries and energies for a given energy level. 

    Inputs:
        df  = Full data set         : pd.DataFrame 
        J   = J quantum number      : value (float) 
        v   = v quantum number      : value (int) 
        ef  = e/f quantum number    : value (str)

    Outputs:
        L   = Geometries            : np.1darray (float)
        E   = Energies              : np.1darray (float)    
    '''

    if state is None:
        df = df[(df["J"]==J)&(df["v"]==v)&(df["e/f"]==ef)].sort_values("L")
        L,E = df["L"].to_numpy(),df["E"].to_numpy()
    elif state is not None:
        df = df[(df["J"]==J)&(df["v"]==v)&(df["e/f"]==ef)&(df["State"]==state)].sort_values("L")
        L,E = df["L"].to_numpy(),df["E"].to_numpy()
    
    return L,E

def cutoff(L: list,E: list,NSigma: float=np.inf, moct=None):
    '''
    Returns Geometries and Energies for Lower < Energy < Upper. 
    Lower = Mean(E)-Std(E)*NSigma 
    Upper = Mean(E)+Std(E)*NSigma 

    Inputs: 
        L       = Geometries            : np.1darray (float)
        E       = Energies              : np.1darray (float)  
        NSigma  = Number of Std         : value      (float)
    
    Outputs:
        L   = Geometries (Adjusted)     : np.1darray (float)
        E   = Energies   (Adjusted)     : np.1darray (float)  
    '''
    moct = "mean" if moct is None else moct

    if moct == "mean":
        moct = np.mean(E)
    elif moct =="median":
        moct = np.median(E)

    std  = np.std(E)

    upper = moct+std*NSigma
    lower = moct-std*NSigma

    data = [[L[num], e] for num, e in enumerate(E) if lower<= e <=upper]
    data = np.transpose(data)
    return data

def allowed_quantum_numbers(df):
    '''
    Returns all quantum number subsets in the total DataFrame

    Inputs:
        df  = Full data set             : pd.DataFrame 
    Outputs:
        QN  = Quantum number subsets    : pd.DataFrame
    '''

    QN = df.groupby(["J","v","e/f"], as_index=False).agg({"L":"count"})[["J","v","e/f","L"]]
    
    if len(QN[QN["L"]!=max(QN["L"])]) != 0:
        print("Inconsistent quantum number representation over geometries. Check your data.")
        #raise ValueError("Inconsistent quantum number representation over geometries. Check your data.")
    
    QN = QN[["J","v","e/f"]]

    return QN

def plot_LE(L,E,line = False,J = None, v = None, ef = None, scale_e=False):
    plt.figure(figsize = (9,9)) 
    E = E*1e-3 if scale_e == True else E

    fmt = "k." if line == False else "k"
    plt.plot(L,E, fmt)

    if v is not None and J is not None and ef is not None:
        title = f"v = {v}, J = {J}, e/f = {ef}"
        plt.title(title)

    xlabel = r"Box Length / $\AA$"
    ylabel = r"Energy / $10^{3}$ cm$^{-1}$" if scale_e == True else r"Energy / cm$^{-1}$"
    
    plt.xlabel(xlabel, fontsize = 20)
    plt.ylabel(ylabel, fontsize = 20)

    plt.grid(which = "both")
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)

def duo_filter(Duo, J, v, ef, manifold="A2Sigma+"):
    Duo=Duo[(Duo["J"]==J)&(Duo["v"]==v)&(Duo["e/f"]==ef)&(Duo["Manifold"]==manifold)]
    return Duo["E"].to_numpy()[0]

def le_wrt_duo(df, Duo, J, v, ef,lock=None):
    # print("getting lock")
    lock = 10 if lock is None else lock

    # print("getting LE")
    L,E = filter(df, J, v, ef)
    # print("got LE")
    Center = duo_filter(Duo, J, v, ef)
    Lower  = Center - lock
    Upper  = Center + lock

    data_cutoff = [[L[num], e] for num, e in enumerate(E) if Lower<= e <=Upper]
    # print("transposing")
    L,E = np.transpose(data_cutoff)
    # print("successful")
    return L,E