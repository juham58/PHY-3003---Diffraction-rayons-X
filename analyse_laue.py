import pandas as pd
import numpy as np

def round_to_half(number):
    return round(number*2)/2

def check_if_half(number):
    if number.is_integer():
        return False
    else:
        return True

def laue_import(distance_cristal, filename, a):
    R = distance_cristal
    data = pd.read_csv(filename)
    data.set_index(" ", inplace=True)
    data.drop(["Area", "Mean", "Min", "Max"], axis=1, inplace=True)
    data["r"] = np.sqrt(data["X"]**2+data["Y"]**2)

    data["u"] = -1/(np.tan(0.5*np.arctan(data["r"]/R)))*(data["X"]/data["r"])
    data["v"] = -1/(np.tan(0.5*np.arctan(data["r"]/R)))*(data["Y"]/data["r"])

    data["h"] = np.nan
    data["k"] = np.nan
    data["l"] = np.nan

    for index, u in enumerate(data["u"]):
        u_rounded = round_to_half(u)
        if check_if_half(u_rounded):
            l = 2
        else:
            l = 1
        h = u_rounded*l
        data.loc[index+1, "h"] = h
        data.loc[index+1, "l"] = l

    for index, v in enumerate(data["v"]):
        v_rounded = round_to_half(v)
        if check_if_half(v_rounded):
            l = 2
        else:
            l = 1
        k = v_rounded*l
        data.loc[index+1, "k"] = k
    
    d_hkl = a/np.sqrt(data["h"]**2+data["k"]**2+data["l"]**2)
    data["lambda_exp"] = 2*d_hkl*np.sin(0.5*np.arctan(data["r"]/R))
    data["lambda_the"] = 2*d_hkl*data["l"]/np.sqrt(data["h"]**2+data["k"]**2+data["l"]**2)
    
    return data

#laue_import(15e-3, 'Results.csv', 562e-12) # test
