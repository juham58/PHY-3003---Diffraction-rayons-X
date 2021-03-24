"""
@author: Justin Hamel

https://github.com/juham58/PHY-3003---Diffraction-rayons-X
"""

from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})

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
    data["X"] = (data["X"]-data.loc[1, "X"])/1e6
    data["Y"] = (data["Y"]-data.loc[1, "Y"])/1e6
    data["r"] = np.sqrt(data["X"]**2+data["Y"]**2)

    data["u"] = -1/(np.tan(0.5*np.arctan(data["r"]/R)))*(data["X"]/data["r"])
    data["v"] = -1/(np.tan(0.5*np.arctan(data["r"]/R)))*(data["Y"]/data["r"])

    data["h"] = np.nan
    data["k"] = np.nan
    data["l"] = np.nan

    for index, u in enumerate(data["u"]):
        if index != 0:
            u_rounded = round_to_half(u)
            data.loc[index+1, "u_rounded"] = u_rounded
            if check_if_half(u_rounded):
                l = 2
            else:
                l = 1
            h = u_rounded*l
            data.loc[index+1, "h"] = h
            data.loc[index+1, "l"] = l
        else:
            data.loc[index+1, "h"] = 0
            data.loc[index+1, "l"] = 0

        for index, v in enumerate(data["v"]):
            if index != 0:
                v_rounded = round_to_half(v)
                data.loc[index+1, "v_rounded"] = v_rounded
                if check_if_half(v_rounded):
                    l = 2
                else:
                    l = 1
                k = v_rounded*l
                data.loc[index+1, "k"] = k
            else:
                data.loc[index+1, "k"] = 0

    
    d_hkl = a/np.sqrt(data["h"]**2+data["k"]**2+data["l"]**2)
    n = data["h"]**2+data["k"]**2+data["l"]**2
    #data["d_hkl"] = d_hkl
    data["n"] = n
    data["lambda_exp"] = 2*d_hkl*np.sin(0.5*np.arctan(data["r"]/R))
    data["lambda_the"] = 2*d_hkl*data["l"]/np.sqrt(data["h"]**2+data["k"]**2+data["l"]**2)
    data["lambda_error"] = np.abs((data["lambda_exp"]-data["lambda_the"])/data["lambda_the"])*100

    return data

#print(laue_import(15e-3, 'Results.csv', 562e-12)) # test
#print(laue_import(15e-3, 'nacl_14mm_35kvp_1images.csv', 562e-12)) # test

def laue_mass_import(crystal):
    if crystal == "si":
        dir = Path.cwd()/"si"
        a = 542e-12
    if crystal == "nacl":
        dir = Path.cwd()/"nacl"
        a = 562e-12
    if crystal == "lif":
        dir = Path.cwd()/"lif"
        a = 402e-12
    files = dir.glob("*.csv")
    for result in files:
        #print(result)
        #print(str(result).split("_")[1][0:2])
        distance = float(str(result).split("_")[1][0:2])*1e-3
        print("-------------\n",result, ":\n")
        print(laue_import(distance, result, a))
        results_path = str(result).split("\\")[-1]+"_results.csv"
        laue_import(distance, result, a).to_csv(Path.cwd()/"Résultats"/results_path)


def laue_graph():
    dir = Path.cwd()/"Résultats"
    files = dir.glob("*.csv")
    for result in files:
        data = pd.read_csv(result)
        crystal_name = str(result).split("\\")[-1].split("_")[0]
        if crystal_name == "nacl":
            crystal_name = "NaCl"
        if crystal_name == "lif":
            crystal_name = "LiF"
        if  crystal_name == "si":
            crystal_name = "Si"
        distance_value = str(result).split("\\")[-1].split("_")[1][0:2]
        voltage_at_peak = str(result).split("\\")[-1].split("_")[2][0:2]
        number_of_images = str(result).split("\\")[-1].split("_")[3][0:2]
        title = "Coordonnées u et v acquises avec le cristal de {}, à une distance de {} mm,\nune tension au pic de {} kV et un nombre d'images moyennées de {}.".format(crystal_name, distance_value, voltage_at_peak, number_of_images)
        print(data.to_latex(caption=title))
        
        plt.figure(figsize=(16,8))
        plt.scatter(data["u"], data["v"])
        plt.title(title)
        plt.xlabel(r"$u=\frac{h}{l}$")
        plt.ylabel(r"$u=\frac{k}{l}$")
        plt.grid()
        plt.savefig(str(result)+".png")


#laue_mass_import("lif")
#laue_mass_import("nacl")
#laue_mass_import("si")
laue_graph()


plt.show()

