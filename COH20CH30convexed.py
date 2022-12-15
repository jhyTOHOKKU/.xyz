import pandas as pd
import numpy as np
import random, math
from math import sqrt, pow

def z_convex_downward(df):
    for i in range(len(df)):
        if (df.at[i,'z'] < -1.4): #& (df.iloc[i].atom == 'C'):
            row_in_cyl = into_cyl(df.iloc[i])
            #converted_radius = ((-1. / radius_max) * row_in_cyl[0]) + 1 
            converted_z = (df.at[i,'z'] - 7.1 * np.cos( (math.pi / 2.) * (row_in_cyl[0] / radius_max) )) + 0.4
            df.at[i,'z'] = converted_z            
        else:
            continue
            
def into_cyl(a): # change rectangular coordinate into cylindrical coordinate
    r = sqrt(pow(a[-3], 2.) + pow(a[-2], 2.))
    if a[-3] == 0.: # if x = 0 
        if a[-2] > 0:
            theta = math.pi / 2. # if (0,positive) then theta = pi / 2
        else:
            theta = -1. * math.pi / 2. # if (0,negative) then theta =  - pi / 2            
    else: # if x is NOT 0 
        if (a[-3] > 0.):
            theta = np.arctan(a[-2] / a[-3]) #first quadrant OR fourth quadrant         
        else:
            theta = np.arctan(a[-2] / a[-3]) + math.pi #second quadrant OR third quadrant  
    cyl_cor = [r, theta, a[-1]] # a[-1] = z
    return cyl_cor

def oxy_bonding(cyl_cor):
    cyl_cor[1] = cyl_cor[1] - 1.428 # CO bond length = 1.428 Angstrom
    return cyl_cor

def hyd_bonding(cyl_cor):
    cyl_cor[1] = cyl_cor[1] - 2.389 # CO bond length + OH bond length = 2.289 Angstrom
    return cyl_cor

def carbohyd_bonding(cyl_cor):
    cyl_cor[1] = cyl_cor[1] - 1.094 # CH bond length  = 1.094 Angstrom
    return cyl_cor

carbontube = pd.read_csv("nt-30-30-30-gra.csv")
ID_number = np.arange(len(carbontube))    # 0 ~ (6320-1)

carbontube.set_index = ID_number
carbontube.insert(0, "ID_number", ID_number)
bool_wall = ( carbontube.z > -1.42 ) & ( carbontube.z < 73.9756 ) # bool_wall : carbon atoms only on tube, except ones on plane 
    # z = 73.9756 : plane above, z= -1.42 : plane below

num_of_functional_COH = int(len(carbontube[bool_wall]) * 0.2)     # 3600 * 20% = 720
num_of_functional_CH = int(len(carbontube[bool_wall]) * 0.3)     # 3600 * 30% = 1080

carbon_on_wall = carbontube.loc[bool_wall]

carbon_on_wall_to_cyl = pd.DataFrame(columns=["atom", "radius", "theta", "z"])

for i in range(len(carbon_on_wall_for_functional_COH20_ID_list)):
    carbon_row = into_cyl(carbon_on_wall.iloc[i])
    carbon_row.insert(0,"C")  
    carbon_on_wall_to_cyl.loc[len(COH_cyl)] = carbon_row
    
radius_max = carbon_on_wall_to_cyl["radius"].max()

carbon_on_wall_for_functional_COH20_ID_list = np.random.choice(len(carbon_on_wall), size=num_of_functional_COH, replace=False)
carbon_on_wall_for_functional_COH20_ID_list.sort()
#replace=False : allow no duplications

carbon_on_wall_randomly_chosen_for_COH20 = carbontube.iloc[carbon_on_wall_for_functional_COH20_ID_list]

#carbontube.drop(carbon_on_wall_for_functional_COH20_ID_list, inplace=True)

carbontube_excluding_COH_functional = carbontube #delete the randomly chosen 720 (20%) carbon atoms from the df carbontube

carbon_on_wall.drop(carbon_on_wall_for_functional_COH20_ID_list, inplace=True)
carbon_on_wall_excluding_COH_functional = carbon_on_wall

#carbon_on_wall_for_functional_CH30_ID_list = np.random.choice(len(carbon_on_wall_excluding_COH_functional), size=num_of_functional_CH, replace=False)

carbon_on_wall_for_functional_CH30_ID_list = random.sample(list(carbon_on_wall_excluding_COH_functional["ID_number"]),num_of_functional_CH)
carbon_on_wall_for_functional_CH30_ID_list.sort()
# random.sample : https://www.geeksforgeeks.org/python-random-sample-function/

carbon_on_wall_for_all_functionals_list = [*carbon_on_wall_for_functional_CH30_ID_list, *carbon_on_wall_for_functional_COH20_ID_list]
carbon_on_wall_for_all_functionals = carbontube.iloc[carbon_on_wall_for_all_functionals_list]

carbon_on_wall_randomly_chosen_for_CH30 = carbontube.iloc[carbon_on_wall_for_functional_CH30_ID_list]
carbontube.drop(carbon_on_wall_for_all_functionals_list, inplace=True)
carbon_excluding_functionals = carbontube

for i in range(len(carbon_on_wall_for_all_functionals_list)):
    a = pd.Series.tolist(carbon_on_wall_for_all_functionals.iloc[i])
    into_cyl(a)

COH_cyl = pd.DataFrame(columns=["atom", "radius", "theta", "z"])
CH_cyl = pd.DataFrame(columns=["atom", "radius", "theta", "z"])

for i in range(len(carbon_on_wall_for_functional_COH20_ID_list)):
    carbon_row = into_cyl(carbon_on_wall_randomly_chosen_for_COH20.iloc[i])
    carbon_row.insert(0,"C")
    carbon_row_copy = carbon_row.copy()  
    COH_cyl.loc[len(COH_cyl)] = carbon_row
    oxy_row = oxy_bonding(carbon_row) # add line for oxygen below carbon
    oxy_row[0] = "O"
    COH_cyl.loc[len(COH_cyl)] = oxy_row    
    hyd_row = hyd_bonding(carbon_row_copy) # add line for hydrogen below oxygen
    hyd_row[0] = "H"
    COH_cyl.loc[len(COH_cyl)] = hyd_row

COH_cartes = pd.DataFrame(columns=["atom", "x", "y", "z"])
COH_cartes["atom"] = COH_cyl["atom"]
COH_cartes["x"] = COH_cyl["radius"] * np.cos(COH_cyl["theta"])
COH_cartes["y"] = COH_cyl["radius"] * np.sin(COH_cyl["theta"])
COH_cartes["z"] = COH_cyl["z"]

for i in range(len(carbon_on_wall_for_functional_CH30_ID_list)):
    carbon_row = into_cyl(carbon_on_wall_randomly_chosen_for_CH30.iloc[i])
    carbon_row.insert(0,"C")
    carbon_row_copy = carbon_row.copy()  
    CH_cyl.loc[len(CH_cyl)] = carbon_row 
    hyd_row = carbohyd_bonding(carbon_row_copy) # add line for hydrogen below carbon
    hyd_row[0] = "H"
    CH_cyl.loc[len(CH_cyl)] = hyd_row

CH_cartes = pd.DataFrame(columns=["atom", "x", "y", "z"])
CH_cartes["atom"] = CH_cyl["atom"]
CH_cartes["x"] = CH_cyl["radius"] * np.cos(CH_cyl["theta"])
CH_cartes["y"] = CH_cyl["radius"] * np.sin(CH_cyl["theta"])
CH_cartes["z"] = CH_cyl["z"]

del carbon_excluding_functionals["ID_number"]

df = pd.concat([carbon_excluding_functionals, COH_cartes, CH_cartes])
ID = np.arange(len(df))
df = df.set_index(ID)
                
z_convex_downward(df)

total_num_of_atoms = str(len(df))
dfAsString = df.to_string(header=False, index=False)
txtfile = open("COH20CH30convex.xyz","w")
txtfile.write(total_num_of_atoms)
txtfile.write("\n\n")
txtfile.write(dfAsString)
txtfile.close()