# This program is for making bonds between carbons and -OH ion or -H ion on the upper plane and the nanotube, but letting the bottom plane be a bare carbon plane.  

import pandas as pd
import numpy as np
import random, math
from math import sqrt, pow

global carbon_carbon_bond_length
global z_coordinate_of_upper_plane

carbon_carbon_bond_length = 1.4200
z_coordinate_of_upper_plane = 73.9756

class Surface(object):
    def __init__(self,df):
        self.df = df
        
    def ramdomly_separate_carbons_into_three_groups(self):
        DF = self.df
        surface_index            = list(DF.index)
        surface_index_shuffle    =  random.sample(surface_index,len(surface_index))
        carbon_to_COH_index      =  surface_index_shuffle[:int(len(surface_index_shuffle) * COH_ratio)]
        carbon_to_CH_index       =  surface_index_shuffle[ int(len(surface_index_shuffle) * COH_ratio):
                                                        int(len(surface_index_shuffle) * (COH_ratio+CH_ratio))]
        carbon_to_no_bond_index  =  surface_index_shuffle[ int(len(surface_index_shuffle) * (COH_ratio+CH_ratio)):]

        carbon_to_COH_index.sort()
        carbon_to_CH_index.sort()
        carbon_to_no_bond_index.sort()

        return carbon_to_COH_index, carbon_to_CH_index, carbon_to_no_bond_index
    
    
def scan_functional_groups_ratio():
    global COH_ratio
    global CH_ratio
    COH_ratio = float( input("Enter C-O-H fuctional group ratio. (between 0 and 1) : ") )
    CH_ratio = float ( input("Enter C-H bond ratio. (between 0 and 1) : ")) 
    return COH_ratio, CH_ratio

def read_tables():
    plane = pd.read_csv("gra.csv")
    carbontube = pd.read_csv("nt-25-25-30.csv")
    return plane, carbontube

def calculate_radius_max(df_cyl):
    radius_max = df_cyl["radius"].max()
    return radius_max

def separete_plane_and_disc(plane_cyl,radius_max):     
    bool_not_greater_than_radius_max = plane_cyl["radius"] <= radius_max
    
    disc_cut_from_plane = plane_cyl.loc[bool_not_greater_than_radius_max] 
    plane_excluding_disc = plane_cyl.loc[~ bool_not_greater_than_radius_max]
    
    disc_cut_from_plane.loc[:,'z'] = -1. * carbon_carbon_bond_length         
    plane_excluding_disc.loc[:,'z'] = z_coordinate_of_upper_plane
    
    disc_cut_from_plane  = into_cartes_df(disc_cut_from_plane)
    plane_excluding_disc = into_cartes_df(plane_excluding_disc)
    
    return disc_cut_from_plane, plane_excluding_disc

def oxy_bonding(cyl_cor):
    cyl_cor[1] = cyl_cor[1] - 1.428 # CO bond length = 1.428 Angstrom (オングストローム [6])
    return cyl_cor

def hyd_bonding(cyl_cor):
    cyl_cor[1] = cyl_cor[1] - 2.389 # CO bond length + OH bond length = 2.289 Angstrom
    return cyl_cor

def carbohyd_bonding(cyl_cor):
    cyl_cor[1] = cyl_cor[1] - 1.094 # CH bond length  = 1.094 Angstrom
    return cyl_cor

def oxy_bonding_plane(car_cor):
    car_cor[-1] = car_cor[-1] + 1.428 # CO bond length = 1.428 Angstrom
    return car_cor

def hyd_bonding_plane(car_cor):
    car_cor[-1] = car_cor[-1] + 2.389 # CO bond length + OH bond length = 2.289 Angstrom
    return car_cor

def carbohyd_bonding_plane(car_cor):
    car_cor[-1] = car_cor[-1] + 1.094 # CH bond length  = 1.094 Angstrom
    return car_cor

def into_cartes(cyl_cor):
    car_cor = ["x","y","z"]
    car_cor[-3] = cyl_cor[-3] * np.cos(cyl_cor[-2]) #x
    car_cor[-2] = cyl_cor[-3] * np.sin(cyl_cor[-2]) #y
    car_cor[-1] = cyl_cor[-1]                       #z  
    return car_cor
        # cyl_cor = [ ......  radius,     theta ,    z]
    
def into_cartes_df(cyl_coord):
    cartes_coord = pd.DataFrame(columns=["atom", "x", "y", "z"])
    cartes_coord["atom"] = cyl_coord["atom"]
    cartes_coord["x"] = cyl_coord["radius"] * np.cos(cyl_coord["theta"])
    cartes_coord["y"] = cyl_coord["radius"] * np.sin(cyl_coord["theta"])
    cartes_coord["z"] = cyl_coord["z"]
    return cartes_coord
    
def into_cyl(a): # change rectangular coordinate into cylindrical coordinate
    r = sqrt(pow(a[-3], 2.) + pow(a[-2], 2.))
    if a[-3] == 0.:  # if x = 0 
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

def into_cyl_df(df):
    cyl = pd.DataFrame(columns = ["atom", "radius", "theta", "z"])
    for i in range(len(df)):
        row = df.loc[i]
        row = into_cyl(row)
        row.insert(0,df.loc[i][0])
        cyl.loc[i] = row
    return cyl

def main():    
    COH_ratio, CH_ratio = scan_functional_groups_ratio()
    plane, carbontube = read_tables()
  
    plane_cyl = into_cyl_df(plane)
    carbontube_cyl = into_cyl_df(carbontube)
    
    radius_max = calculate_radius_max(carbontube_cyl)   
    
    disc_cut_from_plane, plane_excluding_disc = separete_plane_and_disc(plane_cyl,radius_max)    
                                             
    _carbontube = Surface(carbontube)
    _disc_cut_from_plane = Surface(disc_cut_from_plane)
    _plane_excluding_disc = Surface(plane_excluding_disc)
        
    tube_to_COH_index, tube_to_CH_index, tube_no_bond_index =     _carbontube.ramdomly_separate_carbons_into_three_groups()

    disc_to_COH_index, disc_to_CH_index, disc_no_bond_index =     _disc_cut_from_plane.ramdomly_separate_carbons_into_three_groups()
        
    upper_plane_to_COH_index, upper_plane_to_CH_index, upper_plane_no_bond_index =     _plane_excluding_disc.ramdomly_separate_carbons_into_three_groups()
    
    
    COH_cyl = pd.DataFrame(columns=["atom", "radius", "theta", "z"])
    CH_cyl  = pd.DataFrame(columns=["atom", "radius", "theta", "z"])
    
    COH_disc = pd.DataFrame(columns=["atom", "x", "y", "z"])
    CH_disc  = pd.DataFrame(columns=["atom", "x", "y", "z"])
    
    COH_plane_excluding_disc = pd.DataFrame(columns=["atom", "x", "y", "z"])
    CH_plane_excluding_disc  = pd.DataFrame(columns=["atom", "x", "y", "z"])
       
    ################# tube  ###################
    
    for i in tube_to_COH_index:
        carbon_row = into_cyl(carbontube.loc[i])
        carbon_row.insert(0,"C")
        carbon_row_copy = carbon_row.copy()
        COH_cyl.loc[len(COH_cyl)] = carbon_row
              
        oxy_row = oxy_bonding(carbon_row) 
        oxy_row[0] = "O"
        COH_cyl.loc[len(COH_cyl)] = oxy_row
        
        hyd_row = hyd_bonding(carbon_row_copy)
        hyd_row[0] = "H"
        COH_cyl.loc[len(COH_cyl)] = hyd_row
                
    COH_cartes_tube = into_cartes_df(COH_cyl)
  

    for i in tube_to_CH_index:
        carbon_row = into_cyl(carbontube.loc[i])
        carbon_row.insert(0,"C")
        CH_cyl.loc[len(CH_cyl)] = carbon_row
        
        hyd_row = carbohyd_bonding(carbon_row)
        hyd_row[0] = "H"
        CH_cyl.loc[len(CH_cyl)] = hyd_row
    
    CH_cartes_tube = into_cartes_df(CH_cyl)
       
    
 ################# disc #############################   
    
#     for i in disc_to_COH_index:
#         carbon_row = disc_cut_from_plane.loc[i]
#         carbon_row_copied = carbon_row.copy()
#         COH_disc.loc[len(COH_disc)] = carbon_row
# #         carbon_row_copy = carbon_row.copy()
        
#         oxy_row = oxy_bonding_plane(carbon_row) 
#         oxy_row[0] = "O"
#         COH_disc.loc[len(COH_disc)] = oxy_row
        
#         hyd_row = hyd_bonding_plane(carbon_row_copied)
#         hyd_row[0] = "H"
#         COH_disc.loc[len(COH_disc)] = hyd_row


#     for i in disc_to_CH_index:
#         carbon_row = disc_cut_from_plane.loc[i]
#         CH_disc.loc[len(CH_disc)] = carbon_row

#         hyd_row = carbohyd_bonding_plane(carbon_row)
#         hyd_row[0] = "H"
#         CH_disc.loc[len(CH_disc)] = hyd_row

 ##################### upper plane #########################   
    
    for i in upper_plane_to_COH_index:
        carbon_row = plane_excluding_disc.loc[i]
        carbon_row_copied = carbon_row.copy()
        COH_plane_excluding_disc.loc[len(COH_plane_excluding_disc)] = carbon_row
        
        oxy_row = oxy_bonding_plane(carbon_row) 
        oxy_row[0] = "O"
        COH_plane_excluding_disc.loc[len(COH_plane_excluding_disc)] = oxy_row
        
        hyd_row = hyd_bonding_plane(carbon_row_copied)
        hyd_row[0] = "H"
        COH_plane_excluding_disc.loc[len(COH_plane_excluding_disc)] = hyd_row  
    
    for i in upper_plane_to_CH_index:
        carbon_row = plane_excluding_disc.loc[i]
        CH_plane_excluding_disc.loc[len(CH_plane_excluding_disc)] = carbon_row

        hyd_row = carbohyd_bonding_plane(carbon_row)
        hyd_row[0] = "H"
        CH_plane_excluding_disc.loc[len(CH_plane_excluding_disc)] = hyd_row
        
########################## no functionals ####################################     
    tube_without_functionals = carbontube.loc[tube_no_bond_index]
    disc_without_functionals = disc_cut_from_plane.loc[disc_no_bond_index]
    upper_plane_without_functionals = plane_excluding_disc.loc[upper_plane_no_bond_index]


##############################################
    df = pd.concat([upper_plane_without_functionals,
                    tube_without_functionals,
                    disc_cut_from_plane,
                    COH_plane_excluding_disc,
                    COH_cartes_tube,
#                     COH_disc,
                    CH_plane_excluding_disc,
                    CH_cartes_tube])
#                     CH_disc])

    total_num_of_atoms = str(len(df))

    dfAsString = df.to_string(header=False, index=False)

    txtfile = open(f"nt-25-25-30-COH{COH_ratio*100}%_CH{CH_ratio*100}%.xyz","w")
    txtfile.write(f"{total_num_of_atoms}\n\n")
    txtfile.write(dfAsString)
    txtfile.close()
    
    print(f"new CNT = {len(upper_plane_without_functionals)+len(tube_without_functionals)+len(disc_cut_from_plane)}")
    print(f"new COH = {(len(COH_plane_excluding_disc)+len(COH_cartes_tube))//3}")
    print(f"new CH = {(len(CH_plane_excluding_disc)+len(CH_cartes_tube))//2}")
    print(f"radius_max = {radius_max}")
    print(f"depth = {df['z'].max() - df['z'].min()}")
    
if __name__  == "__main__":
    main()

