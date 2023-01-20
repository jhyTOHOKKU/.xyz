import pandas as pd
import numpy as np
import random, math
from math import sqrt, pow

class AppMain:
    
    carbon_carbon_bond_length = 1.42
    z_coordinate_of_upper_plane = 73.9756
    
    def __init__(self):
        pass
    
    def run(self):    
        COH_ratio, CH_ratio = scan_functional_groups_ratio()
        plane, carbontube = read_tables()

        plane_cyl = into_cyl_df(plane)
        carbontube_cyl = into_cyl_df(carbontube)

        radius_max = calculate_radius_max(carbontube_cyl)   

        disc_cut_from_plane, plane_excluding_disc = separete_plane_and_disc(plane_cyl,radius_max)    

        carbon_to_COH_index, carbon_to_CH_index, carbon_no_bond_index = \
        ramdomly_separate_carbons_into_three_groups(carbontube, COH_ratio, CH_ratio)

        COH_cyl = pd.DataFrame(columns=["atom", "radius", "theta", "z"])
        CH_cyl = pd.DataFrame(columns=["atom", "radius", "theta", "z"])

        for i in carbon_to_COH_index:
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

        COH_cartes = into_cartes_df(COH_cyl)

        for i in carbon_to_CH_index:
            carbon_row = into_cyl(carbontube.loc[i])
            carbon_row.insert(0,"C")
            CH_cyl.loc[len(CH_cyl)] = carbon_row

            hyd_row = carbohyd_bonding(carbon_row)
            hyd_row[0] = "H"
            CH_cyl.loc[len(CH_cyl)] = hyd_row

        CH_cartes = into_cartes_df(CH_cyl)

        carbon_without_functionals = carbontube.loc[carbon_no_bond_index]

        df = pd.concat([disc_cut_from_plane, 
                        plane_excluding_disc, 
                        carbon_without_functionals, 
                        COH_cartes, 
                        CH_cartes])

        total_num_of_atoms = str(len(df))

        dfAsString = df.to_string(header=False, index=False)

        txtfile = open(f"nt-20-20-30-COH{COH_ratio*100}%_CH{CH_ratio*100}%.xyz","w")
        txtfile.write(f"{total_num_of_atoms}\n\n")
        txtfile.write(dfAsString)
        txtfile.close()

        print(f"new CNT = {len(carbon_without_functionals)}")
        print(f"new COH = {int(len(COH_cartes)/3)}")
        print(f"new CH = {int(len(CH_cartes)/2)}")
        print(f"radius_max = {radius_max}")
        print(f"depth = {df['z'].max() - df['z'].min()}")
    
    def scan_functional_groups_ratio():
        COH_ratio = input("Enter C-O-H fuctional group ratio. (between 0 and 1) : ") 
        CH_ratio = input("Enter C-H bond ratio. (between 0 and 1) : ")  
        return float(COH_ratio), float(CH_ratio)

    def read_tables():
        plane = pd.read_csv("gra.csv")
        carbontube = pd.read_csv("nt-20-20-30.csv")
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

    def ramdomly_separate_carbons_into_three_groups(carbontube, COH_ratio, CH_ratio):
        carbontube_index = list(carbontube.index)

        carbontube_index_shuffle = random.sample(carbontube_index,len(list(carbontube.index)))
        carbon_to_COH_index  = carbontube_index_shuffle[:int(len(list(carbontube_index_shuffle)) * COH_ratio)]
        carbon_to_CH_index   = carbontube_index_shuffle[ int(len(list(carbontube_index_shuffle)) * COH_ratio):
                                                         int(len(list(carbontube_index_shuffle)) * (COH_ratio+CH_ratio))]
        carbon_no_bond_index = carbontube_index_shuffle[ int(len(list(carbontube_index_shuffle)) * (COH_ratio+CH_ratio)):]

        carbon_to_COH_index.sort()
        carbon_to_CH_index.sort()
        carbon_to_CH_index.sort()

        return carbon_to_COH_index, carbon_to_CH_index, carbon_no_bond_index

    def oxy_bonding(cyl_cor):
        cyl_cor[1] = cyl_cor[1] - 1.428 # CO bond length = 1.428 Angstrom (オングストローム [6])
        return cyl_cor

    def hyd_bonding(cyl_cor):
        cyl_cor[1] = cyl_cor[1] - 2.389 # CO bond length + OH bond length = 2.289 Angstrom
        return cyl_cor

    def carbohyd_bonding(cyl_cor):
        cyl_cor[1] = cyl_cor[1] - 1.094 # CH bond length  = 1.094 Angstrom
        return cyl_cor

    def oxy_bonding_bottom(car_cor):
        car_cor[-1] = car_cor[-1] + 1.428 # CO bond length = 1.428 Angstrom
        return car_cor

    def hyd_bonding_bottom(car_cor):
        car_cor[-1] = car_cor[-1] + 2.389 # CO bond length + OH bond length = 2.289 Angstrom
        return car_cor

    def carbohyd_bonding_bottom(car_cor):
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

    def into_cyl_df(df):
        cyl = pd.DataFrame(columns = ["atom", "radius", "theta", "z"])
        for i in range(len(df)):
            row = df.loc[i]
            row = into_cyl(row)
            row.insert(0,df.loc[i][0])
            cyl.loc[i] = row
        return cyl
    
if __name__  == "__main__":
    app = AppMain()
    app.run()