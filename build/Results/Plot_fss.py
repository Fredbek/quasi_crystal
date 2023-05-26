import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
import os
import sys
import csv

Lengths = ["239", "1393", "8119", "47321"]

scaling_factor = 1 + np.sqrt(2)

sc_exp = [3, 4, 5, 6]

colors = ['red', 'yellow', 'green', 'blue']

Data_Folder, ending = "Data/", ".csv"

Tc = 2.385

def main():

    k = 1.75
    alpha = -0.13
    beta = 0.6
    v = 1

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize = (16, 7))

    for i in range(len(Lengths)):

        Length = Lengths[i]

        L = scaling_factor**(sc_exp[i])

        t, m, s, c, Sm, Ss, Sc =  [], [], [], [], [], [], []

        Data_path = Data_Folder + Length + ending

        Arch = open(Data_path, "r")
        Data = list(csv.reader(Arch, delimiter = ','))
        Arch.close()

        for j in range(len(Data)):
            temperature = float(Data[j][0])/Tc - 1

            if np.abs(temperature*L**(1/v)) <= 0.5:
                t.append(temperature*L**(1/v))
                m.append(float(Data[j][3])*L**(-alpha/v))
                Sm.append(float(Data[j][4])*L**(-alpha/v))
                c.append(float(Data[j][5])*L**(-beta/v))
                Sc.append(float(Data[j][6])*L**(-beta/v))
                s.append(float(Data[j][7])*L**(-k/v))
                Ss.append(float(Data[j][8])*L**(-k/v))
        
        a1.scatter(t, m, label = Length, color = colors[i])
        a1.errorbar(t, m, Sm, color = colors[i])

        a2.scatter(t, c, label = Length, color = colors[i])
        a2.errorbar(t, c, Sc, color = colors[i])

        a3.scatter(t, s, label = Length, color = colors[i])
        a3.errorbar(t, s, Ss, color = colors[i])

    a1.set_xlabel = r"$t\cdot L^{1/\nu}$"
    a1.set_ylabel = r"$m\cdot L^{-\kappa/\nu}$"

    a2.set_xlabel = r"$t\cdot L^{1/\nu}$"
    a2.set_ylabel = r"$C_v\cdot L^{-\kappa/\nu}$"

    a3.set_xlabel = r"$t\cdot L^{1/\nu}$"
    a3.set_ylabel = r"$\chi\cdot L^{-\kappa/\nu}$"

    a1.set_xlim(-0.5, 0.5)
    a2.set_xlim(-0.5, 0.5)
    a3.set_xlim(-0.5, 0.5)

    a1.legend()
    a2.legend()
    a3.legend()
    
    plt.show()

main()





