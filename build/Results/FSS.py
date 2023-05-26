import numpy as np
from scipy.optimize import curve_fit
import csv
import sys 
import os

# We'll have to make it three times so, one for each parameter

Parameters = [7, 3, 5]

Lengths = ["239", "1393", "8119", "47321"]

#Lengths = ["239", "1393", "8119"]

scaling_factor = 1 + np.sqrt(2)

sc_exp = [3, 4, 5, 6]

Data_Folder, ending = "Data/", ".csv"

Tc = 2.388

fit_window = 0.5

better_parameters = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

optimal_parameters = [[], [], [], []]

h = 0.01
Search_Win = [[0.9, 1.1], [1.7, 1.8], [-0.2, -0.1], [-0.5, 0.5]] #v, k, a, b

v_space = np.arange(Search_Win[0][0], Search_Win[0][1] + h, h)

def Media(array):

    media = 0

    for i in range(len(array)):
        media += array[i]/len(array)

    return media

def DesvPadMed(array):

    media = Media(array)

    desvpad = 0

    for i in range(len(array)):
        desvpad += (media - array[i])**2/len(array)

    return np.sqrt(desvpad)

def Polin(x, A, B, C, D, E, F, G):
    return A + x*B + C*x**2 + D*x**3 + E*x**4 + F*x**5 + G*x**6

def Polin2(x, P):
    return P[0] + x*P[1] + P[2]*x**2 + P[3]*x**3 + P[4]*x**4 + P[5]*x**5 + P[6]*x**6

def chi_sq(x_data, y_data, sy_data):

    chi = 0

    P, covariance = curve_fit(Polin, x_data, y_data, sigma = sy_data)

    print(len(y_data))

    if len(y_data) == 0:
        print("Não há dados dentro da região")
        return 10000
    
    else:
        for m in range(len(x_data)):

            chi += (y_data[m] - Polin2(x_data[m], P))**2/sy_data[m]**2
        #print(chi)
        chi = chi/(len(y_data) - len(P) - 1)
        return np.sqrt(chi)#/(len(x_data) - 7)




def main():
    n = 0

    while n < 3:

        ex_space = np.arange(Search_Win[n+1][0], Search_Win[n+1][1] + h, h)

        for v in v_space:
            if v == 0:
                v == 1
            
            for ex in ex_space:
                
                x, y, sy = [], [], []

                for i in range(len(Lengths)):

                    Length = Lengths[i]

                    L = scaling_factor**sc_exp[i]

                    Data_Path = Data_Folder + Length + ending

                    Arch = open(Data_Path, 'r')
                    Data = list(csv.reader(Arch, delimiter = ','))
                    Arch.close()

                    for j in range(len(Data)):

                        temperature = float(Data[j][0])/Tc - 1

                        if np.abs(temperature*L**(1/v)) <= fit_window:
                            x.append(temperature*L**(1/v))
                            y.append(float(Data[j][Parameters[n]])*L**(-ex/v))
                            sy.append(float(Data[j][Parameters[n]+1])*L**(-ex/v))

                X = np.asarray(x)
                Y = np.asarray(y)
                SY = np.asarray(sy)

                chi = chi_sq(X, Y, SY)

                if chi >= 0.9 and chi <= 1.1:
                    optimal_parameters[0].append(v)
                    optimal_parameters[n+1].append(ex)
        print("end")
        n += 1

    print("Resultados:")
    print("v = ", Media(optimal_parameters[0]), " +- ", DesvPadMed(optimal_parameters[0]))
    print("k = ", Media(optimal_parameters[1]), " +- ", DesvPadMed(optimal_parameters[1]))
    print("beta = ", Media(optimal_parameters[2]), " +- ", DesvPadMed(optimal_parameters[2]))
    print("alpha = ", Media(optimal_parameters[3]), " +- ", DesvPadMed(optimal_parameters[3]))
    
main()
    
