import matplotlib.pyplot as plt
import csv
import numpy as np
import sys
from scipy.optimize import curve_fit
import os

h = float(sys.argv[4]) # Tipo a precisão

BuscaMin = float(sys.argv[2])
BuscaMax = float(sys.argv[3]) # Modulo do limite da busca

Lengths = [17, 18, 19, 20, 21, 22, 23, 24, 34, 36, 38, 40, 42, 44, 46, 48]

starting, ending = "Results/Data/DadosRede", ".csv"

SpaceSearch = np.arange(BuscaMin, BuscaMax + h/10, h)

vSpaceSearch = np.arange(1, 1.1 + 1/10, 0.1)

# SpaceSearch = [1.75, 1]

OptimalParameters = [0, 0]

bt_fc_pm = []

Tc = 2.268

fit_length = 1

degree = 4

WithinErrorv, WithinErrork = [], []

def Media(A):
    soma = 0
    for i in range(np.size(A)):
        soma += A[i]/np.size(A)

    return soma

def DesvPadMed(A):
    soma = 0
    Med = Media(A)
    
    for i in range(np.size(A)):
        soma += (A[i] - Med)**2/np.size(A)
        
    return np.sqrt(soma)

def Polin(x, P, deg):
    
    return P[4] + P[3]*x + P[2]*x**2 + P[1]*x**3 + P[0]*x**4
    # return P[1] + P[0]*x

def main():

    LeastChiSq = -1    

    Parameter = 7 #Default is the susceptibility
    if len(sys.argv) > 1:
        if(sys.argv[1] == "e"):
            Parameter = 1
        elif(sys.argv[1] == "m"):
            Parameter = 3
        elif(sys.argv[1] == "c"):
            Parameter = 5
        elif(sys.argv[1] == "s"):
            Parameter = 7
        elif(sys.argv[1] == "b"):
            Parameter = 9
            
            
    for v in vSpaceSearch:
        if v == 0:
            v = 1
        
        # v = 1
        
        for k in SpaceSearch:
            
            # k = 1.75
            
            f_parameter, s_parameter = k, v
                
            X, Y, SY = [], [], []
                
                
            for Length in Lengths:
                    
                meio = str(Length)
                Path = starting + meio + ending
                    
                Arch = open(Path, 'r')
                Dados = list(csv.reader(Arch, delimiter = ','))
                Arch.close()
                    
                for i in range(len(Dados)):
                    temperature = (float(Dados[i][0])/Tc - 1)*Length**(1/s_parameter)
                    if abs(temperature) <= fit_length:
                        X.append(temperature)
                        Y.append(float(Dados[i][Parameter])*Length**(-f_parameter/s_parameter))
                        SY.append(float(Dados[i][Parameter + 1])*Length**(-f_parameter/s_parameter))
                
            Xdata = np.asarray(X)
            Ydata = np.asarray(Y)
            SYdata = np.asarray(SY)
            
            if np.any(Ydata) == True:
                
                parameters = np.polyfit(Xdata, Ydata, degree)
                    
                newChiSq = 0
                    
                for i in range(len(Xdata)):
                    newChiSq += (Ydata[i] - Polin(Xdata[i], parameters, degree))**2/SYdata[i]**2
                    # newChiSq += (Ydata[i] - Gauss(i, parameters[0], parameters[1], parameters[2]))**2
                
                newChiSq = newChiSq/(np.size(Ydata) - degree - 1)
                
                if newChiSq < LeastChiSq or LeastChiSq == -1:
                    OptimalParameters[0] = k
                    OptimalParameters[1] = v
                    LeastChiSq = newChiSq
                    bt_fc_pm = parameters
                    
                
            # if LeastChiSq <= 1 and LeastChiSq != -1:
            #     break
    
    #-----------------------------------------------------------------------
    
    for v in vSpaceSearch:
        if v == 0:
            v = 1
        
        # v = 1
        
        for k in SpaceSearch:
            
            # k = 1.75
            
            f_parameter, s_parameter = k, v
                
            X, Y, SY = [], [], []
                
                
            for Length in Lengths:
                    
                meio = str(Length)
                Path = starting + meio + ending
                    
                Arch = open(Path, 'r')
                Dados = list(csv.reader(Arch, delimiter = ','))
                Arch.close()
                    
                for i in range(len(Dados)):
                    temperature = (float(Dados[i][0])/Tc - 1)*Length**(1/s_parameter)
                    if abs(temperature) <= fit_length:
                        X.append(temperature)
                        Y.append(float(Dados[i][Parameter])*Length**(-f_parameter/s_parameter))
                        SY.append(float(Dados[i][Parameter + 1])*Length**(-f_parameter/s_parameter))
                
            Xdata = np.asarray(X)
            Ydata = np.asarray(Y)
            SYdata = np.asarray(SY)
            
            if np.any(Ydata) == True:
                
                parameters = np.polyfit(Xdata, Ydata, degree)
                    
                newChiSq = 0
                    
                for i in range(len(Xdata)):
                    newChiSq += (Ydata[i] - Polin(Xdata[i], parameters, degree))**2/SYdata[i]**2
                    # newChiSq += (Ydata[i] - Gauss(i, parameters[0], parameters[1], parameters[2]))**2
                
                newChiSq = newChiSq/(np.size(Ydata) - degree - 1)
                
                if newChiSq <= LeastChiSq*1.1 or newChiSq >= LeastChiSq*0.9:
                    WithinErrork.append(k)
                    WithinErrorv.append(v)
    
    #------------------------------------------------------------------------
    
    WEk, WEv = np.asarray(WithinErrork), np.asarray(WithinErrorv)
    
    # f_parameter, s_parameter = OptimalParameters[0], OptimalParameters[1]
    
    f_parameter, s_parameter = Media(WEk), Media(WEv)
    sf_parameter, ss_parameter = DesvPadMed(WEk), DesvPadMed(WEv)
    
    print(WEk)
    
    print(f_parameter)
    print(sf_parameter)
    print(s_parameter)
    print(ss_parameter)
    
    #print(np.size(Ydata))
    
    for Length in Lengths:               
        meio = str(Length)
        Path = starting + meio + ending
                    
        Arch = open(Path, 'r')
        Dados = list(csv.reader(Arch, delimiter = ','))
        Arch.close()
                    
        for i in range(len(Dados)):
           
            X_plt, Y_plt = [], []
           
            temperature = (float(Dados[i][0])/Tc - 1)*Length**(1/s_parameter)
            if abs(temperature) <= fit_length:
                X.append(temperature)
                X_plt.append(temperature)
                Y.append(float(Dados[i][Parameter])*Length**(-f_parameter/s_parameter))
                Y_plt.append(float(Dados[i][Parameter])*Length**(-f_parameter/s_parameter))
    
    # ax.plot(X_plt, Y_plt, label = "Dados Rede: " + meio)
    
    
    Yplt = []
    
    for i in range(len(X)):
        Yplt.append(Polin(X[i], bt_fc_pm, degree))
    
    os.system("python3 PlotaBonito.py" + " " + sys.argv[1] + " " + str(OptimalParameters[0]) + " " + str(OptimalParameters[1]) + " &")
    
    fig, ax = plt.subplots()
    
    ax.scatter(X, Y, label = "Dados")
    ax.scatter(X, Yplt, label = "Função Chute")
    
    plt.legend()
    plt.show()
    
    
    
main()
