from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import csv

Lengths = [239, 1393, 8119, 47321]

# Lengths = [17, 34]

duplas = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]

# duplas = [[0, 1]]
# duplas = [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [10, 11], [12, 13], [14, 15]]

startingBnd, ending ="Data/", ".csv"

range = [2.3, 2.4]

Pms1, Pms2, Covariancia = [], [], []

resultados = []

def DesvPadMed(data, cont):
    
    media = 0
    i = 0
    
    while i < cont:
        media += data[i]
        i += 1
    media = media/cont
    
    i = 0
    desvpad = 0
    
    while i < cont:
        desvpad += (data[i] - media)**2
        i += 1
        
    desvpad = np.sqrt(desvpad/cont)
    
    return desvpad

def Polin(x, a, b, c, d, e, f, g):
    return a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5 + g*x**6

def PolinP(x, P):
    return P[6] + P[5]*x + P[4]*x**2 + P[3]*x**3 + P[2]*x**4 + P[1]*x**5 + P[0]*x**6

def dPolinP(x, P):
    return P[5] + 2*P[4]*x + 3*P[3]*x**2 + 4*P[2]*x**3 + 5*P[1]*x**4 + 6*P[0]*x**5

def difPolin(x, P1, P2):
    return PolinP(x, P1) - PolinP(x, P2)

def g(x, P1, P2):
    return x - difPolin(x, P1, P2)/(dPolinP(x, P1) - dPolinP(x, P2))

def main():
    
    fig, ax = plt.subplots(figsize = (16, 7))
    
    for Length in Lengths:
        meio = str(Length)

        binder_path = startingBnd + meio + ending
        
        Arch_bd = open(binder_path, "r")
        
        bd_data = list(csv.reader(Arch_bd, delimiter = ','))
        
        Arch_bd.close()
        
        X_bd, Y_bd, Yer_bd = [], [], []
        
        for i in bd_data:
            if float(i[0]) >= range[0] and float(i[0]) <= range[1]:
                X_bd.append(float(i[0]))
                Y_bd.append(float(i[9]))
                Yer_bd.append(float(i[10]))
                
        ax.scatter(X_bd, Y_bd, label = "Length: " + meio)
        ax.errorbar(X_bd, Y_bd, Yer_bd)
        
        Xbd = np.asarray(X_bd)
        Ybd = np.asarray(Y_bd)
        
        Param1 = np.polyfit(Xbd, Ybd, 6)
        
        Pms1.append(Param1)
        # Covariancia.append(Cov)
    
    
    for dp in duplas:
        
        contador = 0
        
        a, b, p = range[0], range[1], 0
        
        while contador < 10:
            p = (a+b)/2
            
            if difPolin(a, Pms1[dp[0]], Pms1[dp[1]])*difPolin(p, Pms1[dp[0]], Pms1[dp[1]]) < 0:
                b = p
            else:
                a = p
            
            contador += 1
        
        pbd = p
        
        while difPolin(pbd, Pms1[dp[0]], Pms1[dp[1]]) >= 0.0000001:
            pbd = g(pbd, Pms1[dp[0]], Pms1[dp[1]])
        
        resultados.append(pbd)
        
    soma, ctr = 0, 0
    
    for i in resultados:
        soma += i
        ctr += 1
        
    Tc = soma/ctr
    
    print(resultados)
    print(Tc)
    print(DesvPadMed(resultados, ctr))
    
    plt.legend()
    plt.show()
    
main()
