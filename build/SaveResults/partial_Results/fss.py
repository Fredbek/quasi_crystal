import matplotlib.pyplot as plt
import csv
import numpy as np
import sys

def f(x, a, b, c):
    return a + b*x + c*x**2

Lengths = [239, 1393, 8119, 47321]

starting, ending = "Data/", ".csv"

f_parameter, s_parameter, Tc = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) #The quantity dependent critical exponents and the v critical exponent


def main():
    
    fig, ax = plt.subplots()
    
    Parameter = 7 #Default is the susceptibility
    

    if(sys.argv[1] == "m"):
        Parameter = 3
    elif(sys.argv[1] == "c"):
        Parameter = 5
    elif(sys.argv[1] == "s"):
        Parameter = 7

    
    for i in Lengths:
        
        meio = str(i)
        NetLenght = i
        
        lbl = "Length: " + meio
        
        Path = starting + meio + ending
        
        data = open(Path, 'r')
        dados = list(csv.reader(data, delimiter = ','))
        data.close()
        
        T, Dt = [], []
        
        for i in range(len(dados)):
            temperature = float(dados[i][0])/Tc - 1
            
            if abs(NetLenght**(1/s_parameter)*temperature) <= 2:
                T.append(NetLenght**(1/s_parameter)*temperature)
                Dt.append(NetLenght**(-f_parameter/s_parameter)*float(dados[i][Parameter]))
            
        ax.scatter(T, Dt, label = lbl)
        ax.set_xlabel(r"$L^{1/\nu}t$", size = 10)
        ax.set_ylabel(r"$\chi L^{-\gamma/\nu}$", size = 10)
    plt.show()
        
        
    
    
    
    
main()
