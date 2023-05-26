import matplotlib.pyplot as plt
import csv
import sys

dirpath, ending = "Data/", ".csv"


def main():
    
    NetLenght = sys.argv[1]
    
    filepath = dirpath + NetLenght + ending 
    
    dt8 = open(filepath, "r")
    dados8 = list(csv.reader(dt8, delimiter = ","))
    dt8.close()
    
    T, X8, sX8 = [], [], []
    
    Parameter = 1 ##Default is the energy
    
    if(sys.argv[2] == "e"):
        Parameter = 1
    elif(sys.argv[2] == "m"):
        Parameter = 3
    elif(sys.argv[2] == "c"):
        Parameter = 5
    elif(sys.argv[2] == "s"):
        Parameter = 7
    elif(sys.argv[2] == "b"):
        Parameter = 9
    
    for i in range(len(dados8)):
        T.append(float(dados8[i][0]))
        X8.append(float(dados8[i][Parameter]))
        sX8.append(float(dados8[i][Parameter + 1]))
        
    fig, ax = plt.subplots()
    

    ax.plot(T, X8, color = "black", label = "L = 32", linestyle = "dashed")
    ax.scatter(T, X8, color = "black")
    ax.errorbar(T, X8, sX8, color = "black")
    
    plt.show()
    
main()
