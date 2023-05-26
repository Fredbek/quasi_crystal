import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

#Lengths = [239]
Lengths = [239, 1393, 8119, 47321]
colors = ['red', 'blue', 'violet', 'green']

path, ext = 'Data/', '.csv'

def main():

    p = 7
    lbl = r'$\xi$' #default is susceptibility

    fig, ax = plt.subplots(figsize = (16, 7))

    if(len(sys.argv)) > 1:
        if sys.argv[1] == 'e':
            p = 1
            lbl = r'E/N'

        if sys.argv[1] == 'm':
            p = 3
            lbl = r'm'

        if sys.argv[1] == 'c':
            p = 5
            lbl = r'$c_V$'

        if sys.argv[1] == 'b':
            p = 9
            lbl = r'$B_4$'

    for i in range(len(Lengths)):
        file_path = path + str(Lengths[i]) + ext

        Arch = open(file_path, 'r')
        Data = list(csv.reader(Arch, delimiter = ','))
        Arch.close()

        X, Y, SY = [], [], []

        for j in range(len(Data)):
            X.append(float(Data[j][0]))
            Y.append(float(Data[j][p]))
            SY.append(float(Data[j][p+1]))

        ax.plot(X, Y, color = colors[i], label = 'N: ' + str(Lengths[i]))
        ax.errorbar(X, Y, SY, color = colors[i])

    ax.set_xlim(0, 5)

    ax.set_xlabel("T")
    ax.set_ylabel(lbl)



    ax.grid()
    ax.legend()


    #plt.yscale('log')
    plt.show()
        
main()

