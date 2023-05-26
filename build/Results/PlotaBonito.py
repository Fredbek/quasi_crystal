import sys
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.optimize import curve_fit

Lengths = [239, 1393, 8119, 47321]
WS = [239, 1393, 8119, 47321]

scal_factor = 1 + np.sqrt(2)
coef = [3, 4, 5, 6]

markers = ['v', 'v', 'o', 'o']
colors = ['black', 'red', 'black', 'red']

starting, ending = "Data/", ".csv"

Tc = 2.387

limin = 0.5

temp = np.arange(-limin, limin + 1/100, 1/100)

def Polin1(x, a, b, c, d, f, g,h):
    return a + b*x + c*x**2 + d*x**3 + f*x**4 + g*x**5 + h*x**6

def Polin(x, P):
    
    return P[0] + P[1]*x + P[2]*x**2 + P[3]*x**3 + P[4]*x**4 + P[5]*x**5 + P[6]*x**6

def main():
    
    chi = 0
    
    if len(sys.argv) > 2:
        f_param, s_param = float(sys.argv[2]), float(sys.argv[3])
    else:
        f_param, s_param = 1.75, 1
        
    
    Parameter = 7 #Default is the susceptibility
    exponent = "$\gamma$: "
    
    if len(sys.argv) > 1:
        if(sys.argv[1] == "m"):
            Parameter = 3
            exponent = "$\beta:$ -"    
        elif(sys.argv[1] == "c"):
            Parameter = 5
            exponent = "$\alpha$: "

    X, Y, Yer = [], [], []
    messages = []
    
    
    for i in range(len(Lengths)):
        Length = Lengths[i]
        meio = str(Length)
        Datapath = starting + meio + ending
        
        L = scal_factor**coef[i]

        Arch = open(Datapath, "r")
        data = list(csv.reader(Arch, delimiter = ','))
        Arch.close()
        
        for i in range(len(data)):
            temperature = (float(data[i][0])/Tc - 1)*L**(1/s_param)
            if abs(temperature) <= limin:
                X.append(temperature)
                Y.append(float(data[i][Parameter])*L**(-f_param/s_param))
                Yer.append(float(data[i][Parameter + 1])*L**(-f_param/s_param))
            
    Xdt, Ydt, Ydter = np.asarray(X), np.asarray(Y), np.asarray(Yer)
    
    Param, Covariance = curve_fit(Polin1, Xdt, Ydt, sigma = Ydter)
    
    for i in range(len(X)):
        chi += ((Y[i] - Polin(X[i], Param))**2)/Yer[i]**2
    
    fig, ax = plt.subplots(1, 2, figsize = (16, 7))
    
    for j in range(len(WS)):
        Length = WS[j]
        meio = str(WS[j])
        Datapath = starting + meio + ending
        
        L = scal_factor**coef[j]

        Arch = open(Datapath, "r")
        data = list(csv.reader(Arch, delimiter = ','))
        Arch.close()
        
        X, Y, Yer = [], [], []
        
        for i in range(len(data)):
            temperature = (float(data[i][0])/Tc - 1)*L**(1/s_param)
            if abs(temperature) <= limin:
                X.append(temperature)
                Y.append(float(data[i][Parameter])*L**(-f_param/s_param))
                Yer.append(float(data[i][Parameter + 1])*L**(-f_param/s_param))
        
        ax[0].scatter(X, Y, color = colors[j], marker = markers[j], label = 'Length: ' + meio)
        ax[0].errorbar(X, Y, Yer, color = colors[j], ls = 'none')
        
        Ys = []
        
        for k in range(len(X)):
            y = Y[k] - Polin(X[k], Param)
            Ys.append(y)
            
        ax[1].scatter(X, Ys, color = colors[j], marker = markers[j])
        ax[1].errorbar(X, Ys, Yer, color = colors[j], ls = 'none')
        
    P = []
    one = []
    
    for i in range(len(temp)):
        P.append(Polin(temp[i], Param))
        one.append(float(0))
        
    ax[0].plot(temp, P, color = 'black', linestyle = 'dashed', label = 'Fitted function')
    ax[1].plot(temp, one, color = 'red')
    
    ax[1].plot(messages, messages, label = exponent + str(f_param))
    ax[1].plot(messages, messages, label = r"$\nu: $" + str(s_param))
    
    ax[0].set_xlabel(r"$L^{1/\nu}t$", size = 20)
    ax[1].set_xlabel(r"$L^{1/\nu}t$", size = 20)
    
    ax[0].set_ylabel(r"$\chi L^{-\gamma/\nu}$", size = 20)
    ax[1].set_ylabel(r"$\chi L^{-\gamma/\nu} - g(x)$", size = 15)
    
    ax[0].legend()
    ax[1].legend()
    
    print(chi/(np.size(Xdt) - 5))
    
    plt.show()
    
main()