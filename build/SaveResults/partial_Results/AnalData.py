import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import os

Lengths = [239, 1393, 8119, 47321]
exponents = [3, 4, 5, 6]
scal_factor = (1 + np.sqrt(2))


bt_fc_pm = []

path, ext = "Data/", ".csv"

OptimalParameters = [0, 0, 0]

fit_window = 0.5 #The window for the data adjust

degree = 6 #Degree of the polynomial

ac_err = 0.1

WithinErrorv, WithinErrork, WithinErrorTc = [], [], [] #The coefficients who fall inside the error

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
    
    return P[6] + P[5]*x + P[4]*x**2 + P[3]*x**3 + P[2]*x**4 + P[1]*x**5 + P[0]*x**6 
    # return P[1] + P[0]*x


def main():

    p = int(sys.argv[1]) # 3 for mag, 5 for cv, 7 for susc

    LeastChiSq = -1

    kspace = np.arange(float(sys.argv[1+1]), float(sys.argv[2+1]) + float(sys.argv[3+1]), float(sys.argv[3+1]))
    vspace = np.arange(float(sys.argv[4+1]), float(sys.argv[5+1]) + float(sys.argv[6+1]), float(sys.argv[6+1]))
    Tspace = np.arange(float(sys.argv[7+1]), float(sys.argv[8+1]) + float(sys.argv[9+1]), float(sys.argv[9+1]))

    for T in Tspace:
        for k in kspace:
            for v in vspace:

                if v == 0:
                    v = 1

            f_parameter, s_parameter = k, v

            X, Y ,SY = [], [], []

            for i in range(len(Lengths)):

                name = str(Lengths[i])
                file_path = path + name + ext

                Arch = open(file_path, 'r')
                Data = list(csv.reader(Arch, delimiter = ','))
                Arch.close()

                for j in range(len(Data)):
                    temperature = (float(Data[j][0])/T - 1)*scal_factor**(exponents[i]/s_parameter)

                    if abs(temperature) < fit_window:
                        X.append(temperature)
                        Y.append(float(Data[j][p])*(scal_factor**exponents[i])**(-f_parameter/s_parameter))
                        SY.append(float(Data[j][p])*(scal_factor**exponents[i])**(-f_parameter/s_parameter))

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
                    OptimalParameters = [k, v, T]
                    LeastChiSq = newChiSq
                    bt_fc_pm = parameters

    for T in Tspace:
        for k in kspace:
            for v in vspace:

                if v == 0:
                    v = 1

            f_parameter, s_parameter = k, v

            X, Y ,SY = [], [], []

            for i in range(len(Lengths)):

                name = str(Lengths[i])
                file_path = path + name + ext

                Arch = open(file_path, 'r')
                Data = list(csv.reader(Arch, delimiter = ','))
                Arch.close()

                for j in range(len(Data)):
                    temperature = (float(Data[j][0])/T - 1)*(scal_factor**exponents[i])**(1/s_parameter)

                    if abs(temperature) < fit_window:
                        X.append(temperature)
                        Y.append(float(Data[j][p])*(scal_factor**exponents[i])**(-f_parameter/s_parameter))
                        SY.append(float(Data[j][p])*(scal_factor**exponents[i])**(-f_parameter/s_parameter))

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
                
                #LeastChiSq = 1

                if newChiSq <= LeastChiSq*1.1 and newChiSq >= LeastChiSq*0.9:
                    WithinErrork.append(k)
                    WithinErrorv.append(v)
                    WithinErrorTc.append(T)


    Wek, Wev, WeT = np.asarray(WithinErrork), np.asarray(WithinErrorv), np.asarray(WithinErrorTc)

    print("For k:")
    print("Value: ", Media(Wek))
    print("Uncer: ", DesvPadMed(Wek))
    print("For v:")
    print("Value: ", Media(Wev))
    print("Uncer: ", DesvPadMed(Wev))
    print(r"For $T_c$:")
    print("Value: ", Media(WeT))
    print("Uncer: ", DesvPadMed(WeT))

    print("SQR: ", LeastChiSq)

main()
