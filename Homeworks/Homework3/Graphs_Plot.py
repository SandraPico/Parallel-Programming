#Python code used to generate all the graphs in the report.
#Authors: Sandra Picó and Miquel Larsson


%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def PlotSortedData(name):
    sorted_data = pd.read_csv(name, sep="/", header=None)
    np.asarray(sorted_data)
    time = np.arange(sorted_data.size)
    plt.scatter(time,sorted_data,s = 4,c = 'blue')
    plt.title("Sorted data (P =3)")
    plt.xlabel('N = 10')
    plt.show()

def PlotCPUTime():
    #3.2 Mesure efficiency of your program
    #N = 10^7
    #P = 1,2,3,4,8
    Process = np.array([1,2,3,4,8])
    #Values achieved in the c code
    cpu_time = np.array([6.837924,3.505298,3.326780,3.225967,2.036636]) 

    #Diagrama de barras
    plt.axes(((0.1, 0.3, 0.8, 0.6))) 
    plt.bar(np.arange(5), cpu_time) 
    plt.ylim(0,7)
    P = ["P = 1","P = 2","P = 3","P = 4","P = 8"]
    plt.title('CPU Time')
    plt.xlabel('Number of processes')
    plt.ylabel('Cpu time (seconds)')
    plt.xticks(np.arange(5), P, rotation = 45)

def PlotSpeedUp():
    cpu_time = np.array([3.505298,3.326780,3.225967,2.036636]) 
    time_p1 = 6.837924
    speed_up = np.array([1.95073971,2.05541815,2.11965094,3.35746005])
    P = np.array([2,3,4,8])
    plt.plot(P,speed_up,c = 'green')
    plt.title("Speed-Up")
    plt.xlabel('P values')
    plt.ylabel('cpu_time(1)/cpu_time(p)')
    plt.show()


if __name__ == "__main__":
    print("Start plotting")
    names_fixedP = ["Result_P3_N10.txt","Result_P3_N100.txt","Result_P3_N1000.txt"]
    names_fixedN = ["Result_P1_N1000.txt","Result_P10_N1000.txt","Result_P20_N1000.txt"]
    #Plot graphs with fixed P
    for i in range(0,3):
        PlotSortedData(names_fixedP[i])
    #Plot graphs with fixed N
    for i in range(0,3):
        PlotSortedData(names_fixedN[i])
    #Plot CPU time
    PlotCPUTime()
    #Plot Speed-Up 
    PlotSpeedUp()

