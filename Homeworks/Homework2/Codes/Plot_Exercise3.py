import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def r_function(x):
    return (-(10*(x-0.5)**2))

def u_function(x):
    return (40*(x-0.5)**2-10)

def f_function(x):
    return 80 - 10*(40*(x - 0.5)**2-10)*(x - 0.5)**2


def Plot_rx():
	#Plot the real r(x)
	x = np.arange(-12,12,0.1)
	r = []
	for i in range(0,(x.size)):
    	r.append(r_function(x[i]))

	plt.title("r(x) Function")
	plt.plot(x, r,'blue')
	plt.ylabel('r(x) = -(10(x-0.5)^2)')
	plt.xlabel('x')
	plt.show()

def Plot_ux_real():
	#Plot the real u(x)
	x = np.arange(-100,100,0.1)
	u = []
	for i in range(0,(x.size)):
    	u.append(u_function(x[i]))
	plt.title("Real u(x) Function")
	plt.plot(x, u,'green')
	plt.ylabel('u(x) = (40*(x-0.5)^2)-10')
	plt.xlabel('x')
	plt.show()


def Plot_fx():
	#Plot the real f(x)
	x = np.arange(-100,100,0.1)
	f = []
	for i in range(0,(x.size)):
    	f.append(f_function(x[i]))
	plt.title("f(x) Function")
	plt.plot(x,f,'magenta')
	plt.ylabel('f(x) = 80 - 10*(40*((x - 0.5)^2)-10) *((x - 0.5)^2')
	plt.xlabel('x')
	plt.show()


def Plot_ux_uxreal():
	#We are going to plot the approximation:
	#Compare u(x) with u_approx(x):
	#1) Compute the u_approx(x):
	x = []
	u_approx = []
	for i in range(0,30):
    	name = str(i)+".txt"
    	file = open(name,"r")
    	for line in file: 
        	u_aux, x_aux = line.split(",")
        	u_approx.append(u_aux)
        	x.append(x_aux)

	#Compute the real u(x) using the same x(n).
	u = []
	x_aux = []
	for element in x:
    	e = float(element)
    	u.append(u_function(e))

	plt.title("Comparison between u(x) and u_approx(x)")
	plt.plot(x,u_approx,'magenta')
	plt.plot(x,u,'blue')
	plt.ylabel('Functions')
	plt.xlabel('x')
	plt.show()


def Plot_error():
	#Compute the mean error regarding the different iterations.
	iterations = np.arange(500000)
	file = open("error.txt","r")
	errors = []
	for line in file: 
    	errors.append(line)
    
	plt.title("Error obtained")
	plt.plot(iterations,errors,'red')
	plt.ylabel('Mean Squared Error')
	plt.xlabel('Iterations')
	plt.show()

def main():
	Plot_rx()
	Plot_fx()
	Plot_ux_real()
	Plot_ux_uxreal()
	Plot_error()


if __name__ == "__main__":
    main()

    

