import numpy as np
import matplotlib.pyplot as plt

#########################################

#reading numerical data

time=[]
rho11=[]
rho22=[]

for i in range(0,200):
	data = np.loadtxt("output/step_"+str(i)+".txt", unpack = True, skiprows=1)
	
	#$x $xd $rho11 $rho12 $rho22 $rho_bar11 $rho_bar12 $rho_bar22 $rho11_im $rho12_im $rho22_im $rho_bar11_im $rho_bar12_im $rho_bar22_im $N $N_bar $t
	# 0   1    2     3       4       5           6        7           8        9           10         11         12             13        14   15   16

	time.append(data[16][1])
	rho11.append(data[2][1])
	rho22.append(data[4][1])
	
#########################################
	
rho11_theory=[]
rho22_theory=[]

theta_12=39.0*np.pi/180

mass_1=0.01 #eV/c2
mass_2=0.03 #eV/c2
E=50.0e+6 #eV

hbar=1.054571817e-34 #J/s
eV_to_J=1.60217663e-19 

for t in time:
	rho11_theory.append(1-np.sin(2*theta_12)**2*np.sin(((mass_1**2-mass_2**2)/(4*E*hbar))*eV_to_J*t)**2)
	rho22_theory.append(np.sin(2*theta_12)**2*np.sin(((mass_1**2-mass_2**2)/(4*E*hbar))*eV_to_J*t)**2)
	
#########################################

plt.plot(time,rho11)	
plt.plot(time,rho11_theory,linestyle="dashed")		
plt.plot(time,rho22)	
plt.plot(time,rho22_theory,linestyle="dashed")		
plt.xlabel("Time (s)")
plt.ylabel("Probability")
plt.legend([r"$\rho_{11}^{num}$",r"$\rho_{11}^{theo}$",r"$\rho_{22}^{num}$",r"$\rho_{22}^{theo}$"])
plt.show()
plt.savefig("vacuum_test.pdf",bbox_inches='tight')

