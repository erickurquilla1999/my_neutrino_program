import numpy as np
import matplotlib.pyplot as plt

#########################################

#reading numerical data

time=[]
rho11=[]
rho22=[]
rho11_bar=[]
rho22_bar=[]

for i in range(0,200):
	data = np.loadtxt("output/step_"+str(i)+".txt", unpack = True, skiprows=1)
	
	#$x $xd $rho11 $rho12 $rho22 $rho_bar11 $rho_bar12 $rho_bar22 $rho11_im $rho12_im $rho22_im $rho_bar11_im $rho_bar12_im $rho_bar22_im $N $N_bar $t
	# 0   1    2     3       4       5           6        7           8        9           10         11         12             13        14   15   16

	time.append(data[16][1])
	rho11.append(data[2][1])
	rho22.append(data[4][1])
	rho11_bar.append(data[5][1])
	rho22_bar.append(data[7][1])
	
#########################################

hbar=1.054571817e-34 #J*s
eV_to_J=1.60217663e-19 # (J/eV)
c=299792458 #m/s
c2=89875517873681764 #(m/s)^2
c4=8077608713062490229263800746151696 #(m/s)^4
Gf=4.5437957e14 #J^(-2)

theta_12 = 39.0  * np.pi/180

mass_1=0.01 * eV_to_J / c2 # kg 
mass_2=0.03 * eV_to_J / c2 # kg 
mass_3=0.02 * eV_to_J / c2 # kg

E=50.0e+6 * eV_to_J # J

w=(mass_1**2-mass_2**2)*c4/(4*E*hbar)

background_electron_number_density=1e32 # m^(-3)
background_antielectron_number_density=0.0 # m^(-3)

VW=np.sqrt(2)*Gf*(hbar*c)**3*(background_electron_number_density-background_antielectron_number_density)
vw=2*E*VW/((mass_2**2-mass_1**2)*c4)
k=np.sqrt((np.cos(2*theta_12)-vw)**2+np.sin(2*theta_12)**2)
k_bar=np.sqrt((np.cos(2*theta_12)+vw)**2+np.sin(2*theta_12)**2)
sin2_2theta_matter=np.sin(2*theta_12)**2/k**2
sin2_2theta_matter_bar=np.sin(2*theta_12)**2/k_bar**2

rho11_theory=[]
rho22_theory=[]
rho11_bar_theory=[]
rho22_bar_theory=[]

for t in time:
	rho11_theory.append(1-sin2_2theta_matter*np.sin(w*k*t)**2)
	rho22_theory.append(sin2_2theta_matter*np.sin(w*k*t)**2)
	rho11_bar_theory.append(1-sin2_2theta_matter_bar*np.sin(w*k_bar*t)**2)
	rho22_bar_theory.append(sin2_2theta_matter_bar*np.sin(w*k_bar*t)**2)
	
#########################################

plt.plot(time,rho11)	
plt.plot(time,rho22)
plt.plot(time,rho11_bar)	
plt.plot(time,rho22_bar)	
plt.plot(time,rho11_theory,linestyle="dotted",color="black")		
plt.plot(time,rho22_theory,linestyle="dotted",color="black")
plt.plot(time,rho11_bar_theory,linestyle="dotted",color="black")		
plt.plot(time,rho22_bar_theory,linestyle="dotted",color="black")		
plt.xlabel("Time (s)")
plt.ylabel("Probability")
plt.legend([r"$\rho_{ee}$",r"$\rho_{xx}$",r"$\bar{\rho}_{ee}$",r"$\bar{\rho}_{xx}$","Analytic"])
plt.show()
plt.savefig("matter_vacuum_test.pdf",bbox_inches='tight')
