import numpy as np
import matplotlib.pyplot as plt

#0. x 
#1. x_dir 
#2. rho11 
#3. rho12 
#4. rho13 
#5. rho22 
#6. rho23 
#7. rho33 
#8. rho_bar11 
#9. rho_bar12 
#10. rho_bar13
#11. rho_bar22 
#12. rho_bar23 
#13. rho_bar33
#14. rho11_im 
#15. rho12_im 
#16. rho13_im 
#17. rho22_im 
#18. rho23_im 
#19. rho33_im 
#20. rho_bar11_im 
#21. rho_bar12_im 
#22. rho_bar13_im 
#23. rho_bar22_im 
#24. rho_bar23_im 
#25. rho_bar33_im 
#26. N 
#27. N_bar 
#28. t


time=[]
rho11_av=[]


for i in range(0,10000):
	data = np.loadtxt("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/output/step_"+str(i)+".txt", unpack = True, skiprows=1)
	
	time.append(data[28][30])
	rho11_av.append(data[2][30])
	
#constant

c=299792458 #m/s
c2=c*c 
eV_to_J=1.60217663e-19 

theta_12=1e-6*np.pi/180
theta_13=48.3*np.pi/180
theta_23=8.61*np.pi/180

E=50.0e+6

m1=0.0
m2=0.0
m3=0.049487372

s12=np.sin(theta_12)
c12=np.cos(theta_12)
s13=np.sin(theta_13)
c13=np.cos(theta_13)
s23=np.sin(theta_23)
c23=np.cos(theta_23)

phi32=(m3**2-m2**2)*eV_to_J*c**4/(4*E)
phi31=(m3**2-m1**2)*eV_to_J*c**4/(4*E)
phi21=(m2**2-m1**2)*eV_to_J*c**4/(4*E)


#analityc data

Pee=[]

for t in time:
    Pee.append(1-np.sin(2*theta_13)**2*c12**2*np.sin(phi31*t)**2-np.sin(2*theta_13)**2*s12**2*(np.sin(phi32*t))**2-(c13**4)*np.sin(2*theta_12)**2*(np.sin(phi21*t))**2)

#plots		

plt.plot(time,rho11_av)	
plt.plot(time,Pee)	
plt.show()



