c=299792458 #m/s
c2=c*c 
eV_to_J=1.60217663e-19

import numpy as np
import matplotlib.pyplot as plt

#0. x 
#1. y 
#2. z 
#3. xd 
#4. yd 
#5. zd 
#6. rho11 
#7. rho12 
#8. rho13 
#9. rho22 
#10. rho23 
#11. rho33 
#12. rho_bar11 
#13. rho_bar12 
#14. rho_bar13
#15. rho_bar22 
#16. rho_bar23 
#17. rho_bar33
#18. rho11_im 
#19. rho12_im 
#20. rho13_im 
#21. rho22_im 
#22. rho23_im 
#23. rho33_im 
#24. rho_bar11_im 
#25. rho_bar12_im 
#26. rho_bar13_im 
#27. rho_bar22_im 
#28. rho_bar23_im 
#29. rho_bar33_im 
#30. N 
#31. N_bar 
#32. t


time=[]
rho11_av=[]


for i in range(0,10000):
	data = np.loadtxt("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/output/step_"+str(i)+".txt", unpack = True, skiprows=1)
	
	time.append(data[32][30])
	rho11_av.append(sum(data[6])/len(data[6]))
	
#constant

theta_12=1e-6*np.pi/180
theta_13=48.3*np.pi/180
theta_23=8.61*np.pi/180

E=50.0e6

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

print(phi32)
print(phi31)
print(phi21)

print(np.sin(2*theta_13)**2*c12**2)
print(np.sin(2*theta_13)**2*s12**2)
print((c13**4)*np.sin(2*theta_12)**2)

#analityc data

Pee=[]

print()
print(2*np.pi/1.224499993733192e-11)
print()


for t in time:
    Pee.append(1-np.sin(2*theta_13)**2*c12**2*np.sin(phi31*t)**2-np.sin(2*theta_13)**2*s12**2*(np.sin(phi32*t))**2-(c13**4)*np.sin(2*theta_12)**2*(np.sin(phi21*t))**2)

		

plt.plot(time,rho11_av)	
plt.plot(time,Pee)	
plt.show()



