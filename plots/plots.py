

import numpy as np
import matplotlib.pyplot as plt

#0. x 
#1. xd 
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
rho22_av=[]
rho33_av=[]
rho12_av=[]
rho13_av=[]
rho23_av=[]

rho11=[]
rho22=[]

for i in range(0,1000):
	data = np.loadtxt("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/output/step_"+str(i)+".txt", unpack = True, skiprows=1)
	
	time.append(data[28][30])

	rho11_av.append(sum(data[2])/len(data[2]))
	rho22_av.append(sum(data[5])/len(data[5]))
	rho33_av.append(sum(data[7])/len(data[7]))
	
	rho12_amp=[]
	for i in range(0,len(data[3])):
		rho12_amp.append((data[3][i]**2+data[15][i]**2)**0.5)	
	rho12_av.append(sum(rho12_amp)/len(rho12_amp))	
		
	rho13_amp=[]
	for i in range(0,len(data[4])):
		rho13_amp.append((data[4][i]**2+data[16][i]**2)**0.5)	
	rho13_av.append(sum(rho13_amp)/len(rho13_amp))	

	rho23_amp=[]
	for i in range(0,len(data[6])):
		rho23_amp.append((data[6][i]**2+data[18][i]**2)**0.5)	
	rho23_av.append(sum(rho23_amp)/len(rho23_amp))	

plt.plot(time,rho11_av)	
plt.plot(time,rho12_av)
plt.plot(time,rho13_av)
plt.plot(time,rho22_av)
plt.plot(time,rho23_av)
plt.plot(time,rho33_av)
plt.legend(["rho11_av","rho12_av","rho13_av","rho22_av","rho23_av","rho33_av"])
#plt.yscale("log")  
plt.show()



