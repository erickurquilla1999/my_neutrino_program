
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
rho22_av=[]
rho12_av=[]

rho11=[]
rho22=[]

for i in range(0,5000):
	data = np.loadtxt("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/output/step_"+str(i)+".txt", unpack = True, skiprows=1)
	
	time.append(data[32][30])
	
	rho11.append(data[6][30])
	rho22.append(data[9][30])	
	
	rho11_av.append(sum(data[6])/len(data[6]))
	rho22_av.append(sum(data[9])/len(data[9]))
	
	rho12_amp=[]
	for i in range(0,len(data[7])):
		rho12_amp.append((data[7][i]**2+data[19][i]**2)**0.5)	
	
	rho12_av.append(sum(rho12_amp)/len(rho12_amp))	
	
	
plt.plot(time,rho11_av)
plt.plot(time,rho11)
plt.yscale("log")  
plt.show()
plt.plot(time,rho22_av)
plt.plot(time,rho22)
plt.yscale("log")  
plt.show()
plt.plot(time,rho12_av)
plt.yscale("log")  
plt.show()



