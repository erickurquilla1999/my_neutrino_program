import numpy as np
import matplotlib.pyplot as plt

#########################################

#reading numerical data

time=[]
rho11=[]
rho22=[]
rho33=[]

for i in range(0,200):
	data = np.loadtxt("../source/output/step_"+str(i)+".txt", unpack = True, skiprows=1)
	
	#data indices
	#0. x, 1. x_dir
	#2. rho11, 3. rho12, 4. rho13, 5. rho22, 6. rho23, 7. rho33 
	#14. rho11_im, 15. rho12_im, 16. rho13_im, 17. rho22_im, 18. rho23_im, 19. rho33_im 
	#8. rho_bar11, 9. rho_bar12, 10. rho_bar13, 11. rho_bar22, 12. rho_bar23, 13. rho_bar33
	#20. rho_bar11_im, 21. rho_bar12_im, 22. rho_bar13_im, 23. rho_bar22_im, 24. rho_bar23_im, 25. rho_bar33_im 
	#26. N, 27. N_bar, 28. t

	time.append(data[28][30])
	rho11.append(data[2][30])
	rho22.append(data[5][30])
	rho33.append(data[7][30])
	
#########################################	

#reading theorical data
	
time_teo=[]
rho11_teo=[]
rho22_teo=[]
rho33_teo=[]

data1 = np.loadtxt("vacuum_oscillation_theorical_data.txt", unpack = True, skiprows=1)	
time_teo=data1[0]
rho11_teo=data1[1]
rho22_teo=data1[2]
rho33_teo=data1[3]
	
#########################################
	
#ploting	
	
plt.plot(time,rho11)	
plt.plot(time_teo,rho11_teo)
plt.plot(time,rho22)	
plt.plot(time_teo,rho22_teo)
plt.plot(time,rho33)	
plt.plot(time_teo,rho33_teo)	
plt.legend([r"$\rho_{11}^{num}$",r"$\rho_{11}^{theo}$",r"$\rho_{22}^{num}$",r"$\rho_{22}^{theo}$",r"$\rho_{33}^{num}$",r"$\rho_{33}^{theo}$"])
plt.show()


