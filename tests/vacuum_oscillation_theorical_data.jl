include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/constants.jl")

s12=sin(theta_12)
c12=cos(theta_12)
s13=sin(theta_13)
c13=cos(theta_13)
s23=sin(theta_23)
c23=cos(theta_23)
sdcp=sin(delta_cp)
cdcp=cos(delta_cp)

E=neutrino_energy

mass=[mass_1,mass_2,mass_3]

print("period 1: ")
println(abs(2*pi/(((mass[1]^2-mass[2]^2)/(2*E*hbar))*eV_to_J)))
print("period 2: ")
println(abs(2*pi/(((mass[1]^2-mass[3]^2)/(2*E*hbar))*eV_to_J)))
print("period 3: ")
println(abs(2*pi/(((mass[2]^2-mass[3]^2)/(2*E*hbar))*eV_to_J)))

U=[c12*c13 s12*c13 s13*e^(-im*delta_cp);-s12*c23-c12*s13*s23*e^(im*delta_cp) c12*c23-s12*s13*s23*e^(im*delta_cp) c13*s23;s12*s23-c12*s13*c23*e^(im*delta_cp) -c12*s23-s12*s13*c23*e^(im*delta_cp) c13*c23]	

function compute_prob(a,b,t)
	Prob=0
	for i in collect(1:3)
		for j in collect(1:3)
			Prob+=U[a,i]*conj(U[b,i])*conj(U[a,j])*U[b,j]*(cos(((mass[i]^2-mass[j]^2)/(2*E*hbar))*eV_to_J*t)+im*sin(((mass[i]^2-mass[j]^2)/(2*E*hbar))*eV_to_J*t))			
		end
	end
	return real(Prob)
end

P11=[]
P12=[]
P13=[]

time=time_step*collect(0:number_of_steps-1)

for w in time
	push!(P11,compute_prob(1,1,w))	
	push!(P12,compute_prob(1,2,w))	
	push!(P13,compute_prob(1,3,w))	
end

io = open("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/tests/vacuum_oscillation_theorical_data.txt", "w");
	write(io, "time rho11 rho22 rho33 \n")
	for i in eachindex(P11)
		p11=P11[i]
		p12=P12[i]
		p13=P13[i]
		t=time[i]
		write(io, "$t $p11 $p12 $p13 \n")
	end
close(io);







	


