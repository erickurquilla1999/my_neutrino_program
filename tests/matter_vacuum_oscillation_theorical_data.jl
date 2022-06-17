include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")
include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/constants.jl")

E=neutrino_energy

mass=[mass_1,mass_2,mass_3]

println(sqrt(2)*Gf*(hbar*c)^3)

lambda=sqrt(2)*Gf*(hbar*c)^3*(background_electron_number_density-background_antielectron_number_density)

w_vacuum=abs(((mass[1]^2-mass[2]^2)/(2*E*hbar))*eV_to_J)

w_matter=w_vacuum*sqrt(((lambda/(hbar*w_vacuum))-cos(2*theta_12))^2+sin(2*theta_12)^2)

period_vacuum=2*pi/w_vacuum
period_matter=2*pi/w_matter

theta_12_matter=abs(0.5*atan(sin(2*theta_12)/(cos(2*theta_12)-(lambda/(hbar*w_vacuum)))))

println(" ")
println("lambda: $lambda")
println(" ")
println("omega vacuum: $w_vacuum")
println("omega matter: $w_matter")
println(" ")
println("period vacuum: $period_vacuum")
println("period matter: $period_matter")
println(" ")
println("theta vacuum: $theta_12")
println("theta matter: $theta_12_matter")

P11=[]

time=time_step*collect(0:number_of_steps-1)

for w in time
	push!(P11,1-sin(2*theta_12_matter)^2*sin(w_matter*w)^2)
end

io = open("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/tests/matter_vacuum_oscillation_theorical_data.txt", "w");
	write(io, "time rho11\n")
	for i in eachindex(P11)
		p11=P11[i]
		t=time[i]
		write(io, "$t $p11\n")
	end
close(io);







	


