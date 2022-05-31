include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")

function particles_initial_condition(x__grid_center,y__grid_center,z__grid_center)

	println("Creating particles initial condition")								
	
	#In this function is created the particles initial condition data,
	#two particle per cell, one to right on to left,
	# particles are all in electron and antielectron state including small perturbation
	#the number of neutrinos and antineutrinos is ramdom between a number density
	
	x_particle_position=Array{Float64,1}(undef,0)
	y_particle_position=Array{Float64,1}(undef,0)
	z_particle_position=Array{Float64,1}(undef,0)
	x_particle_direction=Array{Float64,1}(undef,0)
	y_particle_direction=Array{Float64,1}(undef,0)
	z_particle_direction=Array{Float64,1}(undef,0)
	particles_rho=Array{Array{Complex{Float64},2},1}(undef,0)  
	particles_rho_bar=Array{Array{Complex{Float64},2},1}(undef,0)
	number_of_neutrinos=Array{Float64,1}(undef,0)
	number_of_antineutrinos=Array{Float64,1}(undef,0)
	
	for i in eachindex(x__grid_center)
			
			#firs particle to the right
			push!(x_particle_position,x__grid_center[i])
			push!(y_particle_position,y__grid_center[i])
			push!(z_particle_position,z__grid_center[i])			
			push!(x_particle_direction,1)	#x unitary vector velocity					
			push!(y_particle_direction,0)	#y unitary vector velocity					
			push!(z_particle_direction,0)	#z unitary vector velocity					
			rho12=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)					
			rho13=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)
			rho22=perturbation_amplitud*rand(Float64)
			rho33=perturbation_amplitud*rand(Float64)
			push!(particles_rho,[(1.0+0.0im-rho22-rho33) rho12 rho13;conj(rho12) rho22 0.0+0.0im;conj(rho13) 0.0+0.0im rho33])
			rho12=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)					
			rho13=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)
			rho22=perturbation_amplitud*rand(Float64)
			rho33=perturbation_amplitud*rand(Float64)
			push!(particles_rho_bar,[(1.0+0.0im-rho22-rho33) rho12 rho13;conj(rho12) rho22 0.0+0.0im;conj(rho13) 0.0+0.0im rho33])
			cell_vol=cell_x_lenght*cell_y_lenght*cell_z_lenght								
			#push!(number_of_neutrinos,2)						
			#push!(number_of_antineutrinos,2)						
			push!(number_of_neutrinos,cell_vol*particles_number_density+cell_vol*particles_number_density*rand(Float64))						
			push!(number_of_antineutrinos,cell_vol*particles_number_density+cell_vol*particles_number_density*rand(Float64))						
			
			#second particle to the left
			push!(x_particle_position,x__grid_center[i])
			push!(y_particle_position,y__grid_center[i])
			push!(z_particle_position,z__grid_center[i])			
			push!(x_particle_direction,-1)						
			push!(y_particle_direction,0)						
			push!(z_particle_direction,0)	
			rho12=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)					
			rho13=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)
			rho22=perturbation_amplitud*rand(Float64)
			rho33=perturbation_amplitud*rand(Float64)
			push!(particles_rho,[(1.0+0.0im-rho22-rho33) rho12 rho13;conj(rho12) rho22 0.0+0.0im;conj(rho13) 0.0+0.0im rho33])
			rho12=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)					
			rho13=perturbation_amplitud*(rand(Float64)+rand(Float64)*im)
			rho22=perturbation_amplitud*rand(Float64)
			rho33=perturbation_amplitud*rand(Float64)
			push!(particles_rho_bar,[(1.0+0.0im-rho22-rho33) rho12 rho13;conj(rho12) rho22 0.0+0.0im;conj(rho13) 0.0+0.0im rho33])
			cell_vol=cell_x_lenght*cell_y_lenght*cell_z_lenght									
			#push!(number_of_neutrinos,1)						
			#push!(number_of_antineutrinos,1)						
			push!(number_of_neutrinos,cell_vol*particles_number_density+cell_vol*particles_number_density*rand(Float64))						
			push!(number_of_antineutrinos,cell_vol*particles_number_density+cell_vol*particles_number_density*rand(Float64))
	end
	return x_particle_position,y_particle_position,z_particle_position,x_particle_direction,y_particle_direction,z_particle_direction,particles_rho,particles_rho_bar,number_of_neutrinos,number_of_antineutrinos
end



















function rho_perturbation(rho,rho_bar)

end

























function normalization(rho,rho_bar)
	
	gellmann_1=[0.0+0.0im 1.0+0.0im 0.0+0.0im;1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_2=[0.0+0.0im 0.0-1.0im 0.0+0.0im;0.0+1.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_3=[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im -1.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_4=[0.0+0.0im 0.0+0.0im 1.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;1.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_5=[0.0+0.0im 0.0+0.0im 0.0-1.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+1.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_6=[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 1.0+0.0im;0.0+0.0im 1.0+0.0im 0.0+0.0im]
 	gellmann_7=[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0-1.0im;0.0+0.0im 0.0+1.0im 1.0+0.0im]
 	gellmann_8=(1/sqrt(3))*[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 1.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im -2.0+0.0im]

end
















function write_info(name,particle_data,time)
	
	println("Writting $name data")
	
	io = open("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/output/$name", "w");
		write(io, "x_position y_position z_position x_direction y_direction z_direction rho11_re rho12_re rho13_re rho22_re rho23_re rho33_re rho_bar11_re rho_bar12_re rho_bar13_re rho_bar22_re rho_bar23_re rho_bar33_re rho11_im rho12_im rho13_im rho22_im rho23_im rho33_im rho_bar11_im rho_bar12_im rho_bar13_im rho_bar22_im rho_bar23_im rho_bar33_im number_of_neutrinos number_of_antineutrinos time\n")
		for i in eachindex(particle_data[1])
			
			x=particle_data[1][i]
			y=particle_data[2][i]
			z=particle_data[3][i]
			
			xd=particle_data[4][i]
			yd=particle_data[5][i]
			zd=particle_data[6][i]
			
			rho11=real(particle_data[7][i][1,1])
			rho12=real(particle_data[7][i][1,2])
			rho13=real(particle_data[7][i][1,3])
			#rho21=real(particle_data[7][i][2,1])
			rho22=real(particle_data[7][i][2,2])
			rho23=real(particle_data[7][i][2,3])
			#rho31=real(particle_data[7][i][3,1])
			#rho32=real(particle_data[7][i][3,2])
			rho33=real(particle_data[7][i][3,3])
			
			rho_bar11=real(particle_data[8][i][1,1])
			rho_bar12=real(particle_data[8][i][1,2])
			rho_bar13=real(particle_data[8][i][1,3])
			#rho_bar21=real(particle_data[8][i][2,1])
			rho_bar22=real(particle_data[8][i][2,2])
			rho_bar23=real(particle_data[8][i][2,3])
			#rho_bar31=real(particle_data[8][i][3,1])
			#rho_bar32=real(particle_data[8][i][3,2])
			rho_bar33=real(particle_data[8][i][3,3])

			rho11_im=imag(particle_data[7][i][1,1])
			rho12_im=imag(particle_data[7][i][1,2])
			rho13_im=imag(particle_data[7][i][1,3])
			#rho21_im=imag(particle_data[7][i][2,1])
			rho22_im=imag(particle_data[7][i][2,2])
			rho23_im=imag(particle_data[7][i][2,3])
			#rho31_im=imag(particle_data[7][i][3,1])
			#rho32_im=imag(particle_data[7][i][3,2])
			rho33_im=imag(particle_data[7][i][3,3])
			
			rho_bar11_im=imag(particle_data[8][i][1,1])
			rho_bar12_im=imag(particle_data[8][i][1,2])
			rho_bar13_im=imag(particle_data[8][i][1,3])
			#rho_bar21_im=imag(particle_data[8][i][2,1])
			rho_bar22_im=imag(particle_data[8][i][2,2])
			rho_bar23_im=imag(particle_data[8][i][2,3])
			#rho_bar31_im=imag(particle_data[8][i][3,1])
			#rho_bar32_im=imag(particle_data[8][i][3,2])
			rho_bar33_im=imag(particle_data[8][i][3,3])
			
			N=particle_data[9][i]
			N_bar=particle_data[10][i]
			
			t=time

			write(io, "$x $y $z $xd $yd $zd $rho11 $rho12 $rho13 $rho22 $rho23 $rho33 $rho_bar11 $rho_bar12 $rho_bar13 $rho_bar22 $rho_bar23 $rho_bar33 $rho11_im $rho12_im $rho13_im $rho22_im $rho23_im $rho33_im $rho_bar11_im $rho_bar12_im $rho_bar13_im $rho_bar22_im $rho_bar23_im $rho_bar33_im $N $N_bar $t\n")
		end
	close(io);
end












































