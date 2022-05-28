include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")

function particles_initial_condition(x__grid_center,y__grid_center,z__grid_center)
	
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
			push!(x_particle_direction,1)						
			push!(y_particle_direction,0)						
			push!(z_particle_direction,0)						
			push!(particles_rho,[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im])
			push!(particles_rho_bar,[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im])						
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
			push!(particles_rho,[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im])
			push!(particles_rho_bar,[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im])						
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


