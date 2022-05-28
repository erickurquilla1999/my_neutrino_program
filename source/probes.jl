include("/home/erickurquilla1999/Documents/physics/thesis/my_neutrino_program/my_neutrino_program/source/input_parameters.jl")

#=
my_array=[[1,1,1,1,1e-87],[1,1,1,1,1e-87],[1,1,1,1,1e-87],[1,1,1,1,1e-87],[1,1,1,1,1e-87]]

println("dddd")
println(typeof(my_array))
println("dddd")

for i in eachindex(my_array)
	println(my_array[i])	
end

	for j in [1,2,3]
println(j)
end       


if 1<2
println("sssssssssssssss")
end


if 1>=1
println("sssssssssssssss")
end

for i in collect(1:1000)
	println(rand(Float64))
end



q=[1 2;1 3]
w=[1 2;1 3]

println(q*w)

d=1e-3+3e-656im
println(typeof(d))


d=Array{Complex{Float64}}(undef,3,3)
d[1,1]=0.1+0.5im
println(d[1,1])



f=[1.0+5.23im 6.5+6.45im;1.0+5.23im 6.5+6.45im]
println(typeof(f))


	gellmann_1=[0.0+0.0im 1.0+0.0im 0.0+0.0im;1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_2=[0.0+0.0im 0.0-1.0im 0.0+0.0im;0.0+1.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_3=[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im -1.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_4=[0.0+0.0im 0.0+0.0im 1.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;1.0+0.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_5=[0.0+0.0im 0.0+0.0im 0.0-1.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+1.0im 0.0+0.0im 0.0+0.0im]
 	gellmann_6=[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 1.0+0.0im;0.0+0.0im 1.0+0.0im 0.0+0.0im]
 	gellmann_7=[0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0-1.0im;0.0+0.0im 0.0+1.0im 1.0+0.0im]
 	gellmann_8=(1/sqrt(3))*[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 1.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im -2.0+0.0im]
	indentity=[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]
	
	rho=[1.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im;0.0+0.0im 0.0+0.0im 0.0+0.0im]

using LinearAlgebra	

	println(tr(rho))
	uno=0.5^2*(tr(rho*gellmann_1)^2+tr(rho*gellmann_2)^2+tr(rho*gellmann_3)^2+tr(rho*gellmann_4)^2+tr(rho*gellmann_5)^2+tr(rho*gellmann_6)^2+tr(rho*gellmann_7)^2+tr(rho*gellmann_8)^2)
	println(uno)


function my(a,b)
	return a,b
end

println(typeof([1.1 2 2;1.1 2 2;1.1 2 2]))






	cells_limits=Array{Array{Float64,2},2}(undef,2,3)
	println(cells_limits)		
	println("eeeeeeeeeeeeeeeeeeeeeee")
	
	particle_position=Array{Array{Float64,2},2}(undef,2,3)
	println(particle_position)
	
	
	

	

function my(a,b)
	return a,b
end

b,c=my(1,2)
	

println(b)	
println(c)	
	
	
	
	
println(modf(1.2323))
println(typeof(modf(1.2323)[1]))	
println(typeof(modf(1.2323)[2]))	
	
	
#	exit()
	
if 0.5<=1&&0.8<=89
	println("gfgfg")
end

#readline()	

a=[1.1 1.1 1.1;1.1 1.1 1.1;1.1 1.1 1.1]
b=[0.1 0.0 0.0;0.0 0.0 0.0;0.0 0.0 0.0]	
	
a+=b

println(a)	



	n_grid=Array{Array{Float64,2},1}(undef,10)
	n_bar_grid=Array{Array{Float64,2},1}(undef,10)
	println(n_grid)
	for i in eachindex(n_grid)
		n_grid[i]=[0.0 0.0 0.0;0.0 0.0 0.0;0.0 0.0 0.0]
		n_bar_grid[i]=[0.0 0.0 0.0;0.0 0.0 0.0;0.0 0.0 0.0]
	end	
	println(n_grid)








number_of_cells=3
cell_x_lenght=1.0 #cm
cell_y_lenght=1.0 #cm
cell_z_lenght=1.0 #cm

function particles_number_cell(particle_position)

	p_num_cell=Array{Int128,1}(undef,0)

	for i in eachindex(particle_position)
		
		dec,int=modf(particle_position[i][1]/cell_x_lenght)
		
		int+=1.0		
		if (0.0<=dec<=1.0)&&(1<=int<=number_of_cells)&&(0.0<=particle_position[i][2]<=cell_y_lenght)&&(0.0<=particle_position[i][3]<=cell_z_lenght)
			push!(p_num_cell,trunc(Int128,int))
		elseif dec==0.0&&int==number_of_cells+1.0
			push!(p_num_cell,trunc(Int128,int-1.0))
		else
			println("A particle was found out of the grid")
			exit()
		end
	end
	return p_num_cell		
end


parpos=[[0.0 0.5 0.5],[1.0 0.5 0.5],[2.0 0.5 0.5],[3.0 0.5 0.5]]
println(particles_number_cell(parpos))





function shape_function(particle_position,cell_center,num_cell)
	
	particles_shape_function=Array{Array{Float64,2},1}(undef,0)

	for i in eachindex(particle_position)

		#distance from particle to the center of its grid cell
		d_cell=sqrt((particle_position[i][1]-cell_center[num_cell[i]][1])^2+(particle_position[i][2]-cell_center[num_cell[i]][2])^2+(particle_position[i][3]-cell_center[num_cell[i]][3])^2)

		#distance from particle to the center of the next grid cell
		if num_cell[i]==number_of_cells
			d_cell_plus_1=sqrt((particle_position[i][1]-(cell_center[num_cell[i]][1]+cell_x_lenght))^2+(particle_position[i][2]-cell_center[1][2])^2+(particle_position[i][3]-cell_center[1][3])^2)					
		else
			d_cell_plus_1=sqrt((particle_position[i][1]-cell_center[num_cell[i]+1][1])^2+(particle_position[i][2]-cell_center[num_cell[i]+1][2])^2+(particle_position[i][3]-cell_center[num_cell[i]+1][3])^2)				
		end

		#distance from particle to the center of the last grid cell
		if num_cell[i]==1
			d_cell_minus_1=sqrt((particle_position[i][1]-(cell_center[num_cell[i]][1]-cell_x_lenght))^2+(particle_position[i][2]-cell_center[number_of_cells][2])^2+(particle_position[i][3]-cell_center[number_of_cells][3])^2)		
		else
			d_cell_minus_1=sqrt((particle_position[i][1]-cell_center[num_cell[i]-1][1])^2+(particle_position[i][2]-cell_center[num_cell[i]-1][2])^2+(particle_position[i][3]-cell_center[num_cell[i]-1][3])^2)
		end
		
		sf_a=-1.0/cell_x_lenght^2
		sf_c=1.0	
		
		w_cell=sf_a*d_cell^2+sf_c
		w_cell_plus_1=sf_a*d_cell_plus_1^2+sf_c
		w_cell_minus_1=sf_a*d_cell_minus_1^2+sf_c
		
		if w_cell<0.0
			w_cell=0.0
		end
		
		if w_cell_plus_1<0.0
			w_cell_plus_1=0.0
		end
		
		if w_cell_minus_1<0.0
			w_cell_minus_1=0.0
		end

		push!(particles_shape_function,[w_cell_minus_1/(w_cell+w_cell_plus_1+w_cell_minus_1) w_cell/(w_cell+w_cell_plus_1+w_cell_minus_1) w_cell_plus_1/(w_cell+w_cell_plus_1+w_cell_minus_1)])
	end	
	return particles_shape_function
end


cell_center=[[0.5 0.5 0.5],[1.5 0.5 0.5],[2.5 0.5 0.5]]
num_celldd=[3,3,3]
par_pos=[[2.2 0.8 0.6],[2.6 1.0 1.0],[2.4 0.0 0.0]]
number_of_cells=3
cell_x_lenght=1.0 #cm
cell_y_lenght=1.0 #cm
cell_z_lenght=1.0 #cm

println(shape_function(par_pos,cell_center,num_celldd))

=#






















