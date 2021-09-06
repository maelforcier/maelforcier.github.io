using JuMP
using GLPK
using Polyhedra, CDDLib
using LinearAlgebra
using Polymake

include("partition.jl")
include("polyhedra_to_polymake.jl")


#Runner avec julia 1.5 
#../../../../Logiciels_packages/julia-1.5.1/bin/julia



c=[10., 7., 16., 6.]
f=[40 24 4;45 27 4.5; 32 19.2 3.2; 55 33 5.5]
m=12
b=120
ξmin=3
ξmax=7
dmin=[ξmin,3,2]
dmax=[ξmax,3,2]







duallands= Model(GLPK.Optimizer)
@variable(duallands,ν[1:4])
@variable(duallands,𝑈[1:3])
@constraint(duallands, con[i = 1:4,j =1:3], 𝑈[j]-ν[i] <= f[i,j])
@constraint(duallands, con𝑈[i=1:4], ν[i]>=0)
@constraint(duallands, conν[j=1:3], 𝑈[j]>=0)


D = jump_model_to_polymake(duallands) 
RC=polytope.recession_cone(D)

nFanD=fan.normal_fan(D)


normalCones=Array{Polyhedron}(undef,nFanD.N_MAXIMAL_CONES)
println(nFanD.RAYS[nFanD.MAXIMAL_CONES[1,:],:])



# HnFanD=fan.PolyhedralFan(INPUT_RAYS=vcat(hcat(zeros(nFanD.N_RAYS,1),nFanD.RAYS),hcat([1],zeros(1,nFanD.FAN_AMBIENT_DIM)))
# 	,MAXIMAL_CONES=Polymake.IncidenceMatrix(hcat(ones(nFanD.N_MAXIMAL_CONES,1),nFanD.MAXIMAL_CONES)))


HnFanD=hom_fan(nFanD)


Ξ=polytope.Polytope(POINTS=[1 0 0 0 0 ξmin 3 2;1 0 0 0 0 ξmax 3 2])

function dualcostset(x::Vector)
	return polytope.Polytope(POINTS=hcat([1,1],transpose(hcat(vcat(x,-dmin),vcat(x,-dmax)))))
end	



insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x), [x]); v)
#insert and suppress the equal value in a sorted array

function cRlands1dim(x::Vector) #Returns the partition of cR_x on the dimension d_1
	cR=[]
	costset=dualcostset(x)
	for i=1:HnFanD.N_MAXIMAL_CONES
		inter=polytope.intersection(cone(HnFanD,i),costset)
		if polytope.dim(inter)!=-1
			pnts=inter.VERTICES
			for j in 1:length(pnts[:,1])
				insert_and_dedup!(cR,-pnts[j,6])
			end
		end
	end
	return cR
end


x1=[0.833, 3, 4.167,4]
x2=[2.5, 3, 3.5, 3]
x3=[1.833, 4,3.667,2.5]
x4=[2,4.167,3.583,2.250]











function create2stage(d1)
	model= Model(GLPK.Optimizer)
	d=vcat(d1,dmin[2:3])
	@variable(model,x[1:4])
	@variable(model,y[1:4,1:3])
	@objective(model, Min, dot(c,x)+sum(f.*y))
	@constraint(model, capacSnd[i = 1:4], sum(y[i,:]) <= x[i])
	@constraint(model, deman[j = 1:3], sum(y[:,j]) >= d[j])
	@constraint(model, posx[i=1:4], x[i]>=0)
	@constraint(model, posy[i=1:4,j=1:3], y[i,j]>=0)
	@constraint(model,budget,dot(c,x)<=b)
	@constraint(model,cpacFst,sum(x) >=m)
	optimize!(model)
	println(value.(x))
	@show objective_value(model)
	return model
end


create2stage(5);


function Q_lands(x::Vector,d1::Float64)
	model= Model(GLPK.Optimizer)
	d=vcat(d1,dmin[2:3])
	@variable(model,y[1:4,1:3])
	@objective(model, Min, sum(f.*y))
	@constraint(model, capacSnd[i = 1:4], sum(y[i,:]) <= x[i])
	@constraint(model, deman[j = 1:3], sum(y[:,j]) >= d[j])
	@constraint(model, posy[i=1:4,j=1:3], y[i,j]>=0)
	optimize!(model)
	return objective_value(model)
end



function V_part_lands(x::Vector,part::Partition)
	s=part.sze
	return sum(part.prob[k]*Q_lands(x,part.esp[k]) for k=1:s)
end



function solve_part_prob_lands(part::Partition)
	s=part.sze
	model= Model(GLPK.Optimizer)
	@variable(model,x[1:4])
	@variable(model,y[1:4,1:3,1:s])
	@objective(model, Min, dot(c,x)+sum(part.prob[k]*sum(f.*y[:,:,k]) for k=1:s))
	@constraint(model, capacSnd[i = 1:4,k=1:s], sum(y[i,:,k]) <= x[i])
	@constraint(model, deman[j = 2:3,k=1:s], sum(y[:,j,k]) >= dmin[j])
	@constraint(model, demanpart[k=1:s], sum(y[:,1,k]) >= part.esp[k])
	@constraint(model, posx[i=1:4], x[i]>=0)
	@constraint(model, posy[i=1:4,j=1:3,k=1:s], y[i,j,k]>=0)
	@constraint(model,budget,dot(c,x)<=b)
	@constraint(model,cpacFst,sum(x) >=m)
	optimize!(model)
	return(value.(x),objective_value(model))
end

part1=Partition1dim([3,4.167,5,7])


function GAPM_lands(stopindex::Int64)
	i=1
	part=Partition1dim([3,7])
	ub=100000
	lb=-100000
	while i<=stopindex 
	    println(i)
	    (x,lb)=solve_part_prob_lands(part)
	    print("actual point x ")
	    println(round.(x,digits=3))
	    print("lower bound ")
	    println(lb)
	    cR=cRlands1dim(x)
	    insert_in_part!(part,cR)
	    #ub=min(ub,dot(c,x)+V_part_lands(x,part))
	    ub=dot(c,x)+V_part_lands(x,part)
	    print("upper bound ")
	    println(ub)
	    print( "actual adapted partition cR_x ")
	    println(round.(Float64.(cR),digits=3))
	    println( )
	    i=i+1
	end
end

GAPM_lands(6)