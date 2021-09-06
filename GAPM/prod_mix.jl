using JuMP
using GLPK
using Polyhedra, CDDLib
using LinearAlgebra
using Polymake

include("partition.jl")
include("polyhedra_to_polymake.jl")
include("GAPM.jl")

#Runner avec julia 1.5 
#../../../../Logiciels_packages/julia-1.5.1/bin/julia

# include("1d_prod_mix.jl")
include("2d_prod_mix.jl")
# include("4d_prod_mix.jl")


support_model= Model(GLPK.Optimizer)
@variable(support_model,T[1:n*l])
@variable(support_model,h[1:l])
@constraint(support_model,upboundT[i=1:n*l],T[i]<=Txiup[i])
@constraint(support_model,lowboundT[i=1:n*l],T[i]>=Txilow[i])
@constraint(support_model,upboundh[i=1:l],h[i]<=h_up[i])
@constraint(support_model,lowboundh[i=1:l],h[i]>=h_low[i])


support=jump_model_to_polymake(support_model)

vol_supp=prod(vcat(Txiup-Txilow,h_up-h_low))
# vol_supp=26066.784381660804
centr_supp=[1, 4, 1, 9, 1, 7, 3, 10, 40, 6000, 4000]




# Old D without simplification
# D=polytope.Polytope(INEQUALITIES=vcat(hcat(q,hcat(transpose(W),-I(m))),hcat(zeros(l+m),I(l+m))))



dual_prod_mix= Model(GLPK.Optimizer)
@variable(dual_prod_mix,λ[1:2])
@constraint(dual_prod_mix,λ[1]<=q[1])
@constraint(dual_prod_mix,λ[2]<=q[2])
@constraint(dual_prod_mix,λ[1]>=0)
@constraint(dual_prod_mix,λ[2]>=0)

D=jump_model_to_polymake(dual_prod_mix)

nFanDprodmix=fan.normal_fan(D)
# HnFanD=hom_fan(nFanD)
# prod_mix=TSLP_unif_T_gaussian_h(n,m,k,l,c,A,b,q,T_low,T_up,W,h_mean,h_covar,D,nFanD)

nFanDprodmix.MAXIMAL_CONES_FACETS











function Q_prod_mix(x::Vector,T::Matrix,h::Vector)
	model= Model(GLPK.Optimizer)
	@variable(model,y[1:m])
	@objective(model, Min, sum(q.*y))
	@constraint(model, cons[j = 1:m], dot(T[j,:],x) <= h[j]+y[j])
	@constraint(model, posy[i=1:m], y[i]>=0)
	optimize!(model)
	return objective_value(model)
end





function V_part_prod_mix(x::Vector,prob::Vector,T::Array,h::Array)
	s=length(prob)
	return sum(prob[k]*Q_prod_mix(x,T[k][:,:],h[k][:]) for k=1:s)
end

function V_part_prod_mix(x::Vector,part::Partition)
	return V_part_prod_mix(x,part.prob,part.T,part.h)
end




function Q_dual_prod_mix(x::Vector,T::Matrix,h::Vector)
	model= Model(GLPK.Optimizer)
	@variable(model,λ[1:l])
	@objective(model, Max, dot(T*x-h,λ))
	@constraint(model, bound[i = 1:l],λ[i] <= q[i])
	@constraint(model, posλ[i=1:l], λ[i]>=0)
	optimize!(model)
	return objective_value(model)
end





function V_dual_part_prod_mix(x::Vector,prob::Vector,T::Array,h::Array)
	s=length(prob)
	return sum(prob[k]*Q_dual_prod_mix(x,T[k][:,:],h[k][:]) for k=1:s)
end

function V_dual_part_prod_mix(x::Vector,part::Partition)
	return V_dual_part_prod_mix(x,part.prob,part.T,part.h)
end




function solve_part_prob_prod_mix(prob::Vector,T::Array,h::Array)
	s=length(prob)
	model= Model(GLPK.Optimizer)
	@variable(model,x[1:n])
	@variable(model,y[1:m,1:s])
	@objective(model, Min, dot(c,x)+sum(prob[k]*dot(q,y[:,k]) for k=1:s))
	@constraint(model, cons[j = 1:m,k=1:s], dot(T[k][j,:],x) <= h[k][j]+y[j,k])
	@constraint(model, posx[i=1:n], x[i]>=0)
	@constraint(model, posy[i=1:m,k=1:s], y[i,k]>=0)
	optimize!(model)
	return(value.(x),objective_value(model))
end


function solve_part_prob_prod_mix(part::Partition)
	return solve_part_prob_prod_mix(part.prob,part.T,part.h)
end


# function solve_part_prob_prod_mix(polcmplx::Polymake.BigObjectAllocated,vol::Scalar)
# 	(prob,T,h)=polcmplx
# 	return solve_part_prob_prod_mix(part.prob,part.T,part.h)
# end


function probTh(polcmplx::Polymake.BigObjectAllocated,vol::Int64)
	s=polcmplx.N_MAXIMAL_CONES
	prob=zeros(s)
	T=Vector{Matrix{Float64}}(undef,s)
	h=Vector{Vector{Float64}}(undef,s)
	for k in 1:s
		P=poly(polcmplx,k)
		prob[k]=P.VOLUME/vol
		ξ=P.CENTROID
		(T[k],h[k])=ξtoTh(ξ)
	end
	return (prob,T,h)
end



function GAPM_prod_mix(stopindex::Int64)
	i=1
	cPtemp=cP0
	ub=100000
	lb=-100000
	lb=Vector{Float64}(undef,stopindex)
	ub=Vector{Float64}(undef,stopindex)
	cP=Vector{Partition}(undef,stopindex)
	added=Vector{Vector{Int64}}(undef,stopindex)
	x=Vector{Vector{Float64}}(undef,stopindex)
	while i<=stopindex 
	    println(i)
	    (x[i],lb[i])=solve_part_prob_prod_mix(cPtemp)
	    print("actual point x ")
	    println(round.(x[i],digits=3))
	    print("lower bound ")
	    println(lb[i])
	    (added[i],cP[i])=insert_cRx_in_part(cPtemp,x[i],1)
	    cPtemp=cP[i]
	    ub[i]=dot(c,x[i])+V_part_prod_mix(x[i],cPtemp)
	    print("upper bound ")
	    println(ub[i])
	    # print("upper bound dual")
	    # println(dot(c,x[i])+V_dual_part_prod_mix(x[i],cP[i]))
	    print( "size of actual adapted partition ")
	    println(cPtemp.sze)
	    println( )
	    i=i+1
	end
	# tot_add=Vector{Vector{Int64}}(undef,length(added))
	# for i in 1:length(added)
	# 	tot_add[i]=Vector{Int64}(undef,length(added[i]))
	# 	for j in 1:length(added[i])
	# 		tot_add[i][j]=sum(added[i][1:j])
	# 	end
	# end
	return (cP,x,lb,ub,added)
end


# for i in 1:6
#     tot_prob[i]=vec_total(cP[i].prob)
# end


function SAA_prod_mix(s::Int64)
	# Input s number of scenario
	T_rand01=rand(l,m,s)
	h_rand01=rand(l,s)
	T_rand=Vector{Matrix{Float64}}(undef,s)
	h_rand=Vector{Vector{Float64}}(undef,s)
	for k in 1:s
		T_rand[k]=T_low.+T_rand01[:,:,k].*(T_up.-T_low)
		h_rand[k]=h_low.+h_rand01[:,k].*(h_up.-h_low)
	end
	return solve_part_prob_prod_mix(1/s*ones(s),T_rand,h_rand)
end


function MC_SAA_prod_mix(s::Int64,nb_solved::Int64)
	v=Vector{Float64}(undef,nb_solved)
	for i in 1:nb_solved
		println(i)
		(x,v[i])=SAA_prod_mix(s)
	end
	return v
end







T0=zeros(l,n,1)
h0=zeros(l,1)
prob0=ones(1)
T0=[T_mean]
h0=[h_mean]

# cP0=fan.PolyhedralFan(FACET_NORMALS=support.FACETS,MAXIMAL_CONES_FACETS=ones(support.N_FACETS))
# cP0=fan.PolyhedralComplex(FACET_NORMALS=support.FACETS,MAXIMAL_POLYTOPES_FACETS=ones(support.N_FACETS))
cP0=PartitionMultiDim([support],prob0,T0,h0,1,vol_supp)







