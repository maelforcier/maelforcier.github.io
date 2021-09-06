using JuMP
using GLPK
using Polyhedra, CDDLib
using LinearAlgebra
using Polymake

include("partition.jl")
include("polyhedra_to_polymake.jl")


#Runner avec julia 1.5 
#../../../../Logiciels_packages/julia-1.5.1/bin/julia




function ENx(myfan::Polymake.BigObjectAllocated,i::Int,x::Vector) 
	#Compute the set E_{N,x} of GAPM 
	#Input myfan: the normal fan
	 #i: the index of the cone so that N=myfan.MAXIMAL_CONES[i]
	 #x: the point x
	n=length(x)
	p=myfan.N_FACET_NORMALS
	l=myfan.FAN_AMBIENT_DIM
	H=zeros(p,n*l+l)
	M=hrepcone(myfan,i)
	for k in 1:p
		for i in 1:l
			for j in 1:n
				H[k,i+l*(j-1)]=-M[k,i]*x[j]
			end
		end
		for i in 1:l
			H[k,i+l*n]=M[k,i]
		end
	end
	return polytope.Polytope(INEQUALITIES=H)
end



function cRx(myfan::Polymake.BigObjectAllocated,x::Vector,sgn::Int64=1)
	# Returns the polyhedral fan cR_x in the Ξ space
	# H=(-x_1M  -x_2M ... -x_nM M)
	# Or H=(x_1M  x_2M ... x_nM -M) in function of the sign i.e. if the objective of the dual is h-Tx (1) or Tx-h (-1)
	n=length(x)
	my_facet_normals=myfan.FACET_NORMALS
	facet_normal_crx=-sgn*x[1]*my_facet_normals
	for i in 2:n
		facet_normal_crx=hcat(facet_normal_crx,-sgn*x[i]*my_facet_normals)
	end
	facet_normal_crx=hcat(facet_normal_crx,sgn*my_facet_normals)
	return fan.PolyhedralFan(FACET_NORMALS=facet_normal_crx,MAXIMAL_CONES_FACETS=myfan.MAXIMAL_CONES_FACETS)
end

function com_ref_cRx(myfan::Polymake.BigObjectAllocated,x::Vector)
	fan.common_refinement(myfan,cRx(nFanD,x))
end



function ξtoTh(ξ::Vector,n::Int64,l::Int64)
	T=zeros(l,n)
	h=zeros(l)
	for j in 1:n
		T[:,j]=ξ[(j-1)*l+1:j*l]
	end
	h[:]=ξ[n*l+1:n*l+l]
	return (T,h)
end

function ξtoTh(ξ::Polymake.VectorAllocated{Polymake.Rational},n::Int64,l::Int64)
	T=zeros(l,n)
	h=zeros(l)
	for j in 1:n
		T[:,j]=ξ[(j-1)*l+1:j*l]
	end
	h[:]=ξ[n*l+1:n*l+l]
	return (T,h)
end




function insert_cRx_in_part(part::PartitionMultiDim,x::Vector,sgn::Int64=1)
	cR=cRx(nFanDprodmix,x,sgn)
	sze=0
	polylist=Vector{Polymake.BigObjectAllocated}(undef,0)
	problist=Vector{Float64}(undef,0)
	Tlist=Vector{Matrix{Float64}}(undef,0)
	hlist=Vector{Vector{Float64}}(undef,0)
	nbadded=Vector{Int64}(undef,part.sze)
	for k in 1:part.sze
		poly=part.poly[k]
		added_at_k=0
		for j in 1:cR.N_MAXIMAL_CONES
			inter=polytope.Polytope(INEQUALITIES=vcat(poly.INEQUALITIES,hcat(zeros(cR.N_FACET_NORMALS),hrepcone(cR,j))))
			if inter.FULL_DIM
				sze=sze+1
				vol=inter.VOLUME
				ξ=inter.CENTROID
				(T,h)=ξtoTh(ξ[2:(n+1)*l+1],n,l)
				insert!(polylist,sze,inter)
				insert!(problist,sze,vol/part.vol)
				insert!(Tlist,sze,T)
				insert!(hlist,sze,h)
				added_at_k=added_at_k+1
			end
		end
		nbadded[k]=added_at_k
	end
	return (nbadded,PartitionMultiDim(polylist,problist,Tlist,hlist,sze,part.vol))
end