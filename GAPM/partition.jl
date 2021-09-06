using Polymake


abstract type Partition end

mutable struct Partition1dim<:Partition
	low::Float64      			#Lowest value
	up::Float64       			#Upper value
	lgth::Float64				#Length of the interval
	bnds::Vector{Float64}		#List of bounds of the Partition
	prob::Vector{Float64}		#Probabilities of each element
	esp::Vector{Float64}		#Conditional expectation of each element
	sze::Int64					#Size of the partition i.e. length of prob or esp
end




function Partition1dim(bnds::Vector)
	sort!(bnds)
	sze=length(bnds)-1
	low=bnds[1]
	up=bnds[sze+1]
	lgth=up-low
	prob=Array{Float64}(undef,sze)
	esp=Array{Float64}(undef,sze)
	for i in 1:sze
		prob[i]=(bnds[i+1]-bnds[i])/lgth
		esp[i]=(bnds[i+1]+bnds[i])/2
	end
	return Partition1dim(low,up,lgth,bnds,prob,esp,sze)
end






mutable struct PartitionMultiDim<:Partition
	poly::Vector{Polymake.BigObjectAllocated}		#Polyhedral complex representing the partition
	prob::Vector{Float64}		#Probabilities of each element
	T::Vector{Matrix{Float64}}			#Conditional expectation of each element T Matrix
	h::Vector{Vector{Float64}}			#Conditional expectation of each element h Vector
	sze::Int64					#Size of the partition i.e. length of prob or esp
	vol::Float64				#Volume of the whole polyehdral complex
end






function PartitionMultiDim(cmplx::Polymake.BigObjectAllocated,vol::Float64)
	sze=cmplx.N_MAXIMAL_POLYTOPES
	polylist=Array{Polymake.BigObjectAllocated}(undef,sze)
	prob=Array{Float64}(undef,sze)
	T=Vector{Matrix{Float64}}(undef,sze)
	h=Vector{Vector{Float64}}(undef,sze)
	for k in 1:sze
		P=poly(cmplx,k)
		polylist[k]=P
		prob[k]=P.VOLUME/vol
		ξ=P.CENTROID
		(T[k],h[k])=ξtoTh(ξ,n,l)
	end
	return PartitionMultiDim(polylist,prob,T,h,sze,vol)
end



function insert_in_part!(part::Partition1dim,x)
	if x<part.low || x>part.up
		error("Inserted real out of bounds")
	end
	for i in 2:part.sze+1
		if x<part.bnds[i]&& x>part.bnds[i-1]
			insert!(part.esp,i,(part.bnds[i]+x)/2)
			part.esp[i-1]=(x+part.bnds[i-1])/2
			insert!(part.prob,i,(part.bnds[i]-x)/part.lgth)
			part.prob[i-1]=(x-part.bnds[i-1])/part.lgth
			insert!(part.bnds,i,x)
			part.sze=part.sze+1
			break
		end
	end
end


function insert_in_part!(part::Partition1dim,v::Vector)
	for i in 1:length(v)
		insert_in_part!(part,v[i])
	end
end



function refinement(part1::Partition1dim,part2::Partition1dim)
	if part1.sze<part2.sze
		return refinement(part2,part1)
	else
		part=copy(part1)
		insert_in_part!(part,part2.bnds)
		return part
	end
end


function insert_in_part_cmplx!(part::PartitionMultiDim,myfan::Polymake.BigObjectAllocated)
	# TO FINISH
	for k in 1:part.sze
		bool_intersect=false
		added_poly=[]
		print("k")
			println(k)
		for j in 1:myfan.N_MAXIMAL_CONES
			print("j")
			println(j)
			inter=polytope.intersection(part.poly[k],cone(myfan,j))
			println(inter.FEASIBLE)
			if inter.FEASIBLE
				bool_intersect=true
				added_poly=hcat(added_poly,inter)
			end	
		end
		if bool_intersect
			println("You had an intersection")
		end
		println( )
	end
end


