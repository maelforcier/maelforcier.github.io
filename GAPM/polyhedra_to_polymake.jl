
function Hpolyhedra_to_polymake(P::Polyhedron)
	i=1
	T=Polyhedra.coefficient_type(P)
	d=fulldim(P)
	H=halfspaces(hrep(P))
	n=length(H)
	M=Matrix{T}(undef,n,d+1)
	for x in H
		M[i,1]=x.β
		M[i,2:d+1]=-x.a
		i=i+1
	end
	# if T==Float64
	# 	PM= @pm polytope.Polytope{Float}(INEQUALITIES=M)  #Nomal fan does not work with Floats
	# else
		PM= @pm polytope.Polytope{Rational}(INEQUALITIES=M)
	# end
	return PM
end


function jump_model_to_polymake(model::Model)
	P=polyhedron(model,CDDLib.Library())
	removehredundancy!(P)
	return Hpolyhedra_to_polymake(P)
end






function cone(myfan::Polymake.BigObjectAllocated,i::Int)
	return polytope.Cone(RAYS=myfan.RAYS[myfan.MAXIMAL_CONES[i,:],:])
end

function poly_from_fan(myfan::Polymake.BigObjectAllocated,i::Int)
	C=cone(myfan,i)
	R=C.RAYS
	return polytope.Polytope(POINTS=R[:,2:C.CONE_AMBIENT_DIM])
end


function poly(mypolycomplex::Polymake.BigObjectAllocated,i::Int)
	return polytope.Polytope(POINTS=mypolycomplex.VERTICES[mypolycomplex.MAXIMAL_POLYTOPES[i,:],:])
end

function hrepcone(myfan::Polymake.BigObjectAllocated,i::Int64)
	return diagm(myfan.MAXIMAL_CONES_FACETS[i,:])*myfan.FACET_NORMALS
end


function hom_cone(cone::Polymake.BigObjectAllocated)
	#Compute the cone affine cone in homogeneous coordinates
	return polytope.Polytope(POINTS=
		vcat(hcat(zeros(cone.N_RAYS,1),cone.RAYS),hcat([1],zeros(1,cone.CONE_AMBIENT_DIM))) )
end


function hom_fan(myfan::Polymake.BigObjectAllocated)
	return fan.PolyhedralFan(
		RAYS=vcat(hcat([1],zeros(1,myfan.FAN_AMBIENT_DIM)),hcat(zeros(myfan.N_RAYS,1),myfan.RAYS))
		,MAXIMAL_CONES=Polymake.IncidenceMatrix(hcat(ones(myfan.N_MAXIMAL_CONES,1),myfan.MAXIMAL_CONES))
		,LINEALITY_SPACE=hcat(zeros(myfan.LINEALITY_DIM,1),myfan.LINEALITY_SPACE))
end

function fan_to_cmplx(myfan::Polymake.BigObjectAllocated)
	return fan.PolyhedralComplex(POINTS=vcat(hcat([1],zeros(1,myfan.FAN_AMBIENT_DIM)),hcat(zeros(myfan.N_RAYS,1),myfan.RAYS))
		,MAXIMAL_CONES=Polymake.IncidenceMatrix(hcat(ones(myfan.N_MAXIMAL_CONES,1),myfan.MAXIMAL_CONES))
		,LINEALITY_SPACE=hcat(zeros(myfan.LINEALITY_DIM,1),myfan.LINEALITY_SPACE))
end


function hom_cmplx(myfan::Polymake.BigObjectAllocated)
	return fan.PolyhedralComplex(
		POINTS=vcat(hcat([1],zeros(1,myfan.FAN_AMBIENT_DIM)),hcat(zeros(myfan.N_RAYS,1),myfan.RAYS))
		,MAXIMAL_CONES=Polymake.IncidenceMatrix(hcat(ones(myfan.N_MAXIMAL_CONES,1),myfan.MAXIMAL_CONES))
		,LINEALITY_SPACE=hcat(zeros(myfan.LINEALITY_DIM,1),myfan.LINEALITY_SPACE))
end
