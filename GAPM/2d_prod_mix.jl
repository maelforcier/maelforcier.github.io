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



c=[-12,-40]
q=[5,10]

# c=[-12., -40.]
# q=[5.,10.]

# W=[-1. 0.; 0. -1.]

# A=[0. 0.]
# b=[0.]

h_mean=[6000., 4000.]
h_covar=[100. 0.;0. 50.]

dif_to_mean=3*[sqrt(h_covar[1,1]), sqrt(h_covar[2,2])]
h_low=h_mean-dif_to_mean
h_up=h_mean+dif_to_mean #I approximate the gaussian random variable coordinates by uniform random variables

T_low=[3.5 9.; 0.8 36.]

T_up=[4.5 11.; 1.2 44.]

T_mean=(T_up+T_low)/2

Txilow=[3.5, 0.8 , 9., 36.]
Txiup=[4.5,1.2, 11.,44.]
Tximean=(Txilow+Txiup)/2



n=length(c)
m=length(q)
# k=length(b)
l=length(h_mean)


