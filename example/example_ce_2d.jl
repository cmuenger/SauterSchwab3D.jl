using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D

include("singularity_detection.jl")

#=
const pI = point(1,0,0)
const pII = point(0,0,0)
const pIII = point(0.5,0.7,0)
const pIV = point(0.8,-0.3,0)
=#
const pI   = point(0,0,0)
const pII  = point(1,0,0)
const pIII = point(0,1,0)
 
const qIII = point(0,-1,0)


const P = simplex(pI,pIII,pII)
const Q = simplex(pI,qIII,pII)

const sing = singularity_detection(P,Q)


Accuracy2 = 15
ce_ref = CommonEdge4D(sing,SauterSchwab3D._legendre(Accuracy2,0.0,1.0))

function integrand(x,y)
      return ((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
end

function INTEGRAND(u,v)
   n1 = neighborhood(P,u)
   n2 = neighborhood(Q,v)
   x = cartesian(n1)
   y = cartesian(n2)
   output = integrand(x,y)*jacobian(n1)*jacobian(n2)
   return(output)
end

print("Ref: ")
ref = sauterschwab_parameterized(INTEGRAND, ce_ref)
println(ref)
println()

res_tp =ComplexF64[]
res_sp =ComplexF64[]
res_gm =ComplexF64[]
n1 = Int[]
n2 = Int[]
n3 = Int[]
for i in 2:1:14
   Accuracy = i
   ce = CommonEdge4D(sing,SauterSchwab3D._legendre(Accuracy,0.0,1.0))
  
   int_tp = sauterschwab_parameterized(INTEGRAND, ce)
  
   num_pts = 5*length(ce.qps)^4
   push!(n1,num_pts)
   
   push!(res_tp,int_tp)
end

for i in 2:1:7
   Accuracy = i
   ce_s = CommonEdge4D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._shunnham2D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)

   num_pts = length(ce_s.qps[1])^2*length(ce_s.qps[2]) + 4*length(ce_s.qps[2])^2
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

for i in 2:1:12
   Accuracy = i
   ce_gm = CommonEdge4D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._grundmannMoeller2D(2*Accuracy-1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, ce_gm)

   num_pts = length(ce_gm.qps[1])^2*length(ce_gm.qps[2]) + 4*length(ce_gm.qps[2])^2
   push!(n3,num_pts)
   push!(res_gm,int_gm)
end

err_tp = norm.(res_tp.-ref)/norm(ref)
err_sp = norm.(res_sp.-ref)/norm(ref)
err_gm = norm.(res_gm.-ref)/norm(ref)

println(err_tp)
println(err_sp)
println(err_gm)

using Plots
plot( yaxis=:log,xaxis=:log, fontfamily="Times")
plot!(n1,err_tp, label="Gauss Tensor-Product",markershape=:circle)
plot!(n2,err_sp, label="Simplex Tensor-Product",markershape=:rect)
#plot!(n3,err_gm, label="Simplex-Product GM",markershape=:x)
plot!(xlims=(3e1,1e5),ylims=(1e-8,1e-1))
plot!(xlabel="Quad. points/Func. evals.", ylabel="Rel. Error.", title="Common Edge 4D",legend=:bottomleft)


using BenchmarkTools

ref = -0.007051309436221841 + 0.005782096241048667im

ce = CommonEdge4D(sing,SauterSchwab3D._legendre(6,0.0,1.0))
  
int_tp = sauterschwab_parameterized(INTEGRAND, ce)

num_pts = 5*length(ce.qps)^4
println("#Pts: ",num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println("Rel Err: ",err_tp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, ce) end


ce_s = CommonEdge4D_S(sing,(SauterSchwab3D._legendre(6,0.0,1.0),
SauterSchwab3D._shunnham2D(6)))

int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)

num_pts = length(ce_s.qps[1])^2*length(ce_s.qps[2]) + 4*length(ce_s.qps[2])^2

println("#Pts: ",num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println("Rel Err: ",err_sp)
@time   for i in 1:100  sauterschwab_parameterized(INTEGRAND, ce_s) end

#=

ce_gm = CommonEdge4D_S((SauterSchwab3D._legendre(5,0.0,1.0),
SauterSchwab3D._grundmannMoeller2D(2*5-1)))

int_gm = sauterschwab_parameterized(INTEGRAND, ce_gm)

num_pts = length(ce_gm.qps[1])^2*length(ce_gm.qps[2]) + 4*length(ce_gm.qps[2])^2

println("#Pts: ",num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println("Rel Err: ",err_gm)
=#
