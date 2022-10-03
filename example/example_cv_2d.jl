using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D

include("singularity_detection.jl")

#=
const pI = point(1,5,3)
const pII = point(2,5,3)
const pIII = point(7,1,0)
const pIV = point(5,1,-3)
const pV = point(0,0,0)
=#

const pI   = point(0,0,0)
const pII  = point(1,0,0)
const pIII = point(0,1,0)


const qII  = point(-1,0,0)
const qIII = point(0,-1,0)


const P = simplex(pI,pII,pIII)
const Q = simplex(pI,qII,qIII)

const sing = singularity_detection(P,Q)

Accuracy2 = 15
cv_ref = CommonVertex4D(sing,SauterSchwab3D._legendre(Accuracy2,0.0,1.0))

function integrand(x,y)
      return ((x-pII)'*(y-qIII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
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
ref = sauterschwab_parameterized(INTEGRAND, cv_ref)
println(ref)
println()

res_tp =ComplexF64[]
res_sp =ComplexF64[]
res_gm =ComplexF64[]
n1 = Int[]
n2 = Int[]
n3 = Int[]
#Grauss tensor product
for i in 2:1:14
   Accuracy = i
   cv = CommonVertex4D(sing,SauterSchwab3D._legendre(Accuracy,0.0,1.0))

   int_tp = sauterschwab_parameterized(INTEGRAND, cv)

   num_pts = length(cv.qps)^4
   push!(n1,num_pts)

   
   push!(res_tp,int_tp)
end

#ShunnHam simplex product
for i in 2:1:7
   Accuracy = i
   cv_s = CommonVertex4D_S(sing,SauterSchwab3D._shunnham2D(Accuracy))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, cv_s)

   num_pts = length(cv_s.qps)^2
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

#GrundmannMoeller simplex product
for i in 2:1:12
   Accuracy = i
   cv_gm = CommonVertex4D_S(sing,SauterSchwab3D._grundmannMoeller2D(2*Accuracy-1))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, cv_gm)

   num_pts = length(cv_gm.qps)^2
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
plot( yaxis=:log, xaxis=:log, fontfamily="Times")
plot!(n1,err_tp, label="Gauss Tensor-Product",markershape=:circle)
plot!(n2,err_sp, label="Simplex Tensor-Product",markershape=:rect)
#plot!(n3,err_gm, label="Simplex-Product GM",markershape=:x)
plot!(xlims=(5,10^4),ylims=(1e-9,1))
plot!(xlabel="Quad. points/Func. evals.", ylabel="Rel. Error.", title="Common Vertex 4D",legend=:bottomleft)



using BenchmarkTools

ref = 0.004958861812742832 - 0.007274636107126924im

cv = CommonVertex4D(sing,SauterSchwab3D._legendre(6,0.0,1.0))

int_tp = sauterschwab_parameterized(INTEGRAND, cv)

num_pts = length(cv.qps)^4
println(num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println(err_tp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cv) end


cv_s = CommonVertex4D_S(sing,SauterSchwab3D._shunnham2D(6))
   
int_sp = sauterschwab_parameterized(INTEGRAND, cv_s)

num_pts = length(cv_s.qps)^2
println(num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println(err_sp)
 @time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cv_s) end

#=

cv_gm = CommonVertex4D_S(SauterSchwab3D._grundmannMoeller2D(2*6-1))
   
int_gm = sauterschwab_parameterized(INTEGRAND, cv_gm)

num_pts = length(cv_gm.qps)^2
println(num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println(err_gm)
=#