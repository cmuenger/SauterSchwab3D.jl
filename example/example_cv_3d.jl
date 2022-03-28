using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D

const pI   = point(6,5,3)
const pII  = point(5,2,3)
const pIII = point(7,1,0)
const pIV  = point(3,4,0)

const qII  = point(3,8,4)
const qIII = point(7,5,6)
const qIV  = point(9,8,5)

const P = simplex(pI,pII,pIII,pIV)
const Q = simplex(pI,qII,qIII,qIV)

const sing = SauterSchwab3D.singularity_detection(P,Q)

Accuracy2 = 16
cv_ref = CommonVertex6D(sing,SauterSchwab3D._legendre(Accuracy2,0.0,1.0))

function integrand(x,y)
   ((x-pI)'*(y-pIV))*return(exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
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

res_tp =[]
res_sp =[]
res_gm = []
n1 = []
n2 = []
n3 = []

#Grauss tensor product
for i in 2:1:14
   Accuracy = i
   cv = CommonVertex6D(sing,SauterSchwab3D._legendre(Accuracy,0.0,1.0))

   int_tp = sauterschwab_parameterized(INTEGRAND, cv)

   num_pts = 2*length(cv.qps)^6
   push!(n1,num_pts)

   push!(res_tp,int_tp)
end

#ShunnHam simplex product
for i in 2:1:7
   Accuracy = i
   cv_s = CommonVertex6D_S(sing,SauterSchwab3D._shunnham3D(Accuracy))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, cv_s)

   println(length(cv_s.qps))
   num_pts = 2*length(cv_s.qps)^2
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

#GrundmannMoeller simplex product
for i in 2:1:12
   Accuracy = i
   cv_gm = CommonVertex6D_S(sing,SauterSchwab3D._grundmannMoeller3D(2*Accuracy-1))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, cv_gm)

   println(length(cv_gm.qps))
   num_pts = 2*length(cv_gm.qps)^2
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
plot!(xlims=(1e1,10^6),ylims=(1e-7,1))
plot!(xlabel="#Quad. points/Func. evals.", ylabel="Rel. Error.", title="Common Vertex 6D",legend=:bottomleft)



using BenchmarkTools

ref = 0.003477333267758975 - 0.5997621570461863im

cv = CommonVertex6D(sing,SauterSchwab3D._legendre(6,0.0,1.0))

int_tp = sauterschwab_parameterized(INTEGRAND, cv)

num_pts = 2*length(cv.qps)^6
println(num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println(err_tp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cv) end

cv_s = CommonVertex6D_S(sing,SauterSchwab3D._shunnham3D(6))

int_sp = sauterschwab_parameterized(INTEGRAND, cv_s)

num_pts = num_pts = 2*length(cv_s.qps)^2
println(num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println(err_sp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cv_s) end

#=
cv_gm = CommonVertex6D_S(SauterSchwab3D._grundmannMoeller3D(2*6-1))
   
int_gm = sauterschwab_parameterized(INTEGRAND, cv_gm)

num_pts = 2*length(cv_gm.qps)^2
println(num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println(err_gm)
=#