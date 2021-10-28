using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D

const pI   = point(1,5,3)
const pII  = point(2,5,3)
const pIII = point(7,1,0)
const pIV  = point(3,4,0)

const qI  = point(10,11,12)
const qII  = point(10,11,13)
const qIII = point(11,11,12)

const P = simplex(pI,pII,pIII,pIV)
const Q = simplex(qI,qII,qIII)

Accuracy2 = 16
pd_ref = PositiveDistance5D(SauterSchwab3D._legendre(Accuracy2,0.0,1.0))

function integrand(x,y)
			return(exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
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
ref = sauterschwab_parameterized(INTEGRAND, pd_ref)
println(ref)
println()

res_tp =[]
res_sp =[]
res_gm =[]
n1 = []
n2 = []
n3 = []

for i in 2:1:14
   Accuracy = i
   pd = PositiveDistance5D(SauterSchwab3D._legendre(Accuracy,0.0,1.0))
   
   int_tp = sauterschwab_parameterized(INTEGRAND, pd)

   num_pts = length(pd.qps)^5
   push!(n1,num_pts)
   push!(res_tp,int_tp)
end

for i in 2:1:7
   Accuracy = i
   pd_s = PositiveDistance5D_S((SauterSchwab3D._shunnham3D(Accuracy),
                                SauterSchwab3D._shunnham2D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, pd_s)

   num_pts = length(pd_s.qps[1])*length(pd_s.qps[2])
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

for i in 2:1:12
   Accuracy = i
   pd_gm = PositiveDistance5D_S((SauterSchwab3D._grundmannMoeller3D(2*Accuracy-1),
                                SauterSchwab3D._grundmannMoeller2D(2*Accuracy-1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, pd_gm)

   num_pts = length(pd_gm.qps[1])*length(pd_gm.qps[2])
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
plot( yaxis=:log, xaxis=:log)
plot!(n1,err_tp, label="Gauss Tensor-Product",markershape=:x)
plot!(n2,err_sp, label="Simplex-Product SH",markershape=:x)
plot!(n3,err_gm, label="Simplex-Product GM",markershape=:x)
plot!(xlims=(1,10^7),ylims=(1e-15,1))
plot!(xlabel="Quad. points/Func. evals.", ylabel="Rel. Error.", title="Positive Distance 5D")

#=
using BenchmarkTools

ref = -0.0024685313290508833 - 0.0026636388782296638im


pd = PositiveDistance5D(SauterSchwab3D._legendre(4,0.0,1.0))
   
int_tp = sauterschwab_parameterized(INTEGRAND, pd)

num_pts = length(pd.qps)^5
println(num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println(err_tp)

pd_s = PositiveDistance5D_S((SauterSchwab3D._shunnham3D(4),
SauterSchwab3D._shunnham2D(4)))

int_sp = sauterschwab_parameterized(INTEGRAND, pd_s)

num_pts = length(pd_s.qps[1])*length(pd_s.qps[2])
println(num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println(err_sp)

pd_gm = PositiveDistance5D_S((SauterSchwab3D._grundmannMoeller3D(2*4-1),
SauterSchwab3D._grundmannMoeller2D(2*4-1)))

int_gm = sauterschwab_parameterized(INTEGRAND, pd_gm)

num_pts = length(pd_gm.qps[1])*length(pd_gm.qps[2])
println(num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println(err_gm)
=#