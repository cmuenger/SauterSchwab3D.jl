using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D

const pI = point(1,0,0)
const pII = point(0,0,0)
const pIII = point(0.5,0.7,0)
const pIV = point(0.8,-0.3,0)

const P = simplex(pI,pIII,pII)
const Q = simplex(pI,pIV,pII)

Accuracy2 = 30
ce_ref = CommonEdge4D(SauterSchwab3D._legendre(Accuracy2,0.0,1.0))

function integrand(x,y)
      return ((x-pI)'*(y-pII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
end

function INTEGRAND(u,v)
   if u[2]< 0.0 || u[2]>1-u[1] || u[1]<0.0 || u[1]>1.0
      println(u)
   end
   if v[2]< 0.0 || v[2]>1-v[1] || v[1]<0.0 || v[1]>1.0
      println(v)
   end
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

res_tp =[]
res_sp =[]
res_gm = []
n1 = []
n2 = []
n3 = []
for i in 2:1:20
   Accuracy = i
   ce = CommonEdge4D(SauterSchwab3D._legendre(Accuracy,0.0,1.0))
  
   int_tp = sauterschwab_parameterized(INTEGRAND, ce)
  
   num_pts = 5*length(ce.qps)^4
   push!(n1,num_pts)
   
   push!(res_tp,int_tp)
end

for i in 2:1:7
   Accuracy = i
   ce_s = CommonEdge4D_S((SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._shunnham2D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)

   num_pts = length(ce_s.qps[1])^2*length(ce_s.qps[2]) + 4*length(ce_s.qps[2])^2
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

for i in 2:1:12
   Accuracy = i
   ce_gm = CommonEdge4D_S((SauterSchwab3D._legendre(Accuracy,0.0,1.0),
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
plot( yaxis=:log, xaxis=:log)
plot!(n1,err_tp, label="Gauss Tensor-Product",markershape=:x)
plot!(n2,err_sp, label="Simplex-Product SH",markershape=:x)
plot!(n3,err_gm, label="Simplex-Product GM",markershape=:x)
plot!(xlims=(1,10^7),ylims=(1e-15,1))
plot!(xlabel="Quad. points/Func. evals.", ylabel="Rel. Error.", title="Common Edge 4D")

#=
using BenchmarkTools

ref = -0.0028991793711368725 + 0.0012868587156426375im

ce = CommonEdge4D(SauterSchwab3D._legendre(5,0.0,1.0))
  
int_tp = sauterschwab_parameterized(INTEGRAND, ce)

num_pts = 5*length(ce.qps)^4
println("#Pts: ",num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println("Rel Err: ",err_tp)

ce_s = CommonEdge4D_S((SauterSchwab3D._legendre(5,0.0,1.0),
SauterSchwab3D._shunnham2D(5)))

int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)

num_pts = length(ce_s.qps[1])^2*length(ce_s.qps[2]) + 4*length(ce_s.qps[2])^2

println("#Pts: ",num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println("Rel Err: ",err_sp)

ce_gm = CommonEdge4D_S((SauterSchwab3D._legendre(5,0.0,1.0),
SauterSchwab3D._grundmannMoeller2D(2*5-1)))

int_gm = sauterschwab_parameterized(INTEGRAND, ce_gm)

num_pts = length(ce_gm.qps[1])^2*length(ce_gm.qps[2]) + 4*length(ce_gm.qps[2])^2

println("#Pts: ",num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println("Rel Err: ",err_gm)
=#
