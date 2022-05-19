using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D
using StaticArrays

const pI   = point(0,0,0)
const pII  = point(1,0,0)
const pIII = point(0,1,0)
const pIV  = point(0,0,1)

const qIII = point(0,-1,0)
const qIV  = point(0,0,-1)

const P = simplex(pI,pII,pIII,pIV)
const Q = simplex(pI,pII,qIII,qIV)

const sing = SauterSchwab3D.singularity_detection(P,Q)

Accuracy2 = 15
ce_ref = CommonEdge6D(sing,SauterSchwab3D._legendre(Accuracy2,0.0,1.0))

function integrand(x,y)
      return ((x-pI)'*(y-qIII))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
end

function INTEGRAND(u,v)
   j1 = volume(P) * factorial(dimension(P))
   j2 = volume(Q) * factorial(dimension(Q))
   x = barytocart(P,u)
   y = barytocart(Q,v)
   output = integrand(x,y)*j1*j2
   
   return output
end

print("Ref: ")
ref = sauterschwab_parameterized(INTEGRAND, ce_ref)
println(ref)
println()

res_tp =[]
res_sp = []
res_gm = []
n1 = []
n2 = []
n3 = []
for i in 2:1:14
   Accuracy = i
   ce = CommonEdge6D(sing,SauterSchwab3D._legendre(Accuracy,0.0,1.0))
 
   int_tp = sauterschwab_parameterized(INTEGRAND, ce)
 
   num_pts = 5*length(ce.qps)^6
   push!(n1,num_pts)

   push!(res_tp,int_tp)
end

for i in 2:1:6
   Accuracy = i
   ce_s = CommonEdge6D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._shunnham2D(Accuracy),
                          SauterSchwab3D._shunnham3D(Accuracy),
                          SauterSchwab3D._shunnham4D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)

   num_pts = length(ce_s.qps[2])^3+
             2*length(ce_s.qps[1])*length(ce_s.qps[2])*length(ce_s.qps[3])+
             length(ce_s.qps[2])*length(ce_s.qps[4])+
             length(ce_s.qps[1])^2*length(ce_s.qps[2])^2

   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

for i in 2:1:12
   Accuracy = i
   ce_gm = CommonEdge6D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._grundmannMoeller2D(2*Accuracy-1),
                          SauterSchwab3D._grundmannMoeller3D(2*Accuracy-1),
                          SauterSchwab3D._grundmannMoeller4D(2*Accuracy-1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, ce_gm)

   num_pts = length(ce_gm.qps[2])^3+
             2*length(ce_gm.qps[1])*length(ce_gm.qps[2])*length(ce_gm.qps[3])+
             length(ce_gm.qps[2])*length(ce_gm.qps[4])+
             length(ce_gm.qps[1])^2*length(ce_gm.qps[2])^2

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
plot!(xlims=(1e2,5e7),ylims=(1e-11,1))
plot!(xlabel="#Quad. pts/Func. evals", ylabel="Rel. Error.", title="Common Edge 6D",legend=:bottomleft)



using BenchmarkTools

ref =  0.00041274116231997213 - 0.0003417955682495072im

ce = CommonEdge6D(sing,SauterSchwab3D._legendre(5,0.0,1.0))

int_tp = sauterschwab_parameterized(INTEGRAND, ce)

num_pts = 5*length(ce.qps)^6
println("#Pts: ",num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println("Rel Err: ",err_tp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, ce) end

ce_s = CommonEdge6D_S(sing,(SauterSchwab3D._legendre(6,0.0,1.0),
                     SauterSchwab3D._shunnham2D(6),
                     SauterSchwab3D._shunnham3D(6),
                     SauterSchwab3D._shunnham4D(6)))

int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)

num_pts = length(ce_s.qps[2])^3+
2*length(ce_s.qps[1])*length(ce_s.qps[2])*length(ce_s.qps[3])+
length(ce_s.qps[2])*length(ce_s.qps[4])+
length(ce_s.qps[1])^2*length(ce_s.qps[2])^2
println("#Pts: ",num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println("Rel Err: ",err_sp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, ce_s) end

#=
ce_gm = CommonEdge6D_S((SauterSchwab3D._legendre(6,0.0,1.0),
                        SauterSchwab3D._grundmannMoeller2D(2*6-1),
                        SauterSchwab3D._grundmannMoeller3D(2*6-1),
                        SauterSchwab3D._grundmannMoeller4D(2*6-1)))
   
int_gm = sauterschwab_parameterized(INTEGRAND, ce_gm)

num_pts = length(ce_gm.qps[2])^3+
          2*length(ce_gm.qps[1])*length(ce_gm.qps[2])*length(ce_gm.qps[3])+
          length(ce_gm.qps[2])*length(ce_gm.qps[4])+
          length(ce_gm.qps[1])^2*length(ce_gm.qps[2])^2
println("#Pts: ",num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println("Rel Err: ",err_gm)
=#