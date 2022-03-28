using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D
using StaticArrays


const pI   = point(1,0,0)
const pII  = point(0,1,0)
const pIII = point(0,0,0)
const pIV  = point(0,0,1)

const qIV  = point(0,0,-1)

const P = simplex(pI,pII,pIV,pIII)
const Q = simplex(pI,pII,qIV,pIII)

const sing = SauterSchwab3D.singularity_detection(P,Q)

Accuracy2 = 18
cf_ref = CommonFace6D(sing,SauterSchwab3D._legendre(Accuracy2,0.0,1.0))


function integrand(x,y)
      return ((x-pI)'*(y-pIV))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
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
ref = sauterschwab_parameterized(INTEGRAND, cf_ref)
println(ref)
println()

res_tp =[]
res_sp =[]
res_gm =[]
n1 = []
n2 = []
n3 = []
for i in 2:1:15
   Accuracy = i
   cf = CommonFace6D(sing,SauterSchwab3D._legendre(Accuracy,0.0,1.0))

   int_tp = sauterschwab_parameterized(INTEGRAND, cf)

   num_pts = 15*length(cf.qps)^6
   push!(n1,num_pts)

   push!(res_tp,int_tp)

end

for i in 2:1:7
   Accuracy = i
   cf_s = CommonFace6D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._shunnham2D(Accuracy),
                          SauterSchwab3D._shunnham3D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, cf_s)

   num_pts = 3*length(cf_s.qps[1])*length(cf_s.qps[2])*length(cf_s.qps[3])+
            12*length(cf_s.qps[3])^2
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

for i in 2:1:12
   Accuracy = i
   cf_gm = CommonFace6D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._grundmannMoeller2D(2*Accuracy-1),
                          SauterSchwab3D._grundmannMoeller3D(2*Accuracy-1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, cf_gm)

   num_pts = 3*length(cf_gm.qps[1])*length(cf_gm.qps[2])*length(cf_gm.qps[3])+
            12*length(cf_gm.qps[3])^2
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
plot!(xlims=(1e2,1e7),ylims=(1e-8,1))
plot!(xlabel="#Quad. pts/Func. evals", ylabel="Rel. Error.", title="Common Face 6D",legend=:bottomleft)



using BenchmarkTools

ref = 0.0011204365593010864 - 0.000870131729342847im

cf = CommonFace6D(sing,SauterSchwab3D._legendre(4,0.0,1.0))

int_tp = sauterschwab_parameterized(INTEGRAND, cf)

num_pts = 15*length(cf.qps)^6
println("#Pts: ",num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println("Rel. Err: ",err_tp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cf) end

cf_s = CommonFace6D_S(sing,(SauterSchwab3D._legendre(5,0.0,1.0),
                         SauterSchwab3D._shunnham2D(5),
                         SauterSchwab3D._shunnham3D(5)))

int_sp = sauterschwab_parameterized(INTEGRAND, cf_s)

num_pts = 3*length(cf_s.qps[1])*length(cf_s.qps[2])*length(cf_s.qps[3])+
            12*length(cf_s.qps[3])^2
println("#Pts: ",num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println("Rel Err: ",err_sp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cf_s) end
#=
cf_gm = CommonFace6D_S(sing,(SauterSchwab3D._legendre(5,0.0,1.0),
                         SauterSchwab3D._grundmannMoeller2D(2*5-1),
                         SauterSchwab3D._grundmannMoeller3D(2*5-1)))

int_gm = sauterschwab_parameterized(INTEGRAND, cf_gm)

num_pts = 3*length(cf_gm.qps[1])*length(cf_gm.qps[2])*length(cf_gm.qps[3])+
            12*length(cf_gm.qps[3])^2
println("#Pts: ",num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println("Rel Err: ",err_gm)

=#