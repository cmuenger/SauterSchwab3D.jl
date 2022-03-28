using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D

const pI = point(3,4,2)
const pII = point(3,7,4)
const pIII = point(8,3,5)

const P = simplex(pI, pII, pIII)
const Q = simplex(pI, pII, pIII)

const sing = SauterSchwab3D.singularity_detection(P,Q)


Accuracy2 = 24
cf_ref = CommonFace4D(sing,SauterSchwab3D._legendre(Accuracy2,0.0,1.0))


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
ref = sauterschwab_parameterized(INTEGRAND, cf_ref)
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
   cf = CommonFace4D(sing,SauterSchwab3D._legendre(Accuracy,0.0,1.0))
  
   int_tp = sauterschwab_parameterized(INTEGRAND, cf)
  
   num_pts = 6*length(cf.qps)^4
   push!(n1,num_pts)
 
   push!(res_tp,int_tp)
end

for i in 2:1:7
   Accuracy = i
   cf_s = CommonFace4D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._shunnham3D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, cf_s)

   num_pts = 6*length(cf_s.qps[1])*length(cf_s.qps[2]) 
   push!(n2,num_pts)
   push!(res_sp,int_sp)


end

for i in 2:1:12
   Accuracy = i

   cf_gm = CommonFace4D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._grundmannMoeller3D(2*(Accuracy-1)+1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, cf_gm)

   num_pts = 6*length(cf_gm.qps[1])*length(cf_gm.qps[2]) 
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
plot!(xlims=(1e1,1e4),ylims=(1e-6,1))
plot!(xlabel="#Quad. pts/Func. evals", ylabel="Rel. Error.", title="Common Face 4D",legend=:bottomleft)



using BenchmarkTools

ref = 9.199951912100538 - 5.619050336211777im

cf = CommonFace4D(sing,SauterSchwab3D._legendre(5,0.0,1.0))
  
int_tp = sauterschwab_parameterized(INTEGRAND, cf)

num_pts = 6*length(cf.qps)^4
println("#Pts: ",num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println("Rel. Err: ",err_tp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cf) end

cf_s = CommonFace4D_S(sing,(SauterSchwab3D._legendre(6,0.0,1.0),
SauterSchwab3D._shunnham3D(6)))

int_sp = sauterschwab_parameterized(INTEGRAND, cf_s)

num_pts = 6*length(cf_s.qps[1])*length(cf_s.qps[2]) 

println("#Pts: ",num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println("Rel Err: ",err_sp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, cf_s) end

#=
cf_gm = CommonFace4D_S((SauterSchwab3D._legendre(6,0.0,1.0),
SauterSchwab3D._grundmannMoeller3D(2*6-1)))

int_gm = sauterschwab_parameterized(INTEGRAND, cf_gm)

num_pts = 6*length(cf_gm.qps[1])*length(cf_gm.qps[2]) 
println("#Pts: ",num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println("Rel Err: ",err_gm)
=#