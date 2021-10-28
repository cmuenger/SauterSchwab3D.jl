using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D
using StaticArrays


const pI   = point(6,5,3)
const pII  = point(5,2,3)
const pIII = point(7,1,0)
const pIV  = point(3,4,0)

const qIII = point(7,5,6)

const P = simplex(pI,pII,pIV,pIII)
const Q = simplex(pI,qIII,pII)


Accuracy2 = 20
ce_ref = CommonEdge5D(SauterSchwab3D._legendre(Accuracy2,0.0,1.0))


function integrand(x,y)
      return ((x-pI)'*(y-pIV))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
end


function INTEGRAND(u,v)
   #check if point is inside the tetrahedron
   if u[3]< 0.0 || u[3]>1-u[1]-u[2] || u[2]< 0.0 || u[2]>1-u[1] || u[1]<0.0 || u[1]>1.0
      println("u:",u)
   end
   if  v[2]< 0.0 || v[2]>1-v[1] || v[1]<0.0 || v[1]>1.0
      println("v:",v)
   end
   
   
   n1 = neighborhood(P,u)
   n2 = neighborhood(Q,v)
   x = cartesian(n1)
   y = cartesian(n2)
   output = integrand(x,y)*jacobian(n1)*jacobian(n2)
   
   return output
end

#=
print("Ref: ")
ref = sauterschwab_parameterized(INTEGRAND, ce_ref)
println(ref)
#print("Ref pd: ")
#ref_pd = sauterschwab_parameterized(INTEGRAND, pd_ref)
#println(ref_pd)
println()


res_tp =[]
res_sp =[]
res_gm =[]
n1 = []
n2 = []
n3 = []
for i in 2:1:17
   Accuracy = i
   ce = CommonEdge5D(SauterSchwab3D._legendre(Accuracy,0.0,1.0))

   int_tp = sauterschwab_parameterized(INTEGRAND, ce)

  
   num_pts = 5*length(ce.qps)^5
   push!(n1,num_pts)
  
   push!(res_tp,int_tp)
  
end

for i in 2:1:7
   Accuracy = i
   ce_s = CommonEdge5D_S((SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._shunnham2D(Accuracy),
                          SauterSchwab3D._shunnham3D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)

   num_pts = 3*length(ce_s.qps[1])*length(ce_s.qps[2])^2 + 2*length(ce_s.qps[2])*length(ce_s.qps[3])
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

for i in 2:1:12
   Accuracy = i
   ce_gm = CommonEdge5D_S((SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._grundmannMoeller2D(2*Accuracy-1),
                          SauterSchwab3D._grundmannMoeller3D(2*Accuracy-1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, ce_gm)

   num_pts = 3*length(ce_gm.qps[1])*length(ce_gm.qps[2])^2 + 2*length(ce_gm.qps[2])*length(ce_gm.qps[3])
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
plot!(xlims=(1,1e8),ylims=(1e-15,1))
plot!(xlabel="#Quad. pts/Func. evals", ylabel="Rel. Error.", title="Common Edge 5D",legend=:bottomright)

=#

using BenchmarkTools

ref = 4.263868529192273 + 0.1481859134048196im

ce = CommonEdge5D(SauterSchwab3D._legendre(6,0.0,1.0))

int_tp = sauterschwab_parameterized(INTEGRAND, ce)

  
num_pts = 5*length(ce.qps)^5
println("#Pts: ",num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println("Rel Err: ",err_tp)

ce_s = CommonEdge5D_S((SauterSchwab3D._legendre(7,0.0,1.0),
                          SauterSchwab3D._shunnham2D(7),
                          SauterSchwab3D._shunnham3D(7)))
   
int_sp = sauterschwab_parameterized(INTEGRAND, ce_s)
num_pts = 3*length(ce_s.qps[1])*length(ce_s.qps[2])^2 + 2*length(ce_s.qps[2])*length(ce_s.qps[3])

println("#Pts: ",num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println("Rel Err: ",err_sp)

 ce_gm = CommonEdge5D_S((SauterSchwab3D._legendre(7,0.0,1.0),
                          SauterSchwab3D._grundmannMoeller2D(2*7-1),
                          SauterSchwab3D._grundmannMoeller3D(2*7-1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, ce_s)

   num_pts = 3*length(ce_gm.qps[1])*length(ce_gm.qps[2])^2 + 2*length(ce_gm.qps[2])*length(ce_gm.qps[3])
println("#Pts: ",num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println("Rel Err: ",err_gm)
