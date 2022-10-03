using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D

#=
const pI   = point(6,5,3)
const pII  = point(5,2,4)
const pIII = point(7,1,0)
const pIV  = point(3,4,0)
=#

const pI   = point(0,0,0)
const pII  = point(1,0,0)
const pIII = point(0,1,0)
const pIV  = point(0,0,1)

const P = simplex(pI,pII,pIII,pIV)
const Q =  simplex(pI,pII,pIII,pIV)

const sing = SauterSchwab3D.singularity_detection(P,Q)


Accuracy2 = 15
ct_ref = CommonVolume6D(sing,SauterSchwab3D._legendre(Accuracy2,0.0,1.0))

function integrand(x,y)
   return ((x-pI)'*(y-pIII))*(exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
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
ref = sauterschwab_parameterized(INTEGRAND, ct_ref)
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
   ct = CommonVolume6D(sing,SauterSchwab3D._legendre(Accuracy,0.0,1.0))

   int_tp = sauterschwab_parameterized(INTEGRAND, ct)
  
   num_pts = 18*length(ct.qps)^6
   push!(n1,num_pts)

   push!(res_tp,int_tp)
end

for i in 2:1:6
   Accuracy = i
   ct_s = CommonVolume6D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._shunnham2D(Accuracy),
                          SauterSchwab3D._shunnham4D(Accuracy)))
   
   int_sp = sauterschwab_parameterized(INTEGRAND, ct_s)

   num_pts = 2*length(ct_s.qps[1])^2*length(ct_s.qps[3])+
            16*length(ct_s.qps[3])*length(ct_s.qps[2])
   push!(n2,num_pts)
   push!(res_sp,int_sp)
end

for i in 2:1:12
   Accuracy = i
   ct_gm = CommonVolume6D_S(sing,(SauterSchwab3D._legendre(Accuracy,0.0,1.0),
                          SauterSchwab3D._grundmannMoeller2D(2*Accuracy-1),
                          SauterSchwab3D._grundmannMoeller4D(2*Accuracy-1)))
   
   int_gm = sauterschwab_parameterized(INTEGRAND, ct_gm)

   num_pts = 2*length(ct_gm.qps[1])^2*length(ct_gm.qps[3])+
            16*length(ct_gm.qps[3])*length(ct_gm.qps[2])
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
plot( yaxis=:log, xaxis=:log) #, fontfamily="Times" )
plot!(n1,err_tp, label="Gauss Tensor-Product",markershape=:circle)
plot!(n2,err_sp, label="Simplex Tensor-Product",markershape=:rect, color=:green)
#plot!(n3,err_gm, label="GrundmannMoeller Tensor-Product",markershape=:x)
plot!(xlims=(1e2,1e7),ylims=(1e-8,1))
plot!(xlabel="Number of Quadrature points", ylabel="Rel. Error.", title="Identical Tetrahedron 6D",legend=:bottomleft)

savefig("CommonVolume6D.pdf")

using BenchmarkTools

#reference solution
ref= -0.00020240321624519095 + 0.0001297749782688385im

ct = CommonVolume6D(sing,SauterSchwab3D._legendre(5,0.0,1.0))

int_tp = sauterschwab_parameterized(INTEGRAND, ct)

num_pts = 18*length(ct.qps)^6
println("#Pts: ",num_pts)
err_tp = norm.(int_tp.-ref)/norm(ref)
println("Rell Err; ",err_tp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, ct) end


ct_s = CommonVolume6D_S(sing,(SauterSchwab3D._legendre(5,0.0,1.0),
                         SauterSchwab3D._shunnham2D(5),
                         SauterSchwab3D._shunnham4D(5)))

int_sp = sauterschwab_parameterized(INTEGRAND, ct_s)

num_pts = 2*length(ct_s.qps[1])^2*length(ct_s.qps[3])+
          16*length(ct_s.qps[3])*length(ct_s.qps[2])
println("#Pts: ",num_pts)
err_sp = norm.(int_sp.-ref)/norm(ref)
println("Rel Err: ",err_sp)
@time for i in 1:100 sauterschwab_parameterized(INTEGRAND, ct_s) end


#=
ct_gm = CommonVolume6D_S(sing,(SauterSchwab3D._legendre(7,0.0,1.0),
                         SauterSchwab3D._grundmannMoeller2D(2*7-1),
                         SauterSchwab3D._grundmannMoeller4D(2*7-1)))

int_gm = sauterschwab_parameterized(INTEGRAND, ct_gm)

num_pts = 2*length(ct_gm.qps[1])^2*length(ct_gm.qps[3])+
          16*length(ct_gm.qps[3])*length(ct_gm.qps[2])
println("#Pts: ",num_pts)
err_gm = norm.(int_gm.-ref)/norm(ref)
println("Rel Err: ",err_gm)
=#