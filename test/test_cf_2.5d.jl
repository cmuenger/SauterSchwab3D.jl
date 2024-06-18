using LinearAlgebra
using CompScienceMeshes
using SauterSchwab3D
using StaticArrays
using Combinatorics
using Test


include("$(pkgdir(SauterSchwab3D))/example/singularity_detection.jl")

function is_CSM_tet(s)
   CSM_tet = false

   c = cross(s.tangents[1], s.tangents[2])
   dot(c,s.tangents[3]) > 0.0 && (CSM_tet = true)

   return CSM_tet
end




pI   = point(0,0,0)
pII  = point(1,0,0)
pIII = point(0,1,0)
pIV  = point(0,0,1)


P_ = simplex(pI,pII,pIV,pIII)
Q_ = simplex(pI,pIII,pII)

@test is_CSM_tet(P_) == true
c_tri = cartesian(CompScienceMeshes.center(Q_))
d = dot(Q_.normals[1], pI + pII + pIII + pIV - c_tri)
@test d < 0.0

sing_ = singularity_detection(P_,Q_)
cf_ref = CommonFace5D(sing_,SauterSchwab3D._legendre(12,0.0,1.0))


# Tensor Product: Kernel=1 Test
function integrand1(x,y)
   return 1.0
end
function INTEGRAND1(u,v) 
   j1 = volume(P_) * factorial(dimension(P_))
   j2 = volume(Q_) * factorial(dimension(Q_))
   x = barytocart(P_,u)
   y = barytocart(Q_,v)
   output = integrand1(x,y)*j1*j2
   
   return output
end
num1 = sauterschwab_parameterized(INTEGRAND1, cf_ref)
exact1 = Q_.volume*P_.volume
err1 = abs(num1-exact1)/abs(exact1)
@test err1 < 1.0e-12


# Tensor Product: Regular Kernel Test
function integrand2(x,y)
   return 1 +
   2*x[1]-3*x[2]+4*x[3]-5*y[1]+6*y[2]-7*y[3] +
   2*x[1]^2-3*x[2]^3+4*x[3]^4-5*y[1]^5+6*y[2]^6-7*y[3]^7 +
   y[1]*y[2]*y[3]*x[1]*x[2]*x[3]
end
function INTEGRAND2(u,v) 
   j1 = volume(P_) * factorial(dimension(P_))
   j2 = volume(Q_) * factorial(dimension(Q_))
   x = barytocart(P_,u)
   y = barytocart(Q_,v)
   output = integrand2(x,y)*j1*j2
   
   return output
end
num2 = sauterschwab_parameterized(INTEGRAND2, cf_ref)
exact2 = 467/2520
err2 = abs(num2-exact2)/abs(exact2)
@test err2 < 1.0e-12


# Tensor Product: Singular Kernel Reference (for comparison with Simplex Product below)
function integrand(x,y)
      return ((x-pIV)'*(y-pI))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
end
function INTEGRAND(u,v) 
   j1 = volume(P_) * factorial(dimension(P_))
   j2 = volume(Q_) * factorial(dimension(Q_))
   x = barytocart(P_,u)
   y = barytocart(Q_,v)
   output = integrand(x,y)*j1*j2
   
   return output
end
ref = sauterschwab_parameterized(INTEGRAND, cf_ref)




## Simplex Product vs. Tensor Product reference for all permutations
list_tvert = collect(permutations([pI,pII,pIII,pIV]))
list_svert = collect(permutations([pI,pII,pIII]))

for tvert in list_tvert
   for svert in list_svert
      P = simplex(tvert[1],tvert[2],tvert[3],tvert[4])
      Q = simplex(svert[1],svert[2],svert[3])
      
      sing = singularity_detection(P,Q)
      I,J = reorder(sing)

      P_new = simplex(P.vertices[I])
      Q_new = simplex(Q.vertices[J])

      # expected reorder result:
      # Pt = [P1, P2, A1, P3]
      # Ps = [P1, P3, P2]
      @test norm(P_new[1] - Q_new[1]) < 1.0e-14
      @test norm(P_new[2] - Q_new[3]) < 1.0e-14
      @test norm(P_new[4] - Q_new[2]) < 1.0e-14

      if is_CSM_tet(P) == true
         c_tri = cartesian(CompScienceMeshes.center(Q))
         d = dot(Q.normals[1], pI + pII + pIII + pIV - c_tri)

         if d < 0.0 # means triangle is compatible with CSM tetrahedron (direction of n̂ ...)
            @test is_CSM_tet(P_new) == true
            c_new_tri = cartesian(CompScienceMeshes.center(Q_new))
            d_new = dot(Q_new.normals[1], pI + pII + pIII + pIV - c_new_tri)
            @test d_new < 0.0
            @test dot(Q_new.normals[1], Q.normals[1]) > 0.0

            function INTEGRAND(u,v) 
               j1 = volume(P_new) * factorial(dimension(P_new))
               j2 = volume(Q_new) * factorial(dimension(Q_new))
               x = barytocart(P_new,u)
               y = barytocart(Q_new,v)
               output = integrand(x,y)*j1*j2
               
               return output
            end
   
            cf_s = CommonFace5D_S(sing,(SauterSchwab3D._legendre(6,0.0,1.0),
                                       SauterSchwab3D._shunnham2D(6),
                                       SauterSchwab3D._shunnham3D(6)))
            int_sp = sauterschwab_parameterized(INTEGRAND, cf_s)
   
            err_sp = norm(int_sp-ref)/norm(ref)
            @test err_sp < 1.1e-6
         end
      else
         @test is_CSM_tet(P_new) == false
         #@test dot(Q_new.normals[1], Q.normals[1]) > 0.0 # does not work - not needed for CSM
      end
   end
end
print("CommonFace5D_S test finished")