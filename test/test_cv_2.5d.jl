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
pIII = point(1,1,0)
pIV  = point(1,1,1)

qIII = point(0,-1,0)
qII = point(-1,-1,0)


P_ = simplex(pI,pII,pIV,pIII)
Q_ = simplex(pI,qIII,qII)

@test is_CSM_tet(P_) == true
c_tri = cartesian(CompScienceMeshes.center(Q_))
d = dot(Q_.normals[1], pI + pII + pIII + pIV - c_tri)
@test d < 0.0

sing_ = singularity_detection(P_,Q_)
cv_ref = CommonVertex5D(sing_,SauterSchwab3D._legendre(12,0.0,1.0))


# Tensor Product: Kernel=1 Test
function integrand1(x,y)
   return 1.0
end
function INTEGRAND1(u,v)
   n1 = neighborhood(P_,u)
   n2 = neighborhood(Q_,v)
   x = cartesian(n1)
   y = cartesian(n2)
   output = integrand1(x,y)*jacobian(n1)*jacobian(n2)
   return(output)
end
num1 = sauterschwab_parameterized(INTEGRAND1, cv_ref)
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
   n1 = neighborhood(P_,u)
   n2 = neighborhood(Q_,v)
   x = cartesian(n1)
   y = cartesian(n2)
   output = integrand2(x,y)*jacobian(n1)*jacobian(n2)
   return(output)
end
num2 = sauterschwab_parameterized(INTEGRAND2, cv_ref)
exact2 = 89/504
err2 = abs(num2-exact2)/abs(exact2)
@test err2 < 1.0e-12


# Tensor Product: Singular Kernel Reference (for comparison with Simplex Product below)
function integrand(x,y)
    return ((x-pIV)'*(y-pI))*exp(-im*1*norm(x-y))/(4pi*norm(x-y))
end
function INTEGRAND(u,v)
   n1 = neighborhood(P_,u)
   n2 = neighborhood(Q_,v)
   x = cartesian(n1)
   y = cartesian(n2)
   output = integrand(x,y)*jacobian(n1)*jacobian(n2)
   return(output)
end
ref = sauterschwab_parameterized(INTEGRAND, cv_ref)




## Simplex Product vs. Tensor Product reference for all permutations
list_tvert = collect(permutations([pI,pII,pIII,pIV]))
list_svert = collect(permutations([pI,qII,qIII]))

for tvert in list_tvert
   for svert in list_svert
      P = simplex(tvert[1],tvert[2],tvert[3],tvert[4])
      Q = simplex(svert[1],svert[2],svert[3])
      
      sing = singularity_detection(P,Q)
      I,J = reorder(sing)

      P_new = simplex(P.vertices[I])
      Q_new = simplex(Q.vertices[J])


      # expected reorder result:
      # Pt = [P, A1, A2, A3]
      # Ps = [P, B1, B2]
      @test norm(P_new[1] - Q_new[1]) < 1.0e-14

      if is_CSM_tet(P) == true
         @test is_CSM_tet(P_new) == true
         @test dot(Q_new.normals[1], Q.normals[1]) > 0.0

         function INTEGRAND(u,v)
            n1 = neighborhood(P,u)
            n2 = neighborhood(Q,v)
            x = cartesian(n1)
            y = cartesian(n2)
            output = integrand(x,y)*jacobian(n1)*jacobian(n2)
            return(output)
         end
         cv_s = CommonVertex5D_S(sing,(SauterSchwab3D._shunnham3D(6),
                                       SauterSchwab3D._shunnham2D(6)))
   
         int_sp = sauterschwab_parameterized(INTEGRAND, cv_s)

         err_sp = norm(int_sp-ref)/norm(ref)
         @test err_sp < 2.0e-5     
      else
         @test is_CSM_tet(P_new) == false
         @test dot(Q_new.normals[1], Q.normals[1]) > 0.0 
      end
   end
end
print("CommonVertex5D_S test finished")