using SauterSchwab3D: tetrahedron_face_order, tetrahedron_edge_order
using CompScienceMeshes
using SauterSchwab3D
using StaticArrays
using LinearAlgebra


const pI   = point(6,5,3)
const pII  = point(5,2,3)
const pIII = point(7,1,0)
const pIV  = point(3,4,0)

const qII  = point(3,8,4)
const qIII = point(7,5,6)
const qIV  = point(9,8,5)

const P_v = simplex(pI,pII,pIII)
const Q_v = simplex(qII,pI,qIII)


#=
s = singularity_detection(P_v,Q_v)
println(s)

I,J = reorder(s,cv)

K = tetrahedron_face_order(I)
L = tetrahedron_face_order(J)

println(I,J)
println(K,L)
=#

const P_e = simplex(pI,pII,pIII,pIV)
const Q_e = simplex(qIII,pI,pII,qII)

s = singularity_detection(P_e,Q_e)
println(s)

I,J = reorder(s)

K = tetrahedron_face_order(I)
L = tetrahedron_face_order(J)

println(I,J)
println(K,L)
