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
cv = CommonVertex4D(SauterSchwab3D._legendre(1,0.0,1.0))

s = SauterSchwab3D.singularity_detection(P_v,Q_v)
println(s)

I,J = reorder(s,cv)

K = tetrahedron_face_order(I)
L = tetrahedron_face_order(J)

I2,J2,K2,L2 = SauterSchwab3D.reorder(P_v,Q_v,cv)
println(I,J)
println(K,L,"\n")
println(I2,J2)
println(K2,L2)

=#

const P_e = simplex(pI,pII,pIII,pIV)
const Q_e = simplex(qIII,pI,pII,qII)
#ce = CommonEdge4D(SauterSchwab3D._legendre(1,0.0,1.0))

s = SauterSchwab3D.singularity_detection(P_e,Q_e)
println(s)

I,J = reorder(s)

K = tetrahedron_face_order(I)
L = tetrahedron_face_order(J)

println(I,J)
println(K,L)




#println(reverse_face_order(J))

#println(reverse_edge_order(J))

#=
println()



const P_f = simplex(pI,pII,pIII,pIV)
const Q_f = simplex(pIII,qIII,pI,pII)

cf = CommonFace6D(SauterSchwab3D._legendre(1,0.0,1.0))

s = SauterSchwab3D.singularity_detection(P_f,Q_f)
println(s)

I,J = reorder(s,cf)

I2,J2 = SauterSchwab3D.reorder2(s,cf)
println(I,J)
println(I2,J2)



println()



const P_e = simplex(pI,pII,pIII,pIV)
const Q_e = simplex(pII,qIII,qIV,pI)

ce = CommonEdge6D(SauterSchwab3D._legendre(1,0.0,1.0))

s = SauterSchwab3D.singularity_detection(P_e,Q_e)
println(s)

I,J = reorder(s,ce)

I2,J2 = SauterSchwab3D.reorder2(s,ce)
println(I,J)
println(I2,J2)
=#

#=
struct TestFaceOperator end
struct TestEdgeOperator end


struct TestFaceRefSpace end 
struct TestEdgeRefSpace end
const TestRefSpace = Union{TestEdgeRefSpace,TestFaceRefSpace}

function (ϕ::TestFaceRefSpace)(mp)

    d = 1
    return SVector((
        (value=1, divergence=d),
        (value=2, divergence=d),
        (value=3, divergence=d),
        (value=4, divergence=d)
    ))
end

function (ϕ::TestEdgeRefSpace)(mp)

    d = 1
    return SVector((
        (value=1, divergence=d),
        (value=2, divergence=d),
        (value=3, divergence=d),
        (value=4, divergence=d),
        (value=5, divergence=d),
        (value=6, divergence=d)
    ))
end

struct TestFaceIntegrand{C,O,L}
    test_triangular_element::C
    trial_triangular_element::C
    op::O
    test_local_space::L
    trial_local_space::L
end

function (igd::TestFaceIntegrand)(u,v)

    x = neighborhood(igd.test_triangular_element,u)
    y = neighborhood(igd.trial_triangular_element,v)

    f = igd.test_local_space(x)
    g = igd.trial_local_space(y)

    SMatrix{4,4}(
        dot(f[1].value,g[1].value),
        dot(f[1].value,g[2].value),
        dot(f[1].value,g[3].value),
        dot(f[1].value,g[4].value),
        dot(f[2].value,g[1].value),
        dot(f[2].value,g[2].value),
        dot(f[2].value,g[3].value),
        dot(f[2].value,g[4].value),
        dot(f[3].value,g[1].value),
        dot(f[3].value,g[2].value),
        dot(f[3].value,g[3].value),
        dot(f[3].value,g[4].value),
        dot(f[4].value,g[1].value),
        dot(f[4].value,g[2].value),
        dot(f[4].value,g[3].value),
        dot(f[4].value,g[4].value),)
end


struct TestEdgeIntegrand{C,O,L}
    test_triangular_element::C
    trial_triangular_element::C
    op::O
    test_local_space::L
    trial_local_space::L
end


function (igd::TestEdgeIntegrand)(u,v)

    x = neighborhood(igd.test_triangular_element,u)
    y = neighborhood(igd.trial_triangular_element,v)

    f = igd.test_local_space(x)
    g = igd.trial_local_space(y)

    SMatrix{6,6}(
        dot(f[1].value,g[1].value),
        dot(f[1].value,g[2].value),
        dot(f[1].value,g[3].value),
        dot(f[1].value,g[4].value),
        dot(f[1].value,g[5].value),
        dot(f[1].value,g[6].value),
        dot(f[2].value,g[1].value),
        dot(f[2].value,g[2].value),
        dot(f[2].value,g[3].value),
        dot(f[2].value,g[4].value),
        dot(f[2].value,g[5].value),
        dot(f[2].value,g[6].value),
        dot(f[3].value,g[1].value),
        dot(f[3].value,g[2].value),
        dot(f[3].value,g[3].value),
        dot(f[3].value,g[4].value),
        dot(f[3].value,g[5].value),
        dot(f[3].value,g[6].value),
        dot(f[4].value,g[1].value),
        dot(f[4].value,g[2].value),
        dot(f[4].value,g[3].value),
        dot(f[4].value,g[4].value),
        dot(f[4].value,g[5].value),
        dot(f[4].value,g[6].value),
        dot(f[5].value,g[1].value),
        dot(f[5].value,g[2].value),
        dot(f[5].value,g[3].value),
        dot(f[5].value,g[4].value),
        dot(f[5].value,g[5].value),
        dot(f[5].value,g[6].value),
        dot(f[6].value,g[1].value),
        dot(f[6].value,g[2].value),
        dot(f[6].value,g[3].value),
        dot(f[6].value,g[4].value),
        dot(f[6].value,g[5].value),
        dot(f[6].value,g[6].value),)
end


kernel_in_bary(op::TestFaceOperator, test_local_space::TestFaceRefSpace, trial_local_space::TestFaceRefSpace,
    test_chart, trial_chart) = TestFaceIntegrand(
        test_chart, trial_chart,op, test_local_space, trial_local_space)

kernel_in_bary(op::TestEdgeOperator, test_local_space::TestEdgeRefSpace, trial_local_space::TestEdgeRefSpace,
    test_chart, trial_chart) = TestEdgeIntegrand(
            test_chart, trial_chart,op, test_local_space, trial_local_space)


function test_parametrized(igd, strat::SauterSchwabStrategy)
    return igd((0,0,0),(0,0,0))
end

const TestOperator = Union{TestFaceOperator, TestEdgeOperator}
function momintegral(op::TestFaceOperator ,test_local_space::TestRefSpace, trial_local_space::TestRefSpace,
    test_tetrahedral_element, trial_tetrahedral_element, out, sing, strat::SauterSchwabStrategy)
        
    I, J = SauterSchwab3D.reorder2(sing, strat)   
    
    L = tetrahedron_face_order(I)
    K = tetrahedron_face_order(J)
    
    test_tetrahedral_element  = simplex(
        test_tetrahedral_element.vertices[I[1]],
        test_tetrahedral_element.vertices[I[2]],
        test_tetrahedral_element.vertices[I[3]],
        test_tetrahedral_element.vertices[I[4]])
        
    trial_tetrahedral_element = simplex(
        trial_tetrahedral_element.vertices[J[1]],
        trial_tetrahedral_element.vertices[J[2]],
        trial_tetrahedral_element.vertices[J[3]],
        trial_tetrahedral_element.vertices[J[4]])
        
           
    igd = kernel_in_bary(op, test_local_space, trial_local_space,
                test_tetrahedral_element, trial_tetrahedral_element)

    G = test_parametrized(igd, strat)
    
    for j ∈ 1:4, i ∈ 1:4
        out[i,j] += G[K[i],L[j]]
    end
        
    nothing
end


function momintegral(op::TestEdgeOperator ,test_local_space::TestRefSpace, trial_local_space::TestRefSpace,
    test_tetrahedral_element, trial_tetrahedral_element, out, sing, strat::SauterSchwabStrategy)
        
    I, J = SauterSchwab3D.reorder2(sing, strat)   
    
    L,O1 = tetrahedron_edge_order(I)
    K,O2 = tetrahedron_edge_order(J)
    
    test_tetrahedral_element  = simplex(
        test_tetrahedral_element.vertices[I[1]],
        test_tetrahedral_element.vertices[I[2]],
        test_tetrahedral_element.vertices[I[3]],
        test_tetrahedral_element.vertices[I[4]])
        
    trial_tetrahedral_element = simplex(
        trial_tetrahedral_element.vertices[J[1]],
        trial_tetrahedral_element.vertices[J[2]],
        trial_tetrahedral_element.vertices[J[3]],
        trial_tetrahedral_element.vertices[J[4]])
        
           
    igd = kernel_in_bary(op, test_local_space, trial_local_space,
                test_tetrahedral_element, trial_tetrahedral_element)

    G = test_parametrized(igd, strat)
    for j ∈ 1:6, i ∈ 1:6
        out[i,j] += G[K[i],L[j]]
    end
        
    nothing
end
    

cv = CommonVertex6D(SauterSchwab3D._legendre(1,0.0,1.0))
sing = SauterSchwab3D.singularity_detection(P_v,Q_v)
res = zeros(6,6)
ts = TestEdgeRefSpace()
op = TestEdgeOperator()

println(sing)
 

momintegral(op,ts,ts,P_v,Q_v,res,sing,cv)

display(res)
=#
#s = SauterSchwab3D.singularity_detection(P_v,Q_v)
#println(s)


#=
const P_e = simplex(pI,pII,pIII,pIV)
const Q_e = simplex(qIII,pI,qIV,pII)

#ce = CommonEdge6D(SauterSchwab3D._legendre(1,0.0,1.0))
#reorder(P_e,Q_e,ce)

s = SauterSchwab3D.singularity_detection(P_e,Q_e)
println(s)

const P_f = simplex(pI,pII,pIV,pIII)
const Q_f = simplex(pI,qIV,pII,pIII)

#cf = CommonFace6D(SauterSchwab3D._legendre(1,0.0,1.0))
#reorder(P_f,Q_f,cf)

s = SauterSchwab3D.singularity_detection(P_f,Q_f)
println(s)
=#


