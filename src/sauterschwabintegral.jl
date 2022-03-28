


abstract type SauterSchwab3DStrategy end

#= Tetrahedron-Tetrahedron Interaction =#
struct CommonVolume6D{S,Q}      <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonFace6D{S,Q}        <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonEdge6D{S,Q}        <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonVertex6D{S,Q}      <: SauterSchwab3DStrategy sing::S; qps::Q end
struct PositiveDistance6D{S,Q}  <: SauterSchwab3DStrategy sing::S; qps::Q end

struct CommonVolume6D_S{S,Q}     <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonFace6D_S{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonEdge6D_S{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonVertex6D_S{S,Q}     <: SauterSchwab3DStrategy sing::S; qps::Q end
struct PositiveDistance6D_S{S,Q} <: SauterSchwab3DStrategy sing::S; qps::Q end

#= Tetrahedron-Triangle Interaction =#
struct CommonFace5D{S,Q}         <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonEdge5D{S,Q}         <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonVertex5D{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct PositiveDistance5D{S,Q}   <: SauterSchwab3DStrategy sing::S; qps::Q end

struct CommonFace5D_S{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonEdge5D_S{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonVertex5D_S{S,Q}     <: SauterSchwab3DStrategy sing::S; qps::Q end
struct PositiveDistance5D_S{S,Q} <: SauterSchwab3DStrategy sing::S; qps::Q end

#= Triangle-Triangle Interaction =#
struct CommonFace4D{S,Q}         <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonEdge4D{S,Q}         <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonVertex4D{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct PositiveDistance4D{S,Q}   <: SauterSchwab3DStrategy sing::S; qps::Q end

struct CommonFace4D_S{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonEdge4D_S{S,Q}       <: SauterSchwab3DStrategy sing::S; qps::Q end
struct CommonVertex4D_S{S,Q}     <: SauterSchwab3DStrategy sing::S; qps::Q end
struct PositiveDistance4D_S{S,Q} <: SauterSchwab3DStrategy sing::S; qps::Q end


function _legendre(n,a,b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b-a)/2
    x = (x.+1)/2*(b-a).+a
    collect(zip(x,w))
end


@inbounds function calc_vol(vertices::SMatrix{D,N,T}) where {D,N,T}
    @assert N == D + 1
    X = SMatrix{D,D,T}(vertices[i, j + 1] - vertices[i, 1]
                       for i in 1:D, j in 1:D)
    vol = det(X)
    return vol
end

function get_pts_wts(scheme::TnScheme{N,T},  vertices::SMatrix{D,N,U}) where {N,T,D,U}
    @assert N > 0
    @assert N >= D + 1

    ws = scheme.weights
    ps = scheme.points
    @assert length(ws) == length(ps)


  
    xs  = Vector{SVector{D,T}}()
    for i in 1:length(ws)

        p = ps[i]
       
        push!(xs,vertices * p)
    
    end

    vol = calc_vol(vertices) / factorial(N - 1)

    return xs, ws*vol
end

function get_pts_wts(scheme::TnScheme{N,T}, vertices::SVector{N,SVector{D,U}}) where {N,T,D,U}
    return get_pts_wts( scheme,SMatrix{D,N,U}(vertices[n][d] for d in 1:D, n in 1:N))
end

function get_pts_wts( scheme::TnScheme{N}, vertices::AbstractVector) where {N}
    @assert length(vertices) == N
    @assert N > 0
    D = length(vertices[1])
    @assert N >= D + 1
    vertices′ = SVector{N}(map(SVector{D}, vertices))
    return get_pts_wts(scheme, vertices′)
end

function _grundmannMoeller2D(n)
    scheme = GrundmannMoeller.grundmann_moeller(Float64,Val(2),n)
    vertices=[[0,0 ], [1,0], [1,1]]
    x,w = get_pts_wts(scheme,vertices)

    zip(x,w)
end

function _grundmannMoeller3D(n)
    scheme = GrundmannMoeller.grundmann_moeller(Float64,Val(3),n)
    vertices=[[0 0 0], [1,0,0], [1,1,0], [1,1,1] ]
    x,w = get_pts_wts(scheme,vertices)

    zip(x,w)
end

function _grundmannMoeller4D(n)
    scheme = GrundmannMoeller.grundmann_moeller(Float64,Val(4),n)
    vertices=[[0,0,0,0], [1,0,0,0], [1,1,0,0], [1,1,1,0], [1,1,1,1]]
    x,w = get_pts_wts(scheme,vertices)

    zip(x,w)
end


function _shunnham2D(m)
    vertices =  [[0 0], [1 0], [1 1]]
    scheme = ShunnHamQuadrature.shunnham2D_ref(m)

    x,w = ShunnHamQuadrature.get_pts_wts(scheme,vertices)

    zip(x,w)
end

function _shunnham3D(m)
    vertices =  [[0 0 0], [1 0 0], [1 1 0], [1 1 1]]
    scheme = ShunnHamQuadrature.shunnham3D_ref(m)

    x,w = ShunnHamQuadrature.get_pts_wts(scheme,vertices)

    zip(x,w)
end

function _shunnham4D(m)
    vertices =  [[0 0 0 0], [1 0 0 0], [1 1 0 0], [1 1 1 0], [1 1 1 1]]
    scheme = ShunnHamQuadrature.shunnham4D(m)

    x,w = ShunnHamQuadrature.get_pts_wts(scheme,vertices)

    zip(x,w)
end



#=Tetrahedron-tetrahedron interaction =#
#Common Volume Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonVolume6D)

    qps = method.qps
	sum(w1*w2*w3*w4*w5*w6*k6p_ct(integrand, η1, η2, η3, η4, ξ1, ξ2)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (η4, w4) in qps, (ξ1, w5) in qps, (ξ2, w6) in qps)
end
#Common Volume Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonVolume6D_S)

    qps = method.qps
	sum(w1*w2*w3*k6p_ct_A(integrand, x[1], x[2], x[3], x[4], y, z)
		for (x, w1) in qps[3], (y, w2) in qps[1], (z, w3) in qps[1]) +
    sum(w1*w2*k6p_ct_B(integrand, x[1], x[2], x[3], x[4], y[1], y[2])
		for (x, w1) in qps[3], (y, w2) in qps[2])
end

#Common Face Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonFace6D)

	qps = method.qps
	sum(w1*w2*w3*w4*w5*w6*k6p_cf(integrand, η1, η2, η3, η4, ξ1, ξ2)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (η4, w4) in qps, (ξ1, w5) in qps, (ξ2, w6) in qps)
end

#Common Face Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonFace6D_S)

	qps = method.qps
	sum(w1*w2*w3*k6p_cf_A(integrand, x[1], x[2], x[3], y[1], y[2], z)
		for (x, w1) in qps[3], (y, w2) in qps[2], (z, w3) in qps[1])+
    sum(w1*w2*w3*k6p_cf_B(integrand, x[1], x[2], x[3], y[1], y[2], z)
		for (x, w1) in qps[3], (y, w2) in qps[2], (z, w3) in qps[1])+
    sum(w1*w2*k6p_cf_C(integrand, x[1], x[2], x[3], y[1], y[2], y[3])
		for (x, w1) in qps[3], (y, w2) in qps[3])
end

#Common Edge Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonEdge6D)
	qps = method.qps
	sum(w1*w2*w3*w4*w5*w6*k6p_ce(integrand, η1, η2, η3, η4, ξ1, ξ2)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (η4, w4) in qps, (ξ1, w5) in qps, (ξ2, w6) in qps)
end

#Common Simplex Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonEdge6D_S)
	qps = method.qps
	sum(w1*w2*w3*k6p_ce_A(integrand, x[1], x[2], y[1], y[2], z[1], z[2])
		for (x, w1) in qps[2], (y, w2) in qps[2], (z, w3) in qps[2]) +
    sum(w1*w2*w3*k6p_ce_B(integrand, x[1], x[2], y[1], y[2], y[3], z)
		for (x, w1) in qps[2], (y, w2) in qps[3], (z, w3) in qps[1]) + 
    sum(w1*w2*k6p_ce_C(integrand, x[1], x[2], y[1], y[2], y[3], y[4])
		for (x, w1) in qps[2], (y, w2) in qps[4] ) +
    sum(w1*w2*w3*k6p_ce_D(integrand, x[1], x[2], y, z[1], z[2], z[3])
		for (x, w1) in qps[2], (y, w2) in qps[1], (z, w3) in qps[3]) +
    sum(w1*w2*w3*w4*k6p_ce_E(integrand, x[1], x[2], y, z[1], z[2], w)
		for (x, w1) in qps[2], (y, w2) in qps[1], (z, w3) in qps[2], (w, w4) in qps[1])
end

#Common Vertex Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonVertex6D)
	qps = method.qps
    sum(w1*w2*w3*w4*w5*w6*k6p_cv(integrand, η1, η2, η3, η4, ξ1, ξ2)
        for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (η4, w4) in qps, (ξ1, w5) in qps, (ξ2, w6) in qps)
end

#Common Vertex Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonVertex6D_S)
	qps = method.qps
    sum(w1*w2*k6p_cv_A(integrand, x[1], x[2], x[3], y[1], y[2], y[3])
        for (x, w1) in qps, (y, w2) in qps)
end

#Positive Distance Tensor-Product
function sauterschwab_parameterized(integrand, method::PositiveDistance6D)
	qps = method.qps
	sum(w1*w2*w3*w4*w5*w6*k6p_pd(integrand, η1, η2, η3, η4, ξ1, ξ2)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (η4, w4) in qps, (ξ1, w5) in qps, (ξ2, w6) in qps)
end

#Positive Distance Simplex-Product
function sauterschwab_parameterized(integrand, method::PositiveDistance6D_S)

	qps = method.qps
	sum(w1*w2*k6p_pd_A(integrand, x[1], x[2], x[3], y[1], y[2], y[3])
		for (x, w1) in qps, (y, w2) in qps)
end


#=Tetrahedron-triangle interaction =#
#Common Face Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonFace5D)
    qps = method.qps
	sum(w1*w2*w3*w4*w5*k5p_cf(integrand, x1, x2, x3, y1, y2)
		for (x1, w1) in qps, (x2, w2) in qps, (x3, w3) in qps, (y1, w4) in qps, (y2, w5) in qps)
end
#Common Face Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonFace5D_S)
    qps = method.qps
	sum(w1*w2*w3*k5p_cf_A(integrand, x[1], x[2], x[3], y, z)
		for (x, w1) in qps[3], (y, w2) in qps[1], (z, w3) in qps[1]) +
    sum(w1*w2*k5p_cf_B(integrand, x[1], x[2], x[3], y[1], y[2])
        for (x, w1) in qps[3], (y, w2) in qps[2]) 
end

#Common Edge Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonEdge5D)
    qps = method.qps
	sum(w1*w2*w3*w4*w5*k5p_ce(integrand, x1, x2, x3, y1, y2)
		for (x1, w1) in qps, (x2, w2) in qps, (x3, w3) in qps, (y1, w4) in qps, (y2, w5) in qps)
end

#Common Edge Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonEdge5D_S)
    qps = method.qps
	sum(w1*w2*w3*k5p_ce_A(integrand, x[1], x[2], y[1], y[2], z)
		for (x, w1) in qps[2], (y, w2) in qps[2], (z, w3) in qps[1]) +
    sum(w1*w2*k5p_ce_B(integrand, x[1], x[2], y[1], y[2], y[3])
		for (x, w1) in qps[2], (y, w2) in qps[3]) +
    sum(w1*w2*w3*k5p_ce_C(integrand, x[1], x[2], y[1], y[2], z)
		for (x, w1) in qps[2], (y, w2) in qps[2], (z, w3) in qps[1])
end

#Common Vertex Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonVertex5D)
    qps = method.qps
	sum(w1*w2*w3*w4*w5*k5p_cv(integrand, x1, x2, x3, y1, y2)
		for (x1, w1) in qps, (x2, w2) in qps, (x3, w3) in qps, (y1, w4) in qps, (y2, w5) in qps)
end

#Common Vertex Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonVertex5D_S)
    qps = method.qps
	sum(w1*w2*k5p_cv_A(integrand, x[1], x[2], x[3], y[1], y[2])
		for (x, w1) in qps[1], (y, w2) in qps[2])
end

#Positive Distance Tensor-Product
function sauterschwab_parameterized(integrand, method::PositiveDistance5D)

    qps = method.qps
	sum(w1*w2*w3*w4*w5*k5p_pd(integrand, x1, x2, x3, y1, y2)
		for (x1, w1) in qps, (x2, w2) in qps, (x3, w3) in qps, (y1, w4) in qps, (y2, w5) in qps)
end

#Positive Distance Simplex-Product
function sauterschwab_parameterized(integrand, method::PositiveDistance5D_S)

    qps = method.qps
	sum(w1*w2*k5p_pd_A(integrand, x[1], x[2], x[3], y[1], y[2])
		for (x, w1) in qps[1], (y, w2) in qps[2])
end


#=Triangle-triangle interaction =#
#Common Face Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonFace4D)

	qps = method.qps
	sum(w1*w2*w3*w4*k4p_cf(integrand, η1, η2, η3, ξ1)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ1, w4) in qps)
end
#Common Face Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonFace4D_S)

	qps = method.qps
	sum(w1*w2*k4p_cf_A(integrand, x[1], x[2], x[3], y)
		for (x, w1) in qps[2], (y, w2) in qps[1])
end

#Common Edge Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonEdge4D)

	qps = method.qps
	sum(w1*w2*w3*w4*k4p_ce(integrand, η1, η2, η3, ξ1)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ1, w4) in qps)
end

#Common Edge Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonEdge4D_S)

	qps = method.qps
	sum(w1*w2*w3*k4p_ce_A(integrand, x[1], x[2], y, z)
		for (x, w1) in qps[2],  (y, w2) in qps[1], (z, w3) in qps[1]) +
    sum(w1*w2*k4p_ce_B(integrand, x[1], x[2], y[1], y[2])
            for (x, w1) in qps[2], (y, w2) in qps[2])     
end

#Common Vertex Tensor-Product
function sauterschwab_parameterized(integrand, method::CommonVertex4D)

	qps = method.qps
	sum(w1*w2*w3*w4*k4p_cv(integrand, η1, η2, η3, ξ1)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ1, w4) in qps)
end

#Common Vertex Simplex-Product
function sauterschwab_parameterized(integrand, method::CommonVertex4D_S)

	qps = method.qps
	sum(w1*w2*k4p_cv_A(integrand, x[1], x[2], y[1], y[2])
		for (x, w1) in qps, (y, w2) in qps)
end

#Positive Distance Tensor-Product
function sauterschwab_parameterized(integrand, method::PositiveDistance4D)

	qps = method.qps
	sum(w1*w2*w3*w4*k4p_pd(integrand, η1, η2, η3, ξ1)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ1, w4) in qps)
end

#Positive Distance Simplex-Product
function sauterschwab_parameterized(integrand, method::PositiveDistance4D_S)

	qps = method.qps
	sum(w1*w2*k4p_pd_A(integrand, x[1], x[2], y[1], y[2])
		for (x, w1) in qps, (y, w2) in qps)
end



#=
function sauterschwab_parameterized(integrand, method::CommonEdgeS)

    integral = 0.0
    qps_1 = method.qps[1]
    qps_2 = method.qps[2]
    qps_3 = method.qps[3]
    qps_4 = method.qps[4]

	integral += sum(w1*w2*w3*k3p_ce_1(integrand, x[1], x[2], y[1], y[2], z[1], z[2])
		for (x, w1) in qps_2, (y, w2) in qps_2, (z, w3) in qps_2)
    
    integral += sum(w1*w2*w3*k3p_ce_2(integrand, x[1], x[2], y[1], y[2], y[3], z[1])
        for (x, w1) in qps_2, (y, w2) in qps_3, (z, w3) in qps_1)
    
    integral += sum(w1*w2*k3p_ce_3(integrand, x[1], x[2], y[1], y[2], y[3], y[4])
        for (x, w1) in qps_2, (y, w2) in qps_4)
    
    integral += sum(w1*w2*w3*k3p_ce_4(integrand, x[1], x[2], y[1], z[1], z[2], z[3])
        for (x, w1) in qps_2, (y, w2) in qps_1, (z, w3) in qps_3)
    
    integral += sum(w1*w2*w3*w4*k3p_ce_5(integrand, x[1], x[2], y[1], z[1], z[2], w[1])
        for (x, w1) in qps_2, (y, w2) in qps_1, (z, w3) in qps_2, (w, w4) in qps_1)
    
    return integral
end

function sauterschwab_parameterized(integrand, method::CommonVertexS)

    qps = method.qps
	sum(w1*w2*k3p_cv_s( integrand, z[1], z[2], z[3], x[1], x[2], x[3])
		for (x,w1) in qps, (z,w2) in qps)

end

function sauterschwab_parameterized(integrand, method::CommonEdgeS_a)

    integral = 0.0
    qps_1 = method.qps[1]
    qps_2 = method.qps[2]
    qps_3 = method.qps[3]
    qps_4 = method.qps[4]

	integral += sum(w1*w2*w3*k3p_ce_1_slow(integrand, x[1], x[2], y[1], y[2], z[1], z[2])
		for (x, w1) in qps_2, (y, w2) in qps_2, (z, w3) in qps_2)
    #=
    integral += sum(w1*w2*w3*k3p_ce_2(integrand, x[1], x[2], y[1], y[2], y[3], z[1])
        for (x, w1) in qps_2, (y, w2) in qps_3, (z, w3) in qps_1)
    
    integral += sum(w1*w2*k3p_ce_3(integrand, x[1], x[2], y[1], y[2], y[3], y[4])
        for (x, w1) in qps_2, (y, w2) in qps_4)
    
    integral += sum(w1*w2*w3*k3p_ce_4(integrand, x[1], x[2], y[1], z[1], z[2], z[3])
        for (x, w1) in qps_2, (y, w2) in qps_1, (z, w3) in qps_3)
    
    integral += sum(w1*w2*w3*w4*k3p_ce_5(integrand, x[1], x[2], y[1], z[1], z[2], w[1])
        for (x, w1) in qps_2, (y, w2) in qps_1, (z, w3) in qps_2, (w, w4) in qps_1)
     =#
    return integral
end

function sauterschwab_parameterized(integrand, method::CommonVertexGM)

	qps = method.qps
    sum(w1*w2*k3p_cv_gm(integrand, x[1], x[2], x[3], y[1], y[2], y[3])
        for (x, w1) in qps[1], (y, w2) in qps[2])
end

=#

#=
function sauterschwab_parameterized(integrand, method::New)

	qps = method.qps
	sum(w1*w2*w3*w4*k3p_cf_new(integrand, η1, η2, ξ1, ξ2)
		for (η1, w1) in qps, (η2, w2) in qps, (ξ1, w3) in qps, (ξ2, w4) in qps)
end


function sauterschwab_parameterized(integrand, method::GM)

	qps = method.qps
	sum(w1*w2*k3p_cf_gm( integrand, y[1], y[2], y[3], x)
		for (x,w1) in qps[1], (y,w2) in qps[2])
	
end

function sauterschwab_parameterized(integrand, method::Orig2)

	qps = method.qps
	sum(w1*w2*w3*w4*k3p_cf_orig2(integrand, η1, η2, η3, ξ1)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ1, w4) in qps)
end

function sauterschwab_parameterized(integrand, method::NewEdge)

	qps = method.qps
	sum(w1*w2*w3*w4*k3p_ce_new(integrand, η1, η2, η3, ξ1)
		for (η1, w1) in qps, (η2, w2) in qps, (η3, w3) in qps, (ξ1, w4) in qps)
end

function sauterschwab_parameterized(integrand, method::GMEdge)

	qps = method.qps
	sum(w1*w2*w3*k3p_ce_gm( integrand, z[1], z[2], x, y)
		for (x,w1) in qps[1], (y,w2) in qps[1], (z,w3) in qps[2])
	
end
=#
#=
function sauterschwab_parameterized(integrand, method::GMVertex)

	qps = method.qps
	sum(w1*w2*k3p_cv_gm( integrand, z[1], z[2], x[1], x[2])
		for (x,w1) in qps[1], (z,w2) in qps[2])
	
end
=#