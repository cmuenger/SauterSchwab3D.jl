module SauterSchwab3D

using LinearAlgebra
using CompScienceMeshes
using StaticArrays
using FastGaussQuadrature
using GrundmannMoeller
using ShunnHamQuadrature

export sauterschwab_parameterized
export SauterSchwab3DStrategy,Singularity
export PositiveDistance6D, CommonVertex6D, CommonEdge6D, CommonFace6D, CommonVolume6D
export PositiveDistance6D_S, CommonVertex6D_S, CommonEdge6D_S, CommonFace6D_S, CommonVolume6D_S
export PositiveDistance5D, CommonVertex5D, CommonEdge5D, CommonFace5D
export PositiveDistance5D_S, CommonVertex5D_S, CommonEdge5D_S, CommonFace5D_S
export PositiveDistance4D, CommonVertex4D, CommonEdge4D, CommonFace4D
export PositiveDistance4D_S, CommonVertex4D_S, CommonEdge4D_S, CommonFace4D_S
export reorder,reverse_face_order, reverse_edge_order
#=
export CommonEdgeS, CommonEdgeS_a, CommonVertexS
export CommonVertexGM, Orig, Orig2, New, OrigEdge, OrigVertex, NewEdge, GM, GMEdge, GMVertex
=#
function sauterschwab_parameterized end

include("pulled_back_integrals.jl")
include("sauterschwabintegral.jl")
include("parametric_kernel_generator.jl")
include("reorder_vertices.jl")

end # module
