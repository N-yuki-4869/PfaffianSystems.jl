########## References ##########
# GenericTypes.jl: https://github.com/Nemocas/AbstractAlgebra.jl/blob/597ce2fc0006a3ecdc24da2c63fc47ae0e98b7e6/src/generic/GenericTypes.jl
# 
################################

using AbstractAlgebra

const AA = AbstractAlgebra

############################################################
# 
#	WeylAlgebra / WAlgElem 
#
############################################################

# WeylAlgebra{T} = MPolyRing{MPolyRingElem{T}} where T

struct WeylAlgebra{T <: MPolyRing{<:MPolyRingElem}}
	WAlg::T
end

function Base.show(io::IO, D::WeylAlgebra)
	print(io, "WAlg(", D.WAlg, ")")
end

Base.one(D::WeylAlgebra) = WAlgElem(one(D.WAlg))
Base.zero(D::WeylAlgebra) = WAlgElem(zero(D.WAlg))

# TODO: wrapper for constant


struct WAlgElem{T <: MPolyRingElem}
	elem::T
end

function Base.show(io::IO, wae::WAlgElem)
	print(io, wae.elem)
end

Base.:+(x::WAlgElem, y::WAlgElem) = WAlgElem(x.elem + y.elem)
Base.:-(x::WAlgElem, y::WAlgElem) = WAlgElem(x.elem - y.elem)
Base.one(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = WAlgElem(one(wae.elem))
Base.zero(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = WAlgElem(zero(wae.elem))

# TODO: multiplication of WAlgElem
# TODO: multiplication of WAlgElem and constant
# TODO: multiplication of WAlgElem and polynomial

# TODO: coefficients of WAlgElem 
# TODO: monomials of WAlgElem 
# TODO: terms of WAlgElem

############################################################
# 
# weyl_algebra constructor
# 
############################################################

"""
	WeylAlgebra
"""
# weyl_algebra(R::Ring, s::Union{Tuple{Vararg{T}}, AbstractVector{T}}; kw...)
function weyl_algebra(F::Field, s::Vector{Symbol}, ds::Vector{Symbol}; kw...)
	R, gens_R = polynomial_ring(F, s; kw...)
	D, gens_D = polynomial_ring(R, ds; kw...)
	(WeylAlgebra(D), WAlgElem.(D.(gens_R)), WAlgElem.(gens_D))
end

function weyl_algebra(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	weyl_algebra(F, Symbol.(s), Symbol.("d".*s); kw...)
end

function Base.:*(l::WAlgElem, r::WAlgElem)
    # 以下，l = dx, r = x と仮定する
    l_coeffs = l.elem |> coefficients |> collect
    l_mons = l.elem |> monomials |> collect # = [dx]
    r_coeffs = r.elem |> coefficients |> collect    # = [x]
    r_mons = r.elem |> monomials |> collect
    l_coeff = l_coeffs[1]
    l_mon = l_mons[1]   # = dx
    r_coeff = r_coeffs[1]   # = x
    r_mon = r_mons[1]
    ret_dop = l_coeff * (r_coeff * l_mon + derivative(r_coeff, 1)) * r_mon
    return WAlgElem(ret_dop)
end