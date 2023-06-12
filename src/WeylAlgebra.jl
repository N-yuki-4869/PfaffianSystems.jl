import Base: ==, parent

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

Base.one(D::WeylAlgebra) = WAlgElem(one(D.WAlg))
Base.zero(D::WeylAlgebra) = WAlgElem(zero(D.WAlg))

function gens(D::WeylAlgebra)
    gens = D.WAlg |> base_ring |> AA.gens
    gens = D.WAlg.(gens)
    return WAlgElem.(gens)
end
dgens(D::WeylAlgebra) = D.WAlg |> AA.gens .|> WAlgElem

nvars(D::WeylAlgebra) = D.WAlg |> AA.nvars

# function Base.show(io::IO, D::WeylAlgebra)
# 	print(io, "WAlg(", D.WAlg, ")")
# end
function Base.show(io::IO, D::WeylAlgebra)
    print(io, nvars(D),"-dimensional Weyl algebra")
end

# TODO: wrapper for constant


struct WAlgElem{T <: MPolyRingElem{<:MPolyRingElem}} <: AbstractDiffOp
	elem::T
end

Base.parent(wae::WAlgElem) = WeylAlgebra(parent(wae.elem))
gens(wae::WAlgElem) = gens(parent(wae))
dgens(wae::WAlgElem) = dgens(parent(wae))

function Base.show(io::IO, wae::WAlgElem)
	print(io, wae.elem)
end

Base.:+(x::WAlgElem, y::WAlgElem) = WAlgElem(x.elem + y.elem)
Base.:-(x::WAlgElem, y::WAlgElem) = WAlgElem(x.elem - y.elem)
Base.one(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = WAlgElem(one(wae.elem))
Base.zero(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = WAlgElem(zero(wae.elem))

Base.:(==)(x::WAlgElem, y::WAlgElem) = x.elem == y.elem

# TODO: multiplication of WAlgElem
function Base.:*(l::WAlgElem, r::WAlgElem)
    l_coeffs = l.elem |> AA.coefficients |> collect                   
    l_mons = l.elem |> AA.monomials |> collect 
    r_coeffs = r.elem |> AA.coefficients |> collect    
    r_mons = r.elem |> AA.monomials |> collect
    l_size = l_coeffs |> size
    r_size = r_coeffs |> size

    ret_dop = 0
    for i = 1:l_size[1]
        for j = 1:r_size[1]
            ret_dop +=  l_coeffs[i] * Leibniz_rule(l_mons[i],r_coeffs[j]) * r_mons[j]
        end
    end
    return WAlgElem(ret_dop)
end

function Base.:^(x::WAlgElem, y::Integer)
    ret_dop = x
    for _ = 1:y-1
        ret_dop *= x
    end
    return ret_dop
end

function _nth_derivative(f,x,n)
    if n==0
        return f
    else
        for i=1:n
            f = AA.derivative(f,x)
        end
        return f
    end
end


function Leibniz_rule(l_mons,r_coeffs)
    ret_dop = r_coeffs * l_mons
    variable = size(AA.gens(parent(l_mons)))[1]
    for i=1:variable
        coeffs = collect(AA.coefficients(ret_dop))
        mons = collect(AA.monomials(ret_dop))
        coeffs_size = size(coeffs)[1]
        for j=1:coeffs_size
            ret_dop += Leibniz_rule_1(mons[j],coeffs[j],i)
        end
    end
    return ret_dop
end


function Leibniz_rule_1(l_mons,r_coeffs,i)
    ret_dop = 0
    k = 1
    while true
        a = _nth_derivative(r_coeffs,AA.gens(parent(r_coeffs))[i],k) * _nth_derivative(l_mons, AA.gens(parent(l_mons))[i],k) / parent(r_coeffs)(factorial(big(k)))
        a == 0&&break
        ret_dop += a
        k += 1
    end
    return ret_dop
end



# function derivative_roop(l_mons,r_coeffs,n)
#     r_coeff = [r_coeffs]
#     r_mon = [1]
#     r_size = r_coeff |> size
#     if l_mons == 1
#         return r_coeffs * l_mons
#     else
#         ret_dop = 0
#         for i=1:n
#             ret_dop = 0
#             for j=1:r_size[1]
#                 ret_dop += (r_coeff[j] * AA.gens(parent(l_mons))[1] + AA.derivative(r_coeff[j],1)) * r_mon[j] 
#             end
#             r_coeff = collect(coefficients(ret_dop))
#             r_mon = collect(monomials(ret_dop))
#             r_size = r_coeff |> size
#         return ret_dop
#         end
#     end
# end


# TODO: multiplication of WAlgElem and constant
# TODO: multiplication of WAlgElem and polynomial

# TODO: coefficients of WAlgElem 
# TODO: monomials of WAlgElem 
# TODO: terms of WAlgElem

############################################################
# 
# Arithemetic with Rationals and Integers
#
############################################################
Base.:+(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(x + y.elem)
Base.:+(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(x.elem + y)
Base.:*(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(x * y.elem)
Base.:*(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(x.elem * y)

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
    gens_R = D.(gens_R)
    D = WeylAlgebra(D)
	(D, WAlgElem.(gens_R), WAlgElem.(gens_D))
end

function weyl_algebra(F::Field, s::Symbol, ds::Symbol; kw...)
    D, x, dx = weyl_algebra(F, [s], [ds]; kw...)
    return D, x[1], dx[1]
end

function weyl_algebra(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	weyl_algebra(F, Symbol.(s), Symbol.("d".*s); kw...)
end
weyl_algebra(s::AbstractVector{<:AbstractString}; kw...) = weyl_algebra(QQ, s; kw...)

function weyl_algebra(F::Field, s::AbstractString; kw...)
    weyl_algebra(F, Symbol(s), Symbol("d"*s); kw...)
end
weyl_algebra(s::AbstractString; kw...) = weyl_algebra(QQ, s; kw...)




############################################################
# 
# Ideal of Weyl algebra
# 
############################################################

struct Ideal{T <: AbstractDiffOp}
    gens::Vector{T}
    parent::WeylAlgebra

    function Ideal(gens::Vector{T}) where T <: AbstractDiffOp
        par = parent(gens[1])
        if !all([par == parent(g) for g in gens])
            error("All generators must be elements of the same Weyl algebra")
        end
        new{T}(gens, par)
    end
end






