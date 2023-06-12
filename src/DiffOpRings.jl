import Base: ==, parent

const RatFuncElem = AA.Generic.RationalFunctionFieldElem

############################################################
# 
# DiffOpRing / DiffOpRingElem
#	Ring of differential operators over rational functions
#
############################################################

struct DiffOpRing{T <: MPolyRing{<:RatFuncElem}}
	DOR::T
end

struct DORElem{T <: MPolyRingElem{<:RatFuncElem}} <: AbstractDiffOp
	elem::T
end

############################################################
# 
# DiffOpRing constructor
# 
############################################################

"""
	Ring of differential operators over rational functions
"""
function diff_op_ring(F::Field, s::Vector{Symbol}, ds::Vector{Symbol}; kw...)
	R, gens_R = RationalFunctionField(F, s; kw...)
	D, gens_D = polynomial_ring(R, ds; kw...)
    gens_R = D.(gens_R)
    D = DiffOpRing(D)
	(D, DORElem.(gens_R), DORElem.(gens_D))
end

function diff_op_ring(F::Field, s::Symbol, ds::Symbol; kw...)
    D, x, dx = diff_op_ring(F, [s], [ds]; kw...)
    return D, x[1], dx[1]
end

function diff_op_ring(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	diff_op_ring(F, Symbol.(s), Symbol.("d".*s); kw...)
end
diff_op_ring(s::AbstractVector{<:AbstractString}; kw...) = diff_op_ring(QQ, s; kw...)

function diff_op_ring(F::Field, s::AbstractString; kw...)
    diff_op_ring(F, Symbol(s), Symbol("d"*s); kw...)
end
diff_op_ring(s::AbstractString; kw...) = diff_op_ring(QQ, s; kw...)