
############################################################
# 
# DiffOpRing / DiffOpRingElem
#	Ring of differential operators over rational functions
#
############################################################

struct DiffOpRing{T <: MPolyRing{<:RatFuncElem}} <: AbstractDORing
	DOR::T

    function DiffOpRing(F::Field, s::Vector{Symbol}; kw...)
        length(s) != length(unique(s)) && throw(ArgumentError("variables must be unique"))
        ds = Symbol.("d" .* string.(s))
        R, _ = RationalFunctionField(QQ, string.(s))
        raw_D = AbstractAlgebra.polynomial_ring_only(R, ds; kw...)
        new{typeof(raw_D)}(raw_D)
    end

end
unwrap(R::DiffOpRing) = R.DOR

Base.one(R::DiffOpRing) = R |> unwrap |> one |> DORElem
Base.zero(R::DiffOpRing) = R |> unwrap |> zero |> DORElem

base_ring(R::DiffOpRing) = R |> unwrap |> base_ring
function gens(R::DiffOpRing)
	g = R |> unwrap |> base_ring |> gens
	g = unwrap(R).(g)
	return g .|> (s->DORElem(R, s))
end

function dgens(R::DiffOpRing)
    dg = R |> unwrap |> gens
    return dg .|> (s->DORElem(R, s))
end
# dgens(R::DiffOpRing) = R |> unwrap |> gens .|> DORElem

nvars(R::DiffOpRing) = R |> unwrap |> nvars

function Base.show(io::IO, R::DiffOpRing)
	print(io, nvars(R), "-dimensional ring of differential opeartors in [$(join(string.(gens(R)), ","))]")
end

function (R::DiffOpRing)(num::Union{Rational, Integer})
    DORElem(R, unwrap(R)(num))
end

struct DORElem{T <: MPolyRingElem{<:RatFuncElem}} <: AbstractDiffOp
    parent::DiffOpRing
	elem::T
end

Base.:(==)(x::DiffOpRing, y::DiffOpRing) = unwrap(x) == unwrap(y)




Base.://(x::DORElem,y::Union{Rational, Integer}) = x//(parent(x)(y))
Base.://(x::Union{Rational, Integer},y::DORElem) = (parent(y)(x))//y

function Base.://(x::DORElem,y::DORElem)
    ret_dop = 0
    x_coeff = x |> unwrap |> coefficients |> collect                   
    x_mon = x |> unwrap |> monomials |> collect 
    y_coeff = y |> unwrap |> coefficients |> collect    
    y_mon = y |> unwrap |> monomials |> collect
    size(y_mon)[1] != 1 && return throw(DomainError("division by differential operator", string(y)))
    if y_mon[1] != 1
        return throw(DomainError("division by differential operator", string(y)))
    else
        for i=1:size(x_coeff)[1]
            ret_dop += (x_coeff[i]// y_coeff[1]) * x_mon[i]
        end
        return DORElem(parent(x), ret_dop)
    end
end




############################################################
# 
# DiffOpRing constructor
# 
############################################################

"""
	Ring of differential operators over rational functions
"""
function diff_op_ring(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	D = DiffOpRing(F, Symbol.(s))
    return D, gens(D), dgens(D)
end
diff_op_ring(s::AbstractVector{<:AbstractString}; kw...) = diff_op_ring(QQ, s; kw...)

function diff_op_ring(F::Field, s::AbstractString; kw...)
    D = DiffOpRing(F, [Symbol(s)])
    return D, gens(D)[1], dgens(D)[1]
end
diff_op_ring(s::AbstractString; kw...) = diff_op_ring(QQ, s; kw...)





function _coerce_unsafe(x::RatFuncElem, R::Generic.RationalFunctionField, index_map::Dict{<:Integer, <:Integer})
    n = length(index_map)
    MPoly = R |> zero |> numerator |> parent
    #coercion of numerator
    x_nume = numerator(x)
    coerced_nume = _coerce_unsafe(x_nume, MPoly, index_map)
    
    #coercion of denominator
    x_deno = denominator(x)
    coerced_deno = _coerce_unsafe(x_deno, MPoly, index_map)
    
    # return coerced numerator divided ny coerced denominator

    return R(coerced_nume) // R(coerced_deno)
end

function _coerce_unsafe(x::MPolyRingElem, R::Generic.RationalFunctionField, index_map::Dict{<:Integer, <:Integer})
    n = length(index_map)
    MPoly = R |> zero |> numerator |> parent
    cezip = zip(coefficients(x), exponent_vectors(x))
    C = Generic.MPolyBuildCtx(MPoly)
    for (c,e) in cezip
        push_term!(C, c,[get(e, index_map[i], 0) for i in 1:n] )
    end
    return R(finish(C))
end


function coerce(x::DORElem, D::DiffOpRing)
    hom = coercion_homomorphism(parent(x), D)
    n = nvars(D)

    index_map = Dict{Integer, Integer}()
    invhom = inv(hom)
    for i = 1:n
        index_map[i] = get(invhom, i, 0)
    end

    cezip = zip(coefficients(unwrap(x)), exponent_vectors(unwrap(x)))
    M = Generic.MPolyBuildCtx(unwrap(D))

    R = base_ring(D)
    for (c,e) in cezip
        coerced_c = _coerce_unsafe(c, R ,index_map)
        push_term!(M, coerced_c, [get(e, index_map[i], 0) for i in 1:n])
    end
    DORElem(D, finish(M)) 
end
(R::DiffOpRing)(x::DORElem) = coerce(x, R)



function coerce(x::WAlgElem, D::DiffOpRing)
    hom = coercion_homomorphism(parent(x), D) 
    n = nvars(D)

    # define mapping from index of D to index of parent(x)
    index_map = Dict{Integer, Integer}()
    invhom = inv(hom)

    for i = 1:n
        index_map[i] = get(invhom, i, 0)
    end

    cezip = zip(coefficients(unwrap(x)), exponent_vectors(unwrap(x)))
    M = Generic.MPolyBuildCtx(unwrap(D))
    R = base_ring(D)

    for (c, e) in cezip
        coerced_c = _coerce_unsafe(c, R, index_map)
        push_term!(M, coerced_c, [get(e, index_map[i], 0) for i in 1:n])
    end

    DORElem(D, finish(M))
end
(R::DiffOpRing)(x::WAlgElem) = coerce(x, R)