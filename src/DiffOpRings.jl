
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
(R::DiffOpRing)(num::Union{Rational, Integer}) = DORElem(R, unwrap(R)(num))

# Base.one(R::DiffOpRing) = R |> unwrap |> one |> DORElem
# Base.zero(R::DiffOpRing) = R |> unwrap |> zero |> DORElem

# base_ring(R::DiffOpRing) = R |> unwrap |> base_ring
# function gens(R::DiffOpRing)
# 	g = R |> base_ring |> gens
# 	g = unwrap(R).(g)
# 	return g .|> (s->DORElem(R, s))
# end

# function dgens(R::DiffOpRing)
#     dg = R |> unwrap |> gens
#     return dg .|> (s->DORElem(R, s))
# end
# dgens(R::DiffOpRing) = R |> unwrap |> gens .|> DORElem

# nvars(R::DiffOpRing) = R |> unwrap |> nvars

elem_type(D::Union{Type{DiffOpRing{T}}, DiffOpRing{T}}) where {S <: RatFuncElem, T <: MPolyRing{S}} = DORElem{Generic.MPoly{S}}

function Base.show(io::IO, R::DiffOpRing)
	print(io, nvars(R), "-dimensional ring of differential opeartors in [$(join(string.(gens(R)), ","))]")
end

# function (R::DiffOpRing)(num::Union{Rational, Integer})
#     DORElem(R, unwrap(R)(num))
# end

struct DORElem{T <: MPolyRingElem{<:RatFuncElem}} <: AbstractDiffOp
    parent::DiffOpRing
	elem::T
end

# Base.:(==)(x::DiffOpRing, y::DiffOpRing) = unwrap(x) == unwrap(y)


function Base.:^(x::DORElem, y::Integer) 
    if y < 0
        x = 1 // x
        y = - y
    end
    return diff_op_pow(x,y)
end


Base.://(x::DORElem,y::Union{Rational, Integer}) = x//(parent(x)(y))
Base.://(x::Union{Rational, Integer},y::DORElem) = (parent(y)(x))//y

function Base.://(x::DORElem,y::DORElem)
    x_coeff = x |> unwrap |> coefficients                    
    x_mon = x |> unwrap |> monomials  
    y_coeff = y |> unwrap |> coefficients    
    y_mon = y |> unwrap |> monomials

    ret_dop = 0

    # size(y_mon)[1] != 1 && return throw(DomainError("division by differential operator", string(y)))


    # for i=1:size(x_coeff)[1]
    #     ret_dop += (x_coeff[i]// y_coeff[1]) * x_mon[i]
    # end

    isempty(dvars(y)) != true && return throw(DomainError("division by differential operator", string(y)))

    for (xc, xm) in zip(x_coeff, x_mon), (yc, ym) in zip(y_coeff, y_mon)
        ret_dop += (xc // yc) * xm
    end

    return DORElem(parent(x), ret_dop)

end

############################################################
# 
# coersions
# 
############################################################

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

function diff_op_ring(F::Field, s::AbstractString, n::Integer)
    D = DiffOpRing(F, [Symbol(s,i) for i = 1:n])
    return D, gens(D), dgens(D)
end
diff_op_ring(s::AbstractString,n::Integer) = diff_op_ring(QQ, s, n)


############################################################
# 
# Ideal of ring of differential operators
# 
############################################################

function normalform(f::T,g::Vector{T}) where T<: DORElem
    r = zero(parent(f))
    q = Vector{typeof(f)}()
    for _ in g
        push!(q, zero(parent(f)))
    end
    f, r, q = nf_2(f,g,r,q)
    if f == zero(parent(f))
        return r, q
    else
        nf_2(f,g,r,q)
    end
end
function nf_2(f::T, g::Vector{T}, r::T, q::Vector{T}) where T<: DORElem
    r_1, q_1 = wnormalform(f,g)
    r = r + leading_term(r_1)
    for i in 1:size(q)[1]
        q[i] = q[i] + q_1[i]
    end
    f = r_1 - leading_term(r_1)
    return f, r, q
end



function wnormalform(f::T, g::Vector{T}) where T<: DORElem
    r = f
    q = Vector{typeof(f)}()
    for _ in g
        push!(q, zero(parent(f)))
    end
    return wnf_2(f, g, r, q)
    
end

function wnf_2(f::T, g::Vector{T}, r::T, q::Vector{T}) where T<: DORElem
    for i in 1:size(g)[1]
        a = exponent_vectors(leading_term(r))[1] - exponent_vectors(leading_term(g[i]))[1]
        if minimum(a) >= 0
            q_mon = one(parent(f))
            for j in 1:size(dgens(f))[1]
                q_mon *= dgens(f)[j] ^ a[j]
            end
            q[i] = DORElem(parent(f), unwrap(q[i]) + coefficients(leading_term(r))[1] // coefficients(leading_term(g[i]))[1]  * unwrap(q_mon))
            r = DORElem(parent(f), unwrap(r) - coefficients(leading_term(r))[1] // coefficients(leading_term(g[i]))[1]  * unwrap(q_mon) *  unwrap(g[i]))
            wnf_2(f, g, r, q)
        end
    end
    return r, q
end


function leading_term(f::DORElem)
    f_coes = coefficients(f)
    f_mons = monomials(f)
    f_mon = f_mons[1] |> unwrap
    return DORElem(parent(f), f_coes[1] * f_mon)
end