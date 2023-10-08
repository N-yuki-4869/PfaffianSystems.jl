
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

elem_type(D::Union{Type{DiffOpRing{T}}, DiffOpRing{T}}) where {S <: RatFuncElem, T <: MPolyRing{S}} = DORElem{Generic.MPoly{S}}

function Base.show(io::IO, R::DiffOpRing)
	print(io, nvars(R), "-dimensional ring of differential opeartors in [$(join(string.(gens(R)), ","))]")
end

struct DORElem{T <: MPolyRingElem{<:RatFuncElem}} <: AbstractDiffOp
    parent::DiffOpRing
	elem::T
end

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
    f_not0 = true
    while f_not0
        r_1, q_1 = wnormalform(f,g)
        r = r + leading_term(r_1)
        for i in 1:size(q)[1]
            q[i] = q[i] + q_1[i]
        end
        f = r_1 - leading_term(r_1)

        if f == zero(parent(f))
            f_not0 = false
            return r, q
        end
    end
end

function wnormalform(f::T, g::Vector{T}) where T<: DORElem
    r = f
    q = Vector{typeof(f)}()
    for _ in g
        push!(q, zero(parent(f)))
    end

    find_g = true
    while find_g 
        for i in 1:size(g)[1]
            find_g = false
            r == zero(parent(f)) && break
            a = exponent_vectors(leading_term(r))[1] - exponent_vectors(leading_term(g[i]))[1]
            if minimum(a) >= 0
                q_mon = one(parent(f))
                for j in 1:size(dgens(f))[1]
                    q_mon *= dgens(f)[j] ^ a[j]
                end
                q[i] = DORElem(parent(f), unwrap(q[i]) + coefficients(leading_term(r))[1] // coefficients(leading_term(g[i]))[1]  * unwrap(q_mon))
                r = r - DORElem(parent(f), coefficients(leading_term(r))[1] // coefficients(leading_term(g[i]))[1]  * unwrap(q_mon)) * g[i]
                find_g = true
                break
            end
        end
    end
    return r, q 
    
end

function leading_term(f::DORElem, order::Symbol=:lex)
    f == zero(parent(f)) && return zero(parent(f))
    f_coes = coefficients(f)
    f_mons = monomials(f)
    if order == :lex
        f_mon = f_mons[1] |> unwrap
        return DORElem(parent(f), f_coes[1] * f_mon)
    elseif order == :revlex
    elseif order == :grlex
    elseif order == :grevlex
        exp_sum = Vector{Integer}()
        for i in exponent_vectors(f)
            push!(exp_sum, sum(i))
        end
        exp_max = findall(x->x==maximum(exp_sum), exp_sum)
        for i in 1:size(dgens(f))[1]
            for j in exp_max
                # exponent_vectors(f)[j][]


            end

        end
    end
end

function pfaffian_system(G::Vector{T}, S::Vector{T}) where T <: DORElem
    vars = dgens(S[1])
    p = Vector{Vector{Vector{typeof(zero(base_ring(parent(S[1]))))}}}()
    for var in vars
        p_elem = Vector{Vector{typeof(zero(base_ring(parent(S[1]))))}}()
        for s_elem in S
            p_e_elem = Vector{typeof(zero(base_ring(parent(S[1]))))}()
            p_e_elem_coef = coefficients(normalform(var * s_elem, G)[1])
            p_e_elem_mono = monomials(normalform(var * s_elem, G)[1])
            for s in S
                a = true
                for (c,m) in zip(p_e_elem_coef, p_e_elem_mono)
                    if m == s
                        push!(p_e_elem, c)
                        a = false
                    end
                    
                end
                a == true && push!(p_e_elem, zero(base_ring(parent(S[1]))))
                
            end
            push!(p_elem , p_e_elem)
        end
        push!(p, p_elem)
    end
    p_i = Vector{Array{Vector{typeof(zero(base_ring(parent(S[1]))))}}}()
    for i in 1:size(p)[1]
        push!(p_i , vcat(p[i]))
    end
    return p_i
end

