import Base: ==, parent

const RatFuncElem = Generic.RationalFunctionFieldElem

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
# unwrap(wae::DORElem) = wae.elem


# Base.parent(wae::DORElem) = wae |> unwrap |> parent |> DiffOpRing
# Base.parent(wae::DORElem) = wae.parent
# gens(wae::DORElem) = wae |> parent |> gens
# dgens(wae::DORElem) = wae |> parent |> dgens

# function Base.show(io::IO, wae::DORElem)
# 	print(io, unwrap(wae))
# end

# Base.:+(x::DORElem, y::DORElem) = DORElem(parent(x), unwrap(x) + unwrap(y))
# Base.:-(x::DORElem, y::DORElem) = DORElem(parent(x), unwrap(x) - unwrap(y))
# Base.one(wae::Union{Type{DORElem{T}}, DORElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> one |> DORElem
# Base.zero(wae::Union{Type{DORElem{T}}, DORElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> zero |> DORElem

# function vars(wae::DORElem)
#     wae_coeffs = wae |> unwrap |> coefficients |> collect
#     set = Set{typeof(wae)}()
#     for c in wae_coeffs
#         c_nume = numerator(c)
#         c_nume_vars = vars(c_nume)
#         c_nume_vars = unwrap(parent(wae)).(c_nume_vars)
#         c_nume_vars = c_nume_vars .|> (s->DORElem(parent(wae), s))
#         set = union(set, c_nume_vars)

#         c_deno = denominator(c)
#         c_deno_vars = vars(c_deno)
#         c_deno_vars = unwrap(parent(wae)).(c_deno_vars)
#         c_deno_vars = c_deno_vars .|> (s->DORElem(parent(wae), s))
#         set = union(set, c_deno_vars)
#     end
#     wae_vars = collect(set)

#     re_wae_vars = Vector{typeof(wae)}()
#     for i in gens(wae)
#         if i in wae_vars
#             push!(re_wae_vars, i)
#         end
#     end


#     return re_wae_vars
# end

function vars(wae::RatFuncElem)
    set = Set{typeof(wae)}()
    wae_nume = numerator(wae)
    wae_nume_vars = vars(wae_nume)
    set = union(set, wae_nume_vars)

    wae_deno = denominator(wae)
    wae_deno_vars = vars(wae_deno)
    set = union(set, wae_deno_vars)

    return set
end

# function dvars(wae::DORElem)
#     v = vars(unwrap(wae))
#     wae_dvars = collect(Set(v .|> (s->DORElem(parent(wae), s))))

#     re_wae_dvars = Vector{typeof(wae)}()
#     for i in dgens(wae)
#         if i in wae_dvars
#             push!(re_wae_dvars, i)
#         end
#     end

#     return re_wae_dvars
# end

Base.:(==)(x::DiffOpRing, y::DiffOpRing) = unwrap(x) == unwrap(y)
# Base.:(==)(x::DORElem, y::DORElem) = unwrap(x) == unwrap(y)


# function Base.:*(l::DORElem, r::DORElem)
#     l_coeffs = l |> unwrap |> coefficients |> collect                   
#     l_mons = l |> unwrap |> monomials |> collect 
#     r_coeffs = r |> unwrap |> coefficients |> collect    
#     r_mons = r |> unwrap |> monomials |> collect
#     l_size = l_coeffs |> size
#     r_size = r_coeffs |> size

#     ret_dop = 0
#     for i = 1:l_size[1]
#         for j = 1:r_size[1]
#             ret_dop +=  l_coeffs[i] * Leibniz_rule(l_mons[i],r_coeffs[j]) * r_mons[j]
#         end
#     end
#     return DORElem(parent(l), ret_dop)
# end

# function Base.:^(x::DORElem, y::Integer)
#     ret_dop = x
#     for _ = 1:y-1
#         ret_dop *= x
#     end
#     return ret_dop
# end



function derivative(f::RatFuncElem ,x::RatFuncElem)
    f_nume = numerator(f)
    f_deno = denominator(f)
    x_nume = numerator(x)
    nume = parent(x)((derivative(f_nume,x_nume)*f_deno - f_nume*derivative(f_deno,x_nume)))
    deno = parent(x)((f_deno^2))
    #@show nume , deno, typeof(nume) , typeof(deno)
    ret_dop = nume // deno
    return ret_dop
end


# function _nth_derivative(f::T ,x::T ,n::Integer) where T
#     if n==0
#         return f
#     else
#         for i=1:n
#             #@show f,x , typeof(f) , typeof(x)
#             f = derivative(f,x)
#         end
#         return f
#     end
# end


# function Leibniz_rule(l_mons::T ,r_coeffs::U) where {T <: MPolyRingElem{<:RatFuncElem}, U <: RatFuncElem}
#     ret_dop = r_coeffs * l_mons
#     variable = size(gens(parent(l_mons)))[1]
#     for i=1:variable
#         coeffs = collect(coefficients(ret_dop))
#         mons = collect(monomials(ret_dop))
#         coeffs_size = size(coeffs)[1]
#         for j=1:coeffs_size
#             ret_dop += Apply_diff(mons[j],coeffs[j],i)
#         end
#     end
#     return ret_dop
# end


# function Apply_diff(l_mons::T ,r_coeffs:: U, i::Integer) where {T <: MPolyRingElem{<:RatFuncElem}, U <: RatFuncElem}
#     ret_dop = 0
#     k = 1
#     while true
#         a = _nth_derivative(r_coeffs,gens(parent(r_coeffs))[i],k) * _nth_derivative(l_mons, gens(parent(l_mons))[i],k) / parent(r_coeffs)(factorial(big(k)))
#         a == 0 && break
#         ret_dop += a
#         k += 1
#     end
#     return ret_dop
# end


# Base.:+(x::Union{Rational, Integer}, y::DORElem) = DORElem(parent(y), x + unwrap(y))
# Base.:+(x::DORElem, y::Union{Rational, Integer}) = DORElem(parent(x), unwrap(x) + y)
# Base.:-(x::Union{Rational, Integer}, y::DORElem) = DORElem(parent(y), x - unwrap(y))
# Base.:-(x::DORElem, y::Union{Rational, Integer}) = DORElem(parent(x), unwrap(x) - y)
# Base.:*(x::Union{Rational, Integer}, y::DORElem) = DORElem(parent(y), x * unwrap(y))
# Base.:*(x::DORElem, y::Union{Rational, Integer}) = DORElem(parent(x), unwrap(x) * y)

Base.://(x::DORElem,y::Union{Rational, Integer}) = x//(parent(x)(y))
Base.://(x::Union{Rational, Integer},y::DORElem) = (parent(y)(x))//y

# function Base.://(x::DORElem,y::Union{Rational, Integer})
#     ret_dop = 0
#     x_coeff = x |> unwrap |> coefficients |> collect                   
#     x_mon = x |> unwrap |> monomials |> collect 
#     for i=1:size(x_coeff)[1]
#         ret_dop += (x_coeff[i]//y)*x_mon[i]
#     end
#     return DORElem(ret_dop)
# end

# function Base.://(x::Union{Rational, Integer},y::DORElem)
#     y_coeff = y |> unwrap |> coefficients |> collect    
#     y_mon = y |> unwrap |> monomials |> collect
#     size(y_mon)[1] != 1 && return throw(DomainError("division by differential operator", string(y)))
#     if y_mon[1] != 1
#         return throw(DomainError("division by differential operator", string(y)))
#     else
#         ret_dop = parent(x)(x) // y_coeff[1]
#         return DORElem(ret_dop)
#     end
# end


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
# function diff_op_ring(F::Field, s::Vector{Symbol}, ds::Vector{Symbol}; kw...)
# 	R, gens_R = RationalFunctionField(F, s; kw...)
# 	D, gens_D = polynomial_ring(R, ds; kw...)
#     gens_R = D.(gens_R)
#     D = DiffOpRing(D)
# 	(D, DORElem.(gens_R), DORElem.(gens_D))
# end

# function diff_op_ring(F::Field, s::Symbol, ds::Symbol; kw...)
#     D, x, dx = diff_op_ring(F, [s], [ds]; kw...)
#     return D, x[1], dx[1]
# end

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