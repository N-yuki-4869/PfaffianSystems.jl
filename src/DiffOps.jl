############################################################
# 
# Common utilities for differential operators
# 
############################################################

unwrap(wae::T) where T <: AbstractDiffOp = wae.elem

Base.parent(wae::T) where T <: AbstractDiffOp = wae.parent
gens(wae::T) where T <: AbstractDiffOp = wae |> parent |> gens
dgens(wae::T) where T <: AbstractDiffOp = wae |> parent |> dgens

function Base.show(io::IO, wae::T) where T <: AbstractDiffOp
    show(io, unwrap(wae))
end

Base.one(wae::Union{Type{T}, T}) where T <: AbstractDiffOp = T(parent(wae), one(unwrap(wae)))
Base.zero(wae::Union{Type{T}, T}) where T <: AbstractDiffOp = T(parent(wae), zero(unwrap(wae)))

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

"""
    vars(wae::T) where T <: AbstractDiffOp

# Examples 

```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
julia> vars(x+y)
2-element Vector{PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 x
 y

julia> R, (x, y), (dx, dy) = diff_op_ring(["x", "y"])
julia> vars(x+y)
2-element Vector{PfaffianSystems.DORElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}}:
 x
 y

```
"""
function vars(wae::T) where T<: AbstractDiffOp
    wae_coeffs = wae |> unwrap |> coefficients
    set = Set{typeof(wae)}()
    for c in wae_coeffs
        c_vars = vars(c)
        c_vars = unwrap(parent(wae)).(c_vars)
        c_vars = c_vars .|> (s->T(parent(wae), s))
        set = union(set, c_vars)

    end
    wae_vars = collect(set)
    
    re_wae_vars = Vector{typeof(wae)}()
    for i in gens(wae)
        if i in wae_vars
            push!(re_wae_vars, i)
        end
    end

    return re_wae_vars
end


"""
    dvars(wae::T) where T <: AbstractDiffOp

# Examples

```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
julia> dvars(dx+dy)
2-element Vector{PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 dx
 dy

julia> R, (x, y), (dx, dy) = diff_op_ring(["x", "y"])
julia> dvars(dx+dy)
2-element Vector{PfaffianSystems.DORElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}}:
 dx
 dy

```

"""
function dvars(wae::T) where T <: AbstractDiffOp
    v = vars(unwrap(wae))
    wae_dvars = collect(Set(v .|> (s->T(parent(wae), s))))

    re_wae_dvars = Vector{typeof(wae)}()
    for i in dgens(wae)
        if i in wae_dvars
            push!(re_wae_dvars, i)
        end
    end

    return re_wae_dvars 
end


############################################################
# 
# Common arithmetic operations for differential operators
# 
############################################################

Base.:(==)(x::T, y::T) where T <: AbstractDiffOp  = unwrap(x) == unwrap(y)
Base.hash(x::T, h::UInt) where T <: AbstractDiffOp = hash(unwrap(x), h)

Base.:+(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp  = T(parent(y), x + unwrap(y))
Base.:+(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) + y)
Base.:-(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp = T(parent(y), x - unwrap(y))
Base.:-(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) - y)
Base.:-(x::T) where T <: AbstractDiffOp = T(parent(x), -unwrap(x))
Base.:*(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp = T(parent(y), x * unwrap(y))
Base.:*(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) * y)

Base.:+(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) + unwrap(y))
Base.:-(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) - unwrap(y))

function Base.:*(l::T, r::T) where T <: AbstractDiffOp

    l == zero(parent(l)) && return zero(parent(l))
    r == zero(parent(r)) && return zero(parent(r))

    l_coeffs = l |> unwrap |> coefficients
    l_mons = l |> unwrap |> monomials
    r_coeffs = r |> unwrap |> coefficients
    r_mons = r |> unwrap |> monomials

    ret_dop = 0
    for (lc, lm) in zip(l_coeffs, l_mons), (rc, rm) in zip(r_coeffs, r_mons)
        ret_dop +=  lc * leibniz_rule(lm, rc) * rm
    end
    return T(parent(l), ret_dop)
end

function derivative(f::T, x::T) where T <: AbstractDiffOp
    !(isempty(dvars(f)) && isempty(dvars(x))) && error("Must not contain any differential operators")

    f_coeffs = f |> unwrap |> coefficients

    isempty(f_coeffs) && return zero(parent(f))
    x_coeffs = x |> unwrap |> coefficients
    isempty(x_coeffs) && return error("Second argument must not be a variable")
    x_mons = x |> unwrap |> monomials
    diffres = derivative(f_coeffs |> first, x_coeffs |> first)*(x_mons |> first)
    return T(parent(f), diffres)
end

function derivative(f::RatFuncElem ,x::RatFuncElem)

    f_nume = numerator(f)
    f_deno = denominator(f)
    x_nume = numerator(x)

    nume = parent(x)((derivative(f_nume,x_nume)*f_deno - f_nume*derivative(f_deno,x_nume)))
    deno = parent(x)((f_deno^2))

    ret_dop = nume // deno
    return ret_dop
end

function diff_op_pow(x::T, y::Integer) where T <: AbstractDiffOp
	y == 0 && return one(x)

    ret_dop = x
    for _ = 1:y-1
        ret_dop *= x
    end
    return ret_dop
end

Base.:literal_pow(::typeof(^), x::AbstractDiffOp,::Val{y}) where y = x^y


function leibniz_rule(l_mons::T, r_coeffs::U) where {U, T <: MPolyRingElem{<:U}}

    ret_dop = r_coeffs * l_mons
    n = nvars(parent(l_mons))

    for i=1:n, (m, c) in zip(monomials(ret_dop), coefficients(ret_dop))
        ret_dop += apply_diff(m,c,i)
    end

    return ret_dop
end

function apply_diff(l_mons::T ,r_coeffs::U, i::Integer) where {U, T <: MPolyRingElem{<:U}}

    ret_dop = 0
    k = 1
    coeff_var = gens(parent(r_coeffs))[i]
    mon_var = gens(parent(l_mons))[i]

    while true
        a = _nth_derivative(r_coeffs, coeff_var,k) * _nth_derivative(l_mons, mon_var,k) / parent(r_coeffs)(factorial(big(k)))
        a == 0 && break
        ret_dop += a
        k += 1
    end

    return ret_dop
end

function _nth_derivative(f::T, x::T ,n::Integer) where T

    n == 0 && return f 

    for _=1:n
        f = derivative(f,x)
    end

    return f

end

coefficients(f::T) where T <: AbstractDiffOp = f |> unwrap |> coefficients |> collect

function monomials(f::T) where T <: AbstractDiffOp

    f_mons = f |> unwrap |> monomials |> collect
    f_mons = collect(f_mons .|> (s->T(parent(f), s)))

    return f_mons
end

exponent_vectors(f::T) where T <: AbstractDiffOp = f |> unwrap |> exponent_vectors |> collect
        