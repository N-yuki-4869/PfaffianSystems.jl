
############################################################
# 
#	WeylAlgebra / WAlgElem 
#
############################################################

# WeylAlgebra{T} = MPolyRing{MPolyRingElem{T}} where T

# TODO: add inner constructor to enforce the correspondence between variables and derivations and prohibit duplication of variable names
struct WeylAlgebra{T <: MPolyRing{<:MPolyRingElem}} <: AbstractDORing
	WAlg::T

    function WeylAlgebra(F::Field, s::Vector{Symbol}; kw...) 
        length(s) != length(unique(s)) && throw(ArgumentError("variables must be unique"))
        ds = Symbol.("d" .* string.(s))
        polyR = AbstractAlgebra.polynomial_ring_only(F, s; kw...)
        raw_D = AbstractAlgebra.polynomial_ring_only(polyR, ds; kw...)
        new{typeof(raw_D)}(raw_D)
    end
end
unwrap(D::WeylAlgebra) = D.WAlg

Base.one(D::WeylAlgebra) = D |> unwrap |> one |> WAlgElem
Base.zero(D::WeylAlgebra) = D |> unwrap |> zero |> WAlgElem

base_ring(D::WeylAlgebra) = D |> unwrap |> base_ring
function gens(D::WeylAlgebra)
    g = D |> base_ring |> gens
    g = unwrap(D).(g)
    # g .|> WAlgElem
    return g .|> (s->WAlgElem(D, s))
end
# dgens(D::WeylAlgebra) = D |> unwrap |> gens .|> WAlgElem
function dgens(D::WeylAlgebra)
    dg = D |> unwrap |> gens
    return dg .|> (s->WAlgElem(D, s))
end
nvars(D::WeylAlgebra) = D |> unwrap |> nvars

elem_type(D::Union{Type{WeylAlgebra{T}}, WeylAlgebra{T}}) where {S <: MPolyRingElem, T <: MPolyRing{S}} = WAlgElem{Generic.MPoly{S}}

function Base.show(io::IO, D::WeylAlgebra)
    print(io, nvars(D),"-dimensional Weyl algebra in [$(join(string.(gens(D)), ","))]")
end

# TODO: wrapper for constant


struct WAlgElem{T <: MPolyRingElem{<:MPolyRingElem}} <: AbstractDiffOp
    parent::WeylAlgebra
	elem::T
end
# unwrap(wae::WAlgElem) = wae.elem

# unwrap(wae::T) where T <: AbstractDiffOp = wae.elem

# Base.parent(wae::WAlgElem) = wae |> unwrap |> parent |> WeylAlgebra
# Base.parent(wae::WAlgElem) = wae.parent
# Base.parent(wae::T) where T <: AbstractDiffOp = wae.parent
# gens(wae::WAlgElem) = wae |> parent |> gens
# gens(wae::T) where T <: AbstractDiffOp = wae |> parent |> gens
# dgens(wae::WAlgElem) = wae |> parent |> dgens
# dgens(wae::T) where T <: AbstractDiffOp = wae |> parent |> dgens

# function Base.show(io::IO, wae::WAlgElem)
# 	show(io, unwrap(wae))
# end

# function Base.show(io::IO, wae::T) where T <: AbstractDiffOp
#     show(io, unwrap(wae))
# end

# Base.:+(x::WAlgElem{T}, y::WAlgElem{T}) where T <: MPolyRingElem = WAlgElem{T}(parent(x), unwrap(x) + unwrap(y))
# Base.:+(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) + unwrap(y))
# Base.:-(x::WAlgElem{T}, y::WAlgElem{T}) where T <: MPolyRingElem = WAlgElem{T}(parent(x), unwrap(x) - unwrap(y))
# Base.:-(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) - unwrap(y))
# Base.one(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> one |> WAlgElem
# Base.one(wae::Union{Type{T}, T}) where T <: AbstractDiffOp = T(parent(wae), one(unwrap(wae)))
# Base.zero(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> zero |> WAlgElem
# Base.zero(wae::Union{Type{T}, T}) where T <: AbstractDiffOp = T(parent(wae), zero(unwrap(wae)))

# function vars(wae::WAlgElem)
#     wae_coeffs = wae |> unwrap |> coefficients
#     set = Set{typeof(wae)}()
#     for c in wae_coeffs
#         c_vars = vars(c)
#         c_vars = unwrap(parent(wae)).(c_vars)
#         c_vars = c_vars .|> (s->WAlgElem(parent(wae), s))
#         set = union(set, c_vars)

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

# function vars(wae::T) where T<: AbstractDiffOp
#     wae_coeffs = wae |> unwrap |> coefficients
#     set = Set{typeof(wae)}()
#     for c in wae_coeffs
#         c_vars = vars(c)
#         c_vars = unwrap(parent(wae)).(c_vars)
#         c_vars = c_vars .|> (s->T(parent(wae), s))
#         set = union(set, c_vars)

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

# function dvars(wae::WAlgElem)
#     v = vars(unwrap(wae))
#     wae_dvars = collect(Set(v .|> (s->WAlgElem(parent(wae), s))))

#     re_wae_dvars = Vector{typeof(wae)}()
#     for i in dgens(wae)
#         if i in wae_dvars
#             push!(re_wae_dvars, i)
#         end
#     end

#     return re_wae_dvars
# end

# function dvars(wae::T) where T <: AbstractDiffOp
#     v = vars(unwrap(wae))
#     wae_dvars = collect(Set(v .|> (s->T(parent(wae), s))))

#     re_wae_dvars = Vector{typeof(wae)}()
#     for i in dgens(wae)
#         if i in wae_dvars
#             push!(re_wae_dvars, i)
#         end
#     end

#     return re_wae_dvars 
# end

# function dvars(wae::T) where T <: AbstractDiffOp
#     return Set(vars(unwrap(wae)))    
# end


Base.:(==)(x::WeylAlgebra, y::WeylAlgebra) = unwrap(x) == unwrap(y)
# Base.:(==)(x::WAlgElem{T}, y::WAlgElem{T}) where T <: MPolyRingElem = unwrap(x) == unwrap(y)
# Base.:(==)(x::T, y::T) where T <: AbstractDiffOp  = unwrap(x) == unwrap(y)

# TODO: multiplication of WAlgElem
# function Base.:*(l::WAlgElem{T}, r::WAlgElem{T}) where T <: MPolyRingElem
#     # l_coeffs = l.elem |> coefficients |> collect                   
#     # l_mons = l.elem |> monomials |> collect 
#     # r_coeffs = r.elem |> coefficients |> collect    
#     # r_mons = r.elem |> monomials |> collect
#     l_coeffs = l |> unwrap |> coefficients 
#     l_mons = l |> unwrap |> monomials
#     r_coeffs = r |> unwrap |> coefficients
#     r_mons = r |> unwrap |> monomials
#     # l_size = l_coeffs |> size
#     # r_size = r_coeffs |> size

#     ret_dop = 0
#     # for i = 1:l_size[1]
#     #     for j = 1:r_size[1]
#     #         ret_dop +=  l_coeffs[i] * Leibniz_rule(l_mons[i],r_coeffs[j]) * r_mons[j]
#     #     end
#     # end
#     for (lc, lm) in zip(l_coeffs, l_mons), (rc, rm) in zip(r_coeffs, r_mons) 
#         ret_dop +=  lc * Leibniz_rule(lm, rc) * rm
#     end
#     return WAlgElem(parent(l), ret_dop)
# end

# function Base.:*(l::T, r::T) where T <: AbstractDiffOp
#     l_coeffs = l |> unwrap |> coefficients 
#     l_mons = l |> unwrap |> monomials
#     r_coeffs = r |> unwrap |> coefficients
#     r_mons = r |> unwrap |> monomials

#     ret_dop = 0

#     for (lc, lm) in zip(l_coeffs, l_mons), (rc, rm) in zip(r_coeffs, r_mons) 
#         ret_dop +=  lc * Leibniz_rule(lm, rc) * rm
#     end
#     return T(parent(l), ret_dop)
# end

# function Base.:^(x::WAlgElem, y::Integer)
#     ret_dop = x
#     for _ = 1:y-1
#         ret_dop *= x
#     end
#     return ret_dop
# end

# function Base.:^(x::T, y::Integer) where T <: AbstractDiffOp
#     ret_dop = x
#     for _ = 1:y-1
#         ret_dop *= x
#     end
#     return ret_dop
# end



# function _nth_derivative(f::T, x::T ,n::Integer) where T
#     n == 0 && return f 
#     # if n==0
#     #     return f
#     # else
#     for i=1:n
#         f = derivative(f,x)
#     end
#     return f
#     # end
# end


# function Leibniz_rule(l_mons::T, r_coeffs::U) where {T <: MPolyRingElem{<:MPolyRingElem}, U <: MPolyRingElem}
#     ret_dop = r_coeffs * l_mons
#     # variable = size(gens(parent(l_mons)))[1]
#     n = nvars(parent(l_mons))
#     # for i=1:n
#         # coeffs = collect(coefficients(ret_dop))
#         # mons = collect(monomials(ret_dop))
#         # coeffs_size = size(coeffs)[1]
#         # @show coeffs, mons
#         # for j=1:coeffs_size
#             # @show mons[j], coeffs[j], i
#             # ret_dop += Leibniz_rule_1(mons[j],coeffs[j],i)
#             # @show m, c, i
#     # end
#     for i=1:n, (m, c) in zip(monomials(ret_dop), coefficients(ret_dop))
#         ret_dop += Apply_diff(m,c,i)
#     end
#     return ret_dop
# end

# function Leibniz_rule(l_mons::T, r_coeffs::U) where {U, T <: MPolyRingElem{<:U}}
#     ret_dop = r_coeffs * l_mons
#     n = nvars(parent(l_mons))
#     for i=1:n, (m, c) in zip(monomials(ret_dop), coefficients(ret_dop))
#         ret_dop += Apply_diff(m,c,i)
#     end
#     return ret_dop
# end

# function Apply_diff(l_mons::T ,r_coeffs::U, i::Integer) where {T <: MPolyRingElem{<:MPolyRingElem}, U <: MPolyRingElem}
#     ret_dop = 0
#     k = 1
#     coeff_var = gens(parent(r_coeffs))[i]
#     mon_var = gens(parent(l_mons))[i]
#     while true
#         # a = _nth_derivative(r_coeffs,gens(parent(r_coeffs))[i],k) * _nth_derivative(l_mons, gens(parent(l_mons))[i],k) / parent(r_coeffs)(factorial(big(k)))
#         a = _nth_derivative(r_coeffs, coeff_var,k) * _nth_derivative(l_mons, mon_var,k) / parent(r_coeffs)(factorial(big(k)))
#         a == 0 && break
#         ret_dop += a
#         k += 1
#     end
#     return ret_dop
# end

# function Apply_diff(l_mons::T ,r_coeffs::U, i::Integer) where {U, T <: MPolyRingElem{<:U}}
#     ret_dop = 0
#     k = 1
#     coeff_var = gens(parent(r_coeffs))[i]
#     mon_var = gens(parent(l_mons))[i]
#     while true
#         # a = _nth_derivative(r_coeffs,gens(parent(r_coeffs))[i],k) * _nth_derivative(l_mons, gens(parent(l_mons))[i],k) / parent(r_coeffs)(factorial(big(k)))
#         a = _nth_derivative(r_coeffs, coeff_var,k) * _nth_derivative(l_mons, mon_var,k) / parent(r_coeffs)(factorial(big(k)))
#         a == 0 && break
#         ret_dop += a
#         k += 1
#     end
#     return ret_dop
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
# Base.:+(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(parent(y), x + unwrap(y))
# Base.:+(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp  = T(parent(y), x + unwrap(y))
# Base.:+(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(parent(x), unwrap(x) + y)
# Base.:+(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) + y)
# Base.:-(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(parent(y), x - unwrap(y))
# Base.:-(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp = T(parent(y), x - unwrap(y))
# Base.:-(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(parent(x), unwrap(x) - y)
# Base.:-(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) - y)
# Base.:*(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(parent(y), x * unwrap(y))
# Base.:*(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp = T(parent(y), x * unwrap(y))
# Base.:*(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(parent(x), unwrap(x) * y)
# Base.:*(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) * y)

############################################################
# 
# weyl_algebra constructor
# 
############################################################

"""
	WeylAlgebra
"""
function weyl_algebra(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	# weyl_algebra(F, Symbol.(s), Symbol.("d".*s); kw...)
    D = WeylAlgebra(F, Symbol.(s))
    return D, gens(D), dgens(D)
end
weyl_algebra(s::AbstractVector{<:AbstractString}; kw...) = weyl_algebra(QQ, s; kw...)

function weyl_algebra(F::Field, s::AbstractString; kw...)
    # weyl_algebra(F, Symbol(s), Symbol("d"*s); kw...)
    D = WeylAlgebra(F, [Symbol(s)])
    return D, gens(D)[1], dgens(D)[1]
end
weyl_algebra(s::AbstractString; kw...) = weyl_algebra(QQ, s; kw...)




############################################################
# 
# Ideal of Weyl algebra
# 
############################################################

mutable struct DIdeal{T <: AbstractDiffOp} <: Ideal{T}
    base_ring::AbstractDORing
    gens::Vector{T}

    function DIdeal{T}(R::AbstractDORing, gens::Vector) where T <: AbstractDiffOp
        if !all([R == parent(g) for g in gens])
            # error("All generators must be elements of the same Weyl algebra")
            throw(ArgumentError("All generators must be included in", R))
        end
        new{T}(R, gens)
    end 
end

function DIdeal(R::AbstractDORing, gens::Vector{T}) where T <: AbstractDiffOp
    DIdeal{eltype(gens)}(R, gens)
end
DIdeal(gens::Vector{T}) where T <: AbstractDiffOp = DIdeal(parent(gens[1]), gens)

base_ring(I::DIdeal) = I.base_ring
gens(I::DIdeal) = I.gens

function (D::WeylAlgebra)(I::DIdeal)
    DIdeal(D, I |> gens .|> (s->coerce(s, D)))
end

# TOTO: make ideals extensible
# function extension(I::DIdeal, D::AbstractDORing)
#     R = base_ring(D)
#     for g in gens(I)
#         cvzip = zip(coefficients(g.elem), monomials(g.elem))
#         M = Generic.MPolyBuildCtx(parent(g.elem))
#         for (c, v) in cvzip
#             push_term!(M, R(c), D(v))
#         end
#     end
# end


"""
	intersectionIdeal(Is::Ideal...)

Return the intersection of the D-ideals `Is[1]`, `Is[2]`, ...
"""
function intersection(I1::DIdeal{T}, I2::DIdeal{T}) where T <: AbstractDiffOp
    base_ring(I1) != base_ring(I2) && throw(DomainError("Base rings must be the same", string(base_ring(I1)), string(base_ring(I2))))

    var_strs = I1 |> base_ring |> gens .|> string

    Dt, v, dv = weyl_algebra([var_strs; ["t1", "t2"]])
    t = v[end-1:end]
    dt = dv[end-1:end]

    genJ = [reduce(+, t)-1]
    append!(genJ, [t[1]*g for g in Dt.(gens(I1))])
    append!(genJ, [t[2]*g for g in  Dt.(gens(I2))])
end
# function intersection(Is::DIdeal{T}...) where T <: AbstractDiffOp
	# m = length(Is)
    # var_strs = base_ring(Is[1]) |> base_ring |> symbols .|> string
    # n = length(var_strs) 

    # Dt, v, dv = weyl_algebra([var_strs; ["t"*string(i) for i in 1:m]])
    # t = v[n+1:end]
    # dt = dv[n+1:end]

	# t, dt, v2d_t = genVars("t", m)
	# v2d_all = Bijection{Num, Num}()

	# genJ = [reduce(+, t)-1]
    # for i = 1:m
    #     append!(genJ, gens(Is[i]))
    # end
    # @show genJ
	# for i = 1:m
	# 	append!(genJ, Is[i].gens .* t[i])
	# 	for p in Is[i].v2d
	# 		!haskey(v2d_all, p.first) && (v2d_all[p.first] = p.second)
	# 	end
	# end

	# vars = collect(v2d_all.domain)
	# diffops = [v2d_all[v] for v in vars]

	# asir_cmd = 
	# """
	# load("nk_restriction.rr")\$
	# V = [$(vec2str(t, vars, dt, diffops))]\$
	# M = nk_restriction.make_elim_order_mat($m, $(length(vars)))\$
	# J = [$(vec2str(genJ))]\$
	# GB = nd_weyl_gr(J, V, 0, M);
	# """

	# asir_res = asir_cmd |> runAsir |> parseAsir
	# asir_res = filter!((s)->(startswith(s, "[")), asir_res)

	# (length(asir_res) != 1) && return nothing

	# vars_list = cat(t, vars, dt, diffops; dims=1)
	# elimOrdGens = evalAsir(asir_res[1], vars_list)

	# notHasTmpVar(diffOp) = (get_variables(diffOp), [t; dt]) .|> Set |> (s->intersect(s...)) |> isempty

	# return DIdeal(filter(notHasTmpVar, elimOrdGens), Dict((v, d) for (v,d) in zip(vars, diffops)) |> Bijection)
# end

function genNo2str_dict(D::AbstractDORing)
    bj = Bijection{Integer, String}()
    for (i, g) in enumerate(gens(D))
        bj[i] = string(g)
    end
    return bj
end

function coercion_homomorphism(D1::AbstractDORing, D2::AbstractDORing)
    bj1 = genNo2str_dict(D1) 
    bj2 = genNo2str_dict(D2)

    common_vars_str = intersect(bj1.range, bj2.range)
    common_vars_str |> isempty && throw(DomainError("Cannot determine homomorphism from " * string(D1) * " to " * string(D2)))

    hom = Bijection{Integer, Integer}()

    for s in common_vars_str
        hom[bj1(s)] = bj2(s)
    end

    return hom
end

function _coerce_unsafe(x::MPolyRingElem, M::MPolyRing, index_map::Dict{Integer, <:Integer})
    n = length(index_map)
    cezip = zip(coefficients(x), exponent_vectors(x))
    C = Generic.MPolyBuildCtx(M)
    for (c,e) in cezip
        push_term!(C, c,[get(e, index_map[i], 0) for i in 1:n] )
    end
    return finish(C)
end

function coerce(x::WAlgElem, D::WeylAlgebra)
    hom = coercion_homomorphism(parent(x), D) 
    n = nvars(D)

    # define mapping from index of D to index of parent(x)
    index_map = Dict{Integer, Integer}()
    invhom = inv(hom)

    for i = 1:n
        index_map[i] = get(invhom, i, 0)
    end

    cezip = zip(coefficients(unwrap(x)), exponent_vectors(unwrap(x)))
    C = Generic.MPolyBuildCtx(unwrap(D))
    M = base_ring(D)
    for (c, e) in cezip
        coerced_c = _coerce_unsafe(c, M, index_map)
        push_term!(C, coerced_c, [get(e, index_map[i], 0) for i in 1:n])
    end
    WAlgElem(D, finish(C))
end

(D::WeylAlgebra)(x::WAlgElem) = coerce(x, D)


