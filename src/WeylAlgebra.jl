import Base: ==, parent

############################################################
# 
#	WeylAlgebra / WAlgElem 
#
############################################################

# WeylAlgebra{T} = MPolyRing{MPolyRingElem{T}} where T

struct WeylAlgebra{T <: MPolyRing{<:MPolyRingElem}} <: AbstractDORing
	WAlg::T
end
unwrap(D::WeylAlgebra) = D.WAlg

Base.one(D::WeylAlgebra) = D |> unwrap |> one |> WAlgElem
Base.zero(D::WeylAlgebra) = D |> unwrap |> zero |> WAlgElem

base_ring(D::WeylAlgebra) = D |> unwrap |> base_ring
function gens(D::WeylAlgebra)
    g = D |> base_ring |> gens
    g = unwrap(D).(g)
    g .|> WAlgElem
end
dgens(D::WeylAlgebra) = D |> unwrap |> gens .|> WAlgElem
nvars(D::WeylAlgebra) = D |> unwrap |> nvars

function Base.show(io::IO, D::WeylAlgebra)
    print(io, nvars(D),"-dimensional Weyl algebra in [$(join(string.(gens(D)), ","))]")
end

# TODO: wrapper for constant


struct WAlgElem{T <: MPolyRingElem{<:MPolyRingElem}} <: AbstractDiffOp
	elem::T
end
unwrap(wae::WAlgElem) = wae.elem

Base.parent(wae::WAlgElem) = wae |> unwrap |> parent |> WeylAlgebra
gens(wae::WAlgElem) = wae |> parent |> gens
dgens(wae::WAlgElem) = wae |> parent |> dgens

function Base.show(io::IO, wae::WAlgElem)
	show(io, unwrap(wae))
end

Base.:+(x::WAlgElem, y::WAlgElem) = WAlgElem(unwrap(x) + unwrap(y))
Base.:-(x::WAlgElem, y::WAlgElem) = WAlgElem(unwrap(x) - unwrap(y))
Base.one(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> one |> WAlgElem
Base.zero(wae::Union{Type{WAlgElem{T}}, WAlgElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> zero |> WAlgElem

Base.:(==)(x::WAlgElem, y::WAlgElem) = unwrap(x) == unwrap(y)

# TODO: multiplication of WAlgElem
function Base.:*(l::WAlgElem, r::WAlgElem)
    # l_coeffs = l.elem |> coefficients |> collect                   
    # l_mons = l.elem |> monomials |> collect 
    # r_coeffs = r.elem |> coefficients |> collect    
    # r_mons = r.elem |> monomials |> collect
    l_coeffs = l |> unwrap |> coefficients 
    l_mons = l |> unwrap |> monomials
    r_coeffs = r |> unwrap |> coefficients
    r_mons = r |> unwrap |> monomials
    # l_size = l_coeffs |> size
    # r_size = r_coeffs |> size

    ret_dop = 0
    # for i = 1:l_size[1]
    #     for j = 1:r_size[1]
    #         ret_dop +=  l_coeffs[i] * Leibniz_rule(l_mons[i],r_coeffs[j]) * r_mons[j]
    #     end
    # end
    for (lc, lm) in zip(l_coeffs, l_mons), (rc, rm) in zip(r_coeffs, r_mons) 
        ret_dop +=  lc * Leibniz_rule(lm, rc) * rm
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

function _nth_derivative(f::T, x::T ,n::Integer) where T
    n == 0 && return f 
    # if n==0
    #     return f
    # else
    for i=1:n
        f = derivative(f,x)
    end
    return f
    # end
end


function Leibniz_rule(l_mons::T, r_coeffs::U) where {T <: MPolyRingElem{<:MPolyRingElem}, U <: MPolyRingElem}
    ret_dop = r_coeffs * l_mons
    # variable = size(gens(parent(l_mons)))[1]
    n = nvars(parent(l_mons))
    # for i=1:n
        # coeffs = collect(coefficients(ret_dop))
        # mons = collect(monomials(ret_dop))
        # coeffs_size = size(coeffs)[1]
        # @show coeffs, mons
        # for j=1:coeffs_size
            # @show mons[j], coeffs[j], i
            # ret_dop += Leibniz_rule_1(mons[j],coeffs[j],i)
            # @show m, c, i
    # end
    for i=1:n, (m, c) in zip(monomials(ret_dop), coefficients(ret_dop))
        ret_dop += Leibniz_rule_1(m,c,i)
    end
    return ret_dop
end


function Leibniz_rule_1(l_mons::T ,r_coeffs::U, i::Integer) where {T <: MPolyRingElem{<:MPolyRingElem}, U <: MPolyRingElem}
    ret_dop = 0
    k = 1
    coeff_var = gens(parent(r_coeffs))[i]
    mon_var = gens(parent(l_mons))[i]
    while true
        # a = _nth_derivative(r_coeffs,gens(parent(r_coeffs))[i],k) * _nth_derivative(l_mons, gens(parent(l_mons))[i],k) / parent(r_coeffs)(factorial(big(k)))
        a = _nth_derivative(r_coeffs, coeff_var,k) * _nth_derivative(l_mons, mon_var,k) / parent(r_coeffs)(factorial(big(k)))
        a == 0 && break
        ret_dop += a
        k += 1
    end
    return ret_dop
end


# TODO: multiplication of WAlgElem and constant
# TODO: multiplication of WAlgElem and polynomial

# TODO: coefficients of WAlgElem 
# TODO: monomials of WAlgElem 
# TODO: terms of WAlgElem

# TODO: make polynomials coercible
function change_ring(wae::WAlgElem, D::AbstractDORing)
    R = base_ring(D)
    unwrapD = unwrap(D)
    cmzip = zip(coefficients(unwrap(wae)), monomials(unwrap(wae)))
    ret_dop = unwrapD(0)
    for (c, m) in cmzip
        # ret_dop += R(c) * unwrapD(m)
    end
    return D(ret_dop)
end

############################################################
# 
# Arithemetic with Rationals and Integers
#
############################################################
Base.:+(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(x + unwrap(y))
Base.:+(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(unwrap(x) + y)
Base.:-(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(x - unwrap(y))
Base.:-(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(unwrap(x) - y)
Base.:*(x::Union{Rational, Integer}, y::WAlgElem) = WAlgElem(x * unwrap(y))
Base.:*(x::WAlgElem, y::Union{Rational, Integer}) = WAlgElem(unwrap(x) * y)

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

mutable struct DIdeal{T <: AbstractDiffOp} <: Ideal{T}
    base_ring::AbstractDORing
    gens::Vector{T}

    function DIdeal{T}(R::AbstractDORing, gens::Vector) where T <: AbstractDiffOp
        if !all([R == parent(g) for g in gens])
            # error("All generators must be elements of the same Weyl algebra")
            throw(DomainError("All generators must be included in", R))
        end
        new{T}(R, gens)
    end 
end

function DIdeal(R::AbstractDORing, gens::Vector{T}) where T <: AbstractDiffOp
    DIdeal{eltype(gens)}(R, gens)
end

base_ring(I::DIdeal) = I.base_ring
gens(I::DIdeal) = I.gens

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
function intersection(Is::DIdeal{T}...) where T <: AbstractDiffOp
	m = length(Is)
    var_strs = base_ring(Is[1]) |> base_ring |> symbols .|> string
    n = length(var_strs) 

    Dt, v, dv = weyl_algebra([var_strs; ["t"*string(i) for i in 1:m]])
    t = v[n+1:end]
    dt = dv[n+1:end]

	# t, dt, v2d_t = genVars("t", m)
	# v2d_all = Bijection{Num, Num}()

	genJ = [reduce(+, t)-1]
    for i = 1:m
        append!(genJ, gens(Is[i]))
    end
    @show genJ
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
end




