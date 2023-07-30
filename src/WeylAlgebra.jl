
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
(D::WeylAlgebra)(num::Union{Rational, Integer}) = WAlgElem(D, unwrap(D)(num))

# Base.one(D::WeylAlgebra) = D |> unwrap |> one |> WAlgElem
# Base.zero(D::WeylAlgebra) = D |> unwrap |> zero |> WAlgElem
Base.one(D::T) where T <: AbstractDORing = D(1)
Base.zero(D::T) where T <: AbstractDORing = D(0)

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
    print(io, nvars(D),"-d Weyl algebra in [$(join(string.(gens(D)), ","))]")
end

# TODO: wrapper for constant


struct WAlgElem{T <: MPolyRingElem{<:MPolyRingElem}} <: AbstractDiffOp
    parent::WeylAlgebra
	elem::T
end
# unwrap(wae::WAlgElem) = wae.elem


Base.:(==)(x::WeylAlgebra, y::WeylAlgebra) = unwrap(x) == unwrap(y)

############################################################
# 
# weyl_algebra constructor
# 
############################################################

"""
	WeylAlgebra
"""
function weyl_algebra(F::Field, s::AbstractVector{<:AbstractString}; kw...)
    D = WeylAlgebra(F, Symbol.(s))
    return D, gens(D), dgens(D)
end
weyl_algebra(s::AbstractVector{<:AbstractString}; kw...) = weyl_algebra(QQ, s; kw...)

function weyl_algebra(F::Field, s::AbstractString; kw...)
    D = WeylAlgebra(F, [Symbol(s)])
    return D, gens(D)[1], dgens(D)[1]
end
weyl_algebra(s::AbstractString; kw...) = weyl_algebra(QQ, s; kw...)

# TODO: constructor of n-dimensional Weyl algebra
# TODO: make new Weyl algebra with additional variables
function weyl_algebra(F::Field, D::WeylAlgebra, new_vars::AbstractVector{<:AbstractString}; kw...)
    old_vars = gens(D) .|> string
    if !(:new_vars_pos in keys(kw)) || kw[:new_vars_pos] == :append
        vars = [old_vars; new_vars]
    elseif kw[:new_vars_pos] == :prepend
        vars = [new_vars; old_vars]
    else
        throw(ArgumentError("new_vars_pos must be :prepend or :append"))
    end
    De, g, dg = weyl_algebra(F, vars; kw...)
    return De, g, dg
end
weyl_algebra(D::WeylAlgebra, s::AbstractVector{<:AbstractString}; kw...) = weyl_algebra(QQ, D, s; kw...)
# TODO: make new Weyl algebra with a part of variables


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

function Base.show(io::IO, I::DIdeal)
    print(io, "Ideal of ", I.base_ring, " generated by [", join(string.(I.gens), ", "), "]")
end

function (D::WeylAlgebra)(I::DIdeal)
    DIdeal(D, I |> gens .|> (s->coerce(s, D)))
end


"""
	intersectionIdeal(Is::Ideal...)

Return the intersection of the D-ideals `Is[1]`, `Is[2]`, ...
"""
# function intersection(I1::DIdeal{T}, I2::DIdeal{T}) where T <: AbstractDiffOp
function Dintersection(Is::DIdeal{T}...) where T <: AbstractDiffOp
    # base_ring(I1) != base_ring(I2) && throw(DomainError("Base rings must be the same", string(base_ring(I1)), string(base_ring(I2))))
    R = base_ring(Is[1])
    m = length(Is)
    all([R == base_ring(I) for I in Is]) || throw(DomainError("Base rings must be the same"))

    Dt, vall, dvall = weyl_algebra(R, ["t"*string(i) for i in 1:m]; new_vars_pos=:prepend)
    t =  vall[1:end-m]
    v = vall[end-m+1:end]
    # dv = dvall[1:end-m]
    # t = vall[end-m+1:end]
    dt = dvall[1:end-m]

    genJ = [reduce(+, t)-1]
    for (i, I) in enumerate(Is)
        append!(genJ, [t[i]*g for g in Dt.(gens(I))])
    end

	asir_cmd = 
	"""
	load("nk_restriction.rr")\$
	V = [$(vec2str(vall, dvall))]\$
	M = nk_restriction.make_elim_order_mat($m, $(length(v)))\$
	J = [$(vec2str(genJ))]\$
	GB = nd_weyl_gr(J, V, 0, M);
	"""
	asir_res = asir_cmd |> runAsir |> parseAsir
	asir_res = filter!((s)->(startswith(s, "[")), asir_res)

	(length(asir_res) != 1) && throw(DomainError("Invalid result from Asir", asir_res))

	elimOrdGens = evalAsir(asir_res[1], [vall; dvall])

	notHasTmpVar(dop) = (vars(dop), [t; dt]) |> (s-> intersect(s...)) |> isempty

	return R(DIdeal(filter(notHasTmpVar, elimOrdGens)))
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


