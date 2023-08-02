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
julia> vars(x+y)
2-element Vector{PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 x
 y

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
julia> dvars(dx+dy)
2-element Vector{PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 dx
 dy

 julia> dvars(dx+dy)
2-element Vector{PfaffianSystems.DORElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}}:
 dx
 dy

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
Base.:*(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp = T(parent(y), x * unwrap(y))
Base.:*(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) * y)

Base.:+(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) + unwrap(y))
Base.:-(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) - unwrap(y))

function Base.:*(l::T, r::T) where T <: AbstractDiffOp
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

# TODO: case for y < 0
function diff_op_pow(x::T, y::Integer) where T <: AbstractDiffOp
	y == 0 && return one(x)

    ret_dop = x
    for _ = 1:y-1
        ret_dop *= x
    end
    return ret_dop
end

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

    for i=1:n
        f = derivative(f,x)
    end
    return f

end

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

# """
# 	addVars(name::AbstractString[, n::Integer, v2d::Bijection])

# Generate a Symbolics variable `name` and its corresponding differential operator d`name`. When `n` is given, generate vectors [`name`1...`name`n] and [d`name`1...d`name`n]. 
# The correspondence between `name` and d`name` is added to the bijection `v2d`. 

# # Examples

# ```jldoctest
# julia> x, dx, v2d = genVars("x")
# (x, dx, Bijection{Symbolics.Num,Symbolics.Num} (with 1 pairs))

# julia> y, dy, v2d = addVars("y", 2, v2d)
# (Symbolics.Num[y1, y2], Symbolics.Num[dy1, dy2], Bijection{Symbolics.Num,Symbolics.Num} (with 3 pairs))
# ```
# """
# function addVars(name::AbstractString, v2d::Bijection{Num, Num})
# 	var_ex = Symbol(name)
# 	diffop_ex = Symbol("d", var_ex)
# 	var, diffop = @variables $var_ex::Rational, $diffop_ex::Rational
# 	v2d[var] = diffop
# 	return var, diffop, v2d
# end

# function addVars(name::AbstractString, n::Integer, v2d::Bijection{Num, Num})
# 	vars = Vector{Num}(undef, n)
# 	diffops = Vector{Num}(undef, n)

# 	for i = 1:n
# 		vars[i], diffops[i], v2d = addVars(name*string(i), v2d)
# 	end

# 	return vars, diffops, v2d
# end

# """
# 	genVars(name::AbstractString[, n::Integer])

# Generate a Symbolics variable `name` and its corresponding differential operator d`name`. When `n` is given, generate vectors [`name`1...`name`n] and [d`name`1...d`name`n]. 
# The correspondence between `name` and d`name` is returned as a new bijection `v2d`, which will be provided to `addVars`. 

# # Examples

# ```jldoctest
# julia> x, dx, v2d = genVars("x")
# (x, dx, Bijection{Symbolics.Num,Symbolics.Num} (with 1 pairs))

# julia> y, dy, v2d = addVars("y", 2, v2d)
# (Symbolics.Num[y1, y2], Symbolics.Num[dy1, dy2], Bijection{Symbolics.Num,Symbolics.Num} (with 3 pairs))
# ```
# """
# genVars(name::AbstractString) = addVars(name, Bijection{Num, Num}())
# genVars(name::AbstractString, n::Integer) = addVars(name, n, Bijection{Num, Num}())

# """
# 	apply_dmon(DOmon::AbstractTerm, F::Num, p2s::Bijection, v2d::Bijection)

# Apply a term of differential operator `DOterm` to an expression `F`. 
# The bijection `p2s` relates a MulltivariatePolunomials expression to a SymbolicUtils one, and `d2v` does a differential operator to its corresponding variable. 
# """
# function apply_doterm(DOterm::AbstractTerm, F::Num, p2s::Bijection, v2d::Bijection{Num, Num}; use_asir=false)
# 	d2v = inv(v2d)
# 	coef = coefficient(DOterm)
# 	mon = monomial(DOterm)
# 	diffops = variables(mon) .|> (s->p2s[s])
# 	exps = exponents(mon)

# 	retF = F
# 	for (diffop, e) in zip(diffops, exps)
# 		if haskey(d2v, diffop)
# 			for k = 1:e
# 				retF = use_asir ? asir_derivative(retF, d2v[diffop]) : derivative(retF, d2v[diffop])
# 			end
# 		else
# 			coef = diffop^e*coef
# 		end
# 	end
# 	return coef*retF
# end
# function apply_doterm(DOterm::PolyForm, F::Num, v2d::Bijection{Num, Num}; use_asir=false)
# 	@assert length(terms(DOterm.p)) == 1 "Error: Differnetial operator has more than 1 term"
# 	apply_doterm(terms(DOterm.p)[1], F, DOterm.pvar2sym, v2d; use_asir=use_asir)
# end
# apply_doterm(DOterm::Num, F::Num, v2d::Bijection{Num, Num}; use_asir=false) = apply_doterm(DOterm |> value |> PolyForm, F, v2d; use_asir=use_asir)


# function apply_do(DiffOp::BasicSymbolic, F::Num, v2d::Bijection{Num, Num}; use_asir=false)
# 	dp = PolyForm(DiffOp)
# 	p2s = dp.pvar2sym
# 	retF::Num = 0
# 	for DOterm in terms(dp.p)
# 		# retF = retF + apply_doterm(DOterm, F, p2s, v2d; use_asir=use_asir)
# 		retF += apply_doterm(DOterm, F, p2s, v2d; use_asir=use_asir)
# 	end
# 	return retF
# end
# function apply_do(DiffOp::Union{Integer, AbstractFloat}, F::Num, v2d::Bijection{Num, Num}; use_asir=false)
# 	return DiffOp*F
# end
# apply_do(DiffOp::Num, F::Num, v2d::Bijection{Num, Num}; use_asir=false) = apply_do(DiffOp |> value, F, v2d; use_asir=use_asir)

# function dmul(dol::Num, dor::BasicSymbolic, v2d::Bijection{Num, Num}; use_asir=false)
# 	dor_pf = PolyForm(dor)
# 	p2s = dor_pf.pvar2sym
# 	retDO = dor*dol
# 	for dor_term in terms(dor_pf.p)
# 		# retDO += Num(dor_term)*dol
# 		coef = coefficient(dor_term)
# 		mon = monomial(dor_term)
# 		syms = variables(mon) .|> (s->p2s[s])
# 		exps = exponents(mon)
# 		varIdx = map(syms) do s s in v2d.domain end
# 		var_term = reduce(*, [Num(s^e) for (s, e) in zip(syms[varIdx], exps[varIdx])])

# 		retDO += coef*apply_do(dol, Num(var_term), v2d; use_asir=use_asir)*reduce(*, [Num(s^e) for (s, e) in zip(syms[.!varIdx], exps[.!varIdx])])
# 	end
# 	return retDO
# end
# dmul(dol::Num, dor::Num, v2d::Bijection{Num, Num}; use_asir=false) = dmul(dol, value(dor), v2d; use_asir=use_asir)
# dmul(dol, dor, v2d::Bijection{Num, Num}; use_asir=false) = dol*dor