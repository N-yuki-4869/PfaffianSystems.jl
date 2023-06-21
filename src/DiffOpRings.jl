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


Base.one(D::DiffOpRing) = DORElem(one(D.DOR))
Base.zero(D::DiffOpRing) = DORElem(zero(D.DOR))

function gens(D::DiffOpRing)
	gens = D.DOR |> base_ring |> AA.gens
	gens = D.DOR.(gens)
	return DORElem.(gens)
end
dgens(D::DiffOpRing) = D.DOR |> AA.gens .|> DORElem

nvars(D::DiffOpRing) = D.DOR |> AA.nvars

function Base.show(io::IO, D::DiffOpRing)
	print(io, nvars(D), "-dimensional Diff op rings")
end


struct DORElem{T <: MPolyRingElem{<:RatFuncElem}} <: AbstractDiffOp
	elem::T
end


Base.parent(wae::DORElem) = DiffOpRing(parent(wae.elem))
gens(wae::DORElem) = gens(parent(wae))
dgens(wae::DORElem) = dgens(parent(wae))

function Base.show(io::IO, wae::DORElem)
	print(io, wae.elem)
end

Base.:+(x::DORElem, y::DORElem) = DORElem(x.elem + y.elem)
Base.:-(x::DORElem, y::DORElem) = DORElem(x.elem - y.elem)
Base.one(wae::Union{Type{DORElem{T}}, DORElem{T}}) where T <: MPolyRingElem = DiffOpRing(one(wae.elem))
Base.zero(wae::Union{Type{DORElem{T}}, DORElem{T}}) where T <: MPolyRingElem = DiffOpRing(zero(wae.elem))

Base.:(==)(x::DORElem, y::DORElem) = x.elem == y.elem


function Base.:*(l::DORElem, r::DORElem)
    l_coeffs = l.elem |> AA.coefficients |> collect                   
    l_mons = l.elem |> AA.monomials |> collect 
    r_coeffs = r.elem |> AA.coefficients |> collect    
    r_mons = r.elem |> AA.monomials |> collect
    l_size = l_coeffs |> size
    r_size = r_coeffs |> size

    ret_dop = 0
    for i = 1:l_size[1]
        for j = 1:r_size[1]
            ret_dop +=  l_coeffs[i] * Leibniz_rule(l_mons[i],r_coeffs[j]) * r_mons[j]
        end
    end
    return DORElem(ret_dop)
end

function Base.:^(x::DORElem, y::Integer)
    ret_dop = x
    for _ = 1:y-1
        ret_dop *= x
    end
    return ret_dop
end



function AA.derivative(f::RatFuncElem ,x::RatFuncElem)
    f_nume = numerator(f)
    f_deno = denominator(f)
    x_nume = numerator(x)
    nume = parent(x)((AA.derivative(f_nume,x_nume)*f_deno - f_nume*AA.derivative(f_deno,x_nume)))
    deno = parent(x)((f_deno^2))
    #@show nume , deno, typeof(nume) , typeof(deno)
    ret_dop = nume // deno
    return ret_dop
end


function _nth_derivative(f::T ,x::T ,n::Integer) where T
    if n==0
        return f
    else
        for i=1:n
            #@show f,x , typeof(f) , typeof(x)
            f = AA.derivative(f,x)
        end
        return f
    end
end


function Leibniz_rule(l_mons::T ,r_coeffs::U) where {T <: MPolyRingElem{<:RatFuncElem}, U <: RatFuncElem}
    ret_dop = r_coeffs * l_mons
    variable = size(AA.gens(parent(l_mons)))[1]
    for i=1:variable
        coeffs = collect(AA.coefficients(ret_dop))
        mons = collect(AA.monomials(ret_dop))
        coeffs_size = size(coeffs)[1]
        for j=1:coeffs_size
            ret_dop += Leibniz_rule_1(mons[j],coeffs[j],i)
        end
    end
    return ret_dop
end


function Leibniz_rule_1(l_mons::T ,r_coeffs:: U, i::Integer) where {T <: MPolyRingElem{<:RatFuncElem}, U <: RatFuncElem}
    ret_dop = 0
    k = 1
    while true
        a = _nth_derivative(r_coeffs,AA.gens(parent(r_coeffs))[i],k) * _nth_derivative(l_mons, AA.gens(parent(l_mons))[i],k) / parent(r_coeffs)(factorial(big(k)))
        a == 0&&break
        ret_dop += a
        k += 1
    end
    return ret_dop
end


Base.:+(x::Union{Rational, Integer}, y::DORElem) = DORElem(x + y.elem)
Base.:+(x::DORElem, y::Union{Rational, Integer}) = DORElem(x.elem + y)
Base.:*(x::Union{Rational, Integer}, y::DORElem) = DORElem(x * y.elem)
Base.:*(x::DORElem, y::Union{Rational, Integer}) = DORElem(x.elem * y)


function Base.://(x::DORElem,y::DORElem)
    ret_dop = 0
    x_coeff = x.elem |> AA.coefficients |> collect                   
    x_mon = x.elem |> AA.monomials |> collect 
    y_coeff = y.elem |> AA.coefficients |> collect    
    y_mon = y.elem |> AA.monomials |> collect
    if size(y_mon)[1] !== 1
        return "Error"
    else
        for i=1:size(x_coeff)[1]
            ret_dop += (x_coeff[i]// y_coeff[1]) * x_mon[i]
        end
        return DORElem(ret_dop)
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

