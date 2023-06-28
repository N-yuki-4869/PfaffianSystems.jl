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
end
unwrap(R::DiffOpRing) = R.DOR

Base.one(R::DiffOpRing) = R |> unwrap |> one |> DORElem
Base.zero(R::DiffOpRing) = R |> unwrap |> zero |> DORElem

function gens(R::DiffOpRing)
	g = R |> unwrap |> base_ring |> gens
	g = unwrap(R).(g)
	return DORElem.(g)
end
dgens(R::DiffOpRing) = R |> unwrap |> gens .|> DORElem

nvars(R::DiffOpRing) = R |> unwrap |> nvars

function Base.show(io::IO, R::DiffOpRing)
	print(io, nvars(R), "-dimensional ring of differential opeartors in [$(join(string.(gens(R)), ","))]")
end

function (R::DiffOpRing)(num::Union{Rational, Integer})
    unwrap(R)(num) |> DORElem
end

struct DORElem{T <: MPolyRingElem{<:RatFuncElem}} <: AbstractDiffOp
	elem::T
end
unwrap(wae::DORElem) = wae.elem


Base.parent(wae::DORElem) = wae |> unwrap |> parent |> DiffOpRing
gens(wae::DORElem) = wae |> parent |> gens
dgens(wae::DORElem) = wae |> parent |> dgens

function Base.show(io::IO, wae::DORElem)
	print(io, unwrap(wae))
end

Base.:+(x::DORElem, y::DORElem) = DORElem(unwrap(x) + unwrap(y))
Base.:-(x::DORElem, y::DORElem) = DORElem(unwrap(x) - unwrap(y))
Base.one(wae::Union{Type{DORElem{T}}, DORElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> one |> DORElem
Base.zero(wae::Union{Type{DORElem{T}}, DORElem{T}}) where T <: MPolyRingElem = wae |> unwrap |> zero |> DORElem

Base.:(==)(x::DORElem, y::DORElem) = unwrap(x) == unwrap(y)


function Base.:*(l::DORElem, r::DORElem)
    l_coeffs = l |> unwrap |> coefficients |> collect                   
    l_mons = l |> unwrap |> monomials |> collect 
    r_coeffs = r |> unwrap |> coefficients |> collect    
    r_mons = r |> unwrap |> monomials |> collect
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


function _nth_derivative(f::T ,x::T ,n::Integer) where T
    if n==0
        return f
    else
        for i=1:n
            #@show f,x , typeof(f) , typeof(x)
            f = derivative(f,x)
        end
        return f
    end
end


function Leibniz_rule(l_mons::T ,r_coeffs::U) where {T <: MPolyRingElem{<:RatFuncElem}, U <: RatFuncElem}
    ret_dop = r_coeffs * l_mons
    variable = size(gens(parent(l_mons)))[1]
    for i=1:variable
        coeffs = collect(coefficients(ret_dop))
        mons = collect(monomials(ret_dop))
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
        a = _nth_derivative(r_coeffs,gens(parent(r_coeffs))[i],k) * _nth_derivative(l_mons, gens(parent(l_mons))[i],k) / parent(r_coeffs)(factorial(big(k)))
        a == 0 && break
        ret_dop += a
        k += 1
    end
    return ret_dop
end


Base.:+(x::Union{Rational, Integer}, y::DORElem) = DORElem(x + unwrap(y))
Base.:+(x::DORElem, y::Union{Rational, Integer}) = DORElem(unwrap(x) + y)
Base.:-(x::Union{Rational, Integer}, y::DORElem) = DORElem(x - unwrap(y))
Base.:-(x::DORElem, y::Union{Rational, Integer}) = DORElem(unwrap(x) - y)
Base.:*(x::Union{Rational, Integer}, y::DORElem) = DORElem(x * unwrap(y))
Base.:*(x::DORElem, y::Union{Rational, Integer}) = DORElem(unwrap(x) * y)

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
    for (c, e) in cezip
        @show c, typeof(c), e, typeof(e)
        c_nume = numerator(c)
        c_deno = denominator(c)
        ccezip_nume = zip(coefficient(c_nume), exponent_vectors(c_nume))
        ccezip_deno = zip(coefficient(c_deno), exponent_vectors(c_deno))
        ccezip = zip(coefficients(c), exponent_vectors(c))
        C = Generic.MPolyBuildCtx(base_ring(D))
        for (cc, ce) in ccezip
            @show cc, typeof(cc), ce, typeof(ce)
            push_term!(C, cc, [get(ce, index_map[i], 0) for i in 1:n])
        end
        coerced_c = finish(c)

        push_term!(M, coerced_c, [get(e, index_map[i], 0) fpr i in 1:n])
    end
    DORElem(finish(M))
end

(D::DiffOpRing)(x::DORElem) = coerce(x, D)