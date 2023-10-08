module PfaffianSystems

import Base: sort, ==, hash, parent
import AbstractAlgebra: nvars, gens, base_ring, derivative, vars, coefficients, monomials, exponent_vectors

using Bijections
using DataStructures: OrderedSet
using Base: @invokelatest
using AbstractAlgebra

const AA = AbstractAlgebra
const RatFuncElem = Generic.RationalFunctionFieldElem

function Bijection{S, T}(dict::AbstractDict{S, T}) where S where T
	return Bijection(dict)
end

abstract type AbstractDiffOp end
abstract type AbstractDORing end

include("AsirWrapper.jl")
export isAsirAvailable, vec2str, asir_derivative, asir_reduce, asir_fctr

include("DiffOps.jl")
export gens, dgens, base_ring, nvars, vars, dvars 

include("WeylAlgebra.jl")
export weyl_algebra, coerce, elem_type

include("DiffOpRings.jl")
export diff_op_ring, coerce, elem_type, normalform, leading_term

include("DIdeals.jl")
export DIdeal, intersection_DIdeal, stdmon

end
