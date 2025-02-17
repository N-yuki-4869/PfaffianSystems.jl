module PfaffianSystems

import Base: sort, ==, hash, parent
import AbstractAlgebra: nvars, gens, base_ring, derivative, vars, coefficients, monomials, exponent_vectors, elem_type, evaluate

using Bijections
using DataStructures: OrderedSet
using Base: @invokelatest
using AbstractAlgebra


const AA = AbstractAlgebra
const RatFuncElem = Generic.Rat
const MPoly = Generic.MPoly

function Bijection{S, T}(dict::AbstractDict{S, T}) where S where T
	return Bijection(dict)
end

abstract type AbstractDiffOp end
abstract type AbstractDORing end

abstract type AbstractDiffOp end
abstract type AbstractDORing end

include("DiffOps.jl")
export gens, dgens, base_ring, nvars, vars, dvars, isvar, isdvar, evaluate

include("WeylAlgebra.jl")
export weyl_algebra, coerce, elem_type

include("DiffOpRings.jl")
export diff_op_ring, coerce, elem_type, wnormalform, normalform, leading_term, pfaffian_system, pfaffian_system2	

include("DIdeals.jl")
export DIdeal, intersection_DIdeals, restriction_DIdeal, integration_DIdeal, multiplication_DIdeal, stdmon, pfaffian_system2

end
