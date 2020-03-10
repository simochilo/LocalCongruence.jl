module LocalCongruences
	using DataStructures
	using LinearAlgebraicRepresentation
	using NearestNeighbors
	using SparseArrays
	Lar = LinearAlgebraicRepresentation

	include("cea-AA.jl")
	include("cea-SM.jl")
#	include("cea-algorithm.jl")

	export chainCongruence, chainCongruenceAA
end
