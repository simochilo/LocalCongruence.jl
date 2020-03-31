push!(LOAD_PATH,"../src/")

using Documenter, LocalCongruences

makedocs(
	format = :html,
	sitename = "LocalCongruences.jl",
	assets = ["assets/LocalCongruences.css", "assets/logo.png"],
	pages = [
		"1 - Home" => "index.md",
		"2 - Getting Started" => "start.md",
		"3 - Complexes and Chain Operators" => "theory.md",
		"4 - Graph Blas Introduction" => "graph_blas.md",
		"Implementations" => [
			"5.1 - Array of Arrays" => "aa_implementation.md",
			"5.2 - Native Sparse Matrices" => "sm_implementation.md",
			"5.3 - Graph Blas" => "gb_implementation.md",
		],
		"Examples" => [
			"6.1 - Cube Grids" => "example_1.md",
			"6.2 - ..." => "example_2.md",
			"6.3 - ..." => "example_3.md",
		],
		"A - About the Authors" => "authors.md",
		"B - Bibliography" => "bibliography.md"
	]
)
