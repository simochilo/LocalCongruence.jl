push!(LOAD_PATH,"../src/")

using Documenter, LocalCongruences

makedocs(
	format = :html,
	sitename = "LocalCongruences.jl",
	assets = ["assets/LocalCongruences.css", "assets/logo.png"],
	pages = [
		"1 - Home" => "index.md",
		"A - About the Authors" => "authors.md",
		"B - Bibliography" => "bibliography.md"
	]
)
