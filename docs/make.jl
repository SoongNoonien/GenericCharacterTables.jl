using Documenter
using DocumenterCitations
using GenericCharacterTables

include("citation_style.jl")

DocMeta.setdocmeta!(GenericCharacterTables, :DocTestSetup, :(using GenericCharacterTables, Oscar); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=oscar_style,
)

makedocs(
	sitename="GenericCharacterTables.jl",
	format=Documenter.HTML(
	    prettyurls=get(ENV, "CI", "false") == "true",
        assets=[
            "assets/citations.css",
        ],
    ),
	doctest = true,
	pages=[
		"index.md",
		"showinfo.md",
		"calculations.md",
		"modify.md",
		"book.md",
		"unexported.md",
		"references.md",
	],
	plugins=[bib],
)

deploydocs(
   repo   = "github.com/oscar-system/GenericCharacterTables.jl.git",
   target = "build",
   deps = nothing,
   make   = nothing,
)
