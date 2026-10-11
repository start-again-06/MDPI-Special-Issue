using Documenter
using MiniChannelFlow

DocMeta.setdocmeta!(MiniChannelFlow, :DocTestSetup, :(using MiniChannelFlow); recursive = true)

makedocs(;
    modules = [MiniChannelFlow],
    sitename = "MiniChannelFlow.jl",
    authors = "Anjan",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true",
                             canonical = "https://start-again-06/MiniChannelFlow.jl"),
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Validation" => "validation.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

deploydocs(; repo = "github.com/start-again-06/MiniChannelFlow.jl.git", devbranch = "main")
