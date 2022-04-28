using ConvexPolygons
using Documenter

DocMeta.setdocmeta!(ConvexPolygons, :DocTestSetup, :(using ConvexPolygons); recursive=true)

makedocs(;
    modules=[ConvexPolygons],
    authors="Paulo Jabardo <pjabardo@ipt.br>",
    repo="https://github.com/pjsjipt/ConvexPolygons.jl/blob/{commit}{path}#{line}",
    sitename="ConvexPolygons.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
