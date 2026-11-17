using Documenter
using GromovWittenFlags
using PartialFlagVarieties

DocMeta.setdocmeta!(GromovWittenFlags, :DocTestSetup,
  :(using GromovWittenFlags, PartialFlagVarieties); recursive=true)

makedocs(;
  sitename="GromovWittenFlags.jl",
  authors="Pieter Belmans",
  modules=[GromovWittenFlags],
  remotes=nothing,
  checkdocs=:all,
  doctest=true,
  format=Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    repolink="https://github.com/HomogeneousTools/GromovWittenFlags.jl",
    edit_link=nothing,
  ),
  pages=[
    "Home" => "index.md",
    "Getting started" => "getting_started.md",
    "Conventions and scope" => "conventions.md",
    "Certified computation" => "certified.md",
    "Implementation" => "implementation.md",
    "API reference" => "api.md",
    "Internals" => "internals.md",
  ],
)

if get(ENV, "GWF_DEPLOY_DOCS", "false") == "true"
  deploydocs(;
    repo="github.com/HomogeneousTools/GromovWittenFlags.jl.git", devbranch="main"
  )
end
