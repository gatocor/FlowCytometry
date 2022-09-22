push!(LOAD_PATH,"../docs/src/")
push!(LOAD_PATH,"../src")

using Documenter, FlowCytometry

makedocs(sitename="FlowCytometry.jl",
pages = [
    "Home" => "index.md",
    "API.md"
],
format = Documenter.HTML(prettyurls = false)
)

deploydocs(
    repo = "github.com/dsb-lab/FlowCytometry.jl",
)