push!(LOAD_PATH,"../docs/src/")
push!(LOAD_PATH,"../src")

using Documenter, FlowCytometry

makedocs(sitename="FlowCytometry.jl",
pages = [
    "Home" => "index.md",
    "Usage/Usage.md",
    "API.md",
    "References.md",
],
format = Documenter.HTML(prettyurls = false)
)

deploydocs(
    repo = "github.com/gatocor/FlowCytometry.jl",
)