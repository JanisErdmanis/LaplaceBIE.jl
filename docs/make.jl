using LaplaceBIE
using Documenter
using Literate
#using DocumenterLaTeX


#EXAMPLES = joinpath(@__DIR__, "..", "examples")
#OUTPUT = joinpath(@__DIR__, "src/generated")

Literate.markdown(joinpath(@__DIR__, "../examples/homogenous.jl"), joinpath(@__DIR__,"src/"); credit = false, name = "homogenous")

Literate.markdown(joinpath(@__DIR__, "../examples/pointlike.jl"), joinpath(@__DIR__,"src/"); credit = false, name = "pointlike")

#Literate.markdown(joinpath(@__DIR__, "../examples/droplet.jl"), joinpath(@__DIR__,"src/"); credit = false, name = "droplet")


makedocs(sitename="LaplaceBIE.jl",pages = ["index.md","homogenous.md","pointlike.md"])
#makedocs(format = DocumenterLaTeX.LaTeX(),sitename="LaplaceBIE",pages = ["index.md","homogenous.md","pointlike.md"])


#Literate.markdown("../examples/README.jl", "."; documenter=false)

deploydocs(
    repo = "github.com/akels/LaplaceBIE.jl.git",
)
