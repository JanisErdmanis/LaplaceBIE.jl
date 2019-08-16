using LaplaceBIE
using Documenter
using Literate

Literate.markdown(joinpath(@__DIR__, "../examples/homogenous.jl"), joinpath(@__DIR__,"src/"); credit = false, name = "homogenous") #, preprocess = replace_includes)

Literate.markdown(joinpath(@__DIR__, "../examples/pointlike.jl"), joinpath(@__DIR__,"src/"); credit = false, name = "pointlike")

Literate.markdown(joinpath(@__DIR__, "../examples/mdrop.jl"), joinpath(@__DIR__,"src/"); credit = false, name = "mdrop")

cp(joinpath(@__DIR__,"../examples/sphere.jl"),joinpath(@__DIR__,"src/sphere.jl"))
makedocs(sitename="LaplaceBIE.jl",pages = ["index.md","homogenous.md","pointlike.md","mdrop.md"])

deploydocs(
    repo = "github.com/akels/LaplaceBIE.jl.git",
)
