
# to set the Package :
# using Pkg
# Pkg.instantiate("ShastryGS")

# include("src/ShastryGS.jl")
using ShastryGS
using LinearAlgebra
using ITensors

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.enable_threaded_blocksparse()

let 

    Jy_vector::Vector{Float64} = vcat(0:-0.1:-2) 
    Jx_vector::Vector{Float64} = vcat(-2:0.1:0) 
    
    for i = 1:length(Jy_vector)

    if !isempty(ARGS)
        D = parse(Int64, ARGS[1])
        J1 = parse(Float64, ARGS[2])
        J2 = parse(Float64, ARGS[3])
        Delta = parse(Float64, ARGS[4])
        hx = parse(Float64, ARGS[5])
        hz = parse(Float64, ARGS[6])
        N = parse(Int64, ARGS[7])
        Jx = parse(Float64, ARGS[8])
        Jy = parse(Float64, ARGS[9])
        Jz = parse(Float64, ARGS[10])
        model = ARGS[11]
    else
        N::Int64 = 2;
        D::Int64 = 3; 
        J1::Float64 = 1.0;
        J2::Float64 = 3.4;
        Delta::Float64 = 0;
        hx::Float64 = 0.0;
        hz::Float64= 0;
        Jx::Float64= Jy_vector[i];
        Jy::Float64= Jx_vector[i];
        Jz::Float64= 1;
        model::String="XYZ"
    end
    
    dbeta::Float64 = 1e-2
    modit::Int64 = 100
    parameters = Dict("D"=>D, "J1"=>J1, "J2"=>J2, "hz"=>hz, "hx"=>hx, "Delta"=>Delta, "N"=>N, "model"=>model,"Jx"=>Jx,"Jy"=>Jy,"Jz"=>Jz, "dbeta"=>dbeta, "modit"=>modit)
    SU(parameters)

    end

end



