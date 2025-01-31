

using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors

N::Int64 = 2;

# if N == 2
#     include("twoxtwo.jl")
#     using .unit_cell
# elseif N == 6
#     include("sixxsix.jl")
#     using .unit_cell
# end

mutable struct lattice2x1
    A::ITensor
    B::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice2x1()
        name = fieldnames(lattice2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]
        new(ITensor(),ITensor(),g1,g2)
    end
    
    function lattice2x1(a::ITensor,b::ITensor)

        name = fieldnames(lattice2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

        new(a,b,g1,g2)
            
    end

end

mutable struct lattice_ind2x1
    A::Index
    B::Index
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}

    function lattice_ind2x1()
        
        name = fieldnames(lattice_ind2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

        new(ITensor(),ITensor(),g1,g2)
    end

    function lattice_ind2x1(a::Index,b::Index)

        name = fieldnames(lattice2x1)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

        new(a,b,g1,g2)

    end

    function lattice_ind2x1(nature_of_the_legs::String)

        if nature_of_the_legs == "physical"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ia"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ib"))
        elseif nature_of_the_legs == "ancilla"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sa"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sb"))
        elseif nature_of_the_legs == "imaginary"
            new(
            Index([QN(0)=>1]),Index([QN(0)=>1]))
        end

    end

end

println("BEGIN SIMPlE UPDATE")

include("SSMHamiltonian.jl")
include("simple_update_evolution/simpleUpdateJ1.jl")
include("observables/entropy_lambda.jl")
include("observables/correlation_length.jl")
include("observables/magnetisation.jl")
include("observables/energy.jl")
include("observables/PartitionPerSite.jl")

include("get_tens.jl")
include("ctmrg_evolution/ctm.jl")
include("initialisation.jl")

include("simpleupdate.jl")

# BLAS.set_num_threads(10)
# ITensors.Strided.set_num_threads(1)

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

let 
    
    D::Int64 = 5; 
    # J1::Float64 = 0.63;
    # J2::Float64 = 1.;

    J1::Float64 = 1
    J2::Float64 = 5

    h::Float64 = 0.; 

    hs::Float64 = 0; 

    # temperature = LinRange(10.21, 0.05, 100)
    temperature = LinRange(5, 1, 10)
    # sort!(temperature, rev = true) 

    # temperature = LinRange(4.9,3.88, 10)
    # temperature = LinRange(0.065,0.062,10)
    final_temp = temperature[end]-1e-4

    dbetasu = 1e-3;
    Gamma, lambdax, lambday, gt, gg, physical_legs, ancilla_legs, magne, energie, 
    tempe, beta, tmp, PPS = simpleupdate(final_temp,temperature,dbetasu,J1,J2,h,hs,N,D)

    @show real(energie)
    @show magne
    @show tempe
    @show real.(PPS)

    nothing
  
end



println("END SIMPlE UPDATE")



