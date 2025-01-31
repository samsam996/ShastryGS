

module ShastryGS
    
    # Dependencies
    using MKL
    using FileIO
    using JLD2
    using KrylovKit
    using LinearAlgebra
    using ITensors
    using Printf
    using MAT

    # functions
    include("structures.jl")
    include("SSMHamiltonian.jl")
    include("ctmrg_evolution/ctm.jl")
    include("simple_update_evolution/diag_sqrt.jl")

    module movesNabove6
        using ITensors
        using ..ShastryGS
        include("ctmrg_evolution/movesNabove6/OBC_PEPS.jl")
        include("ctmrg_evolution/invert_diag_sqrt.jl")
        include("ctmrg_evolution/find_value.jl")
        include("ctmrg_evolution/movesNabove6/LeftMove.jl")
        include("ctmrg_evolution/movesNabove6/RightMove.jl")
        include("ctmrg_evolution/movesNabove6/UpMove.jl")
        include("ctmrg_evolution/movesNabove6/DownMove.jl")
        include("get_tens.jl")
        export LeftMove!, RightMove!, UpMove!, DownMove!, OBC_PEPS, get_tens
    end

    module movesN2
        using ITensors
        using ..ShastryGS
        include("ctmrg_evolution/movesN2/OBC_PEPS.jl")
        include("ctmrg_evolution/invert_diag_sqrt.jl")
        include("ctmrg_evolution/find_value.jl")
        include("ctmrg_evolution/movesN2/LeftMove.jl")
        include("ctmrg_evolution/movesN2/RightMove.jl")
        include("ctmrg_evolution/movesN2/UpMove.jl")
        include("ctmrg_evolution/movesN2/DownMove.jl")
        include("get_tensN2.jl")
        export LeftMove!, RightMove!, UpMove!, DownMove!, OBC_PEPS, get_tens
    end


    # Dynamically decide which module to load
    function load_moves(N)
        if N == 2
            module_name = :movesN2  
            eval(Meta.parse("using .$(module_name)"))  
        else
            module_name = :movesNabove6  
            eval(Meta.parse("using .$(module_name)"))  
        end
        return nothing
    end

    include("simple_update_evolution/simpleUpdateJ1.jl")
    include("observables/entropy_lambda.jl")
    include("observables/correlation_length.jl")
    include("observables/magnetisation.jl")
    include("observables/energy.jl")
    include("observables/correlation_length.jl")
    include("simpleupdate.jl")


    export SU 

end