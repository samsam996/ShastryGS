

using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors, Printf, MAT


# D = parse(Int64, ARGS[1])
# J1 = parse(Float64, ARGS[2])
# J2 = parse(Float64, ARGS[3])
# Delta = parse(Float64, ARGS[4])
# hx = parse(Float64, ARGS[5])
# hz = parse(Float64, ARGS[6])
# N = parse(Int64, ARGS[7])


function cyclic_matrix(name::NTuple,N)

    g2 = Matrix{typeof(name[1])}(undef, N, N)  
    for i in 1:N
        g2[i, :] = [name[(mod(j - i, N) + 1)] for j in 1:N]
    end
    return g2
end

function countercyclic_matrix(name::NTuple,N)

    g2 = Matrix{typeof(name[1])}(undef, N, N)  
    for i in 1:N
        g2[i, :] = [name[(mod(j - i + 1, N) + 1)] for j in 1:N]
    end
    return g2
end

mutable struct lattice
    A::ITensor
    B::ITensor
    C::ITensor
    D::ITensor
    E::ITensor
    F::ITensor
    G::ITensor
    H::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice(N)

        name = fieldnames(lattice)
        g1 = cyclic_matrix(name, N)
        g2 = countercyclic_matrix(name,N)

        new(
        ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),g1,g2)
    end
    
    # function lattice(a::ITensor,b::ITensor,c::ITensor,d::ITensor,e::ITensor,f::ITensor)

    #     name = fieldnames(lattice)
    #     g1 = [
    #     name[1] name[2] name[3] name[4] name[5] name[6]; 
    #     name[6] name[1] name[2] name[3] name[4] name[5]; 
    #     name[5] name[6] name[1] name[2] name[3] name[4]; 
    #     name[4] name[5] name[6] name[1] name[2] name[3]; 
    #     name[3] name[4] name[5] name[6] name[1] name[2]; 
    #     name[2] name[3] name[4] name[5] name[6] name[1]];

    #     g2 = [
    #     name[1] name[6] name[5] name[4] name[3] name[2]; 
    #     name[2] name[1] name[6] name[5] name[4] name[3]; 
    #     name[3] name[2] name[1] name[6] name[5] name[4]; 
    #     name[4] name[3] name[2] name[1] name[6] name[5]; 
    #     name[5] name[4] name[3] name[2] name[1] name[6]; 
    #     name[6] name[5] name[4] name[3] name[2] name[1]]; 

    #     new(a,b,c,d,e,f,g1,g2)
            
    # end

end



mutable struct lattice_ind
    A::Index
    B::Index
    C::Index
    D::Index
    E::Index
    F::Index
    G::Index
    H::Index
    # gg::Matrix{Symbol}
    # gt::Matrix{Symbol}

    function lattice_ind()
        
        # name = fieldnames(lattice_ind)
        # g1 = [
        # name[1] name[2] name[3] name[4] name[5] name[6]; 
        # name[6] name[1] name[2] name[3] name[4] name[5]; 
        # name[5] name[6] name[1] name[2] name[3] name[4]; 
        # name[4] name[5] name[6] name[1] name[2] name[3]; 
        # name[3] name[4] name[5] name[6] name[1] name[2]; 
        # name[2] name[3] name[4] name[5] name[6] name[1]];

        # g2 = [
        # name[1] name[6] name[5] name[4] name[3] name[2]; 
        # name[2] name[1] name[6] name[5] name[4] name[3]; 
        # name[3] name[2] name[1] name[6] name[5] name[4]; 
        # name[4] name[3] name[2] name[1] name[6] name[5]; 
        # name[5] name[4] name[3] name[2] name[1] name[6]; 
        # name[6] name[5] name[4] name[3] name[2] name[1]]; 

        new(ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor())
    end

    # function lattice_ind(a::Index,b::Index,c::Index,d::Index,e::Index,f::Index)

    #     name = fieldnames(lattice)
    #     g1 = [
    #     name[1] name[2] name[3] name[4] name[5] name[6]; 
    #     name[6] name[1] name[2] name[3] name[4] name[5]; 
    #     name[5] name[6] name[1] name[2] name[3] name[4]; 
    #     name[4] name[5] name[6] name[1] name[2] name[3]; 
    #     name[3] name[4] name[5] name[6] name[1] name[2]; 
    #     name[2] name[3] name[4] name[5] name[6] name[1]];

    #     g2 = [
    #     name[1] name[6] name[5] name[4] name[3] name[2]; 
    #     name[2] name[1] name[6] name[5] name[4] name[3]; 
    #     name[3] name[2] name[1] name[6] name[5] name[4]; 
    #     name[4] name[3] name[2] name[1] name[6] name[5]; 
    #     name[5] name[4] name[3] name[2] name[1] name[6]; 
    #     name[6] name[5] name[4] name[3] name[2] name[1]]; 
          
    #     new(a,b,c,d,e,f,g1,g2)

    # end

    function lattice_ind(nature_of_the_legs::String)

        if nature_of_the_legs == "physical"
            new(
            Index([QN(0)=>4],"ia"),Index([QN(0)=>4],"ib"),
            Index([QN(0)=>4],"ic"),Index([QN(0)=>4],"id"),
            Index([QN(0)=>4],"ie"),Index([QN(0)=>4],"if"),
            Index([QN(0)=>4],"ig"),Index([QN(0)=>4],"ih"))
        elseif nature_of_the_legs == "ancilla"
            new(
            Index([QN(0)=>4],"sa"),Index([QN(0)=>4],"sb"),
            Index([QN(0)=>4],"sc"),Index([QN(0)=>4],"sd"),
            Index([QN(0)=>4],"se"),Index([QN(0)=>4],"sf"),
            Index([QN(0)=>4],"ig"),Index([QN(0)=>4],"ih"))
        elseif nature_of_the_legs == "imaginary"
            new(
            Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]))
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
include("get_tens.jl")
include("ctmrg_evolution/ctm.jl")
include("simpleupdate.jl")


BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.enable_threaded_blocksparse()


function SU(parameters)

    dbetasu = 1e-2;
    J1 = parameters["J1"]
    J2 = parameters["J2"]
    Delta = parameters["Delta"]
    N = parameters["N"]
    hx = parameters["hx"]
    hz = parameters["hz"]
    D = parameters["D"]
    Gamma, lambdax, lambday, gt, gg, physical_legs, magnex, magnez, energie, it, err = simpleupdate(dbetasu,parameters)
    file_name_mat = @sprintf "Results/LocalTensors_N%.0f_J1%.2f_J2%.2f_Delta%.2f_D%.0f_hz%.2f_hx%.2f.mat" N J1 J2 Delta D hz hx

    file = matopen(file_name_mat, "w")
    write(file, "J1", J1)
    write(file, "J2", J2)
    write(file, "hz", hz)
    write(file, "hx", hx)
    write(file, "Delta", Delta)
    write(file, "ener", energie)
    write(file, "magnex", magnex)
    write(file, "magnez", magnez)
    write(file, "err", err)
    write(file, "it", it)
    write(file, "N", N)
    close(file)

end


let 
    
    N::Int64 = 6;
    D::Int64 = 3; 
    J1::Float64 = 1.0;
    J2::Float64 = 2.4;
    Delta::Float64 = 1;
    hx::Float64 = 0;
    hz::Float64= 0.;
    
    parameters = Dict("D"=>D, "J1"=>J1, "J2"=>J2, "hz"=>hz, "hx"=>hx, "Delta"=>Delta, "N"=>N)
    SU(parameters)

    nothing
  
end



println("END SIMPlE UPDATE")



