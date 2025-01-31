

using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors

N::Int64 = 6;

# include("sixxsix.jl")
# using .unit_cell


mutable struct lattice
    A::ITensor
    B::ITensor
    C::ITensor
    D::ITensor
    E::ITensor
    F::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice()
        name = fieldnames(lattice)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 

        new(
        ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),g1,g2)
    end
    
    function lattice(a::ITensor,b::ITensor,c::ITensor,d::ITensor,e::ITensor,f::ITensor)

        name = fieldnames(lattice)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 

        new(a,b,c,d,e,f,g1,g2)
            
    end

end

mutable struct lattice_ind
    A::Index
    B::Index
    C::Index
    D::Index
    E::Index
    F::Index
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}

    function lattice_ind()
        
        name = fieldnames(lattice_ind)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 

        new(ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),g1,g2)
    end

    function lattice_ind(a::Index,b::Index,c::Index,d::Index,e::Index,f::Index)

        name = fieldnames(lattice)
        g1 = [
        name[1] name[2] name[3] name[4] name[5] name[6]; 
        name[6] name[1] name[2] name[3] name[4] name[5]; 
        name[5] name[6] name[1] name[2] name[3] name[4]; 
        name[4] name[5] name[6] name[1] name[2] name[3]; 
        name[3] name[4] name[5] name[6] name[1] name[2]; 
        name[2] name[3] name[4] name[5] name[6] name[1]];

        g2 = [
        name[1] name[6] name[5] name[4] name[3] name[2]; 
        name[2] name[1] name[6] name[5] name[4] name[3]; 
        name[3] name[2] name[1] name[6] name[5] name[4]; 
        name[4] name[3] name[2] name[1] name[6] name[5]; 
        name[5] name[4] name[3] name[2] name[1] name[6]; 
        name[6] name[5] name[4] name[3] name[2] name[1]]; 
          
        new(a,b,c,d,e,f,g1,g2)

    end

    function lattice_ind(nature_of_the_legs::String)

        if nature_of_the_legs == "physical"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ia"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ib"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ic"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"id"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"ie"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"if"))
        elseif nature_of_the_legs == "ancilla"
            new(
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sa"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sb"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sc"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sd"),
            Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"se"),Index([QN(-2)=>1,QN(0)=>2,QN(2)=>1],"sf"))
        elseif nature_of_the_legs == "imaginary"
            new(
            Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]))
        end

    end

end


include("SSMHamiltonian.jl")
include("simple_update_evolution/simpleUpdateJ1.jl")
include("observables/entropy_lambda.jl")
include("observables/correlation_length.jl")
include("fullupdate_evolution/full_update.jl")
include("ctmrg_evolution/ctm.jl")
include("fullupdate_evolution/renormalisation_tens_A.jl")
include("get_tens.jl")
include("observables/magnetisation.jl")
include("observables/energy.jl")
include("observables/z2order_parameter.jl")
include("observables/PartitionPerSite.jl")

include("initialisation.jl")

include("simpleupdate.jl")
include("fullupdate.jl")
include("fastfullupdate.jl")
include("eefullupdate.jl")

println("BEGIN FULL UPDATE")

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

let 
    
    D::Int64 = 4; 
    J1::Float64 = 0.63;
    J2::Float64 = 1.;
    h::Float64 = 1.; 
    hs::Float64 = 0; 

    temperature = LinRange(0.07,0.064,20)
    # temperature = LinRange(10.12,0.09,2)
    temperature = LinRange(4.5, 4, 2)
    temperature = LinRange(10.1,9.12,2)

    temp_su = temperature[end]-1e-4
    dbetasu = 1e-3;
    Gamma, lambdax, lambday, gt, gg, physical_legs, ancilla_legs, magne, energie, 
    tempe, beta, free_energy = simpleupdate(temp_su,temperature,dbetasu,J1,J2,h,hs,N,D)
    temp = 1/beta

    # filename = "LocalTensorsD4temp0.08275J0.63h100.0.jld2"
    # Gamma, lambdax, lambday, gt, gg, physical_legs, ancilla_legs, temp, D, J1, J2, h, hs,nsu = 
    # load(filename,"Gamma","lambdax","lambday","gt","gg","physical_legs",
    # "ancilla_legs","temp","D","J1","J2","h","hs","nsu")
    # beta = 1/temp
    # println("NSU, hs : ",nsu, " ", hs)

    # filename = "U1_ResultsFUD6temp0.05854_nfu_200.0hs0.0h100.0.jld2"    

    # tens_A,cxd,cyd,chi = load(filename,"tens_A","cxd","cyd","chi")
    # physical_legs,ancilla_legs,gt = load(filename,"physical_legs","ancilla_legs","gt")
    # D,J1,J2,h,hs,nfu,free_energy = load(filename,"D","J1","J2","h","hs","nfu","PPS")
    # gx,gy,mu,mu_hs = load(filename,"gx","gy","mu","mu_hs")
    # tens_a = load(filename,"tens_a")
    # f(x) = mod(x-1,6) + 1
    # C,T = load(filename,"C","T")
    # temp = load(filename,"temp")
    # gg = load(filename, "gg")
    # precision_ctm = load(filename, "precision_ctm")

    # beta = 1/temp
    ## FULL UPDATE PART 



    #### BLOCK 1
    dbetafu = 1e-2;
    nfu = 2/dbetafu;
    gx,gy,mu,mu_hs = SSMHamiltonian(gt,physical_legs,nfu,J1,J2,h,hs)
    tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)
   
    ener = [];
    tempera = [];
    order_parameter = [];
    tempprev = 0; it = 0;
    # temperature = LinRange(9., 2, 10)
    temperature = LinRange(0.08,0.05,50)
    temperature = LinRange(9.,5.,5)

    relevant = 1; 
    chi = D*D + 1; 
    precision_ctm = 1e-9;
    @show num_iter_fu = (1/temperature[end] - beta)/dbetafu;
    # # C,T = OBC_PEPS(tens_A,cxd,cyd,gt)
    C,T = OBC_PEPS_ONES(tens_A,cxd,cyd,gt,physical_legs)
    
    #### END OF BLOCK 1 

    C,T,tens_a, tens_A, cxd, cyd, beta, temp = fastfullupdate(temperature,J1,J2,h,hs,dbetafu,C,T,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs,gt,gg,chi,precision_ctm,beta,D,free_energy)


    # C,T,tens_a, tens_A, cxd, cyd, beta, temp = fullupdate(temperature,J1,J2,h,hs,dbetafu,C,T,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs,gt,gg,chi,precision_ctm,beta,D)
    # C,T,tens_a, tens_A, cxd, cyd, beta, temp = eefullupdate(temperature,J1,J2,h,hs,dbetafu,C,T,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs,gt,gg,chi,precision_ctm,beta,D)

    nothing

end