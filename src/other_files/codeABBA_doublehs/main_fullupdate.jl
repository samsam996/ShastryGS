

using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors

N::Int64 = 2;

# include("sixxsix.jl")
# using .unit_cell


mutable struct lattice2x1
    A::ITensor
    B::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice2x1()
        name = fieldnames(lattice2x1)
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]

        new(
        ITensor(),ITensor(),g1,g2)
    end
    
    function lattice2x1(a::ITensor,b::ITensor)

        name = fieldnames(lattice2x1)
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]


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
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]


        new(ITensor(),ITensor(),g1,g2)
    end

    function lattice_ind2x1(a::Index,b::Index)

        name = fieldnames(lattice2x1)
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]

          
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


include("SSMHamiltonian.jl")
include("simple_update_evolution/simpleUpdateJ1.jl")
include("observables/entropy_lambda.jl")
include("observables/correlation_length.jl")
# include("fullupdate_evolution/full_update.jl")
include("ctmrg_evolution/ctm.jl")
include("fullupdate_evolution/renormalisation_tens_A.jl")
include("get_tens.jl")
include("observables/magnetisation.jl")
include("observables/energy.jl")
include("observables/z2order_parameter.jl")

include("initialisation.jl")

include("simpleupdate.jl")
include("fullupdate.jl")
include("fastfullupdate.jl")
# include("eefullupdate.jl")

println("BEGIN FULL UPDATE")

BLAS.set_num_threads(18)
ITensors.Strided.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

let 
    
    D::Int64 = 4; 
    J1::Float64 = 0.63;
    J2::Float64 = 1.;
    h::Float64 = 1.; 
    hs::Float64 = 0; 


    temperature = LinRange(10.1,0.1,2)
    temperature = LinRange(10.1,0.10,2)

    temp_su = temperature[end]-1e-4
    dbetasu = 1e-3;
    Gamma, lambdax, lambday, gt, gg, physical_legs, ancilla_legs, magne, energie, 
    tempe, beta, free_energy, tmp = simpleupdate(temp_su,temperature,dbetasu,J1,J2,h,hs,N,D)
    temp = 1/beta

    # filename = "LocalTensorsD4temp0.08275J0.63h100.0.jld2"
    # Gamma, lambdax, lambday, gt, gg, physical_legs, ancilla_legs, temp, D, J1, J2, h, hs,nsu = 
    # load(filename,"Gamma","lambdax","lambday","gt","gg","physical_legs",
    # "ancilla_legs","temp","D","J1","J2","h","hs","nsu")
    # beta = 1/temp
    # println("NSU, hs : ",nsu, " ", hs)

    # filename = "U1_ResultsFUD4temp0.05627_nfu_200.0hs0.0h100.0.jld2"    
    # tens_A,cxd,cyd,chi = load(filename,"tens_A","cxd","cyd","chi")
    # physical_legs,ancilla_legs,gt = load(filename,"physical_legs","ancilla_legs","gt")
    # D,J1,J2,h,hs,nfu = load(filename,"D","J1","J2","h","hs","nfu")
    # gx,gy,mu,mu_hs = SSMHamiltonian(gt,physical_legs,nfu,J1,J2,h,hs)
    # tens_a = lattice()
    # f(x) = mod(x-1,6) + 1
    # for i = 1:6
    #     for j = 1
    #         ia = getproperty(physical_legs,gt[i,j])
    #         sa = getproperty(ancilla_legs,gt[i,j])
    #         cx1= getproperty(cxd,gt[f(i),f(j)])
    #         cx2 = getproperty(cxd,gt[f(i-1),f(j)])
    #         cy1= getproperty(cxd,gt[f(i),f(j)])
    #         cy2 = getproperty(cxd,gt[f(i),f(j-1)])
    #         A = getproperty(tens_A,gt[i,j])
    #         A_prime = prime(A,noncommoninds(inds(A),[ia,sa]))
    #         AA = (A_prime*dag(A))*cx1*dag(cx2)*cy1*dag(cy2)
    #     end
    # end
    # x = 6 
    # x[4] = 2
    # C,T = load(filename,"C","T")
    # beta = load(filname,"beta")
    ## FULL UPDATE PART 


    @show free_energy

    
    #### BLOCK 1
    dbetafu = 2e-2;
    nfu = 2/dbetafu;
    gx,gy,mu,mu_hs = SSMHamiltonian(gt,physical_legs,nfu,J1,J2,h,hs)
    tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)
   
    ener = [];
    tempera = [];
    order_parameter = [];
    tempprev = 0; it = 0;
    # temperature = LinRange(0.07, 0.055, 100)
    temperature = LinRange(9.,5.,5)

    relevant = 1; 
    chi = D*D + 1; 
    precision_ctm = 1e-9;
    @show num_iter_fu = (1/temperature[end] - beta)/dbetafu;
    C,T = OBC_PEPS(tens_A,cxd,cyd,gt)

    
    #### END OF BLOCK 1 

    C,T,tens_a, tens_A, cxd, cyd, beta, temp = fastfullupdate(temperature,J1,J2,h,hs,dbetafu,C,T,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs,gt,gg,chi,precision_ctm,beta,D,free_energy)

    # C,T,tens_a, tens_A, cxd, cyd, beta, temp = fullupdate(temperature,J1,J2,h,hs,dbetafu,C,T,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs,gt,gg,chi,precision_ctm,beta,D)
    # C,T,tens_a, tens_A, cxd, cyd, beta, temp = eefullupdate(temperature,J1,J2,h,hs,dbetafu,C,T,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs,gt,gg,chi,precision_ctm,beta,D)

    nothing

end