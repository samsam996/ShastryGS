
using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors

mutable struct lattice
    A::ITensor
    B::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice()
        name = fieldnames(lattice)
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]

        new(
        ITensor(),ITensor(),g1,g2)
    end
    
    function lattice(a::ITensor,b::ITensor)

        name = fieldnames(lattice)
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]


        new(a,b,g1,g2)
            
    end

end

mutable struct lattice_ind
    A::Index
    B::Index
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}

    function lattice_ind()
        
        name = fieldnames(lattice_ind)
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]


        new(ITensor(),ITensor(),g1,g2)
    end

    function lattice_ind(a::Index,b::Index)

        name = fieldnames(lattice)
        g1 = [
        name[1] name[2];
        name[2] name[1]]

        g2 = [
            name[1] name[2];
            name[2] name[1]]

          
        new(a,b,g1,g2)

    end

    function lattice_ind(nature_of_the_legs::String)

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
include("observables/z2order_parameter.jl")
include("get_tens.jl")
include("ctmrg_evolution/ctm.jl")
include("initialisation.jl")

include("simpleupdate.jl")

# BLAS.set_num_threads(10)
# ITensors.Strided.set_num_threads(1)

BLAS.set_num_threads(18)
ITensors.Strided.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

let 
    

  
    
    # D::Int64 = 6; 
    # J1::Float64 = 0.63;
    # J2::Float64 = 1.;
    # h::Float64 = 1.; 

    # hs::Float64 = 5e-4; 

    magne = []
    energie = []
    tempe = []
    invcorr1 = []
    invcorr2 = []
    invcorr3 = []
    invcorr4 = []
    wave_vec1 = []
    wave_vec2 = []  
    wave_vec3 = []
    wave_vec4 = []
    eig3 = []
    eig4 = []
    numb_iter = []
    order_parameter_z2 = []

    # name = ["48","4958","5116","5273","5432","5589","5747","5905","6063","6221","6379","6537","6694","6853","701","7168","7326","7484","7642","78"]
    # name = ["59","62","591","619","5921","5931","5941","5952","5962","5972","5983","5993","6003","6014","6014","6024","6034","6044","6055","6065","6076","6086","6096","6107","6117","6127","6138","6148"]
    # name = ["62","64","621","6221","6231","6242","6252","6263","6274","6284","6294","6305","6316","6326","6337","6347","6358","6368","6579","6389"]
    # name = ["55","5579","5658","5737","5816","5894","5973","6053","6131","621","6289","6368","6447","6526","6605","6684","6763","6842","6921","7"]
    name = ["6063","6126","6189","6252","6316","6379","6442","6505","6568","6631","6694","6758","6821","6884","6947","6","701","7074","7137","72"]
    # name = ["6352","63","6405","6458","651","6563","6616","6668","6721","6773","6826","6879","6931","6984","7037","7089","7142","7194","7247","73"]

    # name = ["61","63","611","6121","6131","6142","6152","6163","6174","6184","6195","6205","6215","6226","6237","6247","6258","6268","6279","6289"]
    # name = ["63","65","631","6321","6332","6342","6352","6363","6373","6384","6395","6405","6416","6426","6437","6447","6458","6468","6479","6489"]
    
    name2 = "12424381_D6dhs2e-3"
    name3 = "12424387_D6dhs3e-3"
    name4 = "12424395_D6dhs4e-3"
    name5 = "12424403_D6dhs5e-3"

    simu_name = name2

    name = load(simu_name*"/name.jld2", "name")

    for i = 1:3

 
        # filename = "10778140_D10hs1e-2/LocalTensorsD10temp0.0"*name[i]*"J0.63h100.0.jld2"
        filename = simu_name*"/LocalTensorsD6temp0.0"*name[i]*"J0.63h100.0.jld2"

        D,J1,J2,h,hs,nsu = load(filename,"D","J1","J2","h","hs","nsu")
        Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt = load(filename,"Gamma","lambdax","lambday","physical_legs","ancilla_legs","gt")
        temp = load(filename, "temp")

        chi = D*D + 1; precision_ctm = 1e-7
        tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)
                
        C,T = OBC_PEPS(tens_A,cxd,cyd,gt)

        bc = "open condition"
        C, T, it, err = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,chi,precision_ctm,C,T,bc) 
        mm = magnetisation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs)
        # ener = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,J1,J2,h,hs)

        omega0 = 1;
        omega1 = exp(2*pi*im/3)
        omega2 = exp(4*pi*im/3);
        println([mm[1,1],mm[1,2],mm[1,3],mm[1,4],mm[1,5],mm[1,6]])
     
        order_parameter = 1/6*abs(omega0*mm[1,1]+omega1*mm[1,2]+omega2*mm[1,3]+omega0*mm[1,4]+omega1*mm[1,5]+omega2*mm[1,6])
        
        z2 = z2order_parameter(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs)
        @show push!(order_parameter_z2,z2)
        @show push!(magne,order_parameter)
        @show push!(energie, ener)
        @show push!(tempe, temp)

        xi2, xi3, xi4, dq = correlation_length(C,T,gt)
        # xi2 = zeros(2,2); xi3 = zeros(2,2); xi4 = zeros(2,2); dq = zeros(2,2); 
                    
        push!(invcorr1, xi2[1,1])
        push!(invcorr2, xi2[1,2])
        push!(invcorr3, xi2[2,1])
        push!(invcorr4, xi2[2,2])
        push!(wave_vec1, dq[1,1])
        push!(wave_vec2, dq[1,2])
        push!(wave_vec3, dq[2,1])
        push!(wave_vec4, dq[2,2])
        push!(eig3, xi3[1,1])
        push!(eig4, xi4[1,1])
        push!(numb_iter, it)

        name_data = string("Results/U1_ResultsD",string(D),"temp",string(round(temp,digits = 5)),"_nsu_",string(nsu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
                        save(name_data,
                        "D",D,
                        "magne",order_parameter,
                        "J2",J2,
                        "h",h,
                        "J1",J1,
                        "hs",hs,
                        "nsu",nsu,
                        "energie",ener,
                        "z2",z2,
                        "tempe",temp,
                        "invcorr1",xi2[1,1],
                        "invcorr2",xi2[1,2],
                        "invcorr3",xi2[2,1],
                        "invcorr4",xi2[2,2],
                        "wave_vec1",dq[1,1],
                        "wave_vec2",dq[1,2],
                        "wave_vec3",dq[2,1],
                        "wave_vec4",dq[2,2],
                        "magne_z", mm,
                        "eig3",xi3[1,1],
                        "eig4",xi4[1,1],
                        "numb_iter",it,
                        "chi",chi,
                        "err_ctm",err,
                        "C",C,
                        "T",T,
                        "tens_A",tens_A,
                        "cxd",cxd,
                        "cyd",cyd,
                        "gt",gt)
                

    end

        @show order_parameter_z2
        @show tempe
  
end



println("END SIMPlE UPDATE")



