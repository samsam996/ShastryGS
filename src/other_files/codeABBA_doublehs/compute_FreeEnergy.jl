


using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors


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


include("get_tens.jl")
include("SSMHamiltonian.jl")
include("simple_update_evolution/simpleUpdateJ1.jl")
include("observables/entropy_lambda.jl")
include("observables/correlation_length.jl")
include("observables/magnetisation.jl")
include("observables/energy.jl")
include("observables/z2order_parameter.jl")
include("observables/PartitionPerSite.jl")
include("ctmrg_evolution/ctm.jl")
include("initialisation.jl")
include("simpleupdate.jl")

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)
ITensors.disable_threaded_blocksparse()


let

# file_name = "U1_ResultsD4temp0.0542_nsu_200.0hs0.0h100.0.jld2"
file_name =[] 
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
PPS = []
ener = []



name = "Results"
tt = ["4.0","4.31034","4.65116","5.29101","5.95238","6.25","6.57895","6.89655","7.24638","7.57576","7.87402"]

for i in eachindex(tt)

    @show i
    file_name = name*"/LocalTensorsD4temp"*tt[i]*"J0.63h100.0.jld2"
    D,J1,J2,h,hs,nsu,temp = load(file_name,"D","J1","J2","h","hs","nsu","temp")

    Gamma, lambdax, lambday, gt, FreeEnergy, physical_legs, ancilla_legs = load(file_name,"Gamma","lambdax","lambday","gt","FreeEnergy", "physical_legs","ancilla_legs")


    chi = D*D + 1; precision_ctm = 1e-7
    tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)
    C,T = OBC_PEPS(tens_A,cxd,cyd,gt)

    # can I load C and T right away ? 
    bc = "one wala"
    C, T, it, err = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,chi,precision_ctm,C,T,bc) 
    mm = magnetisation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs)
    ee = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,J1,J2,h,hs)

    PPSf = 1/2*(PartitionPerSite(T,C,tens_a,gt,1,1) ) + FreeEnergy
    omega0 = 1;
    omega1 = exp(2*pi*im/3)
    omega2 = exp(4*pi*im/3);

    order_parameter = 1/2*abs(omega0*mm[1,1] - mm[1,2])

    z2 = z2order_parameter(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs)
    @show push!(ener, real(ee))
    @show push!(order_parameter_z2,z2)
    @show push!(magne,order_parameter)
    @show push!(tempe, temp)
    @show push!(PPS,PPSf)

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
    "PPS",PPSf,
    "gt",gt,
"ener",ener)


end



# name_data = string("Results/U1_ResultsD",string(D),"temp",string(round(temp,digits = 5)),"_nsu_",string(nsu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
#                 save(name_data,
#                 "D",D,
#                 "magne",order_parameter,
#                 "J2",J2,
#                 "h",h,
#                 "J1",J1,
#                 "hs",hs,
#                 "nsu",nsu,
#                 "energie",ener,
#                 "z2",z2,
#                 "tempe",temp,
#                 "invcorr1",xi2[1,1],
#                 "invcorr2",xi2[1,2],
#                 "invcorr3",xi2[2,1],
#                 "invcorr4",xi2[2,2],
#                 "wave_vec1",dq[1,1],
#                 "wave_vec2",dq[1,2],
#                 "wave_vec3",dq[2,1],
#                 "wave_vec4",dq[2,2],
#                 "magne_z", mm,
#                 "eig3",xi3[1,1],
#                 "eig4",xi4[1,1],
#                 "numb_iter",it,
#                 "chi",chi,
#                 "err_ctm",err,
#                 "C",C,
#                 "T",T,
#                 "tens_A",tens_A,
#                 "cxd",cxd,
#                 "cyd",cyd,
#  "gt",gt)
        




 


end