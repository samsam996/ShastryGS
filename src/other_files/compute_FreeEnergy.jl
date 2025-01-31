


using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors


mutable struct lattice2x1
    A::ITensor
    B::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice()
        name = fieldnames(lattice)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]
        new(ITensor(),ITensor(),g1,g2)
    end
    
    function lattice(a::ITensor,b::ITensor)

        name = fieldnames(lattice)
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

    function lattice_ind()
        
        name = fieldnames(lattice_ind)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

        new(ITensor(),ITensor(),g1,g2)
    end

    function lattice_ind(a::Index,b::Index)

        name = fieldnames(lattice)
        g1 = [name[1] name[2]; name[2] name[1]]
        g2 = [name[1] name[2]; name[2] name[1]]

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


include("get_tens2.jl")
include("SSMHamiltonian.jl")
include("simple_update_evolution/simpleUpdateJ1.jl")
include("observables/entropy_lambda.jl")
include("observables/correlation_length.jl")
include("observables/magnetisation.jl")
include("observables/energy_evol_2x2.jl")
include("observables/energy.jl")
include("observables/z2order_parameter.jl")
include("observables/PartitionPerSite.jl")
include("ctmrg_evolution/ctm.jl")
include("initialisation.jl")
include("simpleupdate.jl")

BLAS.set_num_threads(20)
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



name1 = "21045833_scanD7"


tt = ["0.055","0.065","0.0591","0.0631","0.05525","0.05551",
    "0.05577","0.05603","0.05628","0.05654",
    "0.05679","0.05705","0.05731","0.05756","0.05782",
    "0.05808","0.05833","0.05859","0.05884","0.05936",
    "0.05961","0.05987","0.06013","0.06038","0.06064",
    "0.06089","0.06115","0.06141","0.06166","0.06192",
    "0.06218","0.06243","0.06269","0.06295","0.06346",
    "0.06371","0.06397","0.06423","0.06449","0.06474"];




for i in eachindex(tt)

    @show i
    file_name = name1*"/LocalTensorsD7temp"*tt[i]*"J0.63h100.0.jld2"
    D,J1,J2,h,hs,nsu,temp = load(file_name,"D","J1","J2","h","hs","nsu","temp")

    Gamma_old, lambdax_old, lambday_old, gt_old, FreeEnergy, physical_legs_old, ancilla_legs_old = load(file_name,"Gamma","lambdax","lambday","gt","FreeEnergy", "physical_legs","ancilla_legs")
    gt = deepcopy(gt_old)

    Gamma = lattice(getproperty(Gamma_old,gt_old[1,1]),getproperty(Gamma_old,gt_old[1,2]),getproperty(Gamma_old,gt_old[1,1]),
    getproperty(Gamma_old,gt_old[1,2]),getproperty(Gamma_old,gt_old[1,1]),getproperty(Gamma_old,gt_old[1,2]))

    lambdax = lattice(getproperty(lambdax_old,gt_old[1,1]),getproperty(lambdax_old,gt_old[1,2]),getproperty(lambdax_old,gt_old[1,1]),
    getproperty(lambdax_old,gt_old[1,2]),getproperty(lambdax_old,gt_old[1,1]),getproperty(lambdax_old,gt_old[1,2]))

    lambday = lattice(getproperty(lambday_old,gt_old[1,1]),getproperty(lambday_old,gt_old[1,2]),getproperty(lambday_old,gt_old[1,1]),
    getproperty(lambday_old,gt_old[1,2]),getproperty(lambday_old,gt_old[1,1]),getproperty(lambday_old,gt_old[1,2]))

    physical_legs = lattice_ind(getproperty(physical_legs_old,gt[1,1]),getproperty(physical_legs_old,gt[1,2]),
    getproperty(physical_legs_old,gt_old[1,1]),getproperty(physical_legs_old,gt[1,2]),getproperty(physical_legs_old,gt[1,1]),
    getproperty(physical_legs_old,gt[1,2]))

    ancilla_legs = lattice_ind(getproperty(ancilla_legs_old,gt[1,1]),getproperty(ancilla_legs_old,gt[1,2]),
    getproperty(ancilla_legs_old,gt[1,1]),getproperty(ancilla_legs_old,gt[1,2]),getproperty(ancilla_legs_old,gt[1,1]),
    getproperty(ancilla_legs_old,gt[1,2]))


    gt = Gamma.gt
    chi = D*D + 1; 
    chi = 40
    precision_ctm = 1e-7
    tens_a,tens_A,cxd,cyd = get_tens2(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)
    C,T = OBC_PEPS_ONES(tens_A,cxd,cyd,gt,physical_legs)

    # can I load C and T right away ? 
    bc = "one wala"
    C, T, it, err = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,J1,J2,h,chi,precision_ctm,C,T,bc) 
    mm = magnetisation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs)

    PPSf = 1/6*(PartitionPerSite(T,C,tens_a,gt,1,1) + PartitionPerSite(T,C,tens_a,gt,3,1) + PartitionPerSite(T,C,tens_a,gt,5,1)) + 
    2*FreeEnergy

    ener = energy_evol_2x2(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,J1,J2,h,hs)
    # @show vertical = vertical_correlation_evol_2x2(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,1,1)
    # @show gg
    # @show horizontal = horizontal_correlation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,1,1)

    omega0 = 1;
    omega1 = exp(2*pi*im/3)
    omega2 = exp(4*pi*im/3);
    println([mm[1,1],mm[1,2],mm[1,3],mm[1,4],mm[1,5],mm[1,6]])

    order_parameter = 1/6*abs(omega0*mm[1,1]+omega1*mm[1,2]+omega2*mm[1,3]+omega0*mm[1,4]+omega1*mm[1,5]+omega2*mm[1,6])

    z2 = z2order_parameter(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs)
    @show push!(order_parameter_z2,z2)
    @show push!(magne,order_parameter)
    @show push!(tempe, temp)
    @show push!(energie,ener)
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
"gt",gt)

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
#                 "PPS",PPSf,
#  "gt",gt)
        




 


end