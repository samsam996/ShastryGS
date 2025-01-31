# import Pkg; 

# Pkg.add("ITensors")
# # Pkg.add("Plots")
# Pkg.add("FileIO")
# Pkg.add("JLD2")
# Pkg.add("LinearAlgebra")
# Pkg.add("KrylovKit")

using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors

include("SSMHamiltonian.jl")

include("simple_update_evolution/simpleUpdateJ1.jl")
include("observables/entropy_lambda.jl")
include("ctmrg_evolution/ctm.jl")
include("declare_indices.jl")
include("observables/correlation_length.jl")

println("BEGIN SIMPlE UPDATE")

BLAS.set_num_threads(20)
ITensors.Strided.set_num_threads(1)


let 
    
    mutable struct lattice
        A::ITensor
        B::ITensor
        C::ITensor
        D::ITensor
        function lattice()
            new(ITensor(),ITensor(),ITensor(),ITensor())
        end
        function lattice(a::ITensor,b::ITensor,c::ITensor,d::ITensor)
            new(a,b,c,d)
        end
    end

    mutable struct lattice_ind
        A::Index
        B::Index
        C::Index
        D::Index
    end

    d = 4

    nom_fichier = ["6","5786","5571","5357","5143","4929","4714","45","4286","4071","3857","3643","3429","3214"]

    for jj in eachindex(nom_fichier)

        file_name = "9771372D13h185/LocalTensorsD13temp0.0"*nom_fichier[jj]*"J0.63h185.0.jld2"

        # D, temp, J1, J2, nsu, h, hs, Gamma, lambdax, lambday, physical_legs, ancilla_legs, gt, gg = load(file_name,
        # "D", "temp", "J1", "J2", "nsu", "h", "hs", "Gamma", "lambdax", "lambday", "physical_legs", "ancilla_legs", "gt", "gg" )


        filen_name = string("Results/U1_ResultsD",string(D),"temp",string(round(temp,digits = 5)),"_nsu_",string(nsu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
        D,magne,J2,h,J1,hs,nsu,energie,tempe,
        invcorr1,invcorr2,invcorr3,invcorr4,wave_vec1,wave_vec2,wave_vec3,wave_vec4,
        eig3,eig4,numb_iter,chi,C,T,precision = load(file_name,"D",D,"magne",magne,"J2",J2,"h",h,"J1",J1,"hs",hs,"nsu",nsu,"energie",energie,"tempe",tempe,
        "invcorr1",invcorr1,"invcorr2",invcorr2,"invcorr3",invcorr3,"invcorr4",invcorr4,"wave_vec1",wave_vec1,
        "wave_vec2",wave_vec2,"wave_vec3",wave_vec3,"wave_vec4",wave_vec4,
        "eig3",eig3,"eig4",eig4,"numb_iter",numb_iter,"chi",chi,"C",C,"T",T,"precision",precision)
            

        J1 = Float64(J1)
        J2 = Float64(J2)
        check1 = "loadsuccesfully"
        println("load success")
        

        check2 = "Done ctm "
        save("Results/donectm.jld2","check2",check2)

        @show magne = 1/4*abs(mm[1,1]-mm[2,1]+mm[2,2]-mm[1,2])
        @show energie = ener
        tempe = temp
        
      
        invcorr1 = xi2[1,1]; invcorr2 = xi2[1,2]
        invcorr3 = xi2[2,1]; invcorr4 = xi2[2,2]
        wave_vec1 = dq[1,1]; wave_vec2 = dq[1,2]
        wave_vec3 = dq[2,1]; wave_vec4 = dq[2,2]
        eig3 = xi3[1,1]; eig4 = xi4[1,1]
        numb_iter = it;

        name_data = string("Results/U1_ResultsD",string(D),"temp",string(round(temp,digits = 5)),"_nsu_",string(nsu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
        save(name_data,"D",D,"magne",magne,"J2",J2,"h",h,"J1",J1,"hs",hs,"nsu",nsu,"energie",energie,"tempe",tempe,
        "invcorr1",invcorr1,"invcorr2",invcorr2,"invcorr3",invcorr3,"invcorr4",invcorr4,"wave_vec1",wave_vec1,
        "wave_vec2",wave_vec2,"wave_vec3",wave_vec3,"wave_vec4",wave_vec4,
        "eig3",eig3,"eig4",eig4,"numb_iter",numb_iter,"chi",chi,"C",C,"T",T,"precision",precision)
            

    end

end



println("END SIMPlE UPDATE")



