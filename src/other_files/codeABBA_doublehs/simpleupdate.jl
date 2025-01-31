


include("observables/PartitionPerSite.jl")


function simpleupdate(final_temp::Float64,temperature,dbetasu::Float64,J1::Float64,J2::Float64,h::Float64,hs::Float64,N::Int64,D::Int64)

    f(x) = mod(x-1,N) + 1

    Gamma, lambdax, lambday, physical_legs, ancilla_legs, gt, gg = initialisation()

    beta = 0;
    temp =  Inf;

    nsu = 2/dbetasu;
   
    magne = []
    energie = []
    tempe = []
    discarded_weigth = []
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
    PPS = []
    temperature_string = []

    relevant = 1
    FreeEnergy = 0;
    DW = 0
    while final_temp < temp


        tempprev = 1/beta
        beta = beta + dbetasu; 
        temp = 1/beta;
        Gamma,lambdax,lambday, FreeEnergy, DW = simpleupdateJ1(Gamma,lambdax,lambday,physical_legs,gt,nsu,J1,J2,h,hs,D, FreeEnergy)
        # DW = DW + dw;
        # println(DW)
        println(temp)
        # @show FreeEnergy

        if size(temperature)[1] >= relevant && size(temperature)[1] > 0 
            if (temperature[relevant] <= tempprev && temperature[relevant] > temp) 
                
                name_data2 = string("Results/LocalTensorsD",string(D),"temp",string(round(temp,digits = 5)),"J",string(J1),"h",string(round(1e2*h,digits = 3)),".jld2")
                save(name_data2,
                    "D",D,
                    "temp",temp,
                    "J1",J1,
                    "J2",J2,
                    "nsu",nsu,
                    "h",h,
                    "hs",hs,
                    "Gamma",Gamma,
                    "lambdax",lambdax,
                    "lambday",lambday,
                    "FreeEnergy", FreeEnergy,
                    "physical_legs",physical_legs,
                    "ancilla_legs",ancilla_legs,
                    "gt",gt,
                    "DW",DW,
                "gg",gg)

                temp_string = string(round(temp,digits = 5));
                push!(temperature_string, temp_string)
                chi = D*D + 1; precision_ctm = 1e-7
                tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,ancilla_legs,gt)

                C,T = OBC_PEPS(tens_A,cxd,cyd,gt)

                bc = "open condition"
                C, T, it, err = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,chi,precision_ctm,C,T,bc) 
                mm = magnetisation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs)
                ener = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,J1,J2,h,hs)
                pps = 1/2*PartitionPerSite(T,C,tens_a,gt,1,1) + 2*FreeEnergy
            
                # log(κ) + 2* log(λ)
                # 

                omega0 = 1;
                omega1 = exp(2*pi*im/3)
                omega2 = exp(4*pi*im/3);
                order_parameter = 0 #1/6*abs(omega0*mm[1,1]+omega1*mm[1,2]+omega2*mm[1,3]+omega0*mm[1,4]+omega1*mm[1,5]+omega2*mm[1,6])
                
                @show push!(discarded_weigth, DW)
                @show push!(magne,order_parameter)
                @show push!(energie, real(ener))
                @show push!(tempe, temp)
                @show push!(PPS, real(pps))

                xi2, xi3, xi4, dq = correlation_length(C,T,gt)
                xi2 = zeros(2,2); xi3 = zeros(2,2); xi4 = zeros(2,2); dq = zeros(2,2); 
                
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

                relevant = relevant + 1

                # name_data = string("Results/U1_ResultsD",string(D),"temp",string(round(temp,digits = 5)),"_nsu_",string(nsu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
                #     save(name_data,
                #         "D",D,
                #         "magne",order_parameter,
                #         "J2",J2,
                #         "h",h,
                #         "J1",J1,
                #         "hs",hs,
                #         "nsu",nsu,
                #         "energie",ener,
                #         "tempe",temp,
                #         "invcorr1",xi2[1,1],
                #         "invcorr2",xi2[1,2],
                #         "invcorr3",xi2[2,1],
                #         "invcorr4",xi2[2,2],
                #         "wave_vec1",dq[1,1],
                #         "wave_vec2",dq[1,2],
                #         "wave_vec3",dq[2,1],
                #         "wave_vec4",dq[2,2],
                #         "magne_z", mm,
                #         "eig3",xi3[1,1],
                #         "eig4",xi4[1,1],
                #         "numb_iter",it,
                #         "chi",chi,
                #         "err_ctm",err,
                #         "C",C,
                #         "T",T,
                #         "tens_A",tens_A,
                #         "cxd",cxd,
                #         "cyd",cyd,
                #         "FreeEnergy",FreeEnergy,
                #         "physical_legs",physical_legs,
                #         "ancilla_legs",ancilla_legs,
                #         "tens_a",tens_a,
                #     "gt",gt)
            

            end 

        end



        save("Results/temperature.jld2","temp",temperature_string)

    end

    return Gamma, lambdax, lambday, gt, gg, physical_legs, ancilla_legs, magne, energie, tempe, beta, FreeEnergy, PPS  

    

end