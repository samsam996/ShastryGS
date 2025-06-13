

include("initialisation.jl")

function simpleupdate(dbetasu::Float64,parameters::Dict,modit::Int64)


    Gamma, lambdax, lambday, physical_legs, gt, gg = initialisation(parameters)

    nsu = 1/dbetasu;
    J1 = parameters["J1"]
    J2 = parameters["J2"]
    D = parameters["D"]
    hz = parameters["hz"]
    hx = parameters["hx"]
    Delta = parameters["Delta"]
    N = parameters["N"]
    model = parameters["model"]

    load_moves(N)

    f(x) = mod(x-1,N) + 1
    magnex = []
    magney = []
    magnez = []
    energie = []
    numb_iter = []
    Ps = []
    xi = []
    free_energy = 0
    relevant = 1
    ener = 9
    enertmp = 10
    iter = 0
    itctm = 0
    errctm = 0

    while ener < enertmp && abs(ener - enertmp) > 1e-6

        iter = iter + 1
        Gamma,lambdax,lambday,free_energy = simpleupdateJ1(Gamma,lambdax,lambday,physical_legs,gt,nsu,parameters,free_energy)

        if mod(iter,10)==0
            println(iter)
        end

        if mod(iter,modit) == 0
                
            if model=="XY"

                file_name_jld2 = @sprintf "Results/LocalTensors_N%.0f_J1%.2f_J2%.2f_Delta%.2f_D%.0f_hz%.2f_hx%.2f.jld2" N J1 J2 Delta D hz hx
                save(file_name_jld2,
                        "D",D,
                        "Delta",delta,
                        "J1",J1,
                        "J2",J2,
                        "nsu",nsu,
                        "hz",hz,
                        "hx",hx,
                        "Gamma",Gamma,
                        "lambdax",lambdax,
                        "lambday",lambday,
                        "physical_legs",physical_legs,
                        "gt",gt,
                        "N",N,
                    "gg",gg)

            elseif model=="XYZ"

                Jx = parameters["Jx"]
                Jy = parameters["Jy"]
                Jz = parameters["Jz"]

                file_name_jld2 = @sprintf "Results/LocalTensors_N%.0f_J1%.2f_J2%.2f_Jx%.2f_Jy%.2f_Jz%.2f_D%.0f_hz%.2f_hx%.2f.jld2" N J1 J2 Jx Jy Jz D hz hx
                save(file_name_jld2,
                        "D",D,
                        "Delta",delta,
                        "J1",J1,
                        "J2",J2,
                        "nsu",nsu,
                        "hz",hz,
                        "hx",hx,
                        "Gamma",Gamma,
                        "lambdax",lambdax,
                        "lambday",lambday,
                        "physical_legs",physical_legs,
                        "gt",gt,
                        "N",N,
                    "gg",gg)

            else
                println("Please choose a valid model, either XY or XYZ")
                exit()
            end
            
            chi = D*D + 1; 
            precision_ctm = 1e-5
            tens_a,tens_A,cxd,cyd = get_tens(Gamma,lambdax,lambday,physical_legs,gt,N)
            
        
            C,T = OBC_PEPS(tens_A,cxd,cyd,gt,N)
            C, T, itctm, errctm = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,chi,precision_ctm,C,T,N) 

            enertmp = ener
            ener, mmx, mmy, mmz, Ps = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,parameters)
            ener = real(ener)
            # xi2, xi3, xi4, dq = correlation_length(C,T,gt,N)
            xi2 = zeros(2,2)
            xi3 = zeros(2,2)
            xi4 = zeros(2,2)
            dq = 0
            xi = [xi2[1,1], xi3[1,1], xi4[1,1], dq]
            sum_magnex = sum(mmx)
            sum_magney = sum(mmy)
            sum_magnez = sum(mmz)

            @show push!(magnex, sum_magnex)
            @show push!(magney, sum_magney)
            @show push!(magnez, sum_magnez)
            @show push!(energie, real(ener))
        
        end

    end

    return Gamma, lambdax, lambday, gt, gg, physical_legs, magnex, magney, magnez, energie, Ps, itctm, errctm, xi

end


function SU(parameters)

    """
        Performs the imaginary time evolution until either convergence of the energy, or when the energy increases again 

        ARGS: 
            J1,J2, hx, hz..: different coupling constant
            D : bond dimension 
            h : magentic field in the z direction 
            symmetry : either "" (no symmetry) or U1 
            dbeta : time-step of the evolution 
            chi : bond dimension of the corner transfer matrix algorithm 
            modit : will compute the energy at every modit time step

    """

    dbetasu = parameters["dbeta"]
    J1 = parameters["J1"]
    J2 = parameters["J2"]
    N = parameters["N"]
    hx = parameters["hx"]
    hz = parameters["hz"]
    D = parameters["D"]
    model = parameters["model"]
    modit = parameters["modit"]
    N = parameters["N"]
    Jx = parameters["Jx"]
    Jy = parameters["Jy"]
    Jz = parameters["Jz"]
    Delta = parameters["Delta"]

    
    Gamma, lambdax, lambday, gt, gg, physical_legs, magnex, magney, magnez, energie, Ps, it, err, xi = simpleupdate(dbetasu,parameters,modit)

    if model == "XY"
        Delta = parameters["Delta"]
        file_name_mat = @sprintf "Results/LocalTensors_N%.0f_J1%.2f_J2%.2f_Delta%.2f_D%.0f_hz%.2f_hx%.2f.mat" N J1 J2 Delta D hz hx
        file = matopen(file_name_mat, "w")
        write(file, "J1", J1)
        write(file, "J2", J2)
        write(file, "hz", hz)
        write(file, "hx", hx)
        write(file, "Delta", Delta)
        write(file, "ener", energie)
        write(file, "magnex", magnex)
        write(file, "magney", magney)
        write(file, "magnez", magnez)
        write(file, "err", err)
        write(file, "it", it)
        write(file, "N", N)
        write(file, "xi", xi)
        close(file)

    elseif model == "XYZ"
        Jx = parameters["Jx"]
        Jy = parameters["Jy"]
        Jz = parameters["Jz"]
        file_name_mat = @sprintf "Results/LocalTensors_N%.0f_J1%.2f_J2%.2f_Jx%.2f_Jy%.2f_Jz%.2f_D%.0f_hz%.2f_hx%.2f.mat" N J1 J2 Jx Jy Jz D hz hx
        file = matopen(file_name_mat, "w")
        write(file, "J1", J1)
        write(file, "J2", J2)
        write(file, "hz", hz)
        write(file, "hx", hx)
        write(file, "ener", energie)
        write(file, "magnex", magnex)
        write(file, "magney", magney)
        write(file, "magnez", magnez)
        write(file, "err", err)
        write(file, "it", it)
        write(file, "N", N)
        write(file, "Jx", Jx)
        write(file, "Jy", Jy)
        write(file, "Jz", Jz)
        write(file, "Ps", Ps)
        write(file, "xi", xi)
        close(file)
    end
   

end