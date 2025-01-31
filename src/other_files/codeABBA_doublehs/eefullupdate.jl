


include("fullupdate_evolution/full_update_x.jl")
include("fullupdate_evolution/full_update_y.jl")



function eefullupdate(temperature,J1::Float64,J2::Float64,h::Float64,hs::Float64,dbetafu::Float64,C::Vector{lattice},
    T::Vector{lattice},tens_a::lattice,tens_A::lattice,cxd::lattice,cyd::lattice,
    physical_legs::lattice_ind,ancilla_legs::lattice_ind,gt::Matrix{Symbol},gg::Matrix{Symbol},
    chi::Int64,precision_ctm::Float64,beta::Float64,D::Int64)


    nfu = 2/dbetafu;

    gx,gy,mu,mu_hs = SSMHamiltonian(gt,physical_legs,nfu,J1,J2,h,hs)

    ener = [];
    tempera = [];
    order_parameter = [];
    invcorr1 = [];
    eig3 = []
    eig4 = []
    order_parameter_z2 = []

    tempprev = 0; it = 0;    
    relevant = 1; 
    temp = 1/beta;
    
    @show num_iter_fu = (1/temperature[end] - beta)/dbetafu;

    while temp > temperature[end] 

        it = it + 1;
        bc = "no condition"
        C, T, iter,err_ctm = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,J1,J2,h,chi,precision_ctm,C,T,bc) 
       
        @show temperature[relevant]
        if temperature[relevant] <= tempprev && temperature[relevant] > temp
            println("ici")
            mm = magnetisation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs)
            energie = energy(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,ancilla_legs,J1,J2,h,hs)
            xi2, xi3, xi4, dq = correlation_length(C,T,gt)

            relevant = relevant + 1
            omega0 = 1;
            omega1 = exp(2*pi*im/3)
            omega2 = exp(4*pi*im/3);
            magne = 1/6*abs(omega0*mm[1,1]+omega1*mm[1,2]+omega2*mm[1,3]+omega0*mm[1,4]+omega1*mm[1,5]+omega2*mm[1,6])
            z2 = z2order_parameter(C,T,gt,tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs)

            push!(order_parameter_z2,z2)
            push!(order_parameter, magne)
            @show push!(tempera,temp)
            @show push!(ener, real(energie))
            push!(invcorr1, xi2[1,1])
            push!(eig3, xi3[1,1])
            push!(eig4, xi4[1,1])

            name_data = string("Results/U1_ResultsFUD",string(D),"temp",string(round(temp,digits = 5)),"_nfu_",string(nfu),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
                    save(name_data,
                    "C",C,
                    "T",T,
                    "tens_A",tens_A,
                    "cxd",cxd,
                    "cyd",cyd,
                    "physical_legs",physical_legs,
                    "ancilla_legs",ancilla_legs,
                    "D",D,
                    "magne",magne,
                    "magne_z",mm,
                    "z2",z2,
                    "J2",J2,
                    "h",h,
                    "J1",J1,
                    "hs",hs,
                    "nfu",nfu,
                    "energie",energie,
                    "temp",temp,
                    "xi2",xi2,
                    "xi3",xi3,
                    "xi4",xi4,
                    "dq",dq,
                    "numb_iter",iter,
                    "chi",chi,
                    "err_ctm",err_ctm,
                    "tens_a",tens_a,
                    "gt",gt,
                    "gg",gg,
                    "gx",gx,
                    "gy",gy,
                    "mu",mu,
                    "mu_hs",mu_hs,
                    "precision_ctm",precision_ctm)

        end

        @show temp, it
        j = 1
        for i = 1:N
            t1, t2 = full_update_x!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,gx,gy,D,C,T,i,j)
            C, T, iter,err_ctm = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,J1,J2,h,chi,precision_ctm,C,T,bc) 
        end

        i = 1
        for j = 1:N
            t1, t2 = full_update_y!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,gx,gy,D,C,T,i,j)
            C, T, iter,err_ctm = ctm!(tens_a,tens_A,cxd,cyd,gt,physical_legs,ancilla_legs,J1,J2,h,chi,precision_ctm,C,T,bc) 
        end

        for i = 1:N
            for j =1:1
                A = getproperty(tens_A,gt[i,j])
                mu_a = getproperty(mu,gt[i,j])
                A = A*mu_a
                A = noprime(A)
                setproperty!(tens_A,gt[i,j],A)
            end
        end

        for i = [1,4]
            for j =1:1
                A = getproperty(tens_A,gt[i,j])
                mu_a = getproperty(mu_hs,gt[i,j])
                A = A*mu_a
                A = noprime(A)
                setproperty!(tens_A,gt[i,j],A)
            end
        end

        renormalisation_tens_A!(tens_a,tens_A,cxd,cyd,physical_legs,ancilla_legs,gt,N)

        tempprev = 1/beta
        beta = beta + dbetafu
        temp = 1/beta 

    end

    @show ener
    @show order_parameter
    @show tempera

    name_data = string("Results/U1_ResultsD",string(D),"_nfu_",string(nfu),"temp",string(temp),"hs",string(round(1e4*hs,digits = 3)),"h",string(round(1e2*h,digits = 3)),".jld2")
    save(name_data,
        "C",C,
        "T",T,
        "tens_a",tens_a,
        "tens_A",tens_A,
        "cxd",cxd,
        "cyd",cyd,
        "D",D,
        "order_parameter",order_parameter,
        "J2",J2,
        "h",h,
        "J1",J1,
        "hs",hs,
        "nfu",nfu,
        "gt",gt,
        "N",N,
        "physical_legs",physical_legs,
        "ancilla_legs",ancilla_legs,
        "ener",ener,
        "tempera",tempera,
        "precision_ctm",precision_ctm);

    return C, T, tens_a, tens_A, cxd, cyd, beta, temp

    nothing

end