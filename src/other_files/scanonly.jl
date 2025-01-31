
using MKL, FileIO, JLD2, KrylovKit, LinearAlgebra, ITensors

include("SSMHamiltonian.jl")
include("simple_update_evolution/simpleUpdateJ1.jl")
include("observables/entropy_lambda.jl")

include("declare_indices.jl")
include("observables/correlation_length.jl")

println("BEGIN SIMPlE UPDATE")

BLAS.set_num_threads(10)
ITensors.Strided.set_num_threads(1)
# ITensors.enable_threaded_blocksparse(true)

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

    D = 10; d = 4;

    name = fieldnames(lattice);
    @show name
    gg = [name[1] name[2]; name[3] name[4]];
    gt = [name[1] name[3]; name[2] name[4]];    


    N = size(gg)[1];
    f(x) = mod(x-1,N) + 1;

    Gamma, lambdax, lambday, physical_legs, ancilla_legs = declare_indices(gt,d)

    @show eltype(getfield(Gamma, gt[1,1]))


    for i = 1:N
        for j = 1:N
            
            tmp_gamma = ITensor(inds(getfield(Gamma,gt[f(i),f(j)])));
            leg1 = commonind(tmp_gamma,getfield(lambdax,gt[f(i-1),f(j)]))
            leg2 = commonind(tmp_gamma,getfield(lambday,gt[f(i),f(j)]))
            leg3 = commonind(tmp_gamma,getfield(lambdax,gt[f(i),f(j)]))
            leg4 = commonind(tmp_gamma,getfield(lambday,gt[f(i),f(j-1)]))
            leg5 = getproperty(physical_legs,gt[i,j])
            leg6 = getproperty(ancilla_legs,gt[i,j])

            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>1,dag(leg6)=>1] = 1
            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>2,dag(leg6)=>2] = 1
            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>3,dag(leg6)=>3] = 1
            tmp_gamma[leg1=>1,leg2=>1,leg3=>1,leg4=>1,leg5=>4,dag(leg6)=>4] = 1
            
            setproperty!(Gamma,gt[f(i),f(j)],tmp_gamma);

            tmp_lx = ITensor(inds(getfield(lambdax,gt[f(i),f(j)])));
            tmp_ly = ITensor(inds(getfield(lambday,gt[f(i),f(j)])));

            xa11 = commonind(getfield(Gamma,gt[f(i),f(j)]),getfield(lambdax,gt[f(i),f(j)]))
            xa22 = commonind(getfield(Gamma,gt[f(i+1),f(j)]),getfield(lambdax,gt[f(i),f(j)]))

            ya11 = commonind(getfield(Gamma,gt[f(i),f(j)]),getfield(lambday,gt[f(i),f(j)]))
            ya22 = commonind(getfield(Gamma,gt[f(i),f(j+1)]),getfield(lambday,gt[f(i),f(j)]))

            for pp = 1:1                
                tmp_lx[xa11=>pp,xa22=>pp] = 1; 
                tmp_ly[ya11=>pp,ya22=>pp] = 1;
            end

            setproperty!(lambdax,gt[f(i),f(j)],tmp_lx);
            setproperty!(lambday,gt[f(i),f(j)],tmp_ly);

        end
    end

    entrop = entropy_lambda(getfield(lambdax,gt[1,1]));
    beta = 0;
    temp =  Inf;

    dbeta = 1e-3; 
    nsu = 2/dbeta;
    J1 = 0.63; J2 = 1.0; h = 200; hs = 0; #2.4e-4*h;  

    
    # temperature = LinRange(0.04,0.025,15)
    # temperature = LinRange(0.06, 0.03, 15)
    # temperature = LinRange(0.075,0.06,10)
    temperature = LinRange(0.04,0.025,15)
    
    relevant = 1

    nb_iter = 20000; iter = 0
    SS = zeros(nb_iter)
    tempera = zeros(nb_iter)

    while relevant < length(temperature)

        iter = iter + 1
        tempprev = 1/beta
        beta = beta + dbeta; 
        temp = 1/beta;
        tmp = entrop[1]; 
        Gamma,lambdax,lambday = simpleupdateJ1(Gamma,lambdax,lambday,physical_legs,gt,nsu,J1,J2,h,hs,D)
        entrop = entropy_lambda(getfield(lambdax,gt[1,1]));
        dS = tmp - entrop[1];
   
        println(temp)

        if temperature[relevant] < tempprev && temperature[relevant] > temp

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
            "physical_legs",physical_legs,
            "ancilla_legs",ancilla_legs,
            "gt",gt,
            "gg",gg)

            relevant = relevant + 1

        end

    end

end



println("END SIMPlE UPDATE")



