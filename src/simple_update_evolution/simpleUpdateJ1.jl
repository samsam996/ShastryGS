
# include("diag_sqrt.jl")
include("invert_diag.jl")
include("thermal_evolve_x.jl")
include("thermal_evolve_y.jl")

function simpleupdateJ1(Gamma,lambdax,lambday,
    physical_legs,gt::Matrix{Symbol},nsu::Float64,parameters,free_energy)

    N = parameters["N"]
    D = parameters["D"]
    gx,gy,mu = SSMHamiltonian(gt,physical_legs,nsu,parameters)

    Gamma,lambdax,free_energy = thermal_evolve_x(Gamma,lambdax,lambday,physical_legs,gx,gt,D,free_energy,N)
    Gamma,lambday,free_energy = thermal_evolve_y(Gamma,lambdax,lambday,physical_legs,gy,gt,D,free_energy,N)

    for i = 1:N
        for j =1:1
            gamma = getproperty(Gamma,gt[i,j])
            mu_a = getproperty(mu,gt[i,j])
            gamma = gamma*mu_a
            gamma = noprime(gamma)
            setproperty!(Gamma,gt[i,j],gamma)
        end
    end

    return Gamma,lambdax,lambday,free_energy

end