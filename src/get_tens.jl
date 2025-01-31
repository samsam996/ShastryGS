

# include("simple_update_evolution/diag_sqrt.jl")

function get_tens(Gamma,lambdax,lambday,physical_legs,gt,N)

    f(x) = mod(x-1,N) + 1

    tens_a = ShastryGS.lattice(N)
    cx = ShastryGS.lattice(N)
    cy = ShastryGS.lattice(N)

    tens_A = deepcopy(Gamma)
    
    for i = 1:N
        for j = 1:1
    
            Gamma0 = getproperty(tens_A,gt[f(i),f(j)])
            Gamma1 = getproperty(tens_A,gt[f(i+1),f(j)])

            lambdax1 = getproperty(lambdax,gt[f(i),f(j)])
            ind1 = commonind(Gamma0,lambdax1)
            ind2 = commonind(lambdax1,Gamma1)

            Gamma0 = Gamma0*ShastryGS.diag_sqrt(lambdax1)
            Gamma0 = Gamma0*delta((ind1),dag(ind2))
            Gamma1 = ShastryGS.diag_sqrt(lambdax1)*Gamma1
            setproperty!(tens_A,gt[f(i),f(j)],Gamma0)
            setproperty!(tens_A,gt[f(i+1),f(j)],Gamma1)

        end
    end

    for i = 1:N
        for j = 1:1
    
            Gamma0 = getproperty(tens_A,gt[f(i),f(j)])
            Gamma1 = getproperty(tens_A,gt[f(i),f(j+1)])

            lambday1 = getproperty(lambday,gt[f(i),f(j)])
            ind1 = commonind(Gamma0,lambday1)
            ind2 = commonind(lambday1,Gamma1)

            Gamma0 = Gamma0*ShastryGS.diag_sqrt(lambday1)
            Gamma0 =Gamma0*delta((ind1),dag(ind2))
            Gamma1 = ShastryGS.diag_sqrt(lambday1)*Gamma1
            setproperty!(tens_A,gt[f(i),f(j)],Gamma0)
            setproperty!(tens_A,gt[f(i),f(j+1)],Gamma1)

        end
    end

    for i = 1:N
        for j = 1:1

            ia = getproperty(physical_legs,gt[f(i),f(j)])
            A = getproperty(tens_A,gt[f(i),f(j)])
            indA = inds(A)
            ind_dgtb = noncommoninds(indA, [ia])
            A_prime = prime(A,ind_dgtb)
            a = A_prime*dag(A)
            setproperty!(tens_a,gt[f(i),f(j)],a)

        end
    end

    tens_a_tmp = deepcopy(tens_a)

    for i = 1:N
        for j = 1:1
    
            a1 = getproperty(tens_a,gt[f(i),f(j)])
            a2 = getproperty(tens_a,gt[f(i+1),f(j)])
            index = commoninds(a1,a2)
            cx11 = combiner(index[4],index[2],dir = -dir(index[1]))
            a1 = a1*cx11

            a2 = a2*dag(cx11)
            setproperty!(tens_a,gt[f(i),f(j)],a1)
            setproperty!(tens_a,gt[f(i+1),f(j)],a2)
            setproperty!(cx,gt[f(i),f(j)],cx11)
    
        end
    end



    for i = 1:1
        for j = 1:N

                a1 = getproperty(tens_a,gt[f(i),f(j)])
                a2 = getproperty(tens_a,gt[f(i),f(j+1)])
                index = commoninds(a1,a2)
                cy1 = combiner(index[3],index[2], dir = -dir(index[1])) #1,1
                a1 = a1*(cy1)
                a2 = a2*dag(cy1)

                setproperty!(tens_a,gt[f(i),f(j)],a1)
                setproperty!(tens_a,gt[f(i),f(j+1)],a2)
                setproperty!(cy,gt[f(i),f(j)],cy1)
    
        end
    end
    

    return tens_a, tens_A, cx, cy

end