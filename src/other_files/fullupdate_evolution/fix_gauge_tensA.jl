

function fix_gauge_tensA!(tens_A::lattice, gt::Matrix{Symbol}, cxd, cyd, N, D)


    f(x) = mod(x-1,N) + 1

    for j = 1:N
        for i = 1:N

            A = getfield(tens_A, gt[f(i),f(j)])
            C = getfield(tens_A, gt[f(i),f(j+1)])

            cc1 = getfield(cyd,gt[f(i),f(j-1)])
            cc2 = getfield(cyd,gt[f(i),f(j)])

            xa1 = commonind(cc1, A)
            xa2 = commonind(cc2, A)
            
            A = prime(A, xa1)
            M = A*C
            u,s,v = svd(M, noncommoninds(inds(A),[xa2]), maxdim = D)
            s = s/maximum(s)
            u = noprime(u)
            sqrt_s = diag_sqrt(s)
            u_index = commonind(u,s)
            v_index = commonind(s,v)

            A = (u*sqrt_s)
            A = A*delta(dag(v_index),(xa2))

            C = (sqrt_s*v)
            C = C*delta((u_index),dag(xa2))

            # @show A
            # @show C

            setfield!(tens_A, gt[f(i),f(j)], A)            
            setfield!(tens_A, gt[f(i),f(j+1)], C)            

        end
    end


    for i = 1:N
        for j = 1:N

            A = getfield(tens_A, gt[f(i),f(j)])
            B = getfield(tens_A, gt[f(i+1),f(j)])
            # ia = getfield(physical_legs, gt[f(i),f(j)])
            # ib = getfield(physical_legs, gt[f(i+1),f(j)])

            cc1 = getfield(cxd,gt[f(i-1),f(j)])
            cc2 = getfield(cxd,gt[f(i),f(j)])

            xa1 = commonind(cc1, A)
            xa2 = commonind(cc2, A)
            
            A = prime(A, xa1)
            M = A*B
            u,s,v = svd(M, noncommoninds(inds(A),[xa2]), maxdim = D)
            s = s/maximum(s)
            u = noprime(u)
            sqrt_s = diag_sqrt(s)
            u_index = commonind(u,s)
            v_index = commonind(s,v)

            A = (u*sqrt_s)
            A = A*delta(dag(v_index),(xa2))

            B = (sqrt_s*v)
            B = B*delta((u_index),dag(xa2))

            setfield!(tens_A, gt[f(i),f(j)], A)            
            setfield!(tens_A, gt[f(i+1),f(j)], B)            

        end
    end



 

end