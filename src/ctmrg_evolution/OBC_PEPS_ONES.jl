


function OBC_PEPS_ONES(tens_A,cxd,cyd,gt::Matrix{Symbol},physical_legs, N::Int64)

    """
        Set the boundary condition for the iPEPS simulations. 
    """

    f(x) = mod(x-1,N) + 1

    C = [lattice(N) for k = 1:4]
    T = [lattice(N) for k = 1:4]

    cx = lattice(N)
    cy = lattice(N)
    
    for i = 1:N
        for j = 1:1
                        
            A = getproperty(tens_A,gt[i,j])
            cx1 = getproperty(cxd,gt[f(i-1),f(j)])
            cy2 = getproperty(cyd,gt[f(i),f(j)])
            cx3 = getproperty(cxd,gt[f(i),f(j)])
            cy4 = getproperty(cyd,gt[f(i),f(j-1)])

            ind1 = commonind(cx1,A)
            ind2 = commonind(cy2,A)
            ind3 = commonind(cx3,A)
            ind4 = commonind(cy4,A)

            if i == 1 || i == 3 # || i== 6
                ia = getproperty(physical_legs, gt[i,j])
                triplet = ITensor(ia',dag(ia))

                triplet[ia=>1,ia'=>1] = 1.
                triplet[ia=>4,ia'=>4] = 1.
                triplet[ia=>3,ia'=>3] = 1/sqrt(2)
                triplet[ia=>2,ia'=>2] = 1/sqrt(2)
                triplet[ia=>2,ia'=>3] = 1/sqrt(2)
                triplet[ia=>3,ia'=>2] = 1/sqrt(2)
                triplet[ia=>3,ia'=>3] = 1/sqrt(2)
                A = noprime(A*triplet)
                
            else
                ia = getproperty(physical_legs, gt[i,j])
                triplet = ITensor(ia',dag(ia))
                triplet[ia=>3,ia'=>3] = 1/sqrt(2)
                triplet[ia=>2,ia'=>2] = -1/sqrt(2)
                triplet[ia=>2,ia'=>3] = -1/sqrt(2)
                triplet[ia=>3,ia'=>2] = 1/sqrt(2)

                A = noprime(A*triplet)
            end

            C1_prime = prime(A, ind2,ind3)
            c1 = C1_prime*dag(A)
            setproperty!(C[1],gt[i,j],c1)
            T1_prime = prime(A, ind1, ind2, ind3)
            t1 = T1_prime*dag(A)

            # if i == 1 || i == 3 || i==5
            #     t1 = ITensor(inds(t1))
            #     for k = 1:1#dim(ind1)
            #         t1[ind1'=>k, ind2'=>k, ind3'=>k, dag(ind1)=>k, dag(ind2)=>k, dag(ind3)=>k] = 1
            #     end
            # end

            t1 = t1*cy2
            setproperty!(T[1],gt[i,j],t1)

            C2_prime = prime(A, ind1, ind2)
            c2 = C2_prime*dag(A)
            setproperty!(C[2],gt[i,j],c2)
            T2_prime = prime(A, ind4, ind1, ind2)
            t2 = T2_prime*dag(A)

            # if i == 1 || i == 3 || i==5
            #     t2 = ITensor(inds(t2))
            #     for k = 1#:dim(ind1)
            #         t2[ind1'=>k, ind2'=>k, ind4'=>k, dag(ind1)=>k, dag(ind2)=>k, dag(ind4)=>k] = 1
            #     end
            # end

            t2 = t2*dag(cx1)
            setproperty!(T[2],gt[i,j],t2)

            C3_prime = prime(A, ind1, ind4)
            c3 = C3_prime*dag(A)
            setproperty!(C[3],gt[i,j],c3)
            T3_prime = prime(A, ind1, ind4, ind3)
            t3 = T3_prime*dag(A)

            # if i == 1 || i == 4
            #     t3 = ITensor(inds(t3))
            #     for k = 1:dim(ind1)
            #         t3[ind1'=>k, ind4'=>k, ind3'=>k, dag(ind1)=>k, dag(ind4)=>k, dag(ind3)=>k] = 1
            #     end
            # end

            t3 = t3*dag(cy4)
            setproperty!(T[3],gt[i,j],t3)

            C4_prime = prime(A, ind3, ind4)
            c4 = C4_prime*dag(A)
            setproperty!(C[4],gt[i,j],c4)
            T4_prime = prime(A, ind4, ind3, ind2)
            t4 = T4_prime*dag(A)

            # if i == 1 || i == 4
            #     t4 = ITensor(inds(t4))
            #     for k = 1:dim(ind1)
            #         t4[ind4'=>k, ind2'=>k, ind3'=>k, dag(ind4)=>k, dag(ind2)=>k, dag(ind3)=>k] = 1
            #     end
            # end

            t4 = t4*cx3
            setproperty!(T[4],gt[i,j],t4)

        end
    end

    for i = 1:N
        for j = 1:1

            c1 = getproperty(C[1],gt[f(i),f(j)])
            t1 = getproperty(T[1],gt[f(i+1),f(j)])
            t4 = getproperty(T[4],gt[f(i),f(j+1)])
            indx = commoninds(c1,t1)
            indy = commoninds(c1,t4)
            combx = combiner(indx, dir = -dir(indx[1]))
            comby = combiner(indy, dir = -dir(indy[1]))
            setproperty!(cx, gt[i,j], combx)
            setproperty!(cy, gt[i,j], comby)

        end
    end


    
    for i = 1:N
        for j = 1:1

            c1 = getproperty(C[1],gt[f(i-1),f(j-1)])
            t1 = getproperty(T[1],gt[f(i+0),f(j-1)])
            comb = getproperty(cx,gt[f(i-1),f(j-1)])
            c1 = c1*comb;
            t1 = dag(comb)*t1
            setproperty!(C[1],gt[f(i-1),f(j-1)],c1)
            setproperty!(T[1],gt[f(i+0),f(j-1)],t1)

            t1 = getproperty(T[1],gt[f(i+0),f(j-1)])
            c2 = getproperty(C[2],gt[f(i+1),f(j-1)])
            comb = getproperty(cx,gt[f(i+0),f(j-1)])

            t1 = t1*comb;
            c2 = dag(comb)*c2
            setproperty!(T[1],gt[f(i+0),f(j-1)],t1)
            setproperty!(C[2],gt[f(i+1),f(j-1)],c2)

            c2 = getproperty(C[2],gt[f(i+1),f(j-1)])
            t2 = getproperty(T[2],gt[f(i+1),f(j+0)])
            comb = getproperty(cy,gt[f(i+1),f(j-1)])
            c2 = c2*comb;
            t2 = dag(comb)*t2
            setproperty!(C[2],gt[f(i+1),f(j-1)],c2)
            setproperty!(T[2],gt[f(i+1),f(j+0)],t2)

            t2 = getproperty(T[2],gt[f(i+1),f(j+0)])
            c3 = getproperty(C[3],gt[f(i+1),f(j+1)])
            comb = getproperty(cy,gt[f(i+1),f(j+0)])
            t2 = t2*comb;
            c3 = dag(comb)*c3
            setproperty!(T[2],gt[f(i+1),f(j+0)],t2)
            setproperty!(C[3],gt[f(i+1),f(j+1)],c3)



            c4 = getproperty(C[4],gt[f(i-1),f(j+1)])
            t3 = getproperty(T[3],gt[f(i+0),f(j+1)])
            comb = getproperty(cx,gt[f(i-1),f(j+1)])
            c4 = c4*comb;
            t3 = dag(comb)*t3
            setproperty!(C[4],gt[f(i-1),f(j+1)],c4)
            setproperty!(T[3],gt[f(i+0),f(j+1)],t3)

            t3 = getproperty(T[3],gt[f(i+0),f(j+1)])
            c3 = getproperty(C[3],gt[f(i+1),f(j+1)])
            comb = getproperty(cx,gt[f(i+0),f(j+1)])
            t3 = t3*comb;
            c3 = dag(comb)*c3
            setproperty!(T[3],gt[f(i+0),f(j+1)],t3)
            setproperty!(C[3],gt[f(i+1),f(j+1)],c3)



            c1 = getproperty(C[1],gt[f(i-1),f(j-1)])
            t4 = getproperty(T[4],gt[f(i-1),f(j+0)])
            comb = getproperty(cy,gt[f(i-1),f(j-1)])
            c1 = c1*comb;
            t4 = dag(comb)*t4
            setproperty!(C[1],gt[f(i-1),f(j-1)],c1)
            setproperty!(T[4],gt[f(i-1),f(j+0)],t4)

            t4 = getproperty(T[4],gt[f(i-1),f(j+0)])
            c4 = getproperty(C[4],gt[f(i-1),f(j+1)])
            comb = getproperty(cy,gt[f(i-1),f(j+0)])
            t4 = t4*comb;
            c4 = dag(comb)*c4
            setproperty!(T[4],gt[f(i-1),f(j+0)],t4)
            setproperty!(C[4],gt[f(i-1),f(j+1)],c4)


        end
    end




    return C,T

end