


function OBC_PEPS3(tens_A,cxd,cyd,gt)

    N = size(gt)[1]
    f(x) = mod(x-1,N)+1

    DD = size(getproperty(tens_A,gt[1,1]))[1]
    C = [lattice(ITensor(Index(DD),Index(DD)),ITensor(Index(DD),Index(DD)),
    ITensor(Index(DD),Index(DD)),ITensor(Real,Index(DD),Index(DD))) for k=1:4];
    T = [lattice(ITensor(Real,Index(DD),Index(DD),Index(DD)),ITensor(Index(DD),Index(DD),Index(DD)),
    ITensor(Index(DD),Index(DD),Index(DD)),ITensor(Index(DD),Index(DD),Index(DD))) for k=1:4];
    
    cx = lattice()
    cy = lattice()
    
    for i = 1:N
        for j = 1:N
                        
            A = getproperty(tens_A,gt[i,j])
            cx1 = getproperty(cxd,gt[f(i-1),f(j)])
            cy2 = getproperty(cyd,gt[f(i),f(j)])
            cx3 = getproperty(cxd,gt[f(i),f(j)])
            cy4 = getproperty(cyd,gt[f(i),f(j-1)])

            ind1 = commonind(cx1,A)
            ind2 = commonind(cy2,A)
            ind3 = commonind(cx3,A)
            ind4 = commonind(cy4,A)

            c1 = ITensor(ind2',ind3',dag(ind2),dag(ind3))
            for i = 1: max(dim(ind2), dim(ind3))
                c1[ind2'=>i,ind3'=>i,dag(ind2)=>i,dag(ind3)=>i] = 1
            end
            setproperty!(C[1],gt[i,j],c1)
            t1 = ITensor(dag(ind1'),(ind2'),(ind3'),(ind1),dag(ind2),dag(ind3))
            for i = 1: max(dim(ind1),dim(ind2), dim(ind3))
                t1[dag(ind1')=>i,ind2'=>i,ind3'=>i,(ind1)=>i,dag(ind2)=>i,dag(ind3)=>i] = 1
            end
            t1 = t1*cy2
            setproperty!(T[1],gt[i,j],t1)

            c2 = ITensor(dag(ind1'),ind2',ind1,dag(ind2))
            for i = 1:max(dim(ind1),dim(ind2))
                c2[dag(ind1')=>i,ind2'=>i,ind1=>i,dag(ind2)=>i] = 1
            end
            setproperty!(C[2],gt[i,j],c2)
            t2 = ITensor(dag(ind4'),dag(ind1'),ind2',ind4,ind1,dag(ind2))
            for i = 1:max(dim(ind1),dim(ind2),dim(ind4))
                t2[dag(ind4')=>i,dag(ind1')=>i,ind2'=>i,ind4=>i,ind1=>i,dag(ind2)=>i] = 1
            end
            t2 = t2*dag(cx1)
            setproperty!(T[2],gt[i,j],t2)

            c3 = ITensor(dag(ind1'),dag(ind4'),ind1,ind4)
            for i = 1:max(dim(ind1),dim(ind4))
               c3[dag(ind1')=>i,dag(ind4')=>i,ind1=>i,ind4=>i] = 1
            end
            setproperty!(C[3],gt[i,j],c3)
            t3 = ITensor(dag(ind1'),dag(ind4'),ind3',ind1,ind4,dag(ind3))
            for i = 1:max(dim(ind1),dim(ind4),dim(ind3))
                t3[dag(ind1')=>i,dag(ind4')=>i,ind3'=>i,ind1=>i,ind4=>i,dag(ind3)=>i] = 1
            end
            t3 = t3*dag(cy4)
            setproperty!(T[3],gt[i,j],t3)

            c4 = ITensor(dag(ind4'),ind3',ind4,dag(ind3))
            for i = 1:max(dim(ind3),dim(ind4))
                c4[dag(ind4')=>i,ind3'=>i,ind4=>i,dag(ind3)=>i] = 1
            end
            setproperty!(C[4],gt[i,j],c4)
            t4 = ITensor(dag(ind4'),ind3',ind2',ind4,dag(ind3),dag(ind2))
            for i = 1:max(dim(ind4),dim(ind3),dim(ind2))
                t4[dag(ind4')=>i,ind3'=>i,ind2'=>i,ind4=>i,dag(ind3)=>i,dag(ind2)=>i] = 1
            end
            t4 = t4*cx3
            setproperty!(T[4],gt[i,j],t4)

        end
    end

    for i = 1:N
        for j = 1:N

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
        for j = 1:N

            c1 = getproperty(C[1],gt[f(i-1),f(j-1)])
            t1 = getproperty(T[1],gt[f(i+0),f(j-1)])
            comb = getproperty(cx,gt[f(i-1),f(j-1)])
            # @show inds(c1)
            # @show inds(comb)
            c1 = c1*comb;
            t1 = dag(comb)*t1
            setproperty!(C[1],gt[f(i-1),f(j-1)],c1)
            setproperty!(T[1],gt[f(i+0),f(j-1)],t1)

            t1 = getproperty(T[1],gt[f(i+0),f(j-1)])
            c2 = getproperty(C[2],gt[f(i+1),f(j-1)])
            comb = getproperty(cx,gt[f(i+0),f(j-1)])
            # println(i,j)
            # @show inds(comb)
            # @show inds(t1)
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