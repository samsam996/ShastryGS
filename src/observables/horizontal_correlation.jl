


function horizontal_correlation(C,T,tens_a,tens_A,gt,cx,cy,
    physical_legs,i::Int64,j::Int64,parameters::Dict)

    N = parameters["N"]
    model = parameters["model"]

    f(x) = mod(x-1,N) + 1
    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];
    ex = 0

    ss_xyz = zeros(ComplexF64,4,4,4,4)

    if model=="XY"
        Delta=parameters["Delta"]
        if abs(mod(i+j,2) - 0) < 1e-9 #pair
            ss_xyz = Delta*(kron(kron(sx,id),kron(sx,id)) + kron(kron(sy,id),kron(sy,id))) + kron(kron(sz,id),kron(sz,id)) + 
            Delta*(kron(kron(id,sx),kron(sx,id)) + kron(kron(id,sy),kron(sy,id))) + kron(kron(id,sz),kron(sz,id)) 
        end

        if abs(mod(i+j,2) - 1) < 1e-9 # impair
            ss_xyz = Delta*(kron(kron(id,sx),kron(sx,id)) +  kron(kron(id,sy),kron(sy,id))) +  kron(kron(id,sz),kron(sz,id)) +
            Delta*(kron(kron(id,sx),kron(id,sx)) + kron(kron(id,sy),kron(id,sy))) + kron(kron(id,sz),kron(id,sz));
        end
    elseif model=="XYZ"
        Jx = parameters["Jx"]
        Jy = parameters["Jy"]
        Jz = parameters["Jz"]
        if abs(mod(i+j,2) - 0) < 1e-9 #pair
            ss_xyz = Jx*kron(kron(sx,id),kron(sx,id)) + Jy*kron(kron(sy,id),kron(sy,id)) + Jz*kron(kron(sz,id),kron(sz,id)) + 
            Jx*kron(kron(id,sx),kron(sx,id)) + Jy*kron(kron(id,sy),kron(sy,id)) + Jz*kron(kron(id,sz),kron(sz,id)) 
        end

        if abs(mod(i+j,2) - 1) < 1e-9 # impair
            ss_xyz = Jx*kron(kron(id,sx),kron(sx,id)) + Jy*kron(kron(id,sy),kron(sy,id)) + Jz*kron(kron(id,sz),kron(sz,id)) +
            Jx*kron(kron(id,sx),kron(id,sx)) + Jy*kron(kron(id,sy),kron(id,sy)) + Jz*kron(kron(id,sz),kron(id,sz));
        end
    end

    ss_xyz = reshape(ss_xyz,(4,4,4,4))

    Adown = getproperty(tens_A, gt[f(i),f(j)])
    Bdown = getproperty(tens_A, gt[f(i+1),f(j)])

    ia = getproperty(physical_legs, gt[f(i),f(j)])
    ib = getproperty(physical_legs, gt[f(i+1),f(j)])

    # all_inda = noncommoninds(inds(Adown),[sa])
    # all_indb = noncommoninds(inds(Bdown),[sb])
    all_inda = inds(Adown)
    all_indb = inds(Bdown)
 
    Aup = prime(Adown, all_inda)
    Bup = prime(Bdown, all_indb)

    Obs_abx = ITensor(dag(ia),dag(ib),ia',ib')

    c1 = [0,0,0,0]
    for i1 = 1:4
        for i2 = 1:4
            for j1 = 1:4
                for j2 = 1:4
                    if abs(c1[i1] + c1[i2] - c1[j1] - c1[j2]) < 1e-9
                        Obs_abx[ia=>i1,ib=>i2,ia'=>j1,ib'=>j2] = ss_xyz[i1,i2,j1,j2]
                    end
                end
            end
        end
    end

    aa = Aup*dag(Adown)
    bb = Bup*dag(Bdown)
        
    comb1a = getproperty(cx, gt[f(i-1),f(j)])
    comb2a = getproperty(cy, gt[f(i),f(j)])
    comb3a = getproperty(cx, gt[f(i),f(j)])
    comb4a = getproperty(cy, gt[f(i),f(j-1)])

    comb1b = getproperty(cx, gt[f(i),f(j)])
    comb2b = getproperty(cy, gt[f(i+1),f(j)])
    comb3b = getproperty(cx, gt[f(i+1),f(j)])
    comb4b = getproperty(cy, gt[f(i+1),f(j-1)])

    aa = (((aa*dag(comb1a))*comb2a)*comb3a)*dag(comb4a)
    bb = (((bb*dag(comb1b))*comb2b)*comb3b)*dag(comb4b)
   
    #   C1a:i-1,j-1    T1b:i+0,j-1   T1a:i+1,j-1   C2b:i+2,j-1
    #   T4c:i-1,j+0    A1d:i+0,j+0   B1c:i+1,j+0   T2d:i+2,j+0
    #   C4a:i-1,j+1    T3b:i+0,j+1   T3a:i+1,j+1   C3b:i+2,j+1

    c1a = getproperty(C[1],gt[f(i-1),f(j-1)])
    t1b = getproperty(T[1],gt[f(i+0),f(j-1)])
    t1a = getproperty(T[1],gt[f(i+1),f(j-1)])
    c2b = getproperty(C[2],gt[f(i+2),f(j-1)])
    t4c = getproperty(T[4],gt[f(i-1),f(j)])
    a1 = getproperty(tens_a,gt[f(i),f(j)])
    a2 = getproperty(tens_a,gt[f(i+1),f(j)])
    t2d = getproperty(T[2],gt[f(i+2),f(j)])

    c4a = getproperty(C[4],gt[f(i-1),f(j+1)])
    t3b = getproperty(T[3],gt[f(i),f(j+1)])
    t3a = getproperty(T[3],gt[f(i+1),f(j+1)])
    c3b = getproperty(C[3],gt[f(i+2),f(j+1)])

    Left = (c1a*t4c)*c4a;
    Right = (c2b*t2d)*c3b;
    
    num = (((((((((Left)*t1b)*aa)*t3b)*dag(Obs_abx))*t1a)*bb)*t3a)*Right)[]
    denom = ((((Left)*t1b*a1*t3b)*t1a*a2*t3a)*Right)[]

    ex = num/denom

    return ex

end