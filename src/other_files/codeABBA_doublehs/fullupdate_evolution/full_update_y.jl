
include("inverse_matrixM.jl")

function full_update_y!(tens_a::lattice2x1,tens_A::lattice2x1,cxd::lattice2x1,
    cyd::lattice2x1,gt::Matrix{Symbol},physical_legs::lattice_ind2x1,ancilla_legs::lattice_ind2x1,
    gx::lattice2x1,gy::lattice2x1,D::Int64,C::Vector{lattice2x1},T::Vector{lattice2x1},i,j,free_energy)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    precision_full_update = 1e-15


    ggy = getfield(gy, gt[f(i),f(j)])

    C1d = getfield(C[1],gt[f(i-1),f(j-1)])
    T1c = getfield(T[1],gt[f(i),f(j-1)])
    C2d = getfield(C[2],gt[f(i+1),f(j-1)])

    T4b = getfield(T[4],gt[f(i-1),f(j)])
    A = getfield(tens_A,gt[f(i),f(j)])
    T2b = getfield(T[2],gt[f(i+1),f(j)])

    T4d = getfield(T[4],gt[f(i-1),f(j+1)])
    CC = getfield(tens_A,gt[f(i),f(j+1)])
    T2d = getfield(T[2],gt[f(i+1),f(j+1)])

    C4b = getfield(C[4],gt[f(i-1),f(j+2)])
    T3a = getfield(T[3],gt[f(i),f(j+2)])
    C3b = getfield(C[3],gt[f(i+1),f(j+2)])

    Up = C1d*T1c*C2d
    Down = C4b*T3a*C3b
        
    ind_ya = commonind(A,getfield(cyd,gt[f(i),f(j)]))
    ind_yc = commonind(CC,getfield(cyd,gt[f(i),f(j+1)]))

    physical_legs_ia = getfield(physical_legs,gt[f(i),f(j)])
    ancilla_legs_ia = getfield(ancilla_legs,gt[f(i),f(j)])
    physical_legs_ic = getfield(physical_legs,gt[f(i),f(j+1)])
    ancilla_legs_ic = getfield(ancilla_legs,gt[f(i),f(j+1)])

    indsA = inds(A)
    indsC = inds(CC)
    all_but_ph_a = noncommoninds(indsA,[ind_ya, physical_legs_ia])
    all_but_ph_c = noncommoninds(indsC,[ind_yc,physical_legs_ic])

    u,s,v = svd(A,(all_but_ph_a))
    free_energy = free_energy + log(norm(s))
    s = s/norm(s)
    u_index_p = commonind(u,s)
    X = deepcopy(u); p = s*v
    all_but_sx = noncommoninds(inds(X),ancilla_legs_ia)
    X_prime = prime(X,all_but_sx)
            
    u,s,v = svd(CC,[ind_ya,physical_legs_ic])
    q = u*s; Y = deepcopy(v); 
    all_but_sy = noncommoninds(inds(Y),ancilla_legs_ic)
    Y_prime = prime(Y, all_but_sy)

    p0 = deepcopy(p)
    q0 = deepcopy(q)
    p0_prime = prime(p0, noncommoninds(inds(p0), [physical_legs_ia]))
    q0_prime = prime(q0, noncommoninds(inds(q0), [physical_legs_ic]))

    #=      
                | dag(cy22)
    - dag(cx22) - A  - cx11 -
                | cy11
    =#

    c1 = getfield(cxd,gt[f(i-1),f(j)])
    c2 = getfield(cyd,gt[f(i),f(j-1)])
    c3 = getfield(cxd,gt[f(i),f(j)])
    XX2 = X_prime*dag(X)
    XX = XX2*dag(c1)*dag(c2)*(c3)

    c1 = getfield(cxd,gt[f(i-1),f(j+1)])
    c2 = getfield(cyd,gt[f(i),f(j+1)])
    c3 = getfield(cxd,gt[f(i),f(j+1)])
    YY2 = Y_prime*dag(Y);
    YY = YY2*dag(c1)*(c2)*(c3)

    UU = (Up*T4b*XX*T2b)
    DD = (Down*T4d*YY*T2d)
    environment = DD*UU

    ida = ITensor(physical_legs_ia, dag(physical_legs_ia)')
    idc = ITensor(physical_legs_ic, dag(physical_legs_ic)')

    for ii = 1:4
        ida[physical_legs_ia=>ii, dag(physical_legs_ia)' =>ii] = 1;
        idc[physical_legs_ic=>ii, dag(physical_legs_ic)' =>ii] = 1;
    end

    distance = 1
    err_fullupdate = 1;

    p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))
    q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ic]))

    bb = ((environment*p0_prime*q0_prime))*ggy;

    ity = 0;

    while distance > precision_full_update && ity < 50

        ity = ity + 1; 

        p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))
        q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ic]))

        M = ((environment*q_prime)*dag(q))*ida
                
        b = noprime(bb,physical_legs_ic')
        b = noprime(b*dag(q))

        indd1 = commoninds(b,M)
        comb = combiner(indd1, dir = dir(physical_legs_ia))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1],indd2[3],indd2[2]], dir = -dir(physical_legs_ia))

        M_matrix = M*comb*cprime

        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(comb)*(cprime)

        p = inverse_M*b
        p = noprime(p)

        q_prime = prime(q, noncommoninds(inds(q), [physical_legs_ic]))
        p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

        M = (environment*p_prime)*dag(p)*idc
        b = noprime(bb,physical_legs_ia')
        b = noprime(b*dag(p))

        indd1 = commoninds(b,M)
        c = combiner(indd1, dir = dir(physical_legs_ic))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1], indd2[3], indd2[2]])

        M_matrix = M*c*cprime

        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(c)*(cprime)

        q = inverse_M*b
        q = noprime(q)
            
        q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ic]))
        p_prime = prime(p, noncommoninds(inds(p),[physical_legs_ia]))

        middle_term = p0_prime*q0_prime*ggy
        middle_term2 = noprime(middle_term,physical_legs_ia',physical_legs_ic')*ggy

        term1 = noprime(environment*middle_term2)*dag(p0)*dag(q0)
        term2 = environment*(q_prime*dag(q))*(p_prime*dag(p))
        term3 = noprime(environment*middle_term)*dag(p)*dag(q)

        Z = norm(array(((environment*q0_prime)*dag(q0))*p0_prime*dag(p0)))

        previous_err_fulluptade = deepcopy(err_fullupdate)
        err_fullupdate = norm(array(term1) + array(term2) - 2*array(term3))/norm(array(term1))
        @show distance = norm(err_fullupdate - previous_err_fulluptade)

    end

    free_energy = free_energy + (log(norm(q)) + log(norm(p)))

    p = p/norm(((p)))
    q = q/norm(((q)))

    A = X*p; 
    CC = q*Y; 

    setproperty!(tens_A, gt[f(i),f(j)], A)
    setproperty!(tens_A, gt[f(i),f(j+1)], CC)
    Aprime = prime(A, noncommoninds(inds(A), [physical_legs_ia, ancilla_legs_ia]))
    Cprime = prime(CC, noncommoninds(inds(CC), [physical_legs_ic, ancilla_legs_ic]))
            
    aaa = (Aprime*dag(A))*getfield(cxd,gt[f(i),f(j)])*dag(getfield(cxd,gt[f(i-1),f(j)]))*
    getfield(cyd,gt[f(i),f(j)])*dag(getfield(cyd,gt[f(i),f(j-1)]))
    ccc = Cprime*dag(CC)*getfield(cxd,gt[f(i),f(j+1)])*dag(getfield(cxd,gt[f(i-1),f(j+1)]))*
    getfield(cyd,gt[f(i),f(j+1)])*dag(getfield(cyd,gt[f(i),f(j)]))

    tens_a_up = deepcopy(tens_a)
    tens_a_down = deepcopy(tens_a)
    setproperty!(tens_a_up,gt[f(i),f(j)], aaa)
    setproperty!(tens_a_down,gt[f(i),f(j+1)], ccc)
    
    setproperty!(tens_a,gt[f(i),f(j)], aaa)
    setproperty!(tens_a,gt[f(i),f(j+1)], ccc)


    
    return tens_a_up, tens_a_down, free_energy

end


