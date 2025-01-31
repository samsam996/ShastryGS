
include("inverse_matrixM.jl")

function full_update_x!(tens_a::lattice,tens_A::lattice,cxd::lattice,
    cyd::lattice,gt::Matrix{Symbol},physical_legs::lattice_ind,ancilla_legs::lattice_ind,
    gx::lattice,gy::lattice,D::Int64,C::Vector{lattice},T::Vector{lattice},i,j,free_energy)

    N = size(gt)[1]
    f(x) = mod(x-1,N) + 1

    precision_full_update = 1e-15

 
    ggx = getfield(gx, gt[f(i),f(j)])

    C1a = getfield(C[1],gt[f(i-1),f(j-1)])
    T1b = getfield(T[1],gt[f(i),f(j-1)])
    T1a = getfield(T[1],gt[f(i+1),f(j-1)])
    C2b = getfield(C[2],gt[f(i+2),f(j-1)])

    T4c = getfield(T[4],gt[f(i-1),f(j)])
    A = getfield(tens_A,gt[f(i),f(j)])
    B = getfield(tens_A,gt[f(i+1),f(j)])
    T2d = getfield(T[2],gt[f(i+2),f(j)])

    C4a = getfield(C[4],gt[f(i-1),f(j+1)])
    T3b = getfield(T[3],gt[f(i),f(j+1)])
    T3a = getfield(T[3],gt[f(i+1),f(j+1)])
    C3b = getfield(C[3],gt[f(i+2),f(j+1)])

    Left = C1a*T4c*C4a
    Right = C2b*T2d*C3b
        
    ind_xa = commonind(A,getfield(cxd,gt[f(i),f(j)]))
    ind_xb = commonind(B,getfield(cxd,gt[f(i+1),f(j)]))

    physical_legs_ia = getfield(physical_legs,gt[f(i),f(j)])
    ancilla_legs_ia = getfield(ancilla_legs,gt[f(i),f(j)])
    physical_legs_ib = getfield(physical_legs,gt[f(i+1),f(j)])
    ancilla_legs_ib = getfield(ancilla_legs,gt[f(i+1),f(j)])

    indsA = inds(A)
    indsB = inds(B)
    all_but_ph_a = noncommoninds(indsA,[ind_xa, physical_legs_ia])
    all_but_ph_b = noncommoninds(indsB,[ind_xb,physical_legs_ib])

    u,s,v = svd(A,(all_but_ph_a))
    free_energy = free_energy + log(norm(s))
    s = s/norm(s)
    u_index_p = commonind(u,s)
    X = deepcopy(u); p = s*v
    all_but_sx = noncommoninds(inds(X),ancilla_legs_ia)
    X_prime = prime(X,all_but_sx)
            
    u,s,v = svd(B,[ind_xa,physical_legs_ib])
    free_energy = free_energy + log(norm(s))
    s = s/norm(s)
    q = u*s; Y = deepcopy(v); 
    all_but_sy = noncommoninds(inds(Y),ancilla_legs_ib)
    Y_prime = prime(Y, all_but_sy)


    p0 = deepcopy(p)
    q0 = deepcopy(q)
    q0_prime = prime(q0, noncommoninds(inds(q0), [physical_legs_ib]))
    p0_prime = prime(p0, noncommoninds(inds(p0), [physical_legs_ia]))

    #=      
                | dag(cy22)
    - dag(cx22) - A  - cx11 -
                | cy11
    =#

    c1 = getfield(cxd,gt[f(i-1),f(j)])
    c2 = getfield(cyd,gt[f(i),f(j)])
    c3 = getfield(cyd,gt[f(i),f(j-1)])
    XX2 = X_prime*dag(X)
    XX = XX2*dag(c1)*c2*dag(c3)

    c1 = getfield(cyd,gt[f(i+1),f(j)])
    c2 = getfield(cxd,gt[f(i+1),f(j)])
    c3 = getfield(cyd,gt[f(i+1),f(j-1)])
    YY2 = Y_prime*dag(Y);
    YY = YY2*(c1)*c2*dag(c3)

    LL = (Left*T1b*XX*T3b)
    RR = (Right*T1a*YY*T3a)
    environment = LL*RR

    ida = ITensor(physical_legs_ia, dag(physical_legs_ia)')
    idb = ITensor(physical_legs_ib, dag(physical_legs_ib)')

    for ii = 1:4
        ida[physical_legs_ia=>ii, dag(physical_legs_ia)' =>ii] = 1;
        idb[physical_legs_ib=>ii, dag(physical_legs_ib)' =>ii] = 1;
    end

    distance = 1
    err_fullupdate = 1;

    q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ib]))
    p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

    bb = ((environment*p0_prime*q0_prime))*ggx;

    itx = 0; 

    while distance > precision_full_update && itx < 50


        itx = itx + 1

        q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ib]))
        p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

        M = ((environment*q_prime)*dag(q))*ida
                
        b = noprime(bb,physical_legs_ib')
        b = noprime(b*dag(q))

        indd1 = commoninds(b,M)
        c = combiner(indd1, dir = -dir(physical_legs_ia))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1], indd2[3], indd2[2]])

        M_matrix = M*c*cprime
 
        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(c)*(cprime)

        p = inverse_M*b
        p = noprime(p)

        q_prime = prime(q, noncommoninds(inds(q), [physical_legs_ib]))
        p_prime = prime(p, noncommoninds(inds(p), [physical_legs_ia]))

        M = (environment*p_prime)*dag(p)*idb
        b = noprime(bb,physical_legs_ia')
        b = noprime(b*dag(p))

        indd1 = commoninds(b,M)
        c = combiner(indd1, dir = dir(physical_legs_ib))
        indd2 = noncommoninds(inds(M), indd1)
        cprime = combiner([indd2[1],indd2[3],indd2[2]])

        M_matrix = M*c*cprime

        inverse_M_matrix = inverse_matrixM(M_matrix)
        inverse_M = inverse_M_matrix*(c)*(cprime)

        q = inverse_M*b
        q = noprime(q)
            
        q_prime = prime(q, noncommoninds(inds(q),[physical_legs_ib]))
        p_prime = prime(p, noncommoninds(inds(p),[physical_legs_ia]))

        middle_term = p0_prime*q0_prime*ggx
        middle_term2 = noprime(middle_term,physical_legs_ia',physical_legs_ib')*ggx

        term1 = noprime(environment*middle_term2)*dag(p0)*dag(q0)
        term2 = environment*(q_prime*dag(q))*(p_prime*dag(p))
        term3 = noprime(environment*middle_term)*dag(p)*dag(q)

        Z = norm(array(((environment*q0_prime)*dag(q0))*p0_prime*dag(p0)))

        previous_err_fulluptade = deepcopy(err_fullupdate)
        err_fullupdate = norm(array(term1) + array(term2) - 2*array(term3))/Z
        @show distance = norm(err_fullupdate - previous_err_fulluptade)

    end

    free_energy = free_energy + (log(norm(p)) + log(norm(q)))

    p = p/norm(((p)))
    q = q/norm(((q)))

    A = X*p; 
    B = q*Y; 
       
    setproperty!(tens_A, gt[f(i),f(j)], A)
    setproperty!(tens_A, gt[f(i+1),f(j)], B)
    Aprime = prime(A, noncommoninds(inds(A), [physical_legs_ia, ancilla_legs_ia]))
    Bprime = prime(B, noncommoninds(inds(B), [physical_legs_ib, ancilla_legs_ib]))
            
    aaa = (Aprime*dag(A))*getfield(cxd,gt[f(i),f(j)])*dag(getfield(cxd,gt[f(i-1),f(j)]))*
    getfield(cyd,gt[f(i),f(j)])*dag(getfield(cyd,gt[f(i),f(j-1)]))
    bbb = Bprime*dag(B)*getfield(cxd,gt[f(i+1),f(j)])*dag(getfield(cxd,gt[f(i),f(j)]))*
    getfield(cyd,gt[f(i+1),f(j)])*dag(getfield(cyd,gt[f(i+1),f(j-1)]))


    tens_a_left = deepcopy(tens_a)
    tens_a_right = deepcopy(tens_a)
    setproperty!(tens_a_left,  gt[f(i),f(j)], aaa)
    setproperty!(tens_a_right, gt[f(i+1),f(j)], bbb)

    setproperty!(tens_a,gt[f(i),f(j)], aaa)
    setproperty!(tens_a,gt[f(i+1),f(j)], bbb)

   
    
    return tens_a_left, tens_a_right, free_energy

end


