

function PartitionPerSite(T::Vector{lattice2x1},C::Vector{lattice2x1},tens::lattice2x1,gt::Matrix{Symbol},
    i::Int64,j::Int64)

    # N = size(gt)[1]
    # f(x) = mod(x-1,N) + 1


    # t4_down = getfield(T[4],gt[f(i-1),f(j)])
    # c1 = getfield(C[1],gt[f(i-1),f(j-1)])
    # t1_right = getfield(T[1],gt[f(i), f(j-1)])

    # t1 = getfield(T[1],gt[f(i),f(j-1)])
    
    # t1_left = getfield(T[1],gt[f(i),f(j-1)])
    # c2 = getfield(C[2],gt[f(i+1),f(j-1)])
    # t2_down = getfield(T[2],gt[f(i+1),f(j)])

    # t4 = getfield(T[4],gt[f(i-1),f(j)])
    # a = getfield(tens,gt[f(i),f(j)])
    # t2 = getfield(T[2],gt[f(i+1),f(j)])

    # t4_up = getfield(T[4],gt[f(i-1),f(j)])
    # c4 = getfield(C[4],gt[f(i-1),f(j+1)])
    # t3_right = getfield(T[3],gt[f(i),f(j+1)])

    # t3 = getfield(T[3],gt[f(i),f(j+1)])

    # t2_up = getfield(T[2],gt[f(i+1),f(j)])
    # c3 = getfield(C[3],gt[f(i+1),f(j+1)])
    # t3_left = getfield(T[3],gt[f(i),f(j+1)])

    # indc2_left = commonind(c2,t1_left)
    # indc1_right = commonind(c1,t1_right)

    # indc2_down = commonind(c2,t2_down)
    # indc3_up = commonind(c3,t2_up)

    # indc3_left = commonind(c3, t3_left)
    # indc4_right = commonind(c4, t3_right)

    # indc4_up = commonind(c4,t4_up)
    # indc1_down = commonind(c1,t4_down)



    # carre = (c1*delta(dag(indc2_left),dag(indc1_right))*c2)
    # carre = carre*delta(dag(indc2_down),dag(indc3_up))*c3
    # carre = carre*delta(dag(indc4_right),dag(indc3_left))*c4
    # carre = carre*delta(dag(indc4_up),dag(indc1_down))

    # @show carre[]

    # vert = 

    # Z = ((((((((c1*t1)*t4)*a)*c2)*t2)*c4)*t3)*c3)[]


    # return ((Z*carre))#/(vert*hor))

    N = size(gt)[1]
    f(x) = mod(x-1,N)+1;

    c1d = getfield(C[1],gt[f(i-1),f(j-1)])
    t1c = getfield(T[1],gt[f(i),f(j-1)])
    t1d = getfield(T[1],gt[f(i+1),f(j-1)])
    c2c = getfield(C[2],gt[f(i+2),f(j-1)])
    
    t4b = getfield(T[4],gt[f(i-1),f(j)])
    a1 = getfield(tens,gt[f(i),f(j)])
    a2 = getfield(tens,gt[f(i+1),f(j)])
    t2a = getfield(T[2],gt[f(i+2),f(j)])

    t4d = getfield(T[4],gt[f(i-1),f(j+1)])
    a3 = getfield(tens,gt[f(i),f(j+1)])
    a4 = getfield(tens,gt[f(i+1),f(j+1)])
    t2c = getfield(T[2],gt[f(i+2),f(j+1)])

    c4b = getfield(C[4],gt[f(i-1),f(j+2)])
    t3a = getfield(T[3],gt[f(i),f(j+2)])
    t3b = getfield(T[3],gt[f(i+1),f(j+2)])
    c3a = getfield(C[3],gt[f(i+2),f(j+2)])

    # d -- c -- d -- c
    # b -- a -- b -- b
    # d -- c -- d -- c
    # b -- a -- b -- b
    

    C1 = *(c1d,t1c,t4b,a1)
    C2 = *(c2c,t2a,t1d,a2)
    C4 = *(c4b,t3a,t4d,a3)
    C3 = *(c3a,t2c,t3b,a4)
    Z = (C1*C2)*(C3*C4)#*(c1d,t1c,t4b,a1,c2c,t2a,t1d,a2,c4b,t3a,t4d,a3,c3a,t2c,t3b,a4)

    carre = (((c1d*c2c)*c3a)*c4b)[]
    vert = (((((((c1d*c4b)*t1c)*t3a)*t1d)*t3b)*c2c)*c3a)[];
    hor = (((((((c1d*c2c)*t4b)*t2a)*t4d)*t2c)*c4b)*c3a)[];

    return log(array((Z*carre)/(vert*hor))[1])

    
end