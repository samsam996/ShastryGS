

function PartitionPerSite(T::Vector{lattice},C::Vector{lattice},tens::lattice,gt::Matrix{Symbol},i::Int64,j::Int64)

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


    t1a = getfield(T[1], gt[i,j])
    t2a = getfield(T[2], gt[i,j])
    t3a = getfield(T[3], gt[i,j])
    t4a = getfield(T[4], gt[i,j])

    c3a = getfield(C[3],gt[f(i),f(j)])
    c2b = getfield(C[2],gt[f(i+1),f(j)])
    c4f = getfield(C[4],gt[f(i+5),f(j)])

    c1a = getfield(C[1],gt[f(i-1),f(j-1)])
    t1b = getfield(T[1],gt[f(i),f(j-1)])
    t1c = getfield(T[1],gt[f(i+1),f(j-1)])
    c2d = getfield(C[2],gt[f(i+2),f(j-1)])
    
    t4f = getfield(T[4],gt[f(i-1),f(j)])
    a1 = getfield(tens,gt[f(i),f(j)])
    a2 = getfield(tens,gt[f(i+1),f(j)])
    t2c = getfield(T[2],gt[f(i+2),f(j)])

    t4g = getfield(T[4],gt[f(i-1),f(j+1)])
    a3 = getfield(tens,gt[f(i),f(j+1)])
    a4 = getfield(tens,gt[f(i+1),f(j+1)])
    t2b = getfield(T[2],gt[f(i+2),f(j+1)])

    c4d = getfield(C[4],gt[f(i-1),f(j+2)])
    t3g = getfield(T[3],gt[f(i),f(j+2)])
    t3f = getfield(T[3],gt[f(i+1),f(j+2)])
    c3a = getfield(C[3],gt[f(i+2),f(j+2)])


    # a -- b -- c -- d
    # f -- a -- b -- c
    # g -- f -- a -- b
    # d -- g -- f -- a

    carree = (((c1a*c2b)*c3a)*c4f)[]
  
    # carree = carree*delta(dag(commonind(c2b,t2a)),commonind(t2b,c3a))*c3a
    # carree = carree*delta(dag(commonind(c4d,t3a)),(commonind(t3d,c3a)))*c4d
    # carree = carree*delta(dag(commonind(c1a,t4d)), commonind(t4a,c4d))



    C1 = *(c1a,t1b,t4f,a1)
    C2 = *(c2d,t1c,t2c,a2)
    C3 = *(c3a,t2b,t3f,a4)
    C4 = *(c4d,t4g,t3g,a3)
    Z = ((C1*C2)*(C3*C4))[]#*(c1d,t1c,t4b,a1,c2c,t2a,t1d,a2,c4b,t3a,t4d,a3,c3a,t2c,t3b,a4)

    # carre = (((c1d*c2c)*c3a)*c4b)[]
   


    # a -- b -- c -- d
    # f -- a -- b -- c
    # g -- f -- a -- b
    # d -- g -- f -- a

    c4f = getproperty(C[4],gt[f(i+5),j])
    c3c = getproperty(C[3],gt[f(i+2),j])
    t3b = getproperty(T[3],gt[f(i+1),j])
    t3a = getproperty(T[3],gt[i,j])

    c3g = getproperty(C[3],gt[f(i+4),f(j)])
    t2a = getproperty(T[2],gt[f(i),f(j)])
    t2f = getproperty(T[2],gt[f(i+5),f(j)])

    vert = (((((((c1a*c4f)*t1b)*t3a)*t1c)*t3b)*c2d)*c3c)[];
    
    hor = (((((((c1a*c2b)*t4f)*t2a)*t4g)*t2f)*c4d)*c3g)[];

    @show xx = (abs(Z*carree)/abs(vert*hor))

    return log(xx)

    
end