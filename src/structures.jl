


function cyclic_matrix(name::NTuple,N)

    g2 = Matrix{typeof(name[1])}(undef, N, N)  
    for i in 1:N
        g2[i, :] = [name[(mod(j - i, N) + 1)] for j in 1:N]
    end
    return g2
end

function countercyclic_matrix(name::NTuple,N)

    g2 = Matrix{typeof(name[1])}(undef, N, N)  
    for i in 1:N
        g2[i, :] = [name[(mod(j - i + 1, N) + 1)] for j in 1:N]
    end
    return g2
end

mutable struct lattice
    A::ITensor
    B::ITensor
    C::ITensor
    D::ITensor
    E::ITensor
    F::ITensor
    G::ITensor
    H::ITensor
    gg::Matrix{Symbol}
    gt::Matrix{Symbol}
    
    function lattice(N)

        name = fieldnames(lattice)
        g1 = cyclic_matrix(name, N)
        g2 = countercyclic_matrix(name,N)

        new(
        ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),g1,g2)
    end
    
end

mutable struct lattice_ind
    A::Index
    B::Index
    C::Index
    D::Index
    E::Index
    F::Index
    G::Index
    H::Index

    function lattice_ind()
        new(ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor(),ITensor())
    end


    function lattice_ind(nature_of_the_legs::String)

        if nature_of_the_legs == "physical"
            new(
            Index([QN(0)=>4],"ia"),Index([QN(0)=>4],"ib"),
            Index([QN(0)=>4],"ic"),Index([QN(0)=>4],"id"),
            Index([QN(0)=>4],"ie"),Index([QN(0)=>4],"if"),
            Index([QN(0)=>4],"ig"),Index([QN(0)=>4],"ih"))
        elseif nature_of_the_legs == "ancilla"
            new(
            Index([QN(0)=>4],"sa"),Index([QN(0)=>4],"sb"),
            Index([QN(0)=>4],"sc"),Index([QN(0)=>4],"sd"),
            Index([QN(0)=>4],"se"),Index([QN(0)=>4],"sf"),
            Index([QN(0)=>4],"ig"),Index([QN(0)=>4],"ih"))
        elseif nature_of_the_legs == "imaginary"
            new(
            Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]),Index([QN(0)=>1]))
        end

    end

end
