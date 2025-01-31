



function renormalisation_tens_A(tens_a::lattice,tens_A::lattice,cxd::lattice,
    cyd::lattice,physical_legs::lattice_ind,ancilla_legs::lattice_ind,gt::Matrix{Symbol},N::Int64,free_energy)

    f(x) = mod(x-1,N) + 1

    for i = 1:N
        for j = 1:1

            A = getfield(tens_A, gt[i,j])
            free_energy = free_energy + log(norm((A)))
            A = A/norm(((A)))
            setfield!(tens_A, gt[i,j],A)

            physical_legs_ia = getfield(physical_legs,gt[f(i),f(j)])
            ancilla_legs_ia = getfield(ancilla_legs,gt[f(i),f(j)])
  
            Aprime = prime(A, noncommoninds(inds(A), [physical_legs_ia, ancilla_legs_ia]))

            aaa = (Aprime*dag(A))*getfield(cxd,gt[f(i),f(j)])*dag(getfield(cxd,gt[f(i-1),f(j)]))*
            getfield(cyd,gt[f(i),f(j)])*dag(getfield(cyd,gt[f(i),f(j-1)]))
     
            setproperty!(tens_a,gt[f(i),f(j)], aaa)

        end
    end


    return tens_a, tens_A, free_energy

end