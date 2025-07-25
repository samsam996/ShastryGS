
function truncate_ey(C::ITensor,A::ITensor,i::Int,j::Int,gy,
    gt::Matrix{Symbol},physical_legs,ind_new::Index,D::Int64,free_energy,N::Int64)

    f(x) = mod(x-1,N) + 1

    # A : xb1 ya xa2 yc1 ia sa
    # C : xd1 yc2 xc2 ya ic sc

    ia = getproperty(physical_legs,gt[f(i),f(j)])
    ic = getproperty(physical_legs,gt[f(i),f(j+1)])
    exphy = getproperty(gy,gt[f(i),f(j)])

    indsC = inds(C);
    qc = noncommoninds(indsC,[ind_new,ic]);

    uc,sc,vc = svd(C,(qc));
    ua,sa,va = svd(A,(ia,ind_new));
    ucindex = commonind(uc,sc);

    ua = ua*sa;
    vc = vc*sc;

    M = exphy*(vc*ua);
    M = noprime(M)
    um,sm,vm = svd(M,(ucindex,ic),maxdim = D, cutoff = 1e-12);
    free_energy = free_energy + log(norm(sm))
    lambda2 = sm/norm(sm);

    C = uc*um; 
    A = va*vm;

    return C,A,lambda2,free_energy 

end