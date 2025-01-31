

include("horizontal_correlation.jl")
include("vertical_correlation.jl")

function energy(C,T,tens_a,tens_A,gt::Matrix{Symbol},cxd,cyd,physical_legs,parameters::Dict)

    ener_link = 0
    sx = 1/2*[0 1; 1 0]; sy = 1/2*[0 -1im; 1im 0]; sz = 1/2*[1 0; 0 -1]; id = [1 0; 0 1];

    Delta = parameters["Delta"]
    hx = parameters["hx"]
    hz = parameters["hz"]
    J1 = parameters["J1"]
    J2 = parameters["J2"]
    N = parameters["N"]
    model = parameters["model"]

    for i = 1:N
        for j = 1:1
            horizontal = horizontal_correlation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,i,j,parameters)
            vertical = vertical_correlation(C,T,tens_a,tens_A,gt,cxd,cyd,physical_legs,i,j,parameters)
            ener_link = ener_link + horizontal + vertical             
        end 
    end 

    obs1::Matrix{ComplexF64} = zeros(ComplexF64, 4, 4)
    if model=="XY"
        obs1 = kron(sz,sz) + Delta*(kron(sx,sx) + kron(sy,sy)) 
    elseif model=="XYZ"
        Jx = parameters["Jx"]
        Jy = parameters["Jy"]
        Jz = parameters["Jz"]
        obs1 = Jz*kron(sz,sz) + Jx*kron(sx,sx) + Jy*kron(sy,sy) 
    end

    obs2::Matrix{ComplexF64} =  kron(sz,id) + kron(id,sz);
    # obs2::Matrix{ComplexF64} =  kron(id,sy) #+ kron(id,id);
    obs3::Matrix{ComplexF64} =  kron(sx,id) + kron(id,sx);
    obs4::Matrix{ComplexF64} =  kron(sy,id) + kron(id,sy);

    Ps::Matrix{ComplexF64} = 1/2*[0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0]
    Pt1::Matrix{ComplexF64} = 1/2*[1 0 0 -1; 0 0 0 0; 0 0 0 0; -1 0 0 1]
    Pt2::Matrix{ComplexF64} = -1/2*[1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]
    Pt3::Matrix{ComplexF64} = 1/2*[0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0]

    expect_obs1 = 0
    expect_obs2 = 0
    expect_obs3 = 0

    mmx = []
    mmy = []
    mmz = []
    mPs = []
    mPt1 = []
    mPt2 = []
    mPt3 = []

    for i = 1:N
        for j = 1:1

            expect_obs1 = expect_obs1 + 
            local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,obs1)

            tmp2 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,obs2)
            expect_obs2 = expect_obs2 + tmp2
            
            tmp3 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,obs3)
            expect_obs3 = expect_obs3 + tmp3

            tmp4 = local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,obs4)

            push!(mmz,tmp2)
            push!(mmx,tmp3)
            push!(mmy,tmp4)


            push!(mPs, local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,Ps))
            push!(mPt1, local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,Pt1))
            push!(mPt2, local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,Pt2))
            push!(mPt3, local_obs(C,T,tens_a,tens_A,cxd,cyd,gt,physical_legs,i,j,Pt3))

        end
    end

    @show mmx
    @show mmy
    @show mmz

    @show real(mPs)
    @show real(mPt1)
    @show real(mPt2)
    @show real(mPt3)
    Pst = Dict("mpS"=>real(mPs), "mPt1"=>real(mPt1), "mPt2"=>real(mPt2), "mPt3"=>real(mPt3))

    energie = real(J1*ener_link + J2*expect_obs1 - hz*expect_obs2 - hx*expect_obs3)/(N)

    return energie, mmx, mmy, mmz, Pst

end