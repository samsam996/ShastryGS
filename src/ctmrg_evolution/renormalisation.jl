

function renormalisation!(C,T,gt::Matrix{Symbol},N::Int64)

    """
        Renormalise the tensors environement tensors to keep the norm under control. 
    """
    
    for i = 1:N
        for j = 1:N
            for k = 1:4
                tmp = deepcopy(getfield(C[k],gt[i,j]));
                tmp3 = tmp/norm(tmp);
                setproperty!(C[k],gt[i,j],tmp3);

                tmp2 = deepcopy(getfield(T[k],gt[i,j]));
                tmp4 = tmp2/norm(tmp2);
                setproperty!(T[k],gt[i,j],tmp4);
            end
        end
    end


    nothing
    
end