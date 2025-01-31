

function inverse_matrixM(A::ITensor)

    i1 = inds(A)[1];
    i2 = inds(A)[2];
    A_matrix = matrix(A, i1, i2);

    A_inverse = ITensor(dag(inds(A)));
    
    nb_block_i1 = length(i1.space);
    nb_block_i2 = length(i2.space);
    
    dimension_i1 = [];
    dimension_i2 = [];
    
    for i = 1:nb_block_i1
        push!(dimension_i1,blockdim(i1.space,i));
    end
    for i = 1:nb_block_i2
        push!(dimension_i2,blockdim(i2.space,i));
    end
    
    for i = 1:nb_block_i1
        for j = 1:nb_block_i2
            if hassameinds(qn(i1.space,i),qn(i2.space,j))
    
                Block = [];
                
                if i == 1 && j != 1

                    block = A_matrix[1:sum(dimension_i1[1]),
                    sum(dimension_i2[1:j-1])+1:sum(dimension_i2[1:j])] 
                    block = 1/2*(block+transpose(block))
                    inv_block = inv(block)

                    for ii = 1:sum(dimension_i1[1])
                        for jj = sum(dimension_i2[1:j-1])+1:sum(dimension_i2[1:j])
                            A_inverse[i1=>ii,  i2=>jj] = inv_block[jj - sum(dimension_i2[1:j-1]),ii]
                        end
                    end
    
                elseif j == 1 && i !=1

                    block = A_matrix[sum(dimension_i1[1:i-1])+1:sum(dimension_i1[1:i]),
                    1:sum(dimension_i2[1])] 
                    block = 1/2*(block+transpose(block))
                    inv_block = inv(block)

                    for ii = sum(dimension_i1[1:i-1])+1:sum(dimension_i1[1:i]) 
                        for jj = 1:sum(dimension_i2[1])
                            A_inverse[i1=>ii,  i2=>jj] = inv_block[jj, ii - sum(dimension_i1[1:i-1])];
                        end
                    end
    
                elseif i == 1 && j == 1

                    block = A_matrix[1:sum(dimension_i1[1]),1:sum(dimension_i2[1])] 
                    block = 1/2*(block+transpose(block))
                    inv_block = inv(block);
 
                    for ii = 1:sum(dimension_i1[1])
                        for jj = 1:sum(dimension_i2[1])
                            A_inverse[i1=>ii,  i2=>jj] = inv_block[jj, ii];
                        end
                    end
    
                elseif i != 1 && j != 1 

                    block = A_matrix[sum(dimension_i1[1:i-1])+1:sum(dimension_i1[1:i]),
                    sum(dimension_i2[1:j-1])+1:sum(dimension_i2[1:j])] 
                    block = Symmetric(1/2*(block+transpose(block)))
                    inv_block = inv(real(block));

                    for ii = sum(dimension_i1[1:i-1])+1:sum(dimension_i1[1:i])
                        for jj = sum(dimension_i2[1:j-1])+1:sum(dimension_i2[1:j])
                            
                            A_inverse[i1=>ii,  i2=>jj] = inv_block[jj - sum(dimension_i2[1:j-1]),
                            ii - sum(dimension_i1[1:i-1])];

                        end
                    end
    
                end
    
    
            end        
        end
    end

    return A_inverse;

end


