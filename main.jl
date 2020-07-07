

function read_Species_From_Text_File(file_name)
    return reaction_matrix
end



function transformation_matrix_ij(i,j,decay_rates)
    s_ij = ones(eltype(decay_rates[1]))
    # Loop over generations
    for l in 1: (i-j)
        # Loop over species idnex of species ancestor
        for m in j:(i-1)
            # println(l,m)
            # println(decay_rates[m], decay_rates[i])
            # println(decay_rates[m] / (decay_rates[m] - decay_rates[i]))
            s_ij *= decay_rates[m] / (decay_rates[m] - decay_rates[i])
        end
    end
    return s_ij
end



function transformation_matrix_inverse_ij(i,j,decay_rates)
     s_ij = ones(eltype(decay_rates[1]))
     s_ij = decay_rates[j]/ (decay_rates[i] - decay_rates[j])
    # Loop over generations
    for l in 1: (i-j-1)
        # Loop over species idnex of species ancestor
        for m in j:(i-1)
            println(l,m)
            println(decay_rates[m], decay_rates[i])
            println(decay_rates[m] / (decay_rates[m] - decay_rates[ji]))
            s_ij *= decay_rates[m] / (decay_rates[m] - decay_rates[j])
            println(s_ij)
        end
    end
    return s_ij    
end



function assemble_transformation_matrix(reaction_matrix)
    return transformation_matrix
end



function assemble_transformation_matrix_inverse(reaction_matrix)
    return transformation_matrix_inverse
end



function retrieve_generation_number(i,j) 
    return (i-j)
end


