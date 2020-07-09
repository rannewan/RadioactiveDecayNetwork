

function read_Species_From_Text_File(file_name)
    return reaction_matrix
end



function transformation_matrix_ij(i,j,decay_rates)
    # s_ij = ones(eltype(decay_rates[1]))
    s_ij = 1.0
    # Loop over species ancestor index l for species i with (i-j) generations
    for l in j: (i-1)
        # println(l)
        # println(decay_rates[l])
        # println( (decay_rates[l] - decay_rates[i]))
        # println(decay_rates[l] / (decay_rates[l] - decay_rates[i]))
        s_ij *= decay_rates[l] / (decay_rates[l] - decay_rates[i])
    end
    return s_ij
end



function transformation_matrix_inverse_ij(i,j,decay_rates)
     # s_ij = ones(eltype(decay_rates[1]))
     s_ij = decay_rates[j]/ (decay_rates[i] - decay_rates[j])
     # Loop over species ancestor index l for species i with (i-j) generations
    for l in (j+1): (i-1)
        # println(i,j)
        # println(l)
        # println(decay_rates[l])
        # println(decay_rates[j])
        # println( (decay_rates[l] - decay_rates[j]))
        # println(decay_rates[l] / (decay_rates[l] - decay_rates[j]))
        s_ij *= decay_rates[l] / (decay_rates[l] - decay_rates[j])
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

function transformation_matrix_inverse_5x5()
    λ = [1.0; 2.0; 3.0; 4.0; 5.0; ]

    s = zeros(5,5)
    s[1,1] = 1
    s[2,1] = λ[1]/(λ[2]-λ[1])
    s[3,1] = λ[1]/(λ[3]-λ[1]) * λ[2]/(λ[2]-λ[1])
    s[4,1] = λ[1]/(λ[4]-λ[1]) * λ[2]/(λ[2]-λ[1]) * λ[3]/(λ[3]-λ[1])
    s[5,1] = λ[1]/(λ[5]-λ[1]) * λ[2]/(λ[2]-λ[1]) * λ[3]/(λ[3]-λ[1]) * λ[4]/(λ[4]-λ[1])

    s[2,2] = 1
    s[3,2] = λ[2]/(λ[3]-λ[2])
    s[4,2] = λ[2]/(λ[4]-λ[2]) * λ[3]/(λ[3]-λ[2])
    s[5,2] = λ[2]/(λ[5]-λ[2]) * λ[3]/(λ[3]-λ[2]) * λ[4]/(λ[4]-λ[2])

    s[3,3] = 1
    s[4,3] = λ[3]/(λ[4]-λ[3])
    s[5,3] = λ[3]/(λ[5]-λ[3]) * λ[4]/(λ[4]-λ[3])

    s[4,4] = 1
    s[5,4] = λ[4]/(λ[5]-λ[4])

    s[5,5] = 1

    return s
end


function transformation_matrix_5x5()
    λ = [1.0; 2.0; 3.0; 4.0; 5.0; ]

    s = zeros(5,5)
    s[1,1] = 1
    s[2,1] = λ[1]/(λ[1]-λ[2])
    s[3,1] = λ[1]/(λ[1]-λ[3]) * λ[2]/(λ[2]-λ[3])
    s[4,1] = λ[1]/(λ[1]-λ[4]) * λ[2]/(λ[2]-λ[4]) * λ[3]/(λ[3]-λ[4])
    s[5,1] = λ[1]/(λ[1]-λ[5]) * λ[2]/(λ[2]-λ[5]) * λ[3]/(λ[3]-λ[5]) * λ[4]/(λ[4]-λ[5])

    s[2,2] = 1
    s[3,2] = λ[2]/(λ[2]-λ[3])
    s[4,2] = λ[2]/(λ[2]-λ[4]) * λ[3]/(λ[3]-λ[4])
    s[5,2] = λ[2]/(λ[2]-λ[5]) * λ[3]/(λ[3]-λ[5]) * λ[4]/(λ[4]-λ[5])

    s[3,3] = 1
    s[4,3] = λ[3]/(λ[3]-λ[4])
    s[5,3] = λ[3]/(λ[3]-λ[5]) * λ[4]/(λ[4]-λ[5])

    s[4,4] = 1
    s[5,4] = λ[4]/(λ[4]-λ[5])

    s[5,5] = 1

    return s
end
