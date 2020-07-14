

function read_Species_From_Text_File(file_name)
    return reaction_matrix
end



function calculate_concentration(c_0, λ, Δt)
    Λ = Diagonal(λ)
    S_inv = assemble_transformation_matrix_inverse(λ)
    S = assemble_transformation_matrix(λ)
    A = S_inv * Λ * S * Δt
    return exp(A) * c_0
end



function transformation_matrix_inverse_ij(i, j, λ)
    s_ij::eltype(λ) = 1.0
    # Loop over species ancestor index l for species i with (i-j) generations
    for l = j:(i-1)
        s_ij *= λ[l] / (λ[l] - λ[i])
    end
    return s_ij
end



function transformation_matrix_ij(i, j, λ)
    s_ij = λ[j] / (λ[i] - λ[j])
    # Loop over species ancestor index l for species i with (i-j) generations
    for l = (j+1):(i-1)
        s_ij *= λ[l] / (λ[l] - λ[j])
    end
    return s_ij
end



function assemble_transformation_matrix(λ)
    n = length(λ)
    S = zeros(eltype(λ), length(λ), length(λ))
    S = UnitLowerTriangular(S)
    for i = 2:n
        for j = 1:(i-1)
            S[i, j] = transformation_matrix_ij(i, j, λ)
        end
    end
    return S
end



function assemble_transformation_matrix_inverse(λ)
    n = length(λ)
    S_inv = zeros(eltype(λ), length(λ), length(λ))
    S_inv = UnitLowerTriangular(S_inv)
    for i = 2:n
        for j = 1:i
            S_inv[i, j] = transformation_matrix_inverse_ij(i, j, λ)
        end
    end
    return S_inv
end



function retrieve_generation_number(i, j)
    return (i - j)
end



function transformation_matrix_5x5(λ = [1.0; 2.0; 3.0; 4.0; 5.0])

    s = zeros(5, 5)
    s[1, 1] = 1

    s[2, 1] = λ[1] / (λ[2] - λ[1])
    s[3, 1] = λ[1] / (λ[3] - λ[1]) * λ[2] / (λ[2] - λ[1])
    s[4, 1] = λ[1] / (λ[4] - λ[1]) * λ[2] / (λ[2] - λ[1]) * λ[3] / (λ[3] - λ[1])
    s[5, 1] =
        λ[1] / (λ[5] - λ[1]) * λ[2] / (λ[2] - λ[1]) * λ[3] / (λ[3] - λ[1]) *
        λ[4] / (λ[4] - λ[1])

    s[2, 2] = 1
    s[3, 2] = λ[2] / (λ[3] - λ[2])
    s[4, 2] = λ[2] / (λ[4] - λ[2]) * λ[3] / (λ[3] - λ[2])
    s[5, 2] = λ[2] / (λ[5] - λ[2]) * λ[3] / (λ[3] - λ[2]) * λ[4] / (λ[4] - λ[2])

    s[3, 3] = 1
    s[4, 3] = λ[3] / (λ[4] - λ[3])
    s[5, 3] = λ[3] / (λ[5] - λ[3]) * λ[4] / (λ[4] - λ[3])

    s[4, 4] = 1
    s[5, 4] = λ[4] / (λ[5] - λ[4])

    s[5, 5] = 1

    return s
end



function transformation_matrix_inverse_5x5(λ = [1.0; 2.0; 3.0; 4.0; 5.0])

    s = zeros(5, 5)
    s[1, 1] = 1

    s[2, 1] = λ[1] / (λ[1] - λ[2])
    s[3, 1] = λ[1] / (λ[1] - λ[3]) * λ[2] / (λ[2] - λ[3])
    s[4, 1] = λ[1] / (λ[1] - λ[4]) * λ[2] / (λ[2] - λ[4]) * λ[3] / (λ[3] - λ[4])
    s[5, 1] =
        λ[1] / (λ[1] - λ[5]) * λ[2] / (λ[2] - λ[5]) * λ[3] / (λ[3] - λ[5]) *
        λ[4] / (λ[4] - λ[5])

    s[2, 2] = 1
    s[3, 2] = λ[2] / (λ[2] - λ[3])
    s[4, 2] = λ[2] / (λ[2] - λ[4]) * λ[3] / (λ[3] - λ[4])
    s[5, 2] = λ[2] / (λ[2] - λ[5]) * λ[3] / (λ[3] - λ[5]) * λ[4] / (λ[4] - λ[5])

    s[3, 3] = 1
    s[4, 3] = λ[3] / (λ[3] - λ[4])
    s[5, 3] = λ[3] / (λ[3] - λ[5]) * λ[4] / (λ[4] - λ[5])

    s[4, 4] = 1
    s[5, 4] = λ[4] / (λ[4] - λ[5])

    s[5, 5] = 1

    return s
end
