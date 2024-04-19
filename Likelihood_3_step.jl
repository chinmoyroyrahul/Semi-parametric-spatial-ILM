
### Likelihood code for three-step constant piecewise kernel with fixed delta

### packages
using StatsBase # to use statistical tools, like sample command
using Distributions # to generate random number 
using DataFrames

Random.seed!()

### Column used from the data
infec_time = epidata[:,:infec_time]
x = epidata[:,:x]
y = epidata[:,:y]
Id = epidata[:,:Id]
max_time= maximum(infec_time) ## Maximum time in infection time column


#### DIstance matrix based on the X, Y coordinates

dist_matrix_1 = zeros(Float64, size(Id)[1], size(Id)[1])
for i in 1:size(Id)[1]
    for j in 1:size(Id)[1]
        dis = sqrt.((x[i] .- x[j]).^2 .+ (y[i] .- y[j]).^2)
        dist_matrix_1[i,j] = dis
    end
end


#### Likelihood function for three-step

function epiloglike3(theta, delta)
    ### likelihood calculation
    
    likelihood = []
    ####
    for T in 1:(max_time-1)
        infectious_Id = Id[(infec_time.<= T) .& (infec_time .!= 0)]
        infected_Id = Id[infec_time .== T+1]
        sus_Id = Id[(infec_time .> T+1) .| (infec_time .== 0)]
        inf_sum_k_ij = []
        for i in (infected_Id)
            k = []
            sum_k_ij = []
            for j in infectious_Id
                if (dist_matrix_1[i,j] .<delta[1])
                    alpha_1 = theta[1]
                    push!(k,alpha_1)
                elseif ((dist_matrix_1[i,j] .>=delta[1]) .& (dist_matrix_1[i,j] .<delta[2]))
                    alpha_2 = theta[2]
                    push!(k, alpha_2)
                elseif dist_matrix_1[i,j] .>= delta[2]
                    alpha_3 = theta[3]
                    push!(k,alpha_3)
                end
                
            end
            sum_k_ij = sum(k)
            push!(inf_sum_k_ij, sum_k_ij)
        end

        infec_probability = [] ## For infection probability
        if !(T+1 in infec_time)
            prob = 1
            push!(infec_probability, prob)
        else
            prob = 1 .- exp.(- inf_sum_k_ij)
            append!(infec_probability, prob)
        end

        sus_sum_sus_k_ij = [] ## For susceptible individual
        for k in (sus_Id)
            sus_k = []
            sum_sus_k_ij = []
            for l in infectious_Id
                if (dist_matrix_1[k,l] .<delta[1])
                    alpha_1 = theta[1]
                    push!(sus_k,alpha_1)
                elseif ((dist_matrix_1[k,l] .>= delta[1]) .& (dist_matrix_1[k,l] .< delta[2]))
                    alpha_2 = theta[2]
                    push!(sus_k, alpha_2)
                elseif dist_matrix_1[k,l] .>= delta[2]
                    alpha_3 = theta[3]
                    push!(sus_k,alpha_3)
                end
            end
            sum_sus_k_ij = sum(sus_k)
            push!(sus_sum_sus_k_ij,sum_sus_k_ij)
        end

        sus_probability = []
        if (sus_Id == [])
            sus_prob = 0
            push!(sus_probability, sus_prob)
        else
            sus_prob = 1 .- exp.(-sus_sum_sus_k_ij)
            append!(sus_probability, sus_prob)
        end

        one_minus_sus_prob = 1 .- sus_probability
        log_one_minus_sus_prob = log.(one_minus_sus_prob)
        sum_log_one_minus_sus_prob = sum(log_one_minus_sus_prob)
        log_infec_probability = log.(infec_probability)
        sum_log_infec_probability = sum(log_infec_probability)
        indi_likelihood = sum_log_infec_probability .+ sum_log_one_minus_sus_prob
        push!(likelihood,indi_likelihood)

    end

return loglike = sum(likelihood)

end

