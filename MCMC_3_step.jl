## Two step epidemic
## Update simulation data geration code, more time efficient

using StatsBase # to use statistical tools, like sample command
using Distributions # to generate random number 
using DataFrames
using Random
import Random

# Random.seed!(1234)


function epigen(population,size_x,size_y,initial_infection,epi_time,delta_1,alpha_1,alpha_2)
    data = DataFrame()
#     population = 400
    Id = 1: population
    data.Id = Id
#     data
    x = size_x*(1 .- rand(population))
    # println("x: ",x)
    y = size_y*(1 .- rand(population))
    # println("y: ", y)
    data.x = x
    data.y = y
#     data
    data.infec_time = zeros(size(x)[1])
    # data.remove_time = zeros(size(x)[1])
#     data

    # period = 3
    infected = sample(Id,initial_infection,replace = false)
    # println("initial infection at starting time: ", infected)

    data.infec_time[infected] .= 1
    # data.remove_time[infected] .= 1 + period


    t = 1
    distance = zeros(population, population)
    # infection_mat = zeros(population,2)# time and infected id
    # infection_mat[infected,:] = [infected data.infec_time[infected]]
    # infection_mat = convert(DataFrame,infection_mat) 
    # rename!(infection_mat,[:infection,:time])
    # removed = []
#     data

    while t<epi_time
    #     println("########################")
    #     println("Time: ", t)
        infected_in_a_day = []
        for i in 1:size(Id)[1]
    #         println("------id: ", i)
            k = []
            #.& infected != 0
            if !(i in infected) 
                distance[i,infected] = sqrt.((data.x[i] .- data.x[infected]).^2 + (data.y[i] .- data.y[infected]).^2)
    #             println("distance[i,infected]: ",distance[i,infected])

                for j in 1:size(distance)[1]
    #                     println("j: ", j)
    #                 if !(j in removed)
                        if ((distance[i,j] .> 0) .& (distance[i,j] .< delta_1))  
                            alpha_1 = alpha_1
                            push!(k, alpha_1)
    #                         println("k_", k)
                        elseif (distance[i,j] .>= delta_1) 
                            alpha_2 = alpha_2
                            push!(k,alpha_2)
    #                         println("k_", k)
                        end
    #                 end

                end
    #             println("k: ", k)
                sum_k = sum(k)
    #             println("sum k: ", sum_k)
                prob = 1 .- exp(-sum_k)
                u = rand(1)
    #             println("u: ", u)
                if any(prob .> u)
                    data.infec_time[i] = t+1
    #                 data.remove_time[i] = t + 1 + period
                    push!(infected_in_a_day,i)
    #                 println("infected_in_a_day: ",infected_in_a_day)
    #                 infection_mat[i,:]=[i t+1]
                end

            end
        end
        infected = append!(infected,infected_in_a_day)
    #     println("infected: ", infected)

        t = t+1



    end
    return data
end

Random.seed!(1190)
print("Building two step parametric fixed delta: RN: 1190: ")
@time epidata = epigen(400,10,10,1,20,2,0.10,0.0004)

Random.seed!()
## two step

### 2 step likelihood 



### Updated likelihood code
using StatsBase # to use statistical tools, like sample command
using Distributions # to generate random number 
using DataFrames

infec_time = epidata[:,:infec_time]
x = epidata[:,:x]
# println("x = ", x[3])
y = epidata[:,:y]
Id = epidata[:,:Id]
max_time= maximum(infec_time)

# dist_matrix = zeros((size(Id)[1], size(Id)[1]))
# println("dist_mat ", dist_matrix)
# println(typeof(dist_matrix))

#### DIstance matrix

dist_matrix_1 = zeros(Float64, size(Id)[1], size(Id)[1])
for i in 1:size(Id)[1]
#     println("i = ", i)
#     distance = dist_matrix_1
#     println("distance = ",distance)
    for j in 1:size(Id)[1]
#         println("j ", j)
        dis = sqrt.((x[i] .- x[j]).^2 .+ (y[i] .- y[j]).^2)
        dist_matrix_1[i,j] = dis
#         println("dis = ", dis)
    end
#     push!(distance, dis[i])
end

function epiloglike2(theta, delta)
    

    ### likelihood calculation
    
    likelihood = []
    ####
    for T in 1:(max_time-1)
#     for T = 3
#         println("--- T:", T)
#         T = 2
        infectious_Id = Id[(infec_time.<= T) .& (infec_time .!= 0)]
#         println("infectious_Id: ", infectious_Id)
        infected_Id = Id[infec_time .== T+1]
#         println("infected_Id: ", infected_Id)
        sus_Id = Id[(infec_time .> T+1) .| (infec_time .== 0)]
#         println(" Sus_ID: ", sus_Id)
#         likelihood = []
#         k = []
#         sum_k_ij = []
        inf_sum_k_ij = []
        for i in (infected_Id)
            k = []
            sum_k_ij = []
#             println("----infected ID : ", i)
            for j in infectious_Id
#                 println("---infectious ID : ", j)
#                 println("matrix: ", dist_matrix_1[i,j])
                if (dist_matrix_1[i,j] .<delta[1])
                    alpha_1 = theta[1]
                    push!(k,alpha_1)
                elseif dist_matrix_1[i,j] .>= delta[1]
                    alpha_2 = theta[2]
                    push!(k,alpha_2)
                end
                
            end
#             println("k: ", k)
            sum_k_ij = sum(k)
#             println("sum k: ", sum_k_ij)
            push!(inf_sum_k_ij, sum_k_ij)
        end

        infec_probability = []
        if !(T+1 in infec_time)
            prob = 1
            push!(infec_probability, prob)
        else
#             prob = 1 .- exp.(- sum_k_ij)
            prob = 1 .- exp.(- inf_sum_k_ij)
#             prob = 1 .- exp.(- k)
            append!(infec_probability, prob)
        end

#         println("infec_probability: ", infec_probability)


        sus_sum_sus_k_ij = []
#         sus_k = []
#         sum_sus_k_ij = []
        for k in (sus_Id)
            sus_k = []
            sum_sus_k_ij = []
#             println("---sus id : ", k)
            for l in infectious_Id
#                 println("--infectious ID : ", l)
#                 println("matrix: ", dist_matrix_1[k,l])
                if (dist_matrix_1[k,l] .<delta[1])
                    alpha_1 = theta[1]
                    push!(sus_k,alpha_1)
                elseif dist_matrix_1[k,l] .>= delta[1]
                    alpha_2 = theta[2]
                    push!(sus_k,alpha_2)
                end
#                 sum_sus_k_ij = sum(sus_k)
#                 push!(sum_sus_k_ij, sum_sus_k)
#                 println("sus_k: ", sus_k)
            end
            sum_sus_k_ij = sum(sus_k)
#             println("sum sus_k: ", sum_sus_k_ij)
            push!(sus_sum_sus_k_ij,sum_sus_k_ij)

        end
#         println("sus_k: ", sus_k)
#         println("sum sus_k: ", sum_sus_k_ij)

        sus_probability = []
        if (sus_Id == [])
            sus_prob = 0
            push!(sus_probability, sus_prob)
        else
#             sus_prob = 1 .- exp.(-sum_sus_k_ij)
            sus_prob = 1 .- exp.(-sus_sum_sus_k_ij)
#             sus_prob = 1 .- exp.(-sus_k)
            append!(sus_probability, sus_prob)
        end
#         println("sus_probability: ", sus_probability)
        one_minus_sus_prob = 1 .- sus_probability
#         println("one_minus_sus_prob: ", one_minus_sus_prob)
        log_one_minus_sus_prob = log.(one_minus_sus_prob)
#         println("log_one_minus_sus_prob: ", log_one_minus_sus_prob)
        sum_log_one_minus_sus_prob = sum(log_one_minus_sus_prob)
#         println("sum_log_one_minus_sus_prob: ", sum_log_one_minus_sus_prob)
        log_infec_probability = log.(infec_probability)
#         println("log_infec_probability: ",log_infec_probability)
        sum_log_infec_probability = sum(log_infec_probability)
#         println("sum_log_infec_probability: ", sum_log_infec_probability)
        indi_likelihood = sum_log_infec_probability .+ sum_log_one_minus_sus_prob
#         println("indi_like ", T,": " , indi_likelihood )
        push!(likelihood,indi_likelihood)
#         println("likelihood: ", likelihood)

    end

return loglike = sum(likelihood)

end

###MCMC
#### 2 step mcmc (nested function) for fixed delta

### single parameter update, efficient code

## MCMC simulation without constraint ###working # need to show # perfect
using Distributions
using Random
using StatsPlots
using MCMCChains
import Random
# Random.seed!(1234)


# function posterior(alpha_1,alpha_2)
#     prior_alpha_1_distr = truncated(Normal(0, 100000), 0, Inf) # half normal distribution with mean, 
#     #variance, lower and upper values
# #     prior_alpha_1_distr = Uniform(0,10) # continuous uniform(a,b) # here I use delta_1 and delta_2 values
# #     rand(truncated(Normal(0, 0.15), 0,Inf))
#     prior_alpha_1(p) = pdf(prior_alpha_1_distr,p)

# #     prior_alpha_2_distr = Uniform(0,10)
#     prior_alpha_2_distr = truncated(Normal(0, 100000), 0, Inf)
#     prior_alpha_2(p) = pdf(prior_alpha_2_distr,p)
    
#     out = epiloglike([alpha_1,alpha_2],[1.5]) + log(prior_alpha_1(alpha_1)) + log(prior_alpha_2(alpha_2)) 

#     return out
# end

function mhrw(N,burnedinsamples, delta, alpha_1_init,alpha_2_init, alpha_1_propsd,alpha_2_propsd)
    
    function posterior(alpha_1,alpha_2,delta)
        prior_alpha_1_distr = truncated(Normal(0, 100000), 0, Inf) # half normal distribution with mean, 
        #variance, lower and upper values
    #     prior_alpha_1_distr = Uniform(0,10) # continuous uniform(a,b) # here I use delta_1 and delta_2 values
    #     rand(truncated(Normal(0, 0.15), 0,Inf))
        prior_alpha_1(p) = pdf(prior_alpha_1_distr,p)

    #     prior_alpha_2_distr = Uniform(0,10)
        prior_alpha_2_distr = truncated(Normal(0, 100000), 0, Inf)
        prior_alpha_2(p) = pdf(prior_alpha_2_distr,p)

        out = epiloglike2([alpha_1,alpha_2],[delta]) + log(prior_alpha_1(alpha_1)) + log(prior_alpha_2(alpha_2)) 

        return out
    end


    alpha_1_out =  zeros(N)
    alpha_2_out =  zeros(N)


    alpha_1_out[1] = alpha_1_init
    alpha_2_out[1] = alpha_2_init


    posterior_mem_a1 = zeros(N)
    posterior_mem_a1[1] = posterior(alpha_1_out[1],alpha_2_out[1],delta)
#     println("posterior_mem_alpha_1[1]: ", posterior_mem_a1[1])
    
    posterior_mem_a2 = zeros(N)
    posterior_mem_a2[1] = posterior(alpha_1_out[1],alpha_2_out[1],delta)
#     println("posterior_mem_alpha_2[1]: ", posterior_mem_a2[1])


    acceptance1 = 0
    acceptance2 = 0

    #
    for i in 2:N
        
        while i <=N
            
            ## block 1
            alpha_1_prop = (rand(Normal(alpha_1_out[i-1], alpha_1_propsd),1))[1]
#             println("alpha_1_prop: ", alpha_1_prop)
            if alpha_1_prop .< 0
#                 posterior_mem_a1[i] = posterior_mem_a1[i-1]
                posterior_mem_a1[i] = posterior_mem_a2[i-1]
#                 posterior_mem[i] = posterior(alpha_1_out[i-1],alpha_2_out[i-1])
#                 println("posterior_mem_a1[i]: ", posterior_mem_a1[i])
                alpha_1_out[i] = alpha_1_out[i-1]
#                 println("alpha_1_out[i]: ", alpha_1_out[i])
                alpha_2_out[i] = alpha_2_out[i-1]
#                 println("alpha_2_out[i]: ", alpha_2_out[i])
#                 break
            elseif alpha_1_prop .>=0
#                 posterior_mem[i] = posterior(alpha_1_prop,alpha_2_out[i-1])
                posterior_mem_a1[i] = posterior(alpha_1_prop,alpha_2_out[i-1],delta)
                aceptprob = posterior_mem_a1[i] - posterior_mem_a2[i-1]
                u = log(rand())
#                 println("u: ", u)
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_prop
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    posterior_mem_a1[i] = posterior_mem_a1[i]
#                     println("posterior_mem_a1[i]: ",posterior_mem_a1[i])
                    acceptance1 = acceptance1 + 1
                    alpha_2_out[i] = alpha_2_out[i-1]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
#                     break
                else
                    alpha_1_out[i] = alpha_1_out[i-1]
                    posterior_mem_a1[i] = posterior_mem_a2[i-1]
#                     println("posterior_mem_a1[i]: ",posterior_mem_a1[i])
                    alpha_2_out[i] = alpha_2_out[i-1]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
#                     break
                end
                
            end
           
            ## block 2
            
            alpha_2_prop = (rand(Normal(alpha_2_out[i-1], alpha_2_propsd),1))[1]
#             println("alpha_2_prop: ", alpha_2_prop)
            if alpha_2_prop .< 0
#                 posterior_mem[i] = posterior(alpha_1_out[i],alpha_2_out[i-1])
                posterior_mem_a2[i] = posterior_mem_a1[i]
                alpha_1_out[i] = alpha_1_out[i]
#                 println("alpha_1_out[i]: ", alpha_1_out[i])
                alpha_2_out[i] = alpha_2_out[i-1]
#                 println("alpha_2_out[i]: ", alpha_2_out[i])
                break
            elseif alpha_2_prop .>= 0
                posterior_mem_a2[i] = posterior(alpha_1_out[i],alpha_2_prop,delta)
                aceptprob = posterior_mem_a2[i] - posterior_mem_a1[i]
#                 aceptprob = posterior_mem_a2[i] - posterior_mem_a2[i-1]
#                 println("aceptprob: ", aceptprob)
                u = log(rand())
#                 println("u: ", u)
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    alpha_2_out[i] = alpha_2_prop
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    posterior_mem_a2[i] = posterior_mem_a2[i]
#                     println("posterior_mem_a2[i]: ", posterior_mem_a2[i])
                    acceptance2 = acceptance2 + 1
                    break
                else
                    alpha_1_out[i] = alpha_1_out[i]
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    alpha_2_out[i] = alpha_2_out[i-1]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    posterior_mem_a2[i] = posterior_mem_a1[i]
#                     posterior_mem_a2[i] = posterior_mem_a2[i-1]
#                     println("posterior_mem_a2[i]: ", posterior_mem_a2[i])
                    break
                end
            end
#             break
            
            
            
            
        end


# #             i +=1
#         end

    end
    alpha_1_final = deleteat!(alpha_1_out,1:burnedinsamples)
    alpha_2_final = deleteat!(alpha_2_out,1:burnedinsamples)
    
#     println("Acceptance rate = ", acceptance/N)
    println("Acceptance rate = ", ((acceptance1/N)+(acceptance2/N))/2)
    return [alpha_1_final,alpha_2_final]

end

@time mcmc15 = mhrw(30000,0,2,0.08,0.001,0.03612,0.000383)# working

y_15 = mcmc15

println("for delta 2: y_15: ", y_15)

# c15 = (y_15[1])[setdiff(1:end, 1: 5000, )]
# c215 = (y_15[2])[setdiff(1:end, 1: 5000, )]

# ### AIC and DIC
# ##AIC & DIC
# ## AIC = 2k - 2 log-likelihood # where k is the number of estimated parameters in the model

# # parameter

# theta_4 = [mean(c15);mean(c215)]

# ## k

# k_4 = size(theta_4)[1]

# ## Log-likelihood

# l_4 =  epiloglike2(theta_4,[10])

# ## AIC
# AIC_4 = (-2 * l_4) + (2*k_4)

# ## DIC

# ## DIC for parametric

# d_24 = [] # empty set of deviance

# for i in 1:size(c15)[1]
#     d_24i = -2 .* epiloglike2([c15[i],c215[i]],[10])
#     push!(d_24,d_24i)
# end

# d_bar_24 = mean(d_24) # overall deviance mean for two steps
# d_est_24 = -2 .* epiloglike2([mean(c15),mean(c215)],[10]) ## D(theta bar)

# # ## effective number of parameters of the model

# pd_24 = d_bar_24 - d_est_24 


# ## DIC calculation

# DIC_24 = d_est_24 + (2* pd_24)
# DIC_24
# println("AIC_15:  ",AIC_4, "  ","DIC_15:  ", DIC_24)



#### Three step epidemic

## Update simulation data geration code, more time efficient

using StatsBase # to use statistical tools, like sample command
using Distributions # to generate random number 
using DataFrames
using Random
import Random

# Random.seed!(1234)


function epigen(population,size_x,size_y,initial_infection,epi_time,delta_1,delta_2,alpha_1,alpha_2,alpha_3)
    data = DataFrame()
#     population = 400
    Id = 1: population
    data.Id = Id
#     data
    x = size_x*(1 .- rand(population))
    # println("x: ",x)
    y = size_y*(1 .- rand(population))
    # println("y: ", y)
    data.x = x
    data.y = y
#     data
    data.infec_time = zeros(size(x)[1])
    # data.remove_time = zeros(size(x)[1])
#     data

    # period = 3
    infected = sample(Id,initial_infection,replace = false)
    # println("initial infection at starting time: ", infected)

    data.infec_time[infected] .= 1
    # data.remove_time[infected] .= 1 + period


    t = 1
    distance = zeros(population, population)
    # infection_mat = zeros(population,2)# time and infected id
    # infection_mat[infected,:] = [infected data.infec_time[infected]]
    # infection_mat = convert(DataFrame,infection_mat) 
    # rename!(infection_mat,[:infection,:time])
    # removed = []
#     data

    while t<epi_time
    #     println("########################")
    #     println("Time: ", t)
        infected_in_a_day = []
        for i in 1:size(Id)[1]
    #         println("------id: ", i)
            k = []
            #.& infected != 0
            if !(i in infected) 
                distance[i,infected] = sqrt.((data.x[i] .- data.x[infected]).^2 + (data.y[i] .- data.y[infected]).^2)
    #             println("distance[i,infected]: ",distance[i,infected])

                for j in 1:size(distance)[1]
    #                     println("j: ", j)
    #                 if !(j in removed)
                        if ((distance[i,j] .> 0) .& (distance[i,j] .< delta_1))  
                            alpha_1 = alpha_1
                            push!(k, alpha_1)
    #                         println("k_", k)
                        elseif ((distance[i,j] .>= delta_1) .& (distance[i,j] .< delta_2))  
                            alpha_2 = alpha_2
                            push!(k, alpha_2)
                        elseif (distance[i,j] .>= delta_2) 
                            alpha_3 = alpha_3
                            push!(k,alpha_3)
    #                         println("k_", k)
                        end
    #                 end

                end
    #             println("k: ", k)
                sum_k = sum(k)
    #             println("sum k: ", sum_k)
                prob = 1 .- exp(-sum_k)
                u = rand(1)
    #             println("u: ", u)
                if any(prob .> u)
                    data.infec_time[i] = t+1
    #                 data.remove_time[i] = t + 1 + period
                    push!(infected_in_a_day,i)
    #                 println("infected_in_a_day: ",infected_in_a_day)
    #                 infection_mat[i,:]=[i t+1]
                end

            end
        end
        infected = append!(infected,infected_in_a_day)
    #     println("infected: ", infected)

        t = t+1



    end
    return data
end 

Random.seed!(1190)
print("Building three step parametric fixed delta: RN: 1190: ")
@time epidata = epigen(400,10,10,1,20,1.5,3.0,0.08,0.01,0.0003)

Random.seed!()

### Likelihood code

## new likelihood code
# new likelihood code
### Updated code

using StatsBase # to use statistical tools, like sample command
using Distributions # to generate random number 
using DataFrames

infec_time = epidata[:,:infec_time]
x = epidata[:,:x]
# println("x = ", x[3])
y = epidata[:,:y]
Id = epidata[:,:Id]
max_time= maximum(infec_time)


#### DIstance matrix

dist_matrix_1 = zeros(Float64, size(Id)[1], size(Id)[1])
for i in 1:size(Id)[1]
#     println("i = ", i)
#     distance = dist_matrix_1
#     println("distance = ",distance)
    for j in 1:size(Id)[1]
#         println("j ", j)
        dis = sqrt.((x[i] .- x[j]).^2 .+ (y[i] .- y[j]).^2)
        dist_matrix_1[i,j] = dis
#         println("dis = ", dis)
    end
#     push!(distance, dis[i])
end


function epiloglike3(theta, delta)
    

    ### likelihood calculation
    
    likelihood = []
    ####
    for T in 1:(max_time-1)
#     for T = 3
#         println("--- T:", T)
#         T = 2
        infectious_Id = Id[(infec_time.<= T) .& (infec_time .!= 0)]
#         println("infectious_Id: ", infectious_Id)
        infected_Id = Id[infec_time .== T+1]
#         println("infected_Id: ", infected_Id)
        sus_Id = Id[(infec_time .> T+1) .| (infec_time .== 0)]
#         println(" Sus_ID: ", sus_Id)
#         likelihood = []
#         k = []
#         sum_k_ij = []
        inf_sum_k_ij = []
        for i in (infected_Id)
            k = []
            sum_k_ij = []
#             println("----infected ID : ", i)
            for j in infectious_Id
#                 println("---infectious ID : ", j)
#                 println("matrix: ", dist_matrix_1[i,j])
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
#             println("k: ", k)
            sum_k_ij = sum(k)
#             println("sum k: ", sum_k_ij)
            push!(inf_sum_k_ij, sum_k_ij)
        end

        infec_probability = []
        if !(T+1 in infec_time)
            prob = 1
            push!(infec_probability, prob)
        else
#             prob = 1 .- exp.(- sum_k_ij)
            prob = 1 .- exp.(- inf_sum_k_ij)
#             prob = 1 .- exp.(- k)
            append!(infec_probability, prob)
        end

#         println("infec_probability: ", infec_probability)


        sus_sum_sus_k_ij = []
#         sus_k = []
#         sum_sus_k_ij = []
        for k in (sus_Id)
            sus_k = []
            sum_sus_k_ij = []
#             println("---sus id : ", k)
            for l in infectious_Id
#                 println("--infectious ID : ", l)
#                 println("matrix: ", dist_matrix_1[k,l])
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
#                 sum_sus_k_ij = sum(sus_k)
#                 push!(sum_sus_k_ij, sum_sus_k)
#                 println("sus_k: ", sus_k)
            end
            sum_sus_k_ij = sum(sus_k)
#             println("sum sus_k: ", sum_sus_k_ij)
            push!(sus_sum_sus_k_ij,sum_sus_k_ij)
        end

        sus_probability = []
        if (sus_Id == [])
            sus_prob = 0
            push!(sus_probability, sus_prob)
        else
#             sus_prob = 1 .- exp.(-sum_sus_k_ij)
            sus_prob = 1 .- exp.(-sus_sum_sus_k_ij)
#             sus_prob = 1 .- exp.(-sus_k)
            append!(sus_probability, sus_prob)
        end
#         println("sus_probability: ", sus_probability)
        one_minus_sus_prob = 1 .- sus_probability
#         println("one_minus_sus_prob: ", one_minus_sus_prob)
        log_one_minus_sus_prob = log.(one_minus_sus_prob)
#         println("log_one_minus_sus_prob: ", log_one_minus_sus_prob)
        sum_log_one_minus_sus_prob = sum(log_one_minus_sus_prob)
#         println("sum_log_one_minus_sus_prob: ", sum_log_one_minus_sus_prob)
        log_infec_probability = log.(infec_probability)
#         println("log_infec_probability: ",log_infec_probability)
        sum_log_infec_probability = sum(log_infec_probability)
#         println("sum_log_infec_probability: ", sum_log_infec_probability)
        indi_likelihood = sum_log_infec_probability .+ sum_log_one_minus_sus_prob
#         println("indi_like ", T,": " , indi_likelihood )
        push!(likelihood,indi_likelihood)
#         println("likelihood: ", likelihood)

    end

return loglike = sum(likelihood)

end

# single parameter update # new 
using Distributions
using Random
using StatsPlots
using MCMCChains
import Random
# Random.seed!(1234)



function mhrw(N,burnedinsamples,delta_1,delta_2,alpha_1_init,alpha_2_init,alpha_3_init,alpha_1_propsd,alpha_2_propsd,
                alpha_3_propsd)
    
    function posterior(alpha_1,alpha_2,alpha_3,delta_1,delta_2)
        
#         truncated(Normal(0, 100000), 0, Inf)

        prior_alpha_1_distr = truncated(Normal(0, 100000), 0, Inf) 
        prior_alpha_1(p) = pdf(prior_alpha_1_distr,p)
        prior_alpha_2_distr = truncated(Normal(0, 100000), 0, Inf)
        prior_alpha_2(p) = pdf(prior_alpha_2_distr,p)
        prior_alpha_3_distr = truncated(Normal(0, 100000), 0, Inf)
        prior_alpha_3(p) = pdf(prior_alpha_3_distr,p)


        out = epiloglike3([alpha_1,alpha_2,alpha_3],[delta_1,delta_2]) + log(prior_alpha_1(alpha_1)) +
                log(prior_alpha_2(alpha_2)) + log(prior_alpha_3(alpha_3))

    #    return exp(out)
        return out
    #    println("-----")
    #    println("posterior out", (out))
    end

    alpha_1_out =  zeros(N)
    alpha_2_out =  zeros(N)
    alpha_3_out =  zeros(N)
    
    alpha_1_out[1] = alpha_1_init
    alpha_2_out[1] = alpha_2_init
    alpha_3_out[1] = alpha_3_init

#     posterior_mem = zeros(N)
#     posterior_mem[1] = posterior(alpha_1_out[1],alpha_2_out[1],delta_1_out[1])
#     println("posterior_mem[1]: ", posterior_mem[1])
    posterior_mem_a1 = zeros(N)
    posterior_mem_a1[1] = posterior(alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1,delta_2)
#     println("posterior_mem_alpha_1[1]: ", posterior_mem_a1[1])
    
    posterior_mem_a2 = zeros(N)
    posterior_mem_a2[1] = posterior(alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1,delta_2)
#     println("posterior_mem_alpha_2[1]: ", posterior_mem_a2[1])
    
    posterior_mem_a3 = zeros(N)
    posterior_mem_a3[1] = posterior(alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1,delta_2)
#     println("posterior_mem_alpha_a3[1]: ", posterior_mem_a3[1])


    acceptance1 = 0
    acceptance2 = 0
    acceptance3 = 0

    #
    for i in 2:N
#         println("---------")
#         println("i: ", i)
#         alpha_1_prop = (rand(Normal(alpha_1_out[i-1], alpha_1_propsd),1))[1]
#         println("alpha_1_prop: ", alpha_1_prop)

#         alpha_2_prop = (rand(Normal(alpha_2_out[i-1], alpha_2_propsd),1))[1]
#         println("alpha_2_prop: ", alpha_2_prop)
        
        while i <=N
            
            ## block 1
            alpha_1_prop = (rand(Normal(alpha_1_out[i-1], alpha_1_propsd),1))[1]
            # println("alpha_1_prop: ", alpha_1_prop)
            if alpha_1_prop .< 0
                posterior_mem_a1[i] = posterior_mem_a3[i-1]
#                 posterior_mem[i] = posterior(alpha_1_out[i-1],alpha_2_out[i-1])
                # println("posterior_mem_a1[i]: ", posterior_mem_a1[i])
                alpha_1_out[i] = alpha_1_out[i-1]
#                 println("alpha_1_out[i]: ", alpha_1_out[i])
                alpha_2_out[i] = alpha_2_out[i-1]
#                 println("alpha_2_out[i]: ", alpha_2_out[i])
                alpha_3_out[i] = alpha_3_out[i-1]
#                 println("delta_1_out[i]: ", delta_1_out[i])
#                 break
            elseif alpha_1_prop .>=0
                posterior_mem_a1[i] = posterior(alpha_1_prop,alpha_2_out[i-1],alpha_3_out[i-1],delta_1,delta_2)
                # println("posterior_mem_a1[i]: ", posterior_mem_a1[i])
#                 aceptprob = posterior_mem[i] - posterior(alpha_1_out[i-1],alpha_2_out[i-1], delta_1_out[i-1])
                aceptprob = posterior_mem_a1[i] - posterior_mem_a3[i-1]
                # println("posterior_mem_a3[i-1]: ", posterior_mem_a3[i-1])
                # println("aceptprob: ", aceptprob)
                u = log(rand())
                # println("u: ", u)
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_prop
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    posterior_mem_a1[i] = posterior_mem_a1[i]
#                     println("posterior_mem_a1[i]: ", posterior_mem_a1[i])
                    acceptance1 = acceptance1 + 1
                    alpha_2_out[i] = alpha_2_out[i-1]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    alpha_3_out[i] = alpha_3_out[i-1]
#                     println("delta_1_out[i]: ", delta_1_out[i])
#                     break
                else
                    alpha_1_out[i] = alpha_1_out[i-1]
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    posterior_mem_a1[i] = posterior_mem_a3[i-1]
#                     println("posterior_mem_a1[i]: ", posterior_mem_a1[i])
                    alpha_2_out[i] = alpha_2_out[i-1]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    alpha_3_out[i] = alpha_3_out[i-1]
#                     println("delta_1_out[i]: ", delta_1_out[i])
#                     break
                end
                
            end
           
            ## block 2
            
            alpha_2_prop = (rand(Normal(alpha_2_out[i-1], alpha_2_propsd),1))[1]
#             println("alpha_2_prop: ", alpha_2_prop)
            if alpha_2_prop .< 0
                posterior_mem_a2[i] = posterior_mem_a1[i] # posterior_mem_a1[i] # 
                # println("posterior_mem_a2[i]: ", posterior_mem_a2[i])
                alpha_1_out[i] = alpha_1_out[i]
#                 println("alpha_1_out[i]: ", alpha_1_out[i])
                alpha_2_out[i] = alpha_2_out[i-1]
#                 println("alpha_2_out[i]: ", alpha_2_out[i])
                alpha_3_out[i] = alpha_3_out[i-1]
#                 println("delta_1_out[i]: ", delta_1_out[i])
#                 break
            elseif alpha_2_prop .>= 0
                posterior_mem_a2[i] = posterior(alpha_1_out[i],alpha_2_prop,alpha_3_out[i-1],delta_1,delta_2)
                # println("posterior_mem_a2[i]: ", posterior_mem_a2[i])
                aceptprob = posterior_mem_a2[i] - posterior_mem_a1[i]#posterior_mem_a1[i]
#                 println("posterior_mem_a1[i]: ", posterior_mem_a1[i])
                # println("aceptprob: ", aceptprob)
                u = log(rand())
                # println("u: ", u)
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    alpha_2_out[i] = alpha_2_prop
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    alpha_3_out[i] = alpha_3_out[i-1]
#                     println("delta_1_out[i]: ", delta_1_out[i])
                    posterior_mem_a2[i] = posterior_mem_a2[i]
                    # println("posterior_mem_a2[i]: ", posterior_mem_a2[i])
                    acceptance2 = acceptance2 + 1
#                     break
                else
                    alpha_1_out[i] = alpha_1_out[i]
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    alpha_2_out[i] = alpha_2_out[i-1]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    alpha_3_out[i] = alpha_3_out[i-1]
#                     println("delta_1_out[i]: ", delta_1_out[i])
                    posterior_mem_a2[i] = posterior_mem_a1[i]#posterior_mem_a1[i]
                    # println("posterior_mem_a2[i]: ", posterior_mem_a2[i])
#                     break
                end
            end
#             break
            
            ### block 3
            alpha_3_prop = (rand(Normal(alpha_3_out[i-1], alpha_3_propsd),1))[1]
#             println("delta_1_prop: ", delta_1_prop)
            
            if alpha_3_prop .< 0
                posterior_mem_a3[i] = posterior_mem_a2[i] # posterior_mem_a1[i] # 
                # println("posterior_mem_a3[i]: ", posterior_mem_a3[i])
                alpha_1_out[i] = alpha_1_out[i]
#                 println("alpha_1_out[i]: ", alpha_1_out[i])
                alpha_2_out[i] = alpha_2_out[i]
#                 println("alpha_2_out[i]: ", alpha_2_out[i])
                alpha_3_out[i] = alpha_3_out[i-1]
#                 println("delta_1_out[i]: ", delta_1_out[i])
                break
            
            elseif  alpha_3_prop .>= 0
                posterior_mem_a3[i] = posterior(alpha_1_out[i],alpha_2_out[i],alpha_3_prop,delta_1,delta_2)
                # println("posterior_mem_a3[i]: ", posterior_mem_a3[i])
                aceptprob = posterior_mem_a3[i] - posterior_mem_a2[i]#posterior_mem_a1[i]
#                 println("posterior_mem_a2[i]: ", posterior_mem_a2[i])
                # println("aceptprob: ", aceptprob)
                u = log(rand())
                # println("u: ", u)
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    alpha_2_out[i] = alpha_2_out[i]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    alpha_3_out[i] = alpha_3_prop
#                     println("delta_1_out[i]: ", delta_1_out[i])
                    posterior_mem_a3[i] = posterior_mem_a3[i]
                    # println("posterior_mem_a3[i]: ", posterior_mem_a3[i])
                    acceptance3 = acceptance3 + 1
                    break
                else
                    alpha_1_out[i] = alpha_1_out[i]
#                     println("alpha_1_out[i]: ", alpha_1_out[i])
                    alpha_2_out[i] = alpha_2_out[i]
#                     println("alpha_2_out[i]: ", alpha_2_out[i])
                    alpha_3_out[i] = alpha_3_out[i-1]
#                     println("delta_1_out[i]: ", delta_1_out[i])
                    posterior_mem_a3[i] = posterior_mem_a2[i]#posterior_mem_a1[i]
#                     println("posterior_mem_a3[i]: ", posterior_mem_a3[i])
                    break
                    
                end
                
            end
            
            
            
        end


# #             i +=1
#         end

    end
    alpha_1_final = deleteat!(alpha_1_out,1:burnedinsamples)
    alpha_2_final = deleteat!(alpha_2_out,1:burnedinsamples)
    alpha_3_final = deleteat!(alpha_3_out,1:burnedinsamples)
    
#     println("Acceptance rate = ", acceptance/N)
    println("Acceptance rate = ", ((acceptance1/N)+(acceptance2/N)+(acceptance3/N))/3)
    return [alpha_1_final,alpha_2_final,alpha_3_final]
end


@time mcmc715 = mhrw(30000,0,1.5,3.0,0.10,0.01,0.001,0.024,0.0028,0.00042) # working



y_715 = mcmc715
println("delta_1.5_3: y_715: ", y_715)

# c1715 = (y_715[1])[setdiff(1:end, 1: 5000, )]
# c2715 = (y_715[2])[setdiff(1:end, 1: 5000, )]
# c3715 = (y_715[3])[setdiff(1:end, 1: 5000, )]

# ##AIC & DIC, 5,12
# ## AIC = 2k - 2 log-likelihood # where k is the number of estimated parameters in the model

# # parameter

# theta_512 = [mean(c1715);mean(c2715);mean(c3715)]

# ## k

# k_512 = size(theta_512)[1]

# ## Log-likelihood

# l_512 =  epiloglike3(theta_512,[10,20])

# ## AIC
# AIC_512 = (-2 * l_512) + (2*k_512)

# ## DIC

# ## DIC for parametric

# d_3_512 = [] # empty set of deviance

# for i in 1:size(c1715)[1]
#     d_3_512i = -2 .* epiloglike3([c1715[i],c2715[i],c3715[i]],[10,20])
#     push!(d_3_512,d_3_512i)
# end

# d_bar_3_512 = mean(d_3_512) # overall deviance mean for three steps 
# d_est_3_512 = -2 .* epiloglike3([mean(c1715),mean(c2715),mean(c3715)],[10,20]) ## D(theta bar)

# # ## effective number of parameters of the model

# pd_512 = d_bar_3_512 - d_est_3_512 


# ## DIC calculation

# DIC_3_512 = d_est_3_512 + (2* pd_512)
# DIC_3_512
# println("AIC_10_20:  ",AIC_512, "  ","DIC_10_20:  ", DIC_3_512)


## results
# results



using MCMCChains
using StatsPlots
using Plots
using Statistics
   

c1 = (mcmc[1])[setdiff(1:end, 1: 300, )]
p1 = plot(c1, title = "Alpha 1", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c1, [0.025, 0.5, 0.975]))
describe(c1)
println("Standard Deviation: ", std(c1))

c2 = (mcmc[2])[setdiff(1:end, 1: 300, )]
p2 = plot(c2, title = "Alpha 2", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c2, [0.025, 0.5, 0.975]))
describe(c2)
println("Standard Deviation: ", std(c2))

c3 = (mcmc[3])[setdiff(1:end, 1: 3000, )]
p3 = plot(c3, title = "Alpha 3", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c3, [0.025, 0.5, 0.975]))
describe(c3)
println("Standard Deviation: ", std(c3))


c4 = (mcmc[4])[setdiff(1:end, 1: 100, )]
p4 = plot(c4, title = "Alpha 4", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c4, [0.025, 0.5, 0.975]))
describe(c4)
println("Standard Deviation: ", std(c4))

c5 = (mcmc[5])[setdiff(1:end, 1: 100, )]
p5 = plot(c5, title = "Alpha 5", xaxis = (font(5)), yaxis = (font(8)))


println(quantile!(c5, [0.025, 0.5, 0.975]))
describe(c5)
println("Standard Deviation: ", std(c5))

c6 = (mcmc[6])[setdiff(1:end, 1: 100, )]
p6 = plot(c6, title = "Epsilon", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c6, [0.025, 0.5, 0.975]))
describe(c6)
println("Standard Deviation: ", std(c6))

c7 = (mcmc[7])[setdiff(1:end, 1: 100, )]
p7 = plot(c7, title = "Lambda", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c7, [0.025, 0.5, 0.975]))
describe(c7)
println("Standard Deviation: ", std(c7))

###

## convergence test

## Geweke'e test
# using Pkg
# Pkg.add("MCMCDiagnostics")
using MCMCDiagnostics



##
# gewekediag(mcmc; first=0.1, last=0.5, etype=:imse)

# using Pkg
# Pkg.add("GewekeTest")

# gewekediag(mcmc[1]; 0.1, 0.5)
gewekediag(mcmc[1], first=0.1, last=0.5)



###### gelman rubin

## multiple chains with different startig points
# @time mcmc1 = mhrw(10000,0,0.08,0.001,0.03612,0.000383)
# @time mcmc2 = mhrw(10000,0,0.001,0.001,0.07,0.00038)

########
 
 
c11 = (mcmc1[1])[setdiff(1:end, 1: 100, )]
c21 = (mcmc2[1])[setdiff(1:end, 1: 100, )]
c31 = (mcmc3[1])[setdiff(1:end, 1: 100, )]
# println(typeof(c2))

# p1 = plot(c1, title = "Alpha 1", xaxis = (font(5)), yaxis = (font(8)))
p11 = plot(c11)
p12 = plot!(p11,c21,title = "Alpha 1", xaxis = (font(5)), yaxis = (font(8)))
p13 = plot!(p11,c31,title = "Alpha 1", xaxis = (font(5)), yaxis = (font(8)))


### Gelman-Rubin Method

## Chain, M= 3
## iterations after burn in is N
## calculate the within and between chain variance

## within chain variance

# variance of each chain (after throwing out the first N draws)
# variance of alpha_1 from chain 1 , chain 2 and chain 3
var_11 = var(c11)
var_21 = var(c21)
var_31 = var(c31)

M = 3 # number of chain
W = (1/M)* (var_11 + var_21 + var_31)
W
# std(c11)^2

## between chain variance

N = size(c11)[1] # length of chain after burn in

## Chain mean_1
mean_11 = mean(c11)
mean_21 = mean(c21)
mean_31 = mean(c31)

mean_array = [mean_11,mean_21,mean_31]
println("mean: ", mean_array)
mean_1 = mean([mean_11,mean_21,mean_31]) # mean of each of the M means
println("mean_overall: ", mean_1)
B = (N/(M-1)) * sum((mean_array .- mean_1).^2)

## astimated variance of alpha_1 as the weighted sum of between and within chain variance

var_alpha_1 = (1 - (1/N)) * W + (1/N) * B

# calculate the potential scale reduction factor

R = sqrt(var_alpha_1/W) # we want this number close to 1(Gelman-Rubin show that when R is greater 
# then 1.1 or 1.2, we need longer burn-in)

# using Pkg
# Pkg.add("MCMCDiagnostics")
# using MCMCDiagnostics

# mcmc = 

# # Set random seed for reproducibility
# Random.seed!(123)

# # Generate synthetic MCMC samples
# n_samples = 20
# n_params = 2
# n_chains = 3
# chains = randn(n_samples, n_params, n_chains)

# psrf = gelmandiag(chain1)

# size(chain)
# chain1 = [mcmc1[1],mcmc2[1],mcmc3[1]]
# length(chain1)
# size(chain1)
# chain2 = chain1[:,:,:]

# psrf = gelmandiag(chain2)


# m1= reshape(mcmc1[1], length(mcmc1[1]), 1)
# # mcmc5 = reshape(mcmc1, length(mcmc), 1)
# m1

# m2 = reshape(mcmc1[2], length(mcmc1[2]), 1)
# m2

# m3 = reshape(mcmc2[1], length(mcmc2[1]), 1)
# m4 = reshape(mcmc2[2], length(mcmc2[2]), 1)

# m5 = reshape(mcmc3[1], length(mcmc3[1]), 1)
# m6 = reshape(mcmc3[2], length(mcmc[2]), 1)



