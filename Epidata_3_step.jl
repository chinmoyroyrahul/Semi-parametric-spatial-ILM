#### Three step epidemic data
## Generated by: Chinmoy Roy Rahul
## date: APR 19, 2024

### Using built in packages

using StatsBase # to use statistical tools, like sample command
using Distributions # to generate random number 
using DataFrames
using Random
import Random


#### FUnction for generating Epidemic Data for Three step constant piecewise kernel with fixed delta


function epigen(population,size_x,size_y,initial_infection,epi_time,delta_1,delta_2,alpha_1,alpha_2,alpha_3)
    data = DataFrame()
    Id = 1: population
    data.Id = Id
    x = size_x*(1 .- rand(population)) ## X coordinate values
    y = size_y*(1 .- rand(population)) ## Y Coordinate values
    data.x = x
    data.y = y
    data.infec_time = zeros(size(x)[1])
    infected = sample(Id,initial_infection,replace = false) ## initial infection sample
    data.infec_time[infected] .= 1  # setting initial infection 1


    t = 1 ## Time 1
    distance = zeros(population, population)

    while t<epi_time

        infected_in_a_day = []
        for i in 1:size(Id)[1]

            k = []
            if !(i in infected) 
                distance[i,infected] = sqrt.((data.x[i] .- data.x[infected]).^2 + (data.y[i] .- data.y[infected]).^2)


                for j in 1:size(distance)[1]

                        if ((distance[i,j] .> 0) .& (distance[i,j] .< delta_1))  
                            alpha_1 = alpha_1
                            push!(k, alpha_1)
                        elseif ((distance[i,j] .>= delta_1) .& (distance[i,j] .< delta_2))  
                            alpha_2 = alpha_2
                            push!(k, alpha_2)
                        elseif (distance[i,j] .>= delta_2) 
                            alpha_3 = alpha_3
                            push!(k,alpha_3)
                        end

                end

                sum_k = sum(k)
                prob = 1 .- exp(-sum_k)
                u = rand(1)
                if any(prob .> u)
                    data.infec_time[i] = t+1
                    push!(infected_in_a_day,i)
                end

            end
        end
        infected = append!(infected,infected_in_a_day)

        t = t+1

    end
    return data
end 



### Example of a simulated data with parameter values

Random.seed!(1190) ## Random number 
print("Building three step parametric fixed delta: RN: 1190: ") ### 
@time epidata = epigen(400,10,10,1,20,1.5,3.0,0.08,0.01,0.0003) ## @time : will show the time of generating the data



