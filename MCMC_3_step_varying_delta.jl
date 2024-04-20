## MCMC simulation for three step contant piecewise kernel with varying change points

using Distributions
using Random
using StatsPlots
using MCMCChains
import Random
# Random.seed!(1234)


function posterior(a,b,alpha_1,alpha_2,alpha_3,delta_1,delta_2)
    prior_alpha_1_distr = truncated(Normal(0, 100000), 0, Inf) # continuous uniform(a,b) # here I use delta_1 and delta_2 values
    prior_alpha_1(p) = pdf(prior_alpha_1_distr,p)
    prior_alpha_2_distr = truncated(Normal(0, 100000), 0, Inf)
    prior_alpha_2(p) = pdf(prior_alpha_2_distr,p)
    prior_alpha_3_distr = truncated(Normal(0, 100000), 0, Inf)
    prior_alpha_3(p) = pdf(prior_alpha_3_distr,p)
    prior_delta_1_distr = Uniform(0,a)
    prior_delta_1(p) = pdf(prior_delta_1_distr,p)
    prior_delta_2_distr = Uniform(a,b)
    prior_delta_2(p) = pdf(prior_delta_2_distr,p)
    prior_delta_diff_distr = Uniform(0,b)
    prior_delta_diff(p) = pdf(prior_delta_diff_distr,p)
    
    
    
    out = epiloglike3([alpha_1,alpha_2,alpha_3],[delta_1,delta_2]) + log(prior_alpha_1(alpha_1)) +
            log(prior_alpha_2(alpha_2)) + log(prior_alpha_3(alpha_3)) + log(prior_delta_1(delta_1))+
            log(prior_delta_2(delta_2)) + log(prior_delta_diff(delta_2-delta_1))
    
    return out
end

function mhrw(N,burnedinsamples, a,b,alpha_1_init,alpha_2_init,alpha_3_init,delta_1_init,delta_2_init,
                alpha_1_propsd,alpha_2_propsd,alpha_3_propsd,delta_1_propsd,delta_2_propsd)
    alpha_1_out =  zeros(N)
    alpha_2_out =  zeros(N)
    alpha_3_out =  zeros(N)
    delta_1_out =  zeros(N)
    delta_2_out =  zeros(N)
    
    
    alpha_1_out[1] = alpha_1_init
    alpha_2_out[1] = alpha_2_init
    alpha_3_out[1] = alpha_3_init
    delta_1_out[1] = delta_1_init
    delta_2_out[1] = delta_2_init

    posterior_mem_a1 = zeros(N)
    posterior_mem_a1[1] = posterior(a,b,alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1_out[1],delta_2_out[1])
    
    posterior_mem_a2 = zeros(N)
    posterior_mem_a2[1] = posterior(a,b,alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1_out[1],delta_2_out[1])
    
    posterior_mem_a3 = zeros(N)
    posterior_mem_a3[1] = posterior(a,b,alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1_out[1],delta_2_out[1])
    
    posterior_mem_a4 = zeros(N)
    posterior_mem_a4[1] = posterior(a,b,alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1_out[1],delta_2_out[1])

    posterior_mem_a5 = zeros(N)
    posterior_mem_a5[1] = posterior(a,b,alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1_out[1],delta_2_out[1])


    acceptance1 = 0
    acceptance2 = 0
    acceptance3 = 0
    acceptance4 = 0
    acceptance5 = 0

    #
    for i in 2:N

        
        while i <=N
            
            ## block 1
            alpha_1_prop = (rand(Normal(alpha_1_out[i-1], alpha_1_propsd),1))[1]
            if alpha_1_prop .< 0
                posterior_mem_a1[i] = posterior_mem_a5[i-1]

                alpha_1_out[i] = alpha_1_out[i-1]
                alpha_2_out[i] = alpha_2_out[i-1]
                alpha_3_out[i] = alpha_3_out[i-1]
                delta_1_out[i] = delta_1_out[i-1]
                delta_2_out[i] = delta_2_out[i-1]

            elseif alpha_1_prop .>=0
                posterior_mem_a1[i] = posterior(a,b,alpha_1_prop,alpha_2_out[i-1],alpha_3_out[i-1],
                    delta_1_out[i-1],delta_2_out[i-1])

                aceptprob = posterior_mem_a1[i] - posterior_mem_a5[i-1]

                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_prop
                    posterior_mem_a1[i] = posterior_mem_a1[i]

                    acceptance1 = acceptance1 + 1
                    alpha_2_out[i] = alpha_2_out[i-1]
                    alpha_3_out[i] = alpha_3_out[i-1]
                    delta_1_out[i] = delta_1_out[i-1]
                    delta_2_out[i] = delta_2_out[i-1]

                else
                    alpha_1_out[i] = alpha_1_out[i-1]
                    posterior_mem_a1[i] = posterior_mem_a5[i-1]

                    alpha_2_out[i] = alpha_2_out[i-1]
                    alpha_3_out[i] = alpha_3_out[i-1]
                    delta_1_out[i] = delta_1_out[i-1]
                    delta_2_out[i] = delta_2_out[i-1]

                end
                
            end
           
            ## block 2
            
            alpha_2_prop = (rand(Normal(alpha_2_out[i-1], alpha_2_propsd),1))[1]
            if alpha_2_prop .< 0
                posterior_mem_a2[i] = posterior_mem_a1[i] # posterior_mem_a1[i] # 

                alpha_1_out[i] = alpha_1_out[i]
                alpha_2_out[i] = alpha_2_out[i-1]
                alpha_3_out[i] = alpha_3_out[i-1]
                delta_1_out[i] = delta_1_out[i-1]
                delta_2_out[i] = delta_2_out[i-1]

            elseif alpha_2_prop .>= 0
                posterior_mem_a2[i] = posterior(a,b,alpha_1_out[i],alpha_2_prop,alpha_3_out[i-1],
                    delta_1_out[i-1],delta_2_out[i-1])

                aceptprob = posterior_mem_a2[i] - posterior_mem_a1[i]#posterior_mem_a1[i]

                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_prop
                    alpha_3_out[i] = alpha_3_out[i-1]
                    delta_1_out[i] = delta_1_out[i-1]
                    delta_2_out[i] = delta_2_out[i-1]
                    posterior_mem_a2[i] = posterior_mem_a2[i]

                    acceptance2 = acceptance2 + 1
                else
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i-1]
                    alpha_3_out[i] = alpha_3_out[i-1]
                    delta_1_out[i] = delta_1_out[i-1]
                    delta_2_out[i] = delta_2_out[i-1]
                    posterior_mem_a2[i] = posterior_mem_a1[i]
                end
            end
            
            ### block 3
            alpha_3_prop = (rand(Normal(alpha_3_out[i-1], alpha_3_propsd),1))[1]
            
            if alpha_3_prop .< 0
                posterior_mem_a3[i] = posterior_mem_a2[i] 
                alpha_1_out[i] = alpha_1_out[i]
                alpha_2_out[i] = alpha_2_out[i]
                alpha_3_out[i] = alpha_3_out[i-1]
                delta_1_out[i] = delta_1_out[i-1]
                delta_2_out[i] = delta_2_out[i-1]
            elseif alpha_3_prop .>= 0
                posterior_mem_a3[i] = posterior(a,b,alpha_1_out[i],alpha_2_out[i],alpha_3_prop,
                    delta_1_out[i-1],delta_2_out[i-1])
                aceptprob = posterior_mem_a3[i] - posterior_mem_a2[i]
                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_prop
                    delta_1_out[i] = delta_1_out[i-1]
                    delta_2_out[i] = delta_2_out[i-1]
                    posterior_mem_a3[i] = posterior_mem_a3[i]

                    acceptance3 = acceptance3 + 1
                else
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_out[i-1]
                    delta_1_out[i] = delta_1_out[i-1]
                    delta_2_out[i] = delta_2_out[i-1]
                    posterior_mem_a3[i] = posterior_mem_a2[i] 
                end
            end
            
            ### block 4
            
            delta_1_prop = (rand(Normal(delta_1_out[i-1], delta_1_propsd),1))[1]
            
            if delta_1_prop .< 0
                posterior_mem_a4[i] = posterior_mem_a3[i] 
                alpha_1_out[i] = alpha_1_out[i]
                alpha_2_out[i] = alpha_2_out[i]
                alpha_3_out[i] = alpha_3_out[i]
                delta_1_out[i] = delta_1_out[i-1]
                delta_2_out[i] = delta_2_out[i-1]

            elseif delta_1_prop .>= 0
                posterior_mem_a4[i] = posterior(a,b,alpha_1_out[i],alpha_2_out[i],alpha_3_out[i],
                    delta_1_prop,delta_2_out[i-1])
                aceptprob = posterior_mem_a4[i] - posterior_mem_a3[i]
                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_out[i] 
                    delta_1_out[i] = delta_1_prop
                    delta_2_out[i] = delta_2_out[i-1]
                    posterior_mem_a4[i] = posterior_mem_a4[i]

                    acceptance4 = acceptance4 + 1
                else
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_out[i]
                    delta_1_out[i] = delta_1_out[i-1]
                    delta_2_out[i] = delta_2_out[i-1]
                    posterior_mem_a4[i] = posterior_mem_a3[i] 
                end
            end
            
            
            
            
            ### block 5
            delta_2_prop = (rand(Normal(delta_2_out[i-1], delta_2_propsd),1))[1]
            
            if delta_2_prop .< 0
                posterior_mem_a5[i] = posterior_mem_a4[i] 
                alpha_1_out[i] = alpha_1_out[i]
                alpha_2_out[i] = alpha_2_out[i]
                alpha_3_out[i] = alpha_3_out[i]
                delta_1_out[i] = delta_1_out[i]
                delta_2_out[i] = delta_2_out[i-1]

            elseif delta_2_prop .>= 0
                posterior_mem_a5[i] = posterior(a,b,alpha_1_out[i],alpha_2_out[i],alpha_3_out[i],
                    delta_1_out[i],delta_2_prop)
                aceptprob = posterior_mem_a5[i] - posterior_mem_a4[i]
                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_out[i] 
                    delta_1_out[i] = delta_1_out[i]
                    delta_2_out[i] = delta_2_prop
                    posterior_mem_a5[i] = posterior_mem_a5[i]

                    acceptance5 = acceptance5 + 1
                    break
                else
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_out[i]
                    delta_1_out[i] = delta_1_out[i]
                    delta_2_out[i] = delta_2_out[i-1]
                    posterior_mem_a5[i] = posterior_mem_a4[i]
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
    delta_1_final = deleteat!(delta_1_out,1:burnedinsamples)
    delta_2_final = deleteat!(delta_2_out,1:burnedinsamples)
    
    println("Acceptance rate = ", ((acceptance1/N)+(acceptance2/N)+(acceptance3/N)+(acceptance4/N)+(acceptance5/N))/5)
    return [alpha_1_final,alpha_2_final,alpha_3_final,delta_1_final,delta_2_final]

    
    
end


### Example of an MCMC 
@time mcmc = mhrw(30000,0,2,4,0.10,0.01,0.001,1.40,2.90,0.0398,0.00511,0.0004328,0.104,0.1351)


## Traceplots and summary statistics

c1 = (mcmc[1])[setdiff(1:end, 1: 3000, )] ## burning samples
p1 = plot(c1, title = "Alpha 1", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c1, [0.025, 0.5, 0.975]))
describe(c1)
println("Standard Deviation: ", std(c1))

c2 = (mcmc[2])[setdiff(1:end, 1: 3000, )]
p2 = plot(c2, title = "Alpha 2", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c2, [0.025, 0.5, 0.975]))
describe(c2)
println("Standard Deviation: ", std(c2))

c3 = (mcmc[3])[setdiff(1:end, 1: 3000, )]
p3 = plot(c3, title = "Alpha 3", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c3, [0.025, 0.5, 0.975]))
describe(c3)
println("Standard Deviation: ", std(c3))


c4 = (mcmc[4])[setdiff(1:end, 1: 30000, )]
p4 = plot(c4, title = "Delta 1", xaxis = (font(5)), yaxis = (font(8)))

println(quantile!(c4, [0.025, 0.5, 0.975]))
describe(c4)
println("Standard Deviation: ", std(c4))

c5 = (mcmc[5])[setdiff(1:end, 1: 3000, )]
p5 = plot(c5, title = "Delta 2", xaxis = (font(5)), yaxis = (font(8)))


println(quantile!(c5, [0.025, 0.5, 0.975]))
describe(c5)
println("Standard Deviation: ", std(c5))
