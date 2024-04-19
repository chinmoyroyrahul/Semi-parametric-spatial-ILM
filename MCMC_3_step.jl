
# single parameter update MH-MCMC# 

using Distributions
using Random
using StatsPlots
using MCMCChains
import Random


function mhrw(N,burnedinsamples,delta_1,delta_2,alpha_1_init,alpha_2_init,alpha_3_init,alpha_1_propsd,alpha_2_propsd,
                alpha_3_propsd)
    
    function posterior(alpha_1,alpha_2,alpha_3,delta_1,delta_2)
        
        prior_alpha_1_distr = truncated(Normal(0, 100000), 0, Inf) ## Trancated normal/ half normal distribution
        prior_alpha_1(p) = pdf(prior_alpha_1_distr,p)
        prior_alpha_2_distr = truncated(Normal(0, 100000), 0, Inf)
        prior_alpha_2(p) = pdf(prior_alpha_2_distr,p)
        prior_alpha_3_distr = truncated(Normal(0, 100000), 0, Inf)
        prior_alpha_3(p) = pdf(prior_alpha_3_distr,p)


        out = epiloglike3([alpha_1,alpha_2,alpha_3],[delta_1,delta_2]) + log(prior_alpha_1(alpha_1)) +
                log(prior_alpha_2(alpha_2)) + log(prior_alpha_3(alpha_3))

        return out

    end

    alpha_1_out =  zeros(N)
    alpha_2_out =  zeros(N)
    alpha_3_out =  zeros(N)
    
    alpha_1_out[1] = alpha_1_init
    alpha_2_out[1] = alpha_2_init
    alpha_3_out[1] = alpha_3_init

    posterior_mem_a1 = zeros(N)
    posterior_mem_a1[1] = posterior(alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1,delta_2)
    
    posterior_mem_a2 = zeros(N)
    posterior_mem_a2[1] = posterior(alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1,delta_2)
    
    posterior_mem_a3 = zeros(N)
    posterior_mem_a3[1] = posterior(alpha_1_out[1],alpha_2_out[1],alpha_3_out[1],delta_1,delta_2)


    acceptance1 = 0
    acceptance2 = 0
    acceptance3 = 0

    #
    for i in 2:N

        
        while i <=N
            
            ## block 1
            alpha_1_prop = (rand(Normal(alpha_1_out[i-1], alpha_1_propsd),1))[1]
            if alpha_1_prop .< 0
                posterior_mem_a1[i] = posterior_mem_a3[i-1]
                alpha_1_out[i] = alpha_1_out[i-1]
                alpha_2_out[i] = alpha_2_out[i-1]
                alpha_3_out[i] = alpha_3_out[i-1]
            elseif alpha_1_prop .>=0
                posterior_mem_a1[i] = posterior(alpha_1_prop,alpha_2_out[i-1],alpha_3_out[i-1],delta_1,delta_2)
                aceptprob = posterior_mem_a1[i] - posterior_mem_a3[i-1]
                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_prop
                    posterior_mem_a1[i] = posterior_mem_a1[i]
                    acceptance1 = acceptance1 + 1
                    alpha_2_out[i] = alpha_2_out[i-1]
                    alpha_3_out[i] = alpha_3_out[i-1]
                else
                    alpha_1_out[i] = alpha_1_out[i-1]
                    posterior_mem_a1[i] = posterior_mem_a3[i-1]
                    alpha_2_out[i] = alpha_2_out[i-1]
                    alpha_3_out[i] = alpha_3_out[i-1]
                end
                
            end
           
            ## block 2
            
            alpha_2_prop = (rand(Normal(alpha_2_out[i-1], alpha_2_propsd),1))[1]
            if alpha_2_prop .< 0
                posterior_mem_a2[i] = posterior_mem_a1[i] 
                alpha_1_out[i] = alpha_1_out[i]
                alpha_2_out[i] = alpha_2_out[i-1]
                alpha_3_out[i] = alpha_3_out[i-1]

            elseif alpha_2_prop .>= 0
                posterior_mem_a2[i] = posterior(alpha_1_out[i],alpha_2_prop,alpha_3_out[i-1],delta_1,delta_2)
                aceptprob = posterior_mem_a2[i] - posterior_mem_a1[i]

                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_prop
                    alpha_3_out[i] = alpha_3_out[i-1]
                    posterior_mem_a2[i] = posterior_mem_a2[i]
                    acceptance2 = acceptance2 + 1
                else
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i-1]
                    alpha_3_out[i] = alpha_3_out[i-1]
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
                break
            
            elseif  alpha_3_prop .>= 0
                posterior_mem_a3[i] = posterior(alpha_1_out[i],alpha_2_out[i],alpha_3_prop,delta_1,delta_2)
                aceptprob = posterior_mem_a3[i] - posterior_mem_a2[i]
                u = log(rand())
                if u .< aceptprob
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_prop
                    posterior_mem_a3[i] = posterior_mem_a3[i]
                    acceptance3 = acceptance3 + 1
                    break
                else
                    alpha_1_out[i] = alpha_1_out[i]
                    alpha_2_out[i] = alpha_2_out[i]
                    alpha_3_out[i] = alpha_3_out[i-1]
                    posterior_mem_a3[i] = posterior_mem_a2[i]
                    break
                    
                end
            end
   
        end

    end
    alpha_1_final = deleteat!(alpha_1_out,1:burnedinsamples)
    alpha_2_final = deleteat!(alpha_2_out,1:burnedinsamples)
    alpha_3_final = deleteat!(alpha_3_out,1:burnedinsamples)
    
    println("Acceptance rate = ", ((acceptance1/N)+(acceptance2/N)+(acceptance3/N))/3)
    return [alpha_1_final,alpha_2_final,alpha_3_final]
end

### Example of a MH MCMC 
@time mcmc = mhrw(30000,0,1.5,3.0,0.10,0.01,0.001,0.024,0.0028,0.00042) 


## Trace plot and summary Statistics

using MCMCChains
using StatsPlots
using Plots
using Statistics
   

c1 = (mcmc[1])[setdiff(1:end, 1: 3000, )]
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

