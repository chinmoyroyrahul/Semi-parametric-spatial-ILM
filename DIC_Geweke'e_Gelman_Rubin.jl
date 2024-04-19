### Example of How to calculate AIC and  DIC 

### For example we have mcmc
@time mcmc715 = mhrw(30000,0,1.5,3.0,0.10,0.01,0.001,0.024,0.0028,0.00042) # working


y_715 = mcmc715
# println("delta_1.5_3: y_715: ", y_715)

### Arrays
c1715 = (y_715[1])[setdiff(1:end, 1: 5000, )]
c2715 = (y_715[2])[setdiff(1:end, 1: 5000, )]
c3715 = (y_715[3])[setdiff(1:end, 1: 5000, )]

##AIC & DIC,
## AIC = 2k - 2 log-likelihood # where k is the number of estimated parameters in the model

# parameter

theta_512 = [mean(c1715);mean(c2715);mean(c3715)]

## k

k_512 = size(theta_512)[1]
delta = [1.5,3]
## Log-likelihood

l_512 =  epiloglike3(theta_512,delta)

## AIC
AIC = (-2 * l_512) + (2*k_512)

## DIC

## DIC for parametric

d_3_512 = [] # empty set of deviance

for i in 1:size(c1715)[1]
    d_3_512i = -2 .* epiloglike3([c1715[i],c2715[i],c3715[i]],delta)
    push!(d_3_512,d_3_512i)
end

d_bar_3_512 = mean(d_3_512) # overall deviance mean for three steps 
d_est_3_512 = -2 .* epiloglike3([mean(c1715),mean(c2715),mean(c3715)],delta) ## D(theta bar)

# ## effective number of parameters of the model

pd_512 = d_bar_3_512 - d_est_3_512 


## DIC calculation

DIC = d_est_3_512 + (2* pd_512)

println("AIC:  ",AIC, "  ","DIC:  ", DIC)





## convergence test

## Geweke'e test
# using Pkg
# Pkg.add("MCMCDiagnostics")
using MCMCDiagnostics

##
# gewekediag(mcmc; first=0.1, last=0.5, etype=:imse) ## builtin code in julia

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
### Traceplots
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



