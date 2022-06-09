library(ashr)
library(mashr)


dat <- readRDS("formatted/FastQTLSumStats.list.mash.rds")
names(dat)
dim(dat$strong.b)
dim(dat$strong.z)
dim(dat$strong.s)
 
dim(dat$random.b)
dim(dat$random.z)
dim(dat$random.s)


data.strong = list(Bhat=dat$strong.b,  Shat=dat$strong.s )
data.random = list(Bhat=dat$random.b,  Shat=dat$random.s )
 
m.1by1 = mash_1by1(mash_set_data(data.strong$Bhat, data.strong$Shat))
strong.subset = get_significant_results(m.1by1,0.05)

random.subset = sample(1:nrow(data.random$Bhat), 10000)

data.temp = mash_set_data(data.random$Bhat[random.subset,],  data.random$Shat[random.subset,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)




data.random2 = mash_set_data(data.random$Bhat[random.subset,],  data.random$Shat[random.subset,], V=Vhat)
data.strong2 = mash_set_data(data.strong$Bhat[strong.subset,],  data.strong$Shat[strong.subset,], V=Vhat)


U.pca = cov_pca(data.strong2, 5)
U.ed = cov_ed(data.strong2, U.pca)


U.c = cov_canonical(data.random2)
m = mash(data.random2, Ulist = c(U.ed,U.c), outputlevel = 1)

m2 = mash(data.strong2, g=get_fitted_g(m), fixg=TRUE)

head(get_lfsr(m2))
head(get_pm(m2))
head(get_psd(m2))
dim(get_lfsr(m2))
dim(get_pm(m2))
dim(get_psd(m2))



results = get_significant_results(m2, thresh = 0.001, conditions = NULL, sig_fn = get_lfsr )

print(  length(results) ) 

pdf("top5.pdf")
mash_plot_meta(m.c,get_significant_results(m.c)[1], ylab = "Brain regions")
mash_plot_meta(m.c,get_significant_results(m.c)[2], ylab = "Brain regions")
mash_plot_meta(m.c,get_significant_results(m.c)[3], ylab = "Brain regions")
mash_plot_meta(m.c,get_significant_results(m.c)[4], ylab = "Brain regions")
mash_plot_meta(m.c,get_significant_results(m.c)[5], ylab = "Brain regions")
dev.off()





m3 = mash_compute_posterior_matrices(m2, data.strong2, algorithm.version = 'R',  posterior_samples = 100)
names(m3)

m2$result = m3


library(corrplot)

x3 = get_pairwise_sharing(m2, factor=0.5)

pdf("Pairwise_sharing.pdf")
corrplot(x3, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
dev.off()




x3a = get_pairwise_sharing(m2, factor=0.5, lfsr_thresh = 1)
pdf("Pairwise_sharing-2.pdf")
corrplot(x3a, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
dev.off()






x3b = get_pairwise_sharing(m2, factor=0.5, lfsr_thresh = 0.05)
pdf("Pairwise_sharing-3.pdf")
corrplot(x3b, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
dev.off()










x3a = get_pairwise_sharing_from_samples(m2, factor=0.5, lfsr_thresh = 1)
pdf("Pairwise_sharing-4.pdf")
corrplot(x3a, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
dev.off()






x3b = get_pairwise_sharing_from_samples(m2, factor=0.5, lfsr_thresh = 0.05)
pdf("Pairwise_sharing-5.pdf")
corrplot(x3b, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
dev.off()











x3b = get_pairwise_sharing(m2, factor=0.5, lfsr_thresh = 0.001)
pdf("Pairwise_sharing-6.pdf")
corrplot(x3b, method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
dev.off()












##############################  Introduction to mashr
library(ashr)
library(mashr)
set.seed(1)

simdata = simple_sims(500,5,1)
names(simdata)
dim(simdata$B)
dim(simdata$Bhat)
dim(simdata$Shat)
simdata$B
simdata$Bhat
simdata$Shat

## Step 1: Read in the data
## To run mash you need data consisting of a matrix of effects (Bhat) and a matrix of standard errors (Shat), for J effects (rows) in R conditions (columns).
data = mash_set_data(simdata$Bhat, simdata$Shat)
data
names(data)

## Step 2: Set up the covariance matrices
## There are two types of covariance matrix you can use in mash: “canonical” and “data-driven”. The canonical ones are very easy to set up and so we use those here for illustration. However, in applications you will likely also want to use data-driven matrices.
U.c = cov_canonical(data)  
print(names(U.c))

## Step 3: fit the model
m.c = mash(data, U.c)
names(m.c)

## Step 4: Extract Posterior Summaries
## You can extract estimates (posterior means and posterior standard deviations) and measures of significance (local false sign rates) using functions like:
## get_pm (posterior mean), get_psd (posteriore standard deviation) and get_lfsr (local false sign rate)
head(get_lfsr(m.c))
head(get_pm(m.c))
head(get_psd(m.c))
dim(get_lfsr(m.c))
dim(get_pm(m.c))
dim(get_psd(m.c))

## Use get_significant_results to find the indices of effects that are “significant”, which here means they have lfsr less than t in at least one condition, where t is a threshold you specify (default 0.05)
get_significant_results(m.c, thresh = 0.05, conditions = NULL, sig_fn = get_lfsr)
head( get_significant_results(m.c)) 
print(length(get_significant_results(m.c)))
print(length(get_significant_results(m.c, conditions=1)))

## Measure of fit (log-likelihood)
## Use get_loglik to find the log-likelihood of the fit (this will only be useful when you have other fits to compare it with!)
print(get_loglik(m.c))

## Estimated mixture proportions
## Use get_estimated_pi to extract the estimates of the mixture proportions for different types of covariance matrix:
print(get_estimated_pi(m.c))
barplot(get_estimated_pi(m.c),las = 2)



## Metaplot
## The following produces a meta-plot based on the posterior means and posterior variances of an effect. 
mash_plot_meta(m.c,get_significant_results(m.c)[1])
mash_plot_meta(m.c,get_significant_results(m.c)[2])
mash_plot_meta(m.c,get_significant_results(m.c)[3])
mash_plot_meta(m.c,get_significant_results(m.c)[4])
mash_plot_meta(m.c,get_significant_results(m.c)[5])











##############################   Introduction to mash: data-driven covariances
data   = mash_set_data(simdata$Bhat, simdata$Shat)
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)

U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))
U.ed = cov_ed(data, U.pca, subset=strong)


## In general we recommend running mash with both data-driven and canonical covariances. You could do this by combining the data-driven and canonical covariances as in this code:
U.c = cov_canonical(data)  
m   = mash(data, c(U.c,U.ed))
print(get_loglik(m),digits = 10)

print(length(get_significant_results(m)))
print(length(get_significant_results(m, conditions=1)))













































