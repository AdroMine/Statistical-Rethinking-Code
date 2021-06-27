library(rethinking)
data(WaffleDivorce) 
d <- WaffleDivorce

# standardize variables 
d$D <- standardize( d$Divorce ) 
d$M <- standardize( d$Marriage ) 
d$A <- standardize( d$MedianAgeMarriage )

# linear model with Divorce ~ MedianAge at marriage

m5.1 <- quap(
    alist(
        D ~ dnorm(mu, sigma), 
        mu <- a + bA*A, 
        a ~ dnorm(0, 0.2), 
        bA ~ dnorm(0, 0.5), 
        sigma ~ dexp(1)
    ), 
    data = d
)

age_seq_s <- seq(-2, 2, by = 1)
mu <- link(m5.1, list(A = age_seq_s))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI, 0.97)

plot(A ~ D, data = d, col = rangi2, pch = 16, xaxt = 'n', yaxt = 'n')
axis(side = 1, at = age_seq_s, labels = round(age_seq_s * attr(d$D, "scaled:scale") + attr(d$D, "scaled:center"), 1))
axis(side = 2, 
     at = seq(min(d$A),
              max(d$A), l = 5),
     labels = round(seq(min(d$MedianAgeMarriage), 
                        max(d$MedianAgeMarriage), l = 5), 0))
lines(age_seq_s, mu_mean)
shade(mu_PI, age_seq_s)

# Simulating from priors
set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post = prior, data = list(A = c(-2, 2)))
plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) )
for ( i in 1:50 ) lines( c(-2,2) , mu[i,] , col=col.alpha("black",0.4) )

# Simulating from the posterior
# compute percentile interval of mean 
A_seq <- seq( from=-3 , to=3.2 , length.out=30 ) 
mu <- link( m5.1 , data=list(A=A_seq) ) 
mu.mean <- apply( mu , 2, mean ) 
mu.PI <- apply( mu , 2 , PI )
# plot it all 
plot( D ~ A , data=d , col=rangi2 ) 
lines( A_seq , mu.mean , lwd=2 ) 
shade( mu.PI , A_seq )


# Other regression
m5.2 <- quap( 
    alist( 
        D ~ dnorm( mu , sigma ) , 
        mu <- a + bM * M , 
        a ~ dnorm( 0 , 0.2 ) , 
        bM ~ dnorm( 0 , 0.5 ) , 
        sigma ~ dexp( 1 )
    ) ,
    data = d )

M_seq <- seq( from=-3 , to=3.2 , length.out=30 ) 
mu <- link( m5.2 , data=list(M=M_seq) ) 
mu.mean <- apply( mu , 2, mean ) 
mu.PI <- apply( mu , 2 , PI )
# plot it all 
plot( D ~ M , data=d , col=rangi2 ) 
lines( M_seq , mu.mean , lwd=2 ) 
shade( mu.PI , M_seq , col = col.alpha('blue', 0.9))

# Drawing DAG
library(dagitty)
dag5.1 <- dagitty("dag{A -> D; A -> M; M -> D }")
coordinates(dag5.1) <- list(x = c(A = 0, D = 1, M = 2), 
                            y = c(A = 0, D = 1, M = 0))
drawdag(dag5.1)

DMA_dag2 <- dagitty('dag{ D <- A -> M }')
impliedConditionalIndependencies(DMA_dag2)
# D _||_ M | A
# above means D is independent of M conditional on A (M tells nothing new about D after accounting for A)

DMA_dag1 <- dagitty('dag{ D <- A -> M -> D }')
impliedConditionalIndependencies(DMA_dag1) # no conditional independency, so no output



# multiple regression -----------------------------------------------------

m5.3 <-  quap( 
    alist( D ~ dnorm( mu , sigma ) ,
           mu <- a + bM*M + bA*A ,
           a ~ dnorm( 0 , 0.2 ) ,
           bM ~ dnorm( 0 , 0.5 ) ,
           bA ~ dnorm( 0 , 0.5 ) ,
           sigma ~ dexp( 1 )
    ) ,
    data = d 
)
precis( m5.3 )

# plot parameters of all three models until now focusing only on the slope parameters
# all models need to have the same name for the parameters
plot( coeftab(m5.1,m5.2,m5.3), par = c("bA","bM") )

# we can say that once we know median age at marriage for a State, there is little or no additional predictive
# power in also knowing the rate of marriage in that State


# Extra exercise, check relationship between age at Marriage A, and marriage Rate, M
m5.3_5 <- quap(
    alist(
        M ~ dnorm(mu, sigma), 
        mu <- a + b*A, 
        a ~ dnorm(0, 0.5), 
        b ~ dnorm(0, 0.5), 
        sigma ~ dexp(1)
    ), 
    data = d
)

precis(m5.3_5) # does seem to be a negative relationship

seq_A <- seq(-3, 3, by = 1)

mu <- link(m5.3_5, list(A = seq_A))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(M ~ A, data = d, pch = 20)
lines(seq_A, mu.mean)
shade(mu.PI, seq_A)

# Simulating the divorce example
# Simulating the DAG (M <- A -> D)

N <- 50 # number of simulated States 
age <- rnorm( N ) # sim A
mar <- rnorm( N , -age ) # sim A -> M 
div <- rnorm( N , age ) # sim A -> D
div2 <- rnorm(N, age + mar) # both A and M influence D

d2 <- data.frame(A = age, 
                 M = mar, 
                 D = div, 
                 D2 = div2)

# Using above variables in model m5.3 (could use in m5.1 and m5.2 as well)
# which means build above models using these variables?
# most likely yes

m5.3_6 <- quap(
    alist(D ~ dnorm(mu, sigma), 
          mu <- a + bA * A, 
          a ~ dnorm(0, 0.5), 
          bA ~ dnorm(0, 0.5), 
          sigma ~ dexp(1)
    ), 
    data = d2
)
precis(m5.3_6)

m5.3_7 <- quap(
    alist(D ~ dnorm(mu, sigma), 
          mu <- a + bM * M, 
          a ~ dnorm(0, 0.5), 
          bM ~ dnorm(0, 0.5), 
          sigma ~ dexp(1)
    ), 
    data = d2
)
precis(m5.3_7)

m5.3_8 <- quap(
    alist(D ~ dnorm(mu, sigma), 
          mu <- a + bA * A + bM * M, 
          a ~ dnorm(0, 0.5), 
          bA ~ dnorm(0, 0.5), 
          bM ~ dnorm(0, 0.5), 
          sigma ~ dexp(1)
    ), 
    data = d2
)
precis(m5.3_8)

plot( coeftab(m5.3_6, m5.3_7, m5.3_8), par = c("bA","bM") )
# similar to actual data, bA gets wider but stays at similar place
# while bM becomes insignifcant



# Plotting Multivariate posteriors ----------------------------------------

# multple graphs that can be made, we will focus on following three

# - Predictor Residual Plots - plot outcome against residual predictor values. Useful for understanding the model, not much else
# - Posterior prediction plots - model based predictions against raw data (display error in predictions). Tools for checking
#                                fit and assessing predictions. Not causal tools 
# - Counterfactual plots - show implied predictions for imaginary experiments. Allows exploring the causal implications of manipulating
#                          one or more variables


# Predictor Residual Plots  -----------------------------------------------


## we use one predictor to model another predictor

# here we model Marriage Rate M using A (median age at marriage)

m5.4 <- quap(
    alist(
        M ~ dnorm( mu , sigma ), 
        mu <- a + bAM * A, 
        a ~ dnorm( 0, 0.2), 
        bAM ~ dnorm(0, 0.5), 
        sigma ~ dexp(1)
    ), 
    data = d
)

# now residuals can be computed by subracting the observed marriage rate in each State from the predicted rate using above model:

mu <- link(m5.4)
mu_mean <- apply(mu,2, mean)
mu_resid <- d$M - mu_mean

# trying to plot Fig 5.4 from book (page 136)

## top left
# scatter plot of M on A, regression line of M on A + residual lines
plot(M ~ A, data = d, col = "darkblue", cex = 1.5)
ord <- order(d$A)
lines(d$A[ord], mu_mean[ord])
for(i in ord){
    lines(c(d$A[i], d$A[i]), 
          c(mu_mean[i], 
            d$M[i]))
}

## bottom left

# scatterplot of D against resids
# + linear regression of above two

m5.4_1 <- quap(
    alist(
        D ~ dnorm(mu, sigma), 
        mu <- a + bRM * RM, 
        a ~ dnorm(0, 0.5), 
        bRM ~ dnorm(0, 0.5), 
        sigma ~ dexp(1)
    ), data = list(D = d$D, RM = mu_resid)
)
precis(m5.4_1) # bRM not significant

ord <- order(mu_resid) # we want them in sorted order
mur <- mu_resid[ord]
mu <- link(m5.4_1, list(RM = mur))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

plot( d$D ~ mu_resid, col = 'darkblue', cex = 1.5, xlab = "Marriage rate residuals", ylab = "Divorce rate (std)")
lines(mur, mu_mean)
shade(mu_PI, mur)

# similarly plots for A against M can be made



# Posterior Prediction Plots ----------------------------------------------

# two uses - 
# a) did the model correctly approx the posterior dist?
# b) how does the model fail? 


# call link without specifying new data so it uses original data 
mu <- link( m5.3 )
# summarize samples across cases 
mu_mean <- apply( mu , 2 , mean ) 
mu_PI <- apply( mu , 2 , PI )
# simulate observations again no new data, so uses original data 
D_sim <- sim( m5.3 , n=1e4 ) 
D_PI <- apply( D_sim , 2 , PI )

plot( mu_mean ~ d$D , col=rangi2 , ylim=range(mu_PI) , 
      xlab="Observed divorce" , ylab="Predicted divorce" )
abline( a=0 , b=1 , lty=2 ) 
for ( i in 1:nrow(d) ) lines( rep(d$D[i],2) , mu_PI[,i] , col=rangi2 )
identify( x=d$D , y=mu_mean , labels=d$Loc ) # click and identify labels on plot



# Counterfactual Plots ----------------------------------------------------

# counterfactual here means that we use model to make inferences beyond posterior distribution


# Generation steps

## - pick a variable to manipulate, the intervention variable
## - define the range of values to set the intervention variable to
## - for each value of the intervention variable, and for each smaple in posterior, use the causal model 
##   to simulate the values of other variables, including the outcome

## for last part, for the divorce model, m5.3 does not capture A's effect on M which would be required as well
## to estimate influence of A on M, we need to regress A on M. 

data(WaffleDivorce)
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)

m5.3_A <- quap(
    alist(
        ## A -> D <- M
        D ~ dnorm(mu, sigma), 
        mu <- a + bM * M + bA * A, 
        a ~ dnorm(0, 0.2), 
        bM ~ dnorm(0, 0.5), 
        bA ~ dnorm(0, 0.5), 
        sigma ~ dexp(1), 
        ## A -> M
        M ~ dnorm(mu_M, sigma_M), 
        mu_M <- aM + bAM * A, 
        aM ~ dnorm(0, 0.2), 
        bAM ~ dnorm(0, 0.5), 
        sigma_M ~ dexp(1)
    ), data = d
)
precis(m5.3_A)

# let's select A as the variable to manipulate

A_seq <- seq(from = -2, to = 2, length.out = 30)

# next use sim to simulate both M and D in that order (we need to simulate M first before simulating joint influence of A and M on D)

# prep data 
sim_dat <- data.frame( A=A_seq )
# simulate M and then D, using A_seq 
s <- sim( m5.3_A , data=sim_dat , vars=c("M","D") ) # vars tells which variables to simulate and in which order

plot( sim_dat$A , colMeans(s$D) , ylim=c(-2,2) , type="l" , 
      xlab="manipulated A" , ylab="counterfactual D" )
shade( apply(s$D,2,PI) , sim_dat$A ) 
mtext( "Total counterfactual effect of A on D" )

plot(sim_dat$A, colMeans(s$M), ylim = c(-2, 2), type = "l", xlab = "manipulated A", ylab = "counterfactual M")
shade( apply(s$M, 2, PI), sim_dat$A)
mtext( "Counterfactual effect A -> M" )

# new data frame, standardized to mean 26.1 and std dev 1.24 
sim2_dat <- data.frame( A = ( c(20, 30) - 26.1 ) / 1.24 ) 
s2 <- sim( m5.3_A , data=sim2_dat , vars=c("M", "D") ) 
mean( s2$D[,2] - s2$D[,1] ) # effect of increasing median age from 20 to 30
# -4.543678                 # implying 4.5 sd change probably impossibly large



# simulating counterfacutal for an average state with A = 0 and see what changing M does
sim_dat <- data.frame( M=seq(from=-2,to=2,length.out=30) , A=0 ) 
s <- sim( m5.3_A , data=sim_dat , vars="D" ) # see vars, we are only simulating D, not A, becuase M doesn't influence it
plot( sim_dat$M , colMeans(s) , ylim=c(-2,2) , type="l" , xlab="manipulated M" , ylab="counterfactual D" )
shade( apply(s,2,PI) , sim_dat$M ) 
mtext( "Total counterfactual effect of M on D" )



# Masked relationship -----------------------------------------------------


data(milk) 
d <- milk 
str(d)


d$K <- standardize( d$kcal.per.g ) 
d$N <- standardize( d$neocortex.perc ) 
d$M <- standardize( log(d$mass) )

m5.5_draft <-  quap( 
    alist( K ~ dnorm( mu , sigma ) ,
           mu <- a + bN*N ,
           a ~ dnorm( 0 , 1 ) ,
           bN ~ dnorm( 0 , 1 ) ,
           sigma ~ dexp( 1 )
    ) ,
    data=d 
)
# above results in errors since N has NA values

dcc <- d[complete.cases(d$K, d$N, d$M), ]
m5.5_draft <-  quap( 
    alist( K ~ dnorm( mu , sigma ) ,
           mu <- a + bN*N ,
           a ~ dnorm( 0 , 1 ) ,
           bN ~ dnorm( 0 , 1 ) ,
           sigma ~ dexp( 1 )
    ) ,
    data=dcc
)

# first let us consider if priors are reasonable

prior <- extract.prior(m5.5_draft)
xseq <- c(-2, 2)
mu <- link(m5.5_draft, post = prior, data = list(N = xseq))
plot(NULL, xlim = xseq, ylim = xseq)
for(i in 1:50) 
    lines(xseq, mu[i, ], col = col.alpha('black', 0.3))

m5.5 <- quap( 
    alist( 
        K ~ dnorm( mu , sigma ) , 
        mu <- a + bN*N , 
        a ~ dnorm( 0 , 0.2 ) , 
        bN ~ dnorm( 0 , 0.5 ) , 
        sigma ~ dexp( 1 )
    ) , 
    data=dcc )

prior <- extract.prior(m5.5)
xseq <- c(-2, 2)
mu <- link(m5.5, post = prior, data = list(N = xseq))
plot(NULL, xlim = xseq, ylim = xseq)
for(i in 1:50) 
    lines(xseq, mu[i, ], col = col.alpha('black', 0.3))

precis(m5.5)

xseq <- seq( from=min(dcc$N)-0.15 , to=max(dcc$N)+0.15 , length.out=30 ) 
mu <- link( m5.5 , data=list(N=xseq) ) 
mu_mean <- apply(mu,2,mean) 
mu_PI <- apply(mu,2,PI) 
plot( K ~ N , data=dcc ) 
lines( xseq , mu_mean , lwd=2 ) 
shade( mu_PI , xseq )


m5.6 <- quap( 
    alist( 
        K ~ dnorm( mu , sigma ) , 
        mu <- a + bM*M ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) ,
    data=dcc 
) 
precis(m5.6)

mu <- link(m5.6, data = list(M = xseq))
plot( K ~ M, data = dcc)
lines(xseq, colMeans(mu), lwd = 2)
shade(apply(mu, 2, PI), xseq)

m5.7 <- quap( 
    alist( K ~ dnorm( mu , sigma ) , 
           mu <- a + bN*N + bM*M ,
           a ~ dnorm( 0 , 0.2 ) ,
           bN ~ dnorm( 0 , 0.5 ) ,
           bM ~ dnorm( 0 , 0.5 ) ,
           sigma ~ dexp( 1 )
    ) ,
    data=dcc 
)
precis(m5.7)

plot( coeftab( m5.5 , m5.6 , m5.7 ) , pars=c("bM","bN") )

pairs(~ K + M + N, dcc)


# Counterfactual plots using m5.7

## Holding N constant at 0
xseq <- seq( from=min(dcc$M)-0.15 , to=max(dcc$M)+0.15 , length.out=30 ) 
mu <- link( m5.7 , data=data.frame( M=xseq , N=0 ) ) 
mu_mean <- apply(mu,2,mean) 
mu_PI <- apply(mu,2,PI) 
plot( NULL , xlim=range(dcc$M) , ylim=range(dcc$K) ) 
lines( xseq , mu_mean , lwd=2 ) 
shade( mu_PI , xseq )

## Holding M constant at 0
xseq <- seq( from=min(dcc$N)-0.15 , to=max(dcc$N)+0.15 , length.out=30 ) 
mu <- link( m5.7 , data=data.frame( M=0, N=xseq  ) ) 
mu_mean <- apply(mu,2,mean) 
mu_PI <- apply(mu,2,PI) 
plot( NULL , xlim=range(dcc$N) , ylim=range(dcc$K) ) 
lines( xseq , mu_mean , lwd=2 ) 
shade( mu_PI , xseq )

# Simulating a masked relationship

## M -> K <- N
## M -> N
n <- 100
M <- rnorm(n)
N <- rnorm(n, M)
K <- rnorm(n, N - M)
d_sim <- data.frame(K = K, N = N, M = M)
# repeat model 5.5-5.7 with d_sim as data

## M -> K <- N
## N -> M
n <- 100
N <- rnorm(n)
M <- rnorm(n, N)
K <- rnorm(n, N - M)
d_sim2 <- data.frame(K = K, N = N, M = M)


## M -> K <- N
## M <- U -> N
n <- 100
U <- rnorm(n)
N <- rnorm(n, U)
M <- rnorm(n, U)
K <- rnorm(n, N - M)
d_sim3 <- data.frame(K = K, N = N, M = M)

library(dagitty)
dag5.7 <- dagitty( "dag{
                   M -> K <- N
                   M -> N }")
coordinates(dag5.7) <- list(x = c(M = 0, K = 1, N = 2), 
                            y = c(M = 0.5, K = 1, N = 0.5))
MElist <- equivalentDAGs(dag5.7)

drawdag(MElist)


# Categorical Variables ---------------------------------------------------


# Binary Variable

data(Howell1) 
d <- Howell1
str(d)


# h ~ Normal(mu, sigma)
# mu = a + bM * male
# a ~ Normal(178, 20)
# bM ~ Normal(0, 10)
# sigma ~ Uniform(0, 50)

## Above, bM represents the expected difference between males and females in height. 
## a is just the average height of females not everyone (this could imply changing the priors accordingly)
##                                                      (but with lots of data, we can go with a weak prior)
## this approach also assumes that there is more uncertainty about one of the categories 'male' (since it's prediction involves two params)

# prior distributions for mu for females and males

mu_female <- rnorm(1e4, 178, 20)
mu_male   <- rnorm(1e4, 178, 20) + rnorm(1e4, 0, 10)
precis(data.frame(mu_female, mu_male))

# prior for males is wider now

# we don't want prior for males to be more uncertain (before seeing the data, the priors for the two categories should be similar)

# Another approach is to use INDEX VARIABLE (label encoding vs dummy encoding)
# encode male/female as 1, 2 and have separate alpha with same prior for both

# h ~ Normal(mu, sigma)
# mu = a_sex[i]
# a ~ Normal(178, 20) for j = 1..2
# sigma ~ Uniform(0, 50)

## We get two different alphas now a1 and a2 and this solves our problem of males height becoming more uncertain

d$sex <- ifelse(d$male == 1, 2, 1)
m5.8 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <-  a[sex], 
    a[sex] <- dnorm(178, 20), 
    sigma ~ dunif(0, 50)
  ), data = d
)
precis( m5.8 , depth=2 )

# gives expected height in each category

# Extract posterior samples and calculate difference in heights between males and females

post <- extract.samples(m5.8)
post$diff_fm <- post$a[,1] - post$a[,2]
precis(post, depth = 2)

# diff_fm gives the expected difference between a female and male in the sample
# this kind of calculation is called a CONTRAST. 



# Multiple Categories

## We could use dummy variable approach or index variable approach
## with dummy variable, indicatory variables explode (too many created)
## Multilevel models depend upon index variables

data(milk) 
d <- milk
levels(d$clade)

d$clade_id <- as.integer(d$clade)

## Build following model

## K ~ Normal(mu, sigma)               # K is standardised kilocalories
## mu = a_clade[i]
## a ~ Normal(0, 0.5) for j = 1..4
## sigma ~ Exp(1)

d$K <- standardize(d$kcal.per.g)

m5.9 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a[clade_id], 
    a[clade_id] ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d
)
precis(m5.9, depth = 2)

labels <- paste( "a[" , 1:4 , "]:" , levels(d$clade) , sep="" )
plot( precis( m5.9 , depth=2 , pars="a" ) , labels=labels , xlab="expected kcal (std)" )


# Randomly assign these primates to Hogwarts houses

set.seed(63) 
d$house <- sample( rep(1:4, each = 8) , size = nrow(d) )
# 1 Gryffindor, 2 - Hufflepuff, 3 - Ravenclaw, 4 - Slytherin
house_labels <- paste0("h[", 1:4, "]:", c("Gryffindor", "HufflePuff", "Ravenclaw", "Slytherin"))


m5.10 <- quap( 
  alist( 
    K ~ dnorm( mu , sigma ),
    mu <- a[clade_id] + h[house],
    a[clade_id] ~ dnorm( 0 , 0.5 ),
    h[house] ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=d 
)
precis(m5.10, depth = 2) 

plot(precis(m5.10, depth = 2, pars = c("a", "h")), labels = c(labels, house_labels))
