library(rethinking)
# gaussian distr ----------------------------------------------------------

mu.list <- seq( from=150, to=170 , length.out=200 ) 
sigma.list <- seq( from=4 , to=20 , length.out=200 ) 
post2 <- expand.grid( mu=mu.list , sigma=sigma.list ) 
post2$LL <- sapply( 1:nrow(post2) , function(i) sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] , log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) + dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) ) 
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE , prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ] 
sample2.sigma <- post2$sigma[ sample2.rows ] 
plot( sample2.mu , sample2.sigma , cex=0.5 , col=col.alpha(rangi2,0.1) , xlab="mu" , ylab="sigma" , pch=16 )



# regression --------------------------------------------------------------

flist <- alist( height ~ dnorm( mu , sigma ) , 
                mu ~ dnorm( 178 , 20 ) ,
                sigma ~ dunif( 0 , 50 ))

m4.1 <- quap( flist , data=d2 )
precis( m4.1 )
post <- extract.samples(m4.1, n = 1e4)
precis(post, hist = FALSE)
plot(post, col = col.alpha(rangi2, 0.3), cex = 0.5, pch = 16)


# Priors ------------------------------------------------------------------

set.seed(2971) 
N <- 100
a <- rnorm( N , 178 , 20 ) 
b <- rnorm( N , 0 , 10 )

plot( NULL , xlim=range(d2$weight) , ylim=c(-100,400) , xlab="weight" , ylab="height" )
abline( h=0 , lty=2 )
abline( h=272 , lty=1 , lwd=0.5 )
mtext( "b ~ dnorm(0,10)" )
xbar <- mean(d2$weight)
for ( i in 1:N ) 
    curve( a[i] + b[i]*(x - xbar) , 
           from = min(d2$weight) , 
           to=max(d2$weight) , 
           add=TRUE , 
           col=col.alpha("black",0.2) )

b <- rlnorm( 1e4 , 0 , 1 ) 
dens( b , xlim=c(0,5) , adj=0.1 )

set.seed(2971) 
N <- 100
a <- rnorm( N , 178 , 20 ) 
b <- rlnorm( N , 0 , 1 )
plot( NULL , xlim=range(d2$weight) , ylim=c(-100,400) , xlab="weight" , ylab="height" )
abline( h=0 , lty=2 )
abline( h=272 , lty=1 , lwd=0.5 )
mtext( "b ~ dnorm(0,10)" )
xbar <- mean(d2$weight)
for ( i in 1:N ) 
    curve( a[i] + b[i]*(x - xbar) , 
           from = min(d2$weight) , 
           to=max(d2$weight) , 
           add=TRUE , 
           col=col.alpha("black",0.2) )



# posterior ---------------------------------------------------------------
data(Howell1); 
d <- Howell1; 
d2 <- d[ d$age >= 18 , ]

xbar <- mean(d2$weight)

m4.3 <- quap( 
    alist(
        height ~ dnorm( mu , sigma ) , 
        mu <- a + b*( weight - xbar ) ,
        a ~ dnorm( 178 , 20 ) ,
        b ~ dlnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 50 )
    ) , 
    data=d2 
)



# plotting model ----------------------------------------------------------

plot( height ~ weight , data=d2 , col=rangi2, pch = 20 ) 
post <- extract.samples( m4.3 ) 
a_map <- mean(post$a) 
b_map <- mean(post$b) 
curve(a_map + b_map*(x - xbar) , add=TRUE, col = 'red')

N <- 300
dN <- d2[ 1:N , ] 
mN <- quap( 
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b*( weight - mean(weight) ) ,
        a ~ dnorm( 178 , 20 ) ,
        b ~ dlnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=dN )

# extract 20 samples from the posterior 
post <- extract.samples( mN , n=20 )
# display raw data and sample size 
plot( dN$weight , dN$height , xlim=range(d2$weight) , ylim=range(d2$height) , col=rangi2 , xlab="weight" , ylab="height", pch = 20 ) 
mtext(concat("N = ",N))
# plot the lines, with transparency 
for ( i in 1:20 ) 
    curve( post$a[i] + post$b[i]*(x-mean(dN$weight)) , 
           col=col.alpha("black",0.3) , add=TRUE )


post <- extract.samples( m4.3 ) 
mu_at_50 <- post$a + post$b * ( 50 - xbar )
PI(mu_at_50)

mu <- link(m4.3)

weight.seq <- seq( from=25 , to=70 , by=1 )
mu <- link( m4.3 , data=data.frame(weight=weight.seq) )


plot( height ~ weight , d2 , type="n" )
for ( i in 1:100 ) 
    points( weight.seq , mu[i,] , pch=16 , col=col.alpha(rangi2,0.1) )

mu.mean <- apply( mu , 2 , mean ) 
mu.PI <- apply( mu , 2 , PI , prob=0.89 )


plot( height ~ weight , data=d2 , col=col.alpha(rangi2,0.5) )
lines( weight.seq , mu.mean )
shade( mu.PI , weight.seq )

# simulate with sigma as well
# prediction interval

sim.height <- sim( m4.3 , data=list(weight=weight.seq), n = 1e4 )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )


plot( height ~ weight , d2 , col=col.alpha(rangi2,0.5) )
lines( weight.seq , mu.mean )
shade( mu.PI , weight.seq )
shade( height.PI , weight.seq )



# polynomial regression ---------------------------------------------------

# standardise variables, imp for polynomial regression
data("Howell1")
d <- Howell1
d$weight_s <- ( d$weight - mean(d$weight) )/sd(d$weight) 
d$weight_s2 <- d$weight_s^2 
m4.5 <-  quap( 
    alist(
        height ~ dnorm( mu , sigma ) ,
        mu <- a + b1*weight_s + b2*weight_s2 ,
        a ~ dnorm( 178 , 20 ) ,
        b1 ~ dlnorm( 0 , 1 ) ,
        b2 ~ dnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )

precis(m4.5)

# Plotting

# generate x axis weights
weight.seq <- seq( from = -2.2 , to = 2 , length.out=30 ) 
pred_dat <- list( weight_s=weight.seq , weight_s2=weight.seq^2 ) 

# use link function to calculate multiple values for each weight using multiple sampels of mu
mu <- link( m4.5 , data=pred_dat ) 

# use above to get the mean value at each weight
mu.mean <- apply( mu , 2 , mean ) 

# get Percentile Interval at each point
mu.PI <- apply( mu , 2 , PI , prob=0.89 ) 

# Take the sigma into account as well for computing each point's y-value
sim.height <- sim( m4.5 , data=pred_dat ) 

# Now again compute PI using these predictions (where we are taking sampling variance into place as well)
# above we were only considering the variance in values of beta (just looking at the mean value of y)
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

# plotting the whole thing above

plot( height ~ weight_s , d , col=col.alpha(rangi2,0.5) )  # plot the real data
lines( weight.seq , mu.mean ) # draw the predicted posterior line (the path of \mu)
shade( mu.PI , weight.seq )   # draw the PI interval around predictions (89% interval of the mean of y)
shade( height.PI , weight.seq ) # draw the prediction interval (much wider) where the actual values can tend to be in


# Cubic poly
d$weight_s3 <- d$weight_s^3 
m4.6 <- quap( 
    alist(
        height ~ dnorm( mu , sigma ) , 
        mu <- a + b1*weight_s + b2*weight_s2 + b3*weight_s3 , 
        a ~ dnorm( 178 , 20 ) , b1 ~ dlnorm( 0 , 1 ) , 
        b2 ~ dnorm( 0 , 10 ) , 
        b3 ~ dnorm( 0 , 10 ) , 
        sigma ~ dunif( 0 , 50 )
    ) ,
    data=d )

pred_data3 <- list( weight_s=weight.seq , weight_s2=weight.seq^2 , weight_s3 = weight.seq^3) 
mu3 <- link(m4.6, data = pred_data3)

mu3.mean <- apply(mu3, 2, mean)
mu3.PI <- apply(mu3, 2, PI, 0.89)
mu3.h <- sim(m4.6, pred_data3)
mu3.he <- apply(mu3.h, 2, PI, 0.89)


# plotting but adding real axis labels (not standardised ones)
plot(height ~ weight_s, data = d, col = col.alpha(rangi2, 0.5), xaxt = 'n')
lines(weight.seq, mu3.mean)
shade(mu3.PI, weight.seq)
shade(mu3.he, weight.seq)
at <- seq(-2, 2)
labels <- at * sd(d$weight) + mean(d$weight)
axis(side = 1, at = at, labels = round(labels, 1))

# splines -----------------------------------------------------------------


data("cherry_blossoms")
d <- cherry_blossoms

d2 <- d[ complete.cases(d$doy) , ] # complete cases on doy 

# create basis splintes with 15 knots
num_knots <- 15 
knot_list <- quantile( d2$year , probs=seq(0,1,length.out=num_knots) )
library(splines) 
B <- bs(d2$year,
        knots=knot_list[-c(1,num_knots)] , 
        degree=3 , intercept=TRUE )

# plot the basis splines (not the resulting model, only the basis splines)
plot( NULL , xlim=range(d2$year) , ylim=c(0,1) , xlab="year" , ylab="basis" )
for ( i in 1:ncol(B) ) lines( d2$year , B[,i] , col = i)

# Linear regression now with these basis splines (each column in B as a variable)

m4.7 <- quap( 
    alist( D ~ dnorm( mu , sigma ) ,
           mu <- a + B %*% w ,
           a ~ dnorm(100,10),
           w ~ dnorm(0,10),
           sigma ~ dexp(1)
    ),
    data = list( D = d2$doy , B = B ) ,
    start = list( w = rep( 0 , ncol(B) ) ) )

# can check the posterior means for all params
precis(m4.7, depth = 2)

# plotting posterior distribution
post <- extract.samples( m4.7 )
w <- apply( post$w , 2 , mean )

# Plot the basis functions * weights calculated using quap
plot( NULL , xlim=range(d2$year) , ylim=c(-6,6) , 
      xlab="year" , ylab="basis * weight" )
for ( i in 1:ncol(B) ) lines( d2$year , w[i]*B[,i] , col = col.alpha(i, 0.7), lwd = 4)
points(knot_list, y = rep(6, num_knots), pch = 3)

# Posterior interval for mu (97%)
mu <- link( m4.7 )
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu,2,PI,0.97)

doy <- sim(m4.7, list( D = d2$doy , B = B ))
doy_h <- apply(doy, 2, PI, 0.97)

plot( d2$year , d2$doy , col=col.alpha(rangi2,0.3) , pch=16 ) 
lines(d2$year, mu_mean,  col = 'darkred')
shade( mu_PI , d2$year , col=col.alpha("black",0.5) )
shade(doy_h, d2$year, col = col.alpha('red', 0.2))



# 4.7 Chapter End of Questions --------------------------------------------




