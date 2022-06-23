setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import dataset

if(!exists("vl_screen"))
  vl_screen <- readRDS("Rdata/final_results_table.rds")
clean <- vl_screen[vllib=="vllib002" & class== "enh./enh."]

#simulate some data
set.seed(20160227)
x <- clean$multiplicative
y <- clean$log2FoldChange
#for simple models nls find good starting values for the parameters even if it throw a warning
# m <- nls(y ~ a*x/(b+x))
m <- nls(y ~ x-a*(1-exp(x)), start = list(a= 1))
#get some estimation of goodness of fit
cor(y, predict(m))

#plot
smoothScatter(x, y)
lines(x,
      predict(m),
      lty=2,
      col="red",
      lwd=3)


# model <- nls(y ~ SSlogis(x, a, b, c), data = data.frame(x=x, y= y))
clean <- clean[log2FoldChange>0]
x <- clean$multiplicative
y <- clean$log2FoldChange
model <- nls(y ~ NLS.asymReg(x, init, m, plateau))
smoothScatter(x,y)
abline(0,1)
lines(seq(0,15, length.out=100),
      predict(model, newdata = data.frame(x= seq(0,15, length.out=100))), col= "red",
      lwd= 2, lty= "11")
smoothScatter(predict(model),y)
abline(0,1)






library(deSolve)
#simulating some population growth from the logistic equation and estimating the parameters using nls
log_growth <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN <- R*N*(1-N/K)
    return(list(c(dN)))
  })
}
#the parameters for the logisitc growth
pars  <- c(R=0.2,K=1000)
#the initial numbers
N_ini  <- c(N=1)
#the time step to evaluate the ODE
times <- seq(0, 50, by = 1)
#the ODE
out   <- ode(N_ini, times, log_growth, pars)
#add some random variation to it
N_obs<-out[,2]+rnorm(51,0,50)
#numbers cannot go lower than 1
N_obs<-ifelse(N_obs<1,1,N_obs)
#plot
plot(times,N_obs)

#find the parameters for the equation
SS <- getInitial(y~SSlogis(x,alpha,xmid,scale),
                 data= data.frame(y= y,x= x))

#we used a different parametrization
K_start<-SS["alpha"]
R_start<-1/SS["scale"]
N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
#the formula for the model
log_formula <- formula(y~K*N0*exp(R*x)/(K+N0*(exp(R*x)-1)))
#fit the model
m <- nls(log_formula,start= list(K= K_start,R= R_start,N0= N0_start))
#estimated parameters
summary(m)
#get some estimation of goodness of fit
cor(y,predict(m))
#plot
smoothScatter(x, y)
# tx <- seq(0,15, length.out= 100)
# ty <- predict(m, newdata = data.frame(x= seq(0,15, length.out= 100)))
# lines(tx, ty, col="red", lty=2, lwd=3)
# predict(m, newdata = data.frame(x= seq(0,15, length.out= 100)))


require(aomisc)
# nls fit
model <- nls(y ~ NLS.asymReg(x, init, m, plateau))

model <- nls(y ~ NLS.expoGrowth(x, a, b),
             data = data.frame(x= x, y= y))
tx <- seq(0,15, length.out= 100)
ty <- predict(m, newdata = data.frame(x= seq(0,15, length.out= 100)))
lines(tx, ty, col="red", lty=2, lwd=3)
 

X <- c(1, 3, 5, 7, 9, 11, 13, 20)
Y <- c(8.22, 14.0, 17.2, 16.9
        
        
         
        
       , 19.2, 19.6, 19.4, 19.6)
# nls fit
model <- nls(Y ~ NLS.asymReg(X, init, m, plateau))
