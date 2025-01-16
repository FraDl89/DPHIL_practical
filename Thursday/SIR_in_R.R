library(deSolve) #To solve ODEs
library(ggplot2) #Library to plot results
library(tidyverse)
library(socialmixr)
#####################################################################
#
# COMPARTMENTAL MODELS: SIR (note, not everything will be done during practicals)
#
####################################################################

# Let’s study the SIR model for a closed population, i.e., 
# one in which we can neglect births and deaths.

closed.sir.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  N <- S+I+R
  ## now code the model equations
  dSdt <- -beta*S*I/N
  dIdt <- beta*S*I/N-gamma*I
  dRdt <- gamma*I
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}

# Note that the order and type of the arguments and output of this function must 
# exactly match ode’s expectations. Thus, for instance, the time variable t 
# must be the first argument even if, as is the case here, nothing 
# in the function depends on time.

#if we’re thinking of a disease such as early Covid, and measuring time in days

R_0 <- 3 #Reproduction number ~ 2-3
gamma <- 1/7 #infectious period ~ a week
beta <- R_0*gamma #from the definition of beta
parms <- c("beta"=beta,"gamma"=gamma)

#define time and initial populations
N<-100000
times <- seq(from=0,to=120,by=1)
xstart <- c(S=N-1,I=1,R=0)

ode(
  func=closed.sir.model,
  y=xstart,
  times=times,
  parms=parms
) %>%
  as.data.frame() -> out
#put results in a dataframe and plot it (ggplot2)
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')



#####################################################################
#
# COMPARTMENTAL MODELS: Growth rate and Attack Rate
#
####################################################################


#exponential growth comes from the approximation of the equation for I(t):
#indeed, when S ~~ N (at the beginning of the epidemic) we can say
# dI/dt = beta*S*I/N - gamma *I ~  (beta - gamma)*I = gamma*(R_0-1)*I -> 
#I(t) ~ I_0*exp(gamma*(R_0-1)*t)

out %>%
  dplyr::select(time,I) %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  geom_line(aes(x=time, y=exp(gamma*(R_0-1)*time), color="Approximation")) +
  theme_classic()+
  ylim(c(0,max(50000)))+
  labs(x='time (days)',y='number of individuals')


#R(\infty) is the final attack rate. We find it from the equations:
# dS/dR = - R_0 S/N
# ln S = -R_0 R(t)/N
# S(\infty) = exp (-R_0 (1-S(infty)/N) ) 

#This is a fixed point equation, so we need a particular library
library(spuRs)

function.fixed.point <- function(x,R_0){
  R_0<-3
  return(exp(-R_0 *(1-x)))
}

#we need to give an initial guess
fixedpoint<-fixedpoint(function.fixed.point, x0=0.05, tol = 1e-4, max.iter = 100)

out %>%
  dplyr::select(time,S) %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  geom_hline(aes(yintercept=fixedpoint*100000, color="R infinity")) +
  theme_classic()+
  labs(x='time (days)',y='number of individuals')


#####################################################################
#
# First exercise: can we model in a SIR setting some form of epidemic
# control?
####################################################################

closed.sir.model.controlled <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  c <- params["control"] #control parameter  
  
  N <- S+I+R
  ## now code the model equations
  
  if(t>20 & t<=60){
  beta <- beta*(1-c) #careful here with the numeric stability
  }
  dSdt <- -beta*S*I/N
  dIdt <- beta*S*I/N-gamma*I
  dRdt <- gamma*I
  
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}
#beta -> beta*(1-c)  for 25<t<50, 0<c<1
control=0.3
parms <- c("beta"=beta,"gamma"=gamma, "control"=control)
times <- seq(from=0,to=240,by=1)


ode(
  func=closed.sir.model.controlled,
  y=xstart,
  times=times,
  parms=parms,
) %>%
  as.data.frame() -> out


out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')




#####################################################################
#
# Second exercise: can we add model uncertainty (stochasticity) to this system?
# 
####################################################################

#the infection rate may come from a distribution rather than a point value

#Plot the distribution
shape <- 2
rate <- 10
grid<- seq(0,1,by=0.01)
y<-dgamma(grid, shape, rate)

ggplot()+
  geom_line(aes(x=grid,y=y))+
  xlab("Value")+
  ylab("pdf")+
  theme_bw()

#take a number of samples from this distribution

betasamples <- rgamma(50,shape=shape, rate=rate)

#we use lapply to take advantage of parallelisation, but we need an auxiliary function

auxiliary.function <- function(x){
    ode(
    func=closed.sir.model,
    y=xstart,
    times=times,
    parms=c("beta"=x,"gamma"=gamma),
  ) %>%
    as.data.frame() -> out
  return(out)
}

results <- lapply(betasamples, auxiliary.function)

#To plot all the dataframes in one single figure, we need to assign to each an id
#and create a single dataframe that can be passed to ggplot
ggplot(bind_rows(results, .id="id"), aes(time, I, color=id)) +
  geom_line() +
  theme_bw()


#####################################################################
#IMPORTANT NOTE: this is not the only form of uncertainty we can input.
#We can make a stochastic version of the SIR model: the Gillespie algorithm
#
# the Idea is that infected people infect the susceptible at a rate beta/N*S *I
# and then each recover at a rate gamma. The simulation consider the two events 
# (infection and recovery) as two separate poisson processes.  This is an 'event-based'
#  simulation:  time moves forward as events happen.
# The time at which the next event happen is distributed as exp(gamma*I+beta/N*S*I) (why?)
# and the type of event is 'infection' with probability  beta*S/(N*gamma) (why?)
#
#try to code this and compare the solution with the deterministic SIR. What do you notice?
#####################################################################
sir_Gillespie <- function(beta, gamma, N, S0, I0, R0, tf) {
  time <- 0
  S <- S0
  I <- I0
  R <- R0
  ta <- numeric(0)
  Sa <- numeric(0)
  Ia <- numeric(0)
  Ra <- numeric(0)
  while (time < tf) {
    ta <- c(ta, time)
    Sa <- c(Sa, S)
    Ia <- c(Ia, I)
    Ra <- c(Ra, R)
    pf1 <- beta * S * I
    pf2 <- gamma * I
    pf <- pf1 + pf2
    dt <- rexp(1, rate = pf)
    time <- time + dt
    if (time > tf) {
      break
    }
    ru <- runif(1)
    if (ru < (pf1/pf)) {
      S <- S - 1
      I <- I + 1
    } else {
      I <- I - 1
      R <- R + 1
    }
    if (I == 0) {
      break
    }
  }
  results <- data.frame(time = ta, S = Sa, I = Ia, R = Ra)
  return(results)
}


sir_out <- sir_Gillespie(beta/N,gamma,N,N-5,5,0,100)

ggplot(sir_out_long,aes(x=time,y=value,colour=variable,group=variable))+
  # Add line
  geom_line(lwd=2)+
  #Add labels
  xlab("Time")+ylab("Number")


#####################################################################
#
# Rest of the first part: we can include more and more features, play with the model.
# For example: can we include waning immunity? (SIRS) or incubation time (SEIR)? 
#
####################################################################








#What happens if there is waning of immunity?
closed.sirs.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  xi <- params["xi"]
  N <- S+I+R
  ## now code the model equations
  dSdt <- -beta*S*I/N + xi*R
  dIdt <- beta*S*I/N-gamma*I
  dRdt <- gamma*I - xi*R
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}


#Assume that immunity wanes after an average of 3 months
newparms <- c("beta"=beta,"gamma"=gamma, "xi"=1/90)
#we need to observe the process for longer to see an effect
newtimes <- seq(from=0,to=360,by=1)

ode(
  func=closed.sirs.model,
  y=xstart,
  times=newtimes,
  parms=newparms
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  ggtitle("The SIRS model")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x='time (days)',y='number of individuals')


#(Q. Does this process have a stationary state? Can you find it on pen and paper?

#We should also consider the incubation period



closed.seirs.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  xi <- params["xi"]
  nu <- params["nu"]
  N <- S+E+I+R
  ## now code the model equations
  dSdt <- -beta*S*I/N + xi*R
  dEdt <- beta*S*I/N - nu*E
  dIdt <- nu*E-gamma*I
  dRdt <- gamma*I - xi*R
  ## combine results into a single vector
  dxdt <- c(dSdt,dEdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}


#Assume that incubation lasts for 3 days on average
newparms <- c("beta"=beta,"gamma"=gamma, "xi"=1/180, "nu"=1/3)
#we need to observe the process for longer to see an effect
newtimes <- seq(from=0,to=480,by=1)
xstart <- c(S=99999,E=0,I=1,R=0)

ode(
  func=closed.seirs.model,
  y=xstart,
  times=newtimes,
  parms=newparms
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')


#We could also consider that 
    #1) not all infected people will recover: some will die
    #2) but before this happens, they may be hospitalised
    #2) But some of the people may be asymptomatic: in this case they will recover

more.complete.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  E <- x[2]
  A <- x[3]
  I <- x[4]
  H <- x[5]
  D <- x[6]
  R <- x[7]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  xi <- params["xi"]
  nu <- params["nu"]
  gamma_H <- params["gamma_H"]
  N <- S+E+A+I+H+D+R
  ## now code the model equations
  dSdt <- -beta*S*I/N + xi*R       
  dEdt <- beta*S*I/N - nu*E
  dAdt <- 1/3*nu*E - gamma*A        #say one in three is asymptomatic
  dIdt <- 2/3*nu*E-gamma*I
  dHdt <- 1/10*gamma*I -gamma_H*H   #say one in ten needs hospital
  dDdt <- 1/5*gamma*H 
  dRdt <- 9/10*gamma*I+4/5*gamma_H*H  +gamma*A - xi*R
  ## combine results into a single vector
  dxdt <- c(dSdt,dEdt,dAdt,dIdt,dHdt,dDdt,dRdt)
  ## return result as a list!
  list(dxdt)
}


newparms <- c("beta"=beta,"gamma"=gamma, "gamma_H"=1/10, "xi"=1/180, "nu"=1/3)
#we need to observe the process for longer to see an effect
newtimes <- seq(from=0,to=720,by=1)
xstart <- c(S=99999,E=0,A=0,I=1,H=0,D=0,R=0)

ode(
  func=more.complete.model,
  y=xstart,
  times=newtimes,
  parms=newparms
) %>%
  as.data.frame() -> out
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#it becomes messy to visualise this. But we are interested only in deaths, 
# hospitalisations, total number of infected (E+A+I)

out %>%
  mutate(I=E+I+A) %>%
  dplyr::select(time, I, H, D) %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  ggtitle("Modified SAEIRS")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x='time (days)',y='number of individuals')




#So far we considered "beta" as the per-capita probability of infecting another
#person in the population. We know that this depends on many covariates 
# (e.g. sex,age, job/school - but also policies such as social distancing, vaccines...)

#say we want to consider the age matrix. This should come from the POLYMOD study

# POLYMOD <- read_excel("contact_matrices_152_countries/MUestimates_all_locations_2.xlsx",
#                         sheet = "United Kingdom of Great Britain")

#POLYMOD <- contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0, 5, 10, 15,20), symmetric = TRUE)


POLYMOD_stud<-contact_matrix(polymod, countries = "United Kingdom",
                             age.limits = c(5,10,15,20), symmetric=T, 
                             filter = list(cnt_school=1))

#This has two components: rate of contacts, population.
#We need the obvious constraint mijNi=mjiNj

metapopulation_model <- function(t,x,params,POLYMOD_stud){
  S<-x[1:4]
  I<-x[5:8]
  R<-x[9:12]
  beta <- params["beta"]
  gamma <- params["gamma"]
  contact_m=POLYMOD_stud$matrix
  demography=POLYMOD_stud$demography$population
  #we use matrix formalism to derive the new equations
  dSdt = - beta*S* contact_m %*%I/demography
  dIdt = beta*S* contact_m %*%I/demography - gamma*I
  dRdt = gamma* matrix(data=I)
  
  
  return(list(c(dSdt,dIdt,dRdt)))
}




newparms <- c("beta"=0.2,"gamma"=gamma)
#we need to observe the process for longer to see an effect
newtimes <- seq(from=0,to=70,by=0.5)
#susceptibles 

S <- POLYMOD_stud$demography$population
#Infect 100 people in the 15-20 age class
S[3] <- S[3]-2
S[1] <- S[1]-10

I <- c(10,0,2,0)
R <- rep(0,4)
xstart <- c(S,I,R)

ode(
  func=metapopulation_model,
  y=xstart,
  times=newtimes,
  parms=newparms,
  POLYMOD_stud=POLYMOD_stud
) %>% as.data.frame() -> solution 
  
  names(solution) <- c("t","S_5","S_10","S_15","S_20","I_5","I_10","I_15","I_20","R_5","R_10","R_15","R_20")
  solution %>%
  dplyr::select(t,I_5,I_10,I_15,I_20) %>%
  gather(variable,value,-t) %>%
  ggplot(aes(x=t,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  ggtitle("Metapopulation SIR")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x='time (days)',y='number of individuals')

  
  
  #Imagine that you wanted to be more precise and add many compartments, and many
  #age groups. Further, you'd like to understand where people get infected 
  #(school/work/random?). Finally, you'd like to consider individual heterogeneity
  
  #at this point, an agent-based model is more useful.

