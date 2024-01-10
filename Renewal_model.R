###############################################################################
#Simulate a renewal process (infection) with fixed reproduction number and serial
#interval, with multiple introductions at different times
##
###############################################################################
library(ggplot2)
library(cowplot) #this is to improve plotting

set.seed(2012)
###########################################################################
#############################################
#### Renewal equation
#############################################
###########################################################################

# the Renewal equation is I(t) =  R(t)∫_0^∞ ω(τ) I(t-τ) dτ (with some assumptions)
# where I(t) is the expected number of infected individuals at time t
# note that it depends on two quantities: R(t), the instantaneous reproduction number, and 
# ω(τ), the generation time distribution. R(t) represents how infectious a typical
# individual is at time t, whereas ω(τ) can be seen as a weighting term, that tells
# how likely is for an infected individual to make an infectious contact after a time τ 

# the equation shows that the number of infected people at time t is a function of
# the number of infected people at earlier times, weighted by a function that describes 
# how individual infectiousness evolves in time, times the population level
# average number of contacts at risk at time t per infected individual.


###############################################################################
###############################################################################



# Start with some realistic data for Covid (discretise first):

# The generation time distribution: the nth element is the fraction of
# transmissions that occur n days after infection.

#it would be a gamma distribution with mean 5.5 and std dev 2.14 ((Ferretti, Ledda et al 2020)
x <- seq(0,14,by=0.05)
#mean <- a*scale
#var <- a*scale**2
var<-2.14**2
mean<-5.5
scale <- var/mean
a<-mean**2/var
plot(x,dgamma(x, shape=a, scale=scale), xlab = "time (days)", ylab="pdf (generation time distribution)")
# This is a discretised version, as
# befits covid (Ferretti, Ledda et al 2020)
#
si <- dgamma(seq(1,15), shape=a, scale=scale)
#then, normalise it
si <- si/sum(si)
#for simplicity, assume that after the last day the generation interval distribution is 0
si <- c(0.0028993756, 0.0424730927, 0.1240492512, 0.1872169056, 0.1967815246,
        0.1645307574, 0.1174711149, 0.0747170618, 0.0435081787, 0.0236309301,
        0.0121317746, 0.0059451849, 0.0028018312, 0.0012772359, 0.0005657808)
max_si <- length(si)

# Number of days of seeding infections, before the renewal process starts:
n0 <- 10
# Number of days of renewal process:
N<- 30
# The mean of the exponential distribution from which the number of seed
# infections is drawn each day:
tau <- 10

#
df <- data.frame(group = "Try",
                 day = seq(from = 1, to = n0 + N),
                 date = seq(from = as.Date(epidemic_day_1),
                            by = 1, length.out = N + n0),
                 delta_I = NA_real_,
                 I = NA_real_)
# Incidence due to seed infections:
df$delta_I[1:n0] <- rexp(n = n0, rate = 1/tau)
df$I <- cumsum(df$delta_I)
               
R_t <- 1.5

df$R_tilde <- R_t

# After the seeding period, calculate incidence on each new day, as a convolution
# between previous incidence and the generation time distribution, weighted by
# today's R:
Renewal.process <- function(df,n0,N){
  for (day in seq(from = n0 + 1, to = n0 + N)) {
    # The time window in the recent past for which the incidence then contributes
    # to incidence now:
    infectiousness_range <- seq(from = max(day - max_si, 1), to = day - 1)
    # The incidence for those days:
    contributing_incidences <- df$delta_I[infectiousness_range]
    # The weights for their contribution to incidence now:
    weights <- si[day - infectiousness_range]
    # R today, it will include some effects:
    R_today <- df$R_tilde[[day]]
    # Put together for incidence now:
    incidence_today <- R_today * sum(contributing_incidences * weights)
    df$delta_I[[day]] <- incidence_today
    df$I[[day]] <- df$I[[day - 1]] + incidence_today
    
  }
  return(df)
}
df <- Renewal.process(df,n0,N)
ggplot(df) +
  geom_line(aes(x=date, y=I))+
  theme_bw()

#################################################################################
# Let's look at some data from Covid
################################################################################


library(openxlsx) #to open the data
setwd("/Users/adminaccount/Infectious Disease Dropbox/Francesco Di Lauro/Teaching") 
covid_data <- read.xlsx("Cases_hospitalisations_deaths_COVID_19_untilJuly2022.xlsx")

#translate date to usual format
covid_data$date = as.Date(covid_data$date, origin='1899-12-30')
plot(covid_data$date, covid_data$new_diagnoses, xlab="date", ylab="new diagnoses")

#focus on the period between May and June, 2021
epidemic_day_1 <- as.Date("2021-05-01")
epidemic_day_N <- as.Date("2021-06-30")
covid_data_restricted <- subset(covid_data,covid_data$date >= epidemic_day_1)

covid_data_restricted <- subset(covid_data_restricted,covid_data_restricted$date <= epidemic_day_N)

plot(covid_data_restricted$date, covid_data_restricted$new_diagnoses, xlab="date", ylab="new diagnoses")

#get days as integers
days <- as.numeric( covid_data_restricted$date-covid_data_restricted$date[1])

naive_estimate <- lm(formula = log(new_diagnoses) ~ days, data=covid_data_restricted )
#plot estimate on top of data
lines(covid_data_restricted$date, exp(naive_estimate$coefficients[2]*days) *exp(naive_estimate$coefficients[1]))

#the growth rate according to this model is
r <- naive_estimate$coefficients[2]

#from which our estimate of R
R_est <- ((r+1/scale)^a) / (1/scale^a) 

################################################################################
# Let's seed the renewal equation using this information, see how it fits
###############################################################################
n0<-10
N= as.numeric(epidemic_day_N-epidemic_day_1)
df <- data.frame(group = "Estimate",
                 day = seq(from = 1, to = n0 + N),
                 date = seq(from = as.Date(epidemic_day_1),
                            by = 1, length.out = N + n0),
                 delta_I = NA_real_,
                 I = NA_real_)

#distribute initially infected people across seeding days
tau <- exp(naive_estimate$coefficients[1])

# Incidence due to seed infections:
df$delta_I[1:n0] <- rexp(n = n0, rate = 1/tau)
df$I <- cumsum(df$delta_I)

df$R_tilde <- R_est


df <- Renewal.process(df,n0,N)
ggplot() +
  geom_line(aes(x=df$date[n0:N+n0], y=df$delta_I[n0:N+n0]))+
  theme_bw()+
  geom_point(aes(x=covid_data_restricted$date, y=covid_data_restricted$new_diagnoses))

#The initial condition seems off (assumption: growing epidemic, we should distribute cases accordingly, further R is not stable)


##############################################################################
# Why predictions fail?
#############################################################################

#what if, based on this, we project further in time? say one month?
N= as.numeric(epidemic_day_N-epidemic_day_1) + 30
df <- data.frame(group = "Prediction",
                 day = seq(from = 1, to = n0 + N),
                 date = seq(from = as.Date(epidemic_day_1),
                            by = 1, length.out = N + n0),
                 delta_I = NA_real_,
                 I = NA_real_)

#distribute initially infected people across seeding days
tau <- exp(naive_estimate$coefficients[1])

# Incidence due to seed infections:
df$delta_I[1:n0] <- rep(tau, n0)
df$I <- cumsum(df$delta_I)

df$R_tilde <- R_est

df <- Renewal.process(df,n0,N)
ggplot() +
  geom_line(aes(x=df$date[n0:N+n0], y=df$delta_I[n0:N+n0]))+
  theme_bw()+
  geom_point(aes(x=covid_data_restricted$date, y=covid_data_restricted$new_diagnoses))


#what happened in reality??
#focus on the period between May and June, 2021
epidemic_day_N <- as.Date("2021-07-30")
covid_data_prolong<- subset(covid_data,covid_data$date >= epidemic_day_1)

covid_data_prolong <- subset(covid_data_prolong,covid_data_prolong$date <= epidemic_day_N)



df <- Renewal.process(df,n0,N)
ggplot() +
  geom_line(aes(x=df$date[n0:N+n0], y=df$delta_I[n0:N+n0]))+
  theme_bw()+
  geom_point(aes(x=covid_data_prolong$date, y=covid_data_prolong$new_diagnoses), color='red')+
  geom_point(aes(x=covid_data_restricted$date, y=covid_data_restricted$new_diagnoses))



#############################################################################
## Add some noise to the process?
############################################################################

library(arm)



# Number of days of seeding infections, before the renewal process starts:
n0 <- 10
# Number of days of renewal process:
N<- 30
# The mean of the exponential distribution from which the number of seed
# infections is drawn each day:
tau <- 10

#
df <- data.frame(group = "Complication",
                 day = seq(from = 1, to = n0 + N),
                 date = seq(from = as.Date(epidemic_day_1),
                            by = 1, length.out = N + n0),
                 delta_I = NA_real_,
                 I = NA_real_)
# Incidence due to seed infections:
df$delta_I[1:n0] <- rexp(n = n0, rate = 1/tau)
df$I <- cumsum(df$delta_I)

R_t <- 1.5

#define a maximum value for R_t
R_max <- 4
#gaussian noise standard deviation
stdev<-0.5

#The classic way of adding noise to a quantity which is bound (e.g. 0<=R_t<=R_max always)
# is as follows: 
# 1) consider the quantity scaled in [0-1]: R/R_max
# 2) apply a logistic transformation  f(x) = x/(1-x): this function is not bound
# 3) add noise to the transformed variable
# 4) use the inverse logistic function to return to the process
df$R_tilde = R_t
df$R_transformed <- logit( df$R_tilde/R_max) + rnorm(n=length(df$R_tilde), mean=0, sd=stdev)

df$R_tilde <- invlogit(df$R_transformed)*R_max


df_noise <- Renewal.process(df,n0,N)
ggplot(df_noise) +
  geom_line(aes(x=date, y=I))+
  theme_bw()

###############################################################################
#In this model, R_t does NOT depend on R_{t-1}. If we want to add a dependence
###############################################################################


df$R_tilde <- R_t
for(day in seq(n0,n0+N)){
  transformed_rt <- logit( df$R_tilde[[day-1]]/R_max) + rnorm(n=1, mean=0, sd=stdev)
  df$R_tilde[[day]] <- invlogit(transformed_rt)*R_max
}

df_noise <- Renewal.process(df,n0,N)
plot<-ggplot(df_noise) +
  geom_line(aes(x=date, y=I))+
  theme_bw()

#add an inset with R_t
inset<-ggplot(df_noise) +
  geom_line(aes(x=date, y=R_tilde), color='red')+
  theme_bw()

#this is used from cowplot to make insets
ggdraw() +
  draw_plot(plot) +
  draw_plot(inset, x = 0.15, y = .6, width = .4, height = .3)

