#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 14:58:21 2022

"""


import numpy as np
import matplotlib.pyplot as plt

# This model simulates a population in which 
# 1) - every individual is either working or at home (Depending on age)
# 2) -  SEIR: when an infected individual shares places with susceptible, it may infect them 


p_inf = 0.05 #probability per part of day of infecting other (this can be improved a lot)
incubation_period_avg=3
recovery_period_avg=6

class individual:
    def __init__(self, Id,  age, home_id, work_id, status):
        self.id = Id  #individual id
        self.age = age 
        self.status = status
        self.home_id = home_id #each individual has a home and a work id
        self.work_id = work_id  
        self.infectiontime= np.inf #each susceptible is not be infected
        self.location=self.home_id
        #takes care of the whereabouts of people...
    def whereabouts(self, time_window):
        if self.age>65 or self.age<18:  #youngs and elder at home... this is a toy model!
            self.location=self.home_id
        else:
            if time_window==1:  # 0 (from 00 to 8am), 1 (from 9am to 5pm), 2(6pm to 12am)
                self.location=self.work_id
            else:
                self.location=self.home_id
                
    def infection(self,d=0):  #this updates the person infection process
        if self.status=="S":
            self.infectiontime = d+1+np.random.exponential(scale=recovery_period_avg)
            self.status="E"   
            self.place_of_infection = self.location
        elif self.status=="E" and d>=self.infectiontime:
            self.status="I"
            self.recovery=d+1+np.random.exponential(scale=incubation_period_avg)
        #recovery   
        elif self.status=="I" and d>=self.recovery:
            self.status="R"
#generate a distribution of individuals, households, workplaces


#People get assigned a unique house id and work id, which will determine their whereabouts
N = 2_500
N_Households = 400
N_offices = 200
number_of_days=120

#age: to be replaced with something more realistic 
age =np.random.choice(range(80), size=N, replace=True)
households = np.random.choice(range(N_Households), size=N, replace=True)
offices = np.random.choice(range(N_Households, N_Households+N_offices), size=N, replace=True)



#Initialise a population
people=[]
for i in range(N):
    people.append(individual(i,age[i],households[i],offices[i],"S"))

#ten people infected
for i in range(10):
    people[i].status="I"
    people[i].infectiontime=0
    people[i].recovery = 1+np.random.exponential(scale=incubation_period_avg)
    
infected_per_day = []
susceptible_per_day = []
exposed_per_day=[]
recovered_per_day=[]

infected_at_work=np.zeros(number_of_days)
infected_at_home=np.zeros(number_of_days)
#now simulation begins
for d in range(number_of_days):
    for t in range(3):
        
        infected_people = [person for person in people if(person.status=="I") ]
        #check location of people
        for person in people:
            person.whereabouts(t)
        
        #list of locations where infected people are
    
        unique_locations = np.unique(np.array([person.location for person in infected_people]))

           
        #count where people get infected
        infections_at_home=0
        infections_at_work=0
        for person in people:         
            #now we need to find all individuals who share place with infected
            if (person.location in unique_locations and person.status=="S"):
                #generate random number. This should be replaced by a more complicated function
                random_numbers = np.random.uniform()<=p_inf 
                
                if (random_numbers): 
                    person.infection(d)
                    if person.location in households: #check where infections happen
                        infections_at_home+=1
                    else:
                        infections_at_work+=1
   
            if (person.status!='S'):
                person.infection(d)
        infected_at_work[d] += infections_at_work
        infected_at_home[d] += infections_at_home
    infected_per_day.append(sum( [person.status=="I" for person in people] ))
    susceptible_per_day.append(sum( [person.status=="S" for person in people] ))    
    exposed_per_day.append(sum( [person.status=="E" for person in people] ))            
    recovered_per_day.append(sum( [person.status=="R" for person in people] ))            
    
plt.scatter(range(number_of_days),infected_per_day, color='r', label='I')
plt.scatter(range(number_of_days),recovered_per_day, color='b', label='R')
plt.legend()
plt.xlabel("day")

plt.show()

plt.scatter(range(number_of_days),infected_at_work, color='r', label="Work")
plt.scatter(range(number_of_days),infected_at_home, color='b', label="Home")
plt.xlabel("day")
plt.ylabel("Number of infections")
plt.legend()






    
        
        
        
        
        
        
        
        
        
        
        