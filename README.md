# MAM---Flexible-Monte-Carlo-Agent-based-Model-for-Modelling-COVID-19-Spread


The purpose of the model is the model the spread of COVID-19. Here we are modeling three major outbreaks in Israel using MAM_deleon.m, where day1= 01 Dec 2020 (before the first vaccination campaign in Israel. e modeling three major outbreaks in Israel, where day1= 01 Dec 2020
(before the first vaccination campaign in Israel). 

First outbreak - Day 1-Day 90
Second outbreak - Day 200-Day 320
Third outbreak - Day 365 - Day 450

For all of the outbreaks, the data for the vaccination rate is given in who_vac.mat where the first column is the date of the first vaccination for each person (assuming that all the people that were vaccinated in the first dose also got the second dose).
R.csv gives the R_t as a function of time from day 1 to day 500
Fit. mat contains the analytical function for the vaccine effectiveness as a function of time. 

For the validation of the model, the following files were used:

Prevalence_of_virants_in_Israel.txt contains the Prevalence of delta and omicron in Israel as a function of time.

vaccinated-per-day-2022-07-26.csv contains data on the vaccination rate in Israel as a function of time, divided into ages.

CC_DATA.csv  contains all the daily confirmed cases in Israel, divided by age group

Confirmed cases from abroad.csv contains data on the daily confirmed cases in Israel from abroad. 
