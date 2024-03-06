
#Script where EpiEstim is used to reconstruct daily incidence from weekly aggregated data.
#Uses Marburg data to illustrate the example, corresponding to case study #2 in the Severiy applications,
#where we have access to weekly reports of cases and only total deaths

#install.packages('EpiEstim', repos = c('https://mrc-ide.r-universe.dev', 'https://cloud.r-project.org'))
library(EpiEstim)
library(incidence2)
library(epiparameter)
library(tidyverse)

#Loading data and aggregating by week
Marburg_EqGuinea_linelist <- read_csv("Data/Marburg_EqGuinea_linelist.csv")
Marburg_EqGuinea_linelist$Onset_week <- as.Date(Marburg_EqGuinea_linelist$Onset_week)
Marburg_aggregated <- Marburg_EqGuinea_linelist %>% incidence(date_index = "Onset_week") %>% complete_dates(by = 7L)

#Getting si for Marburg and formatting it to work with EpiEstim
Marburg_si <- epidist_db(disease = "Marburg", epi_dist = "serial_interval", subset = is_parameterised)
Marburg_si <- discretise(Marburg_si)
wrap_si <- function(si) {
  domain <- seq(1L, to = si$prob_dist$qf(0.999), by = 1L)
  pmf <- si$prob_dist$d(domain)
  pmf[1] <- 0
  pmf <- pmf / sum(pmf)
  pmf
}
si_distr <- wrap_si(Marburg_si)

#Using estimate_R with dt=7 to reconstruct daily incidence
Rt_est <- estimate_R(incid = Marburg_aggregated$count,
                     dt = 7L,
                     dt_out = 7L,
                     recon_opt = "naive",
                     iter = 10L,
                     tol = 1e-6,
                     grid = list(precision = 0.001, min = -1, max = 1),
                     config = make_config(si_distr = si_distr),
                     method = "non_parametric_si")



plot(Rt_est, legend = FALSE)


#Reconstructed daily incidence data:
daily_incidence <- Rt_est$I #Why is the length not equal to the interval of in the original dataset?
#The first 7 days have the same value for the incidence, presumably because, as it says on the vignette, 
#the first week doesn't have previous data to reconstruct the incidence from 
#"...the earliest the incidence reconstruction can start is at least the first day of the second aggregation window."

#Now using daily incidence to estimate CFR

#First we create the data frame for {cfr} with daily incidence and the corresponding dates

daily_incidence <- as.data.frame(daily_incidence[-c(1:6)]) #Here I deleted the first 6 rows as a temporary solution to the above question (I know this isn't right)
names(daily_incidence) <- "cases"

Marburg_daily <- Marburg_EqGuinea_linelist %>% incidence(date_index = "Onset_week") %>% complete_dates()
complete_days <- Marburg_daily$date_index

daily_incidence$date <- as.Date(paste(complete_days))
daily_incidence$deaths <- paste(complete_days) #This column is needed but not actually used (to my understanding)

#Creating epidist with Marburg's onset to death delay

set.seed(1)
extract_param(type = "range", values = c(8, 2, 16), distribution = "gamma", samples = 77) 
marburg_onset_death <- epidist(disease = "marburg", epi_dist = "onset_to_death", prob_distribution = "gamma", prob_distribution_params = c(shape=2.095, scale=4.513))
plot(marburg_onset_death)

#Function to estimate ci

ci_text <- function(x,n){
  bin_out <- binom.test(as.numeric(x), as.numeric(n)); est <- round(100*c(bin_out$estimate, bin_out$conf.int))
  paste(est[1], "% (95 CI: ", est[2], "-", est[3],"%)", sep= "")
}

#Estimating known outcomes

known_outcomes_MVD <- cfr::estimate_outcomes(daily_incidence, delay_density = function(x) density(marburg_onset_death, x))
total_known_outcomes <- sum(known_outcomes_MVD$estimated_outcomes)

#Total number of deaths during the outbreak

total_MVD_deaths <- Marburg_EqGuinea_linelist %>% filter(Status == "dead") %>% summarise(total_deaths = length(Onset_week))

#The CFR adjusted for delays is:

ci_text(35, round(total_known_outcomes))

#The CFR not adjusted for delays is:

ci_text(35, 40)


## Comparison to using aggregated data directly and simply adding all cases to the first day of the week (previous version)

#Formatting data for CFR
Marburg_EqGuinea_linelist$Death_week <- Marburg_EqGuinea_linelist$Onset_week
MVD_cases_deaths <- incidence2::incidence(Marburg_EqGuinea_linelist, c("Onset_week","Death_week")) %>% complete_dates()
MVD_cases_deaths <- cfr::prepare_data(MVD_cases_deaths, cases_variable = "Onset_week", deaths_variable = "Death_week")

#Estimating known outcomes
known_outcomes_aggregated <- mes_MVD <- cfr::estimate_outcomes(MVD_cases_deaths, delay_density = function(x) density(marburg_onset_death, x))
total_known_aggregated <- sum(known_outcomes_aggregated$estimated_outcomes)

#CFR of aggregated weekly data adjusted for delays:

ci_text(35, round(total_known_aggregated))


#The difference isn't so obvious in this case when rounding the data (total known outcomes with EpiEstim- 37.6
#total known outcomes with og weekly aggregated data- 38.1)

35/38.1

35/37.6


