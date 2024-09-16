library(tidyverse)
library(mgcv)

### This provides the regression approach to estimating the force of infection as described in the paper
### Time spent in every month formally fit to allow for the time-varying covariates piece, see methods

##### A. Get the data and tidy ##################################################
df <- read.csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/mnj_rdt_data_foi_20230821.csv")

## Making individual ID code and household ID code a factor for fitting random effects
df$unique_ind_id <- as.factor(df$unique_ind_id)
## First 5 digits of the individual ID code give the site code + household code
df$household <- as.factor(substring(df$unique_ind_id,1,5)) 

##### B. REMIND ONESELF OF BASIC SURVIVAL ANALYSIS ##################################################

## simple binomial likelihood function for sanity check
likelihood.foi <- function(par, time.interval, event) { 
              #event #size #prob                   #log(p)
	u <- dbinom(event, 1, 1-exp(-par*time.interval), log=TRUE)
	return(-sum(u))
}

## calculate likelihood across range of parameters
test.par <- seq(from=0.00001, to=0.002, length=100)
rc <- rep(NA, length(test.par))
for (j in 1:length(test.par)){ 
	rc[j] <- likelihood.foi(par=test.par[j], time.interval=df$interlude.days, event=df$rdt.result)
}

## plot out
plot(test.par, -rc, xlab="force of infection", ylab="log likelihood", type="l")
abline(v=test.par[rc==min(rc)])
test.par[rc==min(rc)]

## show that you can recover this value using logit with cloglog
fit1 <- glm(df$rdt.result~offset(log(df$interlude.days)),family=binomial(link="cloglog"), data=df)
exp(fit1$coeff[1])      ## This should be the same as above
test.par[rc==min(rc)]   ## checks out


##### C. INCLUDING COVARIATES ##################################################

## Fit the force of infection with effects of age and site
## (and then time and household)

## Time: Fitting seasonality as a set of fixed effects using MGCV
## Household: Including random effect of household

#Fit: Smoothed age + Days in month 1-12 + site + offset log (interlude) + smoothed household

fit3 <- gam(rdt.result~s(age_yrs)+
              days.in.MO1+days.in.MO2+days.in.MO3+days.in.MO4+days.in.MO5+days.in.MO6+
              days.in.MO7+days.in.MO8+days.in.MO9+days.in.M10+days.in.M11+days.in.M12+
              site_code+offset(log(interlude.days))+s(household, bs="re"),
            family=binomial(link="cloglog"), data=df)

summary(fit3)
gam.check(fit3)

## look at age pattern
par(mfrow=c(1,1))
plot(fit3)

## Plot month effects (subtracting one, because didn't do -1 in the fit above, so doing contrasts!)
plot(1:12, fit3$coeff[1]+fit3$coeff[2:13], xlab="month", ylab="FOI per day", pch=19)
for (j in 1:12) points(c(j,j),fit3$coeff[1]+fit3$coeff[j+1]+c(-1.96,1.96)*summary(fit3)$se[j+1],type="l",lty=3)


###############################################################################################
## D. WRITING FUNCTIONS TO SIMULATE AND PLOT ##################################################
###############################################################################################

## Coded as a function, so nominally could include any chosen fit - leverages R's 'predict' functions
## Running off the fit above

## First, a function to work with times/seasons using calendar dates
## Input: Start date month; start date day; timing of any interventions
## Create a vector of the number of days in month M.01 to M.12

###############################################################################################
#Function 1
#getting date intervals for a given start date and timing of interventions
f.date.intervals <- function(start.date.mo, start.date.day, timing.interventions){
  
  #Creating a vector of all 365 days + a second year to allow interventions in late months to carry over
  #For standardization, using 2023-2024
  v.sequence.i <- seq(ymd("2023-01-01"), ymd("2024-12-31"), 1)
  
  #Timing interventions is an integer of the days after the start date where there are interventions
  #E.g., if the start date is Mar 01 and the interventions are 20, 40, 60 then there are interventions on
  #Mar01 + 20 = Mar21 ... etc
  #Note: must start with 1
  
  #Using mdy() to get a date format for start date
  start.date <- mdy(paste0(start.date.mo, "-", start.date.day, "-", "2023"))
  #Finding the index of the start date in the 2023-2024 sequence
  start.date.index <- which(v.sequence.i == start.date)
  #Finding the end date (stop at last date of intervention)
  end.date.count <- max(timing.interventions)
  end.date <- ymd(v.sequence.i[start.date.index + end.date.count])
  end.date.index <- which(v.sequence.i == end.date)
  #Pulling out our interval of interest
  my.interval <- v.sequence.i[(start.date.index):end.date.index]
  
  return(my.interval)
}


###############################################################################################
## Second, a function to count the number of days in a given month for an interval
## E.g., for an interval for January 15 to Mar 2: 16 days in Jan, 28 days in Feb, 2 days in Mar
## Input: an interval sequence (e.g., Mar 15 to April 28)
## Output: A list of days per month example (list length = 12)
## List used to allow easy access when predicting

#Function 2
f.date.wrangler.v3 <- function(my.interval){
  
  #counting days per month for a vector that is a sequence of dates
  #adding to a list
  days.in.M.l <- list(NA)
  for(i in 1:12){
    days.in.M.l[[i]] <- length(which(month(my.interval) == i))
  }
  return(days.in.M.l)
}


###############################################################################################
## Third, a function to take a fitted model and predict; 
## with timing.interventions a vector giving what day(s) the interventions happened
##** scale.interventions currently not implemented

#Function 3
#Using a fitted model, intervention start date, intervention timing, and ages of interest
findYearTimeCourse.v3 <- function(fitted.model,           #fitted model from above
                                  start.date.mo,          #Integer from 1:12
                                  start.date.day,         #Integer from 1:31 (or last day of month)
                                  timing.interventions,   #e.g., c(1, 30, 60)
                                  my.ages,                #Integers of age in years
                                  scale.interventions){   #scale.interventions not used yet
  
  #** month a legacy
  
  #Making a vector of site codes
  sites <- c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5")
  ages <- my.ages
  #Passing my start date to the date.interval function to get the interval as a sequence of dates
  my.interval <- f.date.intervals(start.date.mo, start.date.day, timing.interventions)
  
  #time steps evaluated
  tvals <- 1:(max(timing.interventions))
  intervention <- rep(0,length(tvals))
  #noting which days in the sequence have interventions
  intervention[timing.interventions] <- 1
  
  #** make month continuous to smooth out prediction
  month <- rep(seq(1,12,length=30*12),3)[1:length(tvals)]
  
  # storage
  pred <- se.pred <- array(dim=c(length(sites),length(ages),max(timing.interventions)))
  
  #Loop through each site
  for (s in 1:length(sites)) {
    #Loop through each age
    for (a in 1:length(ages)) {
      #Loop through each day of the segments of the intervention sequence
      #E.g., if interventions on day 1, day 50, and day 100 then
      #   run from day 1 to day 50, restart on day 51 and run to day 100, restart on day 101... then bind together
      for (j in 2:length(timing.interventions)){
        #Offsetting by 1 in indices to deal with start day
        interval.days.indices <- tvals[timing.interventions[j-1]:timing.interventions[j]]-tvals[timing.interventions[j-1]]+1
        interval.days <- my.interval[interval.days.indices]
        #Passing the interval segment to date wrangler function to get days per month
        days.in.M.l <- f.date.wrangler.v3(interval.days)
        #Making a temp data frame with the ages, days in months, site code, etc to pass to predict function
        newData <- data.frame(age_yrs = ages[a],
                              ##* Is this mid.month used?
                              mid.month=month[timing.interventions[j-1]:timing.interventions[j]],
                              site_code=sites[s],
                              days.in.MO1=days.in.M.l[[1]],
                              days.in.MO2=days.in.M.l[[2]],
                              days.in.MO3=days.in.M.l[[3]],
                              days.in.MO4=days.in.M.l[[4]],
                              days.in.MO5=days.in.M.l[[5]],
                              days.in.MO6=days.in.M.l[[6]],
                              days.in.MO7=days.in.M.l[[7]],
                              days.in.MO8=days.in.M.l[[8]],
                              days.in.MO9=days.in.M.l[[9]],
                              days.in.M10=days.in.M.l[[10]],
                              days.in.M11=days.in.M.l[[11]],
                              days.in.M12=days.in.M.l[[12]],
                              household=as.factor("S1.37"), #blandest available household, random effect close to 0
                              interlude.days=tvals[timing.interventions[j-1]:timing.interventions[j]]-tvals[timing.interventions[j-1]]+1)
        
        tmp <- predict(fitted.model, newdata=newData, se=TRUE, type="response")
        pred[s,a,timing.interventions[j-1]:timing.interventions[j]] <- tmp$fit
        se.pred[s,a,timing.interventions[j-1]:timing.interventions[j]] <- tmp$se.fit
      }
    }
  }
  
  return(list(pred=pred, se.pred=se.pred, timing.interventions=timing.interventions, month=month, ages=ages))

}

###############################################################################################
## Fourth, a function to take a prediction output and make a plottable data frame

#Function 4
f.plot.tpp <- function(my.start.date.mo,
                       my.start.date.day,
                       my.timing.interventions,
                       sites.to.plot,
                       my.ages){
  #Params
  # my.start.date.mo        : an integer 1:12
  # my.start.date.day       : an integer 1:31 (or last day of that month)
  # my.timing.interventions : a sequence of counting since start date
  #                           Note: Must start with 1
  #                           E.g., c(1, 90) or c(1, 30, 60, 90)
  # my.ages                 : sequence of integers: ages in years for which to simulate
  #                           (e.g, 1:20 or c(2, 5, 15, 25))
  # sites.to.plot           : the subset of sites to plot
  #                           (e.g, c("S1", "S2"))
  
  #Calling the prediction function from above
  tpp <- findYearTimeCourse.v3(
    fitted.model = fit3,
    start.date.mo = my.start.date.mo,
    start.date.day = my.start.date.day,
    timing.interventions = my.timing.interventions,
    my.ages = my.ages,
    scale.interventions=NULL
  )
  
  sites <- c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5")
  sites.indices <- which(sites %in% sites.to.plot)
  
  #loop-loop then rowbind for each site and each age
  
  #tpp structure: tpp[site, age, day]
  
  #creating an empty list to hold temp dfs for each site
  l.df <- list(NA)
  
  for(i in 1:length(sites.indices)){
    
    #creating an empty list to hold temp dfs for each age
    l.df.i <- list(NA)
    
    for(j in 1:length(my.ages)){
      mydf.j <- tibble(day       = 1:length(tpp$pred[sites.indices[i], j, ]),
                       site_code = sites[sites.indices[i]],
                       age       = my.ages[j],
                       prob.inf  = tpp$pred[sites.indices[i], j, ])
      l.df.i[[j]] <- mydf.j
      df.i <- bind_rows(l.df.i)
    }
    l.df[[i]] <- df.i
  }

  df.p <- bind_rows(l.df)
  
  return(df.p)
}



###############################################################################################
## E. FROM REGRESSIONS START SIMULATING INFECTION PROBS UNDER  DISRUPTIONS ####################
###############################################################################################

#Creating data frames, exporting as csv files for downstream scripts that plot

## Figure S2: FOI over age ####################################################################

#FOI for the 12 months
df.1.01 <- f.plot.tpp(my.start.date.mo =  1, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.02 <- f.plot.tpp(my.start.date.mo =  2, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.03 <- f.plot.tpp(my.start.date.mo =  3, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.04 <- f.plot.tpp(my.start.date.mo =  4, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.05 <- f.plot.tpp(my.start.date.mo =  5, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.06 <- f.plot.tpp(my.start.date.mo =  6, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.07 <- f.plot.tpp(my.start.date.mo =  7, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.08 <- f.plot.tpp(my.start.date.mo =  8, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.09 <- f.plot.tpp(my.start.date.mo =  9, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.10 <- f.plot.tpp(my.start.date.mo = 10, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.11 <- f.plot.tpp(my.start.date.mo = 11, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.1.12 <- f.plot.tpp(my.start.date.mo = 12, my.start.date.day = 1, my.timing.interventions = c(1, 1), my.ages = 1:99, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))

#Adding a column for each month and binding the 12 months together
df.1.01 <- df.1.01 %>% mutate(month = "JAN")
df.1.02 <- df.1.02 %>% mutate(month = "FEB")
df.1.03 <- df.1.03 %>% mutate(month = "MAR")
df.1.04 <- df.1.04 %>% mutate(month = "APR")
df.1.05 <- df.1.05 %>% mutate(month = "MAY")
df.1.06 <- df.1.06 %>% mutate(month = "JUN")
df.1.07 <- df.1.07 %>% mutate(month = "JUL")
df.1.08 <- df.1.08 %>% mutate(month = "AUG")
df.1.09 <- df.1.09 %>% mutate(month = "SEP")
df.1.10 <- df.1.10 %>% mutate(month = "OCT")
df.1.11 <- df.1.11 %>% mutate(month = "NOV")
df.1.12 <- df.1.12 %>% mutate(month = "DEC")

df.foi.by.age <- rbind(df.1.01, df.1.02, df.1.03, df.1.04, df.1.05, df.1.06, df.1.07, df.1.08, df.1.09, df.1.10, df.1.11, df.1.12) %>% 
  mutate(month = factor(month, levels = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"))) %>%
  mutate(site_code = factor(site_code, levels = c("S1", "N2", "N3", "N1", "N4", "S2", "S3", "S5", "N5", "S6"))) %>%
  mutate(age_cat = case_when(
    age >   0 & age <  6 ~ "Young children",
    age >=  6 & age < 14 ~ "School aged children",
    age >= 14 & age < 20 ~ "Teenagers",
    age >= 20            ~ "Adults")) %>%
  mutate(age_cat = factor(age_cat, levels = c("Young children", "School aged children", "Teenagers", "Adults"))) %>%
  mutate(code.new = case_when(
    site_code == "S1" ~ "MNJ.06",
    site_code == "S2" ~ "MNJ.07",
    site_code == "S3" ~ "MNJ.08",
    site_code == "S5" ~ "MNJ.09",
    site_code == "S6" ~ "MNJ.10",
    site_code == "N1" ~ "MNJ.05",
    site_code == "N2" ~ "MNJ.04",
    site_code == "N3" ~ "MNJ.03",
    site_code == "N4" ~ "MNJ.02",
    site_code == "N5" ~ "MNJ.01")) %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10")))

#Exporting for use in plotting
# write_csv(df.foi.by.age, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_by_age.csv")



## Figure 2C: FOI per month by age group and site ####################################################################

#FOI for 30 days for the 12 months
df.2.01 <- f.plot.tpp(my.start.date.mo =  1, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.02 <- f.plot.tpp(my.start.date.mo =  2, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.03 <- f.plot.tpp(my.start.date.mo =  3, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.04 <- f.plot.tpp(my.start.date.mo =  4, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.05 <- f.plot.tpp(my.start.date.mo =  5, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.06 <- f.plot.tpp(my.start.date.mo =  6, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.07 <- f.plot.tpp(my.start.date.mo =  7, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.08 <- f.plot.tpp(my.start.date.mo =  8, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.09 <- f.plot.tpp(my.start.date.mo =  9, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.10 <- f.plot.tpp(my.start.date.mo = 10, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.11 <- f.plot.tpp(my.start.date.mo = 11, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.2.12 <- f.plot.tpp(my.start.date.mo = 12, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))

#binding the 12 months together
df.2.01 <- df.2.01 %>% mutate(month = "JAN")
df.2.02 <- df.2.02 %>% mutate(month = "FEB")
df.2.03 <- df.2.03 %>% mutate(month = "MAR")
df.2.04 <- df.2.04 %>% mutate(month = "APR")
df.2.05 <- df.2.05 %>% mutate(month = "MAY")
df.2.06 <- df.2.06 %>% mutate(month = "JUN")
df.2.07 <- df.2.07 %>% mutate(month = "JUL")
df.2.08 <- df.2.08 %>% mutate(month = "AUG")
df.2.09 <- df.2.09 %>% mutate(month = "SEP")
df.2.10 <- df.2.10 %>% mutate(month = "OCT")
df.2.11 <- df.2.11 %>% mutate(month = "NOV")
df.2.12 <- df.2.12 %>% mutate(month = "DEC")

df.monthly.day.foi <- rbind(df.2.01, df.2.02, df.2.03, df.2.04, df.2.05, df.2.06, df.2.07, df.2.08, df.2.09, df.2.10, df.2.11, df.2.12) %>% 
  mutate(month = factor(month, levels = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"))) %>%
  mutate(site_code = factor(site_code, levels = c("S1", "N2", "N3", "N1", "N4", "S2", "S3", "S5", "N5", "S6"))) %>%
  mutate(age_cat = case_when(
    age >   0 & age <  6 ~ "Young children",
    age >=  6 & age < 14 ~ "School aged children",
    age >= 14 & age < 20 ~ "Teenagers",
    age >= 20            ~ "Adults")) %>%
  mutate(age_cat = factor(age_cat, levels = c("Young children", "School aged children", "Teenagers", "Adults"))) %>%
  mutate(code.new = case_when(
    site_code == "S1" ~ "MNJ.06",
    site_code == "S2" ~ "MNJ.07",
    site_code == "S3" ~ "MNJ.08",
    site_code == "S5" ~ "MNJ.09",
    site_code == "S6" ~ "MNJ.10",
    site_code == "N1" ~ "MNJ.05",
    site_code == "N2" ~ "MNJ.04",
    site_code == "N3" ~ "MNJ.03",
    site_code == "N4" ~ "MNJ.02",
    site_code == "N5" ~ "MNJ.01")) %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10")))

#Exporting for use in plotting
# write_csv(df.monthly.day.foi, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_monthly.csv")




## Figure 2D: FOI for 90 days following two storms ####################################################################

#Cyclone Batsirai (BATS):  5 February
#Cyclone Freddy (FRED:    21 February
df.BATS <- f.plot.tpp(my.start.date.mo =  2, my.start.date.day =  5, my.timing.interventions = c(1, 90), my.ages = 1:50, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
df.FRED <- f.plot.tpp(my.start.date.mo =  2, my.start.date.day = 21, my.timing.interventions = c(1, 90), my.ages = 1:50, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))

#binding the 12 months together
df.BATS <- df.BATS %>% mutate(storm = "CYCLONE BATSIRAI")
df.FRED <- df.FRED %>% mutate(storm = "CYCLONE FREDDY")

df.storm.day.foi <- rbind(df.BATS, df.FRED) %>% 
  mutate(storm = factor(storm, levels = c("CYCLONE BATSIRAI", "CYCLONE FREDDY"))) %>%
  mutate(site_code = factor(site_code, levels = c("S1", "N2", "N3", "N1", "N4", "S2", "S3", "S5", "N5", "S6"))) %>%
  mutate(age_cat = case_when(
    age >   0 & age <  6 ~ "Young children",
    age >=  6 & age < 14 ~ "School aged children",
    age >= 14 & age < 20 ~ "Teenagers",
    age >= 20            ~ "Adults")) %>%
  mutate(age_cat = factor(age_cat, levels = c("Young children", "School aged children", "Teenagers", "Adults"))) %>%
  mutate(code.new = case_when(
    site_code == "S1" ~ "MNJ.06",
    site_code == "S2" ~ "MNJ.07",
    site_code == "S3" ~ "MNJ.08",
    site_code == "S5" ~ "MNJ.09",
    site_code == "S6" ~ "MNJ.10",
    site_code == "N1" ~ "MNJ.05",
    site_code == "N2" ~ "MNJ.04",
    site_code == "N3" ~ "MNJ.03",
    site_code == "N4" ~ "MNJ.02",
    site_code == "N5" ~ "MNJ.01")) %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10")))

#Exporting for use in plotting
write_csv(df.storm.day.foi, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_post_storm.csv")






##############################################################################################################################
#FOR DECAY FIGURE
##############################################################################################################################
#FOI for 2 years for an exemplar start date (Jan 1) and exemplar age (age = 1 y)
df.foi.decay <- f.plot.tpp(my.start.date.mo =  1, my.start.date.day = 1, my.timing.interventions = c(1, 365*2), my.ages = 1, sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5"))
#Exporting for use in plotting
#write_csv(df.foi.decay, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/foi_decay_exemplar.csv")




##############################################################################################################################
#FOI for 20 days for March 1st
df.MAR <- f.plot.tpp(my.start.date.mo =  3, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13), sites.to.plot = c("S5", "S6", "N5")) %>%
  filter(day %in% c(15, 20, 30))

df.APR <- f.plot.tpp(my.start.date.mo =  4, my.start.date.day = 1, my.timing.interventions = c(1, 30), my.ages = c(5, 13), sites.to.plot = c("S5", "S6", "N5")) %>%
  filter(day %in% c(15, 20, 30))



## Figure 3C: Plotting return times ####################################################################

## Scenario: March 1st for 254 days (sufficiently long time series)
df.RT.plotter <- f.plot.tpp(my.start.date.mo =  3, my.start.date.day = 1, my.timing.interventions = c(1, 254),  my.ages = c(5, 13, 18, 27), sites.to.plot = c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5")) %>%
  filter(day == 1) %>%
  mutate(code.new = case_when(
    site_code == "S1" ~ "MNJ.06",
    site_code == "S2" ~ "MNJ.07",
    site_code == "S3" ~ "MNJ.08",
    site_code == "S5" ~ "MNJ.09",
    site_code == "S6" ~ "MNJ.10",
    site_code == "N1" ~ "MNJ.05",
    site_code == "N2" ~ "MNJ.04",
    site_code == "N3" ~ "MNJ.03",
    site_code == "N4" ~ "MNJ.02",
    site_code == "N5" ~ "MNJ.01"))

#Exporting for use in plotting
# write_csv(df.RT.plotter, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_for_return_time_plotting.csv")








