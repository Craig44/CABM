library(ibm)

ibm_run = extract.run("../ibm/run.log")
ibm_run = extract.run("../ibm/check1.log")

ibm_run$model_attributes$Recruitment

names(ibm_run)
ibm_run$fisher_age_freq$Values
ibm_run$warnings_encounted

age_length_freq = ibm_run$Trwl_age_length_method$Values

head(age_length_freq)

years = unique(age_length_freq$year)
ages = unique(age_length_freq$age)

#
census_male =  ibm_run$Fishing_Mortality$`age_freq-census-male-FishingTrwl`
census_female =  ibm_run$Fishing_Mortality$`age_freq-census-female-FishingTrwl`
census_total = census_male + census_female
census_total <- sweep(x = census_total, MARGIN = 1, STATS = rowSums(census_total), FUN = "/")

ibm_years = rownames(ibm_run$Fishing_Mortality$`age_freq-male-FishingTrwl`)
ibm_male = ibm_run$Fishing_Mortality$`age_freq-male-FishingTrwl`
ibm_female = ibm_run$Fishing_Mortality$`age_freq-female-FishingTrwl`
## convert to overall props
ibm_total = ibm_female + ibm_male
ibm_total <- sweep(x = ibm_total, MARGIN = 1, STATS = rowSums(ibm_total), FUN = "/")

## compare age-freq from observation class to that from census data stored in the fishing process
## to confirm generate unbiased observations
# - red purple = census info
# - blue is observation from sub-sampling catches and using age-length matrix to scale up to catch.
for(i in 1:length(years)) {
  fish_year_ndx = ibm_years == years[i]
  this_age_year = age_length_freq[age_length_freq$year == years[i],]
  age_length_prop = this_age_year$expected / sum(this_age_year$expected)
  plot(ages, age_length_prop, type = "l", lwd = 2, xlab = "age", ylab = "Proportion", ylim = c(0,0.2), col = "blue", main = years[i])
  lines(ages, ibm_total[fish_year_ndx,], lwd = 3, col = "red")
  lines(ages, census_total[fish_year_ndx,], lwd = 2, col = "purple", lty = 2)
  
  #Sys.sleep(3)
}


#################################
## compare male and female results
##################################
age_length_freq = ibm_run$Trwl_age_length_method$Values
age_length_freq_male = ibm_run$Trwl_age_length_method_male$Values
age_length_freq_female = ibm_run$Trwl_age_length_method_female$Values

head(age_length_freq)
years = unique(age_length_freq$year)
ages = unique(age_length_freq$age)

ibm_years = rownames(ibm_run$Fishing_Mortality$`age_freq-male-FishingTrwl`)
ibm_male = ibm_run$Fishing_Mortality$`age_freq-male-FishingTrwl`
ibm_female = ibm_run$Fishing_Mortality$`age_freq-female-FishingTrwl`
ibm_male <- sweep(x = ibm_male, MARGIN = 1, STATS = rowSums(ibm_male), FUN = "/")
ibm_female <- sweep(x = ibm_female, MARGIN = 1, STATS = rowSums(ibm_female), FUN = "/")

census_female =  ibm_run$Fishing_Mortality$`age_freq-census-female-FishingTrwl`
census_male =  ibm_run$Fishing_Mortality$`age_freq-census-male-FishingTrwl`
census_male <- sweep(x = census_male, MARGIN = 1, STATS = rowSums(census_male), FUN = "/")
census_female <- sweep(x = census_female, MARGIN = 1, STATS = rowSums(census_female), FUN = "/")

## males
for(i in 1:length(years)) {
  fish_year_ndx = ibm_years == years[i]
  this_age_year = age_length_freq_male[age_length_freq_male$year == years[i],]
  age_length_prop = this_age_year$expected / sum(this_age_year$expected)
  plot(ages, age_length_prop, type = "l", lwd = 2, xlab = "age", ylab = "Proportion", ylim = c(0,0.14), col = "blue", main = years[i])
  lines(ages, ibm_male[fish_year_ndx,], lwd = 3, col = "red")
  lines(ages, census_male[fish_year_ndx,], lwd = 2, col = "purple", lty = 2)
  Sys.sleep(1)
}

## females
for(i in 1:length(years)) {
  fish_year_ndx = ibm_years == years[i]
  this_age_year = age_length_freq_female[age_length_freq_female$year == years[i],]
  age_length_prop = this_age_year$expected / sum(this_age_year$expected)
  plot(ages, age_length_prop, type = "l", lwd = 2, xlab = "age", ylab = "Proportion", ylim = c(0,0.14), col = "blue", main = years[i])
  lines(ages, ibm_female[fish_year_ndx,], lwd = 3, col = "red")
  lines(ages, census_female[fish_year_ndx,], lwd = 2, col = "purple", lty = 2)
  Sys.sleep(3)
}

## Male and female should be different
# - different selectivities.
# - different recuritment
# - different growth etc
for(i in 1:length(years)) {
  fish_year_ndx = ibm_years == years[i]
  this_age_year = age_length_freq_female[age_length_freq_female$year == years[i],]
  this_male_age_year = age_length_freq_male[age_length_freq_male$year == years[i],]
  
  age_length_prop = this_age_year$expected / sum(this_age_year$expected)
  male_age_length_prop = this_male_age_year$expected / sum(this_male_age_year$expected)
  
  plot(ages, age_length_prop, type = "l", lwd = 2, xlab = "age", ylab = "Proportion", ylim = c(0,0.14), col = "blue", main = years[i])
  lines(ages, male_age_length_prop, lwd = 3, col = "red")
  legend('topright', legend = c("Male", "Female"), col = c("red","blue"), lwd = 2)
  Sys.sleep(2)
}



