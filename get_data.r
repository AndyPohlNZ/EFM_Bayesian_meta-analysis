#########################################################
## STAT 619 Term Project - Data Processing
## Andy Pohl - 30058005
## University of Calgary
## March 2017
#########################################################


#### Functions to process data ####
RD = function(r1,r0, n1,n0){
  # Given incidence and number calculate RD - the difference
  # in death probability per 1000 births
  return(((r1/n1) - (r0/n0)) * 1000)
}
var_RD = function(r1,r0,n1,n0){
  # Given incidence and number calculate the variance in the 
  # RD statistic with a continuity correction
  if(r1 == 0 | r0 == 0){
    r1 = r1+0.5
    r0 = r0+0.5
    n1 = n1 + 1
    n0 = n0 +1
  }
  
  vrd = ((((r1/n1)*(1-(r1/n1)))/n1) + (((r0/n0)*(1-(r0/n0)))/n0)) * 10^6
  return(vrd)
}

emperical_logit = function(r1, r0, n1,n0, continuity_correction = TRUE){
  if(continuity_correction){
    r1 = r1+0.5; r0 = r0+0.5; n1 = n1+1; n0 = n0+1
  }
  el = log((r1/(n1-r1))) - log((r0/(n0-r0)))
  return(el)
}

var_emperical_logit = function(r1, r0, n1, n0, continutiy_correction = TRUE){
  if(continutiy_correction){
    r1 = r1+0.5; r0 = r0+0.5; n1 = n1+1; n0 = n0+1
  }
  var_el = (1/r1) + (1/(n1-r1)) + (1/r0)  + (1/(n0-r0))
  return(var_el)
}

extend_vec <- function(vec, n, value){
  # extends a vector vec to length n but repeting 'value' n-length(vec) times
  output <- c(vec, rep(value, n-length(vec)))
  return(output)
} 

#### Data ####

efm_data <- data.frame(
  study_type = c(rep(1, 9), rep(2, 7), rep(3, 10)),
  study_id = c(1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,1,2,3,4,5,6,7,8,9,10),
  year = c(1976,1976,1978,1979,1981,1985,1985,1987,1993,1973,1973,1975,
           1977,1978,1979,1982,1975,1975,1975,1976,1977,1978,1980,1980,1984,1986),
  n_treatment = c(175,242,253,463,445,485,6530,122,746, 1162,150,608,4210,554,4978,
                  45880,1024,1080,1950,3529,3852,7312,4503,8174,7911,17586),
  n_control = c(175,241,251,232,482,493,6554,124,682, 5427,6836,6179,2923,692,8634,
                66208,991,1161,11599,4323,4114,15357,4240,6740,7582,17409),
  death_treatment = c(1,2,0,3,1,0,14,17,2,2,0,1,1,1,0,10,0,9,1,1,21,6,2,5,2,5),
  death_control = c(1,1,1,0,0,1,14,18,9,17,15,37,9,3,2,45,4,7,14,15,53,35,19,15,13,7)
)

efm_data$n = efm_data$n_control + efm_data$n_treatment
efm_data$RD <- apply(efm_data, MARGIN = 1, function(x){RD(x[6], x[7], x[4], x[5])})
efm_data$var_RD <- apply(efm_data, MARGIN = 1, function(x){var_RD(x[6], x[7], x[4], x[5])})
efm_data$year_adjusted <- efm_data$year - 1979

