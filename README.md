# boldness-avoidance-foraging
Repository for data and analytical R code for the manuscript entitled "Intricate covariation between exploration and avoidance learning in a generalist predator", currently being considered for publication in Behavioral Ecology  

List of files:

1. exploration.csv: exploration data

**id:** the unique ID of experimental individuals
**sex:** sex of the individual
**population:** population of origin. Shou and Mei: two nonfocal populations in Kaohsiung, KS: the focal population in Kaohsiung, PT: Pingtung
**exp:** exploration trial number, values are 1 (the first trial) or 2 (the second trial)
**1stMove:** moving latency in seconds
**cross:** time of first entering the novel environment, values in seconds
**percent_visible:** precent time being visible 
**moving_bouts:** the number of moving bouts in a trial
**precent_moving:** percent time spent moving
**percent_novel:** percent time in the novel environment
**pc1:** score of the first principal component. See Materials and Methods for more detail
**pc2:** score of the second principal component. See Materials and Methods for more detail
 
2. dff5b.csv: raw data from exploration and avoidance learning trials
**id:** the unique ID of the individual
  **sex:** sex of the individual
  **group:** treatment group in avoidance learing experiments. R: individuals encountered red-bitter and yellow-normal color taste combination. Y:   
             individuals encountered yellow-bitter and red-normal color taste combination
  **exp:** trial number, ranging from 1 to 5. See Materials and Methods for more details
  **exp_date:** date of experiment
  **latency:** foraging latency in seconds
  **1st_prey:** first prey attacked. P: palatable, U: unpalatable
  **unpalatable_50:** how many unpalatable prey attacked among the first 50% of prey attacked
  **total_50:** the number of first 50% prey attacked
  **palatable_50:** how many palatable prey attacked among the first 50% of prey attacked. Equals total_50 - unpalatable_50
  **unpalatable_consumed:** how many unpalatable prey consumed in a trial
  **total_consumed:** total amount of prey consumed in a trial
  **palatable_consumed:** how many palatable prey consumed in a trial. Equals total_consumed - palatable_consumed
  **population:** population of origin. KS: Kaohsiung, PT: Pingtung
  **pc1:** the first principal component score of exploration PCA, a proxy for general exploration tendency. Individuals with smaller numbers are faster              explorers. 
  **1st_prey2:** using 0 and 1 to code the same information as in 1st_prey for statistical purposes
  **percent_unpalatable_50:** the proportion of unpalatable prey among the first 50% prey attacked
  **percent_unpalatable_consumed:** the proportion of unpalatable prey among all prey consumed

3. HighstatLibV10.R: from Zuur et al. (2009) containing function to calculate variance inflation factor values
4. stats.R: R code for all statistical analyses 
