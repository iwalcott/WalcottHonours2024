##############################################################################
#                      Bottleneck output data analysis
##############################################################################

### load packages
library(geohippos)
library(data.table)
library(viridis)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

#==============================================================================#
#              Load demographic inference output data
#==============================================================================#

### set results folders
res_folder1 <- "/data/scratch/isobel/bottle_res1/"
res_folder2 <- "/data/scratch/isobel/bottle_res2/"
res_folder3 <- "/data/scratch/isobel/bottle_res3/"

### create list of files in each folder
fn <- list.files(res_folder1, pattern = "*.rds", full.names = TRUE)
fn2 <- list.files(res_folder2, pattern = "*.rds", full.names = TRUE)
fn3 <- list.files(res_folder3, pattern = "*.rds", full.names = TRUE)

### Load data for sample
dat_list = sapply(fn3, function (x) data.table(readRDS(x)))

### Convert list names to run numbers ##
names(dat_list) <- unlist(lapply(names(dat_list), function(x) substr(x, 42, 54)))

#==============================================================================#
#                      Load logged simulated data
#==============================================================================#

### set logfile folder
log_folder <- "/data/scratch/isobel/bottle_log/"

### Generate list of potential files from folder
sim_list <- list.files(log_folder, pattern = "*.txt", full.names = TRUE)

### Generate list of simulated data
sim_files <- lapply(sim_list, function(x) read.delim2(x, header = T, sep = ","))

### extract model and pop_init for runname
names(sim_files) <- lapply(sim_files, function(x) {
  name_list <- paste0(x$model[1], x$pop_init[1])
  return(name_list)
})

#==============================================================================#
#              Add loci subsampling metadata
#==============================================================================#

### import run_list dataframe from sim script and change runnum name
bottledf <- bottledf %>% rename(runnum = runnumb)

### add subsampled loci data to dataframe
df <-  left_join(bottledf, df[, c(1:2, 15)], by = c("model", "pop_init"), multiple = "first")
df <- df %>% rename(loci_total = loci, runnum = runnumb)

### Add simulated data to main data list
for (i in 1:length(dat_list)) {
  id <- grep(substr(names(dat_list)[i], 1, 9), df$runnum)
  code <- paste0(df$model[id], df$pop_init[id])
  
  dat_list[i][[1]][4] <- sim_files[code]
  dat_list[i][[1]][[4]]$year <- 2301 - dat_list[i][[1]][[4]]$cycle
}

#==============================================================================#
#                              Tidy Data
#==============================================================================#

### Trim data to analysis range
data_trimmed <- trim(dat_list, ncores = 2, trim_pt = 350)

### Convert from list to one long dataframe
data_long <- par_bottle_long(dat_list, ncores = 5)

### add runnum column without loci
data_long$runnum2 <- substr(data_long$runnum, 1, 9)

### add sample size column by runnum
dl3 <- left_join(data_long, df[, c(12:13, 15)], by = c("runnum2" = "runnum"))

### Each set has been run through previous code individually

### Allocate set number
dl1$set <- "set1"
dl2$set <- "set2"
dl3$set <- "set3"

### Combine intoone data frame
data <- rbind(dl1, dl2, dl3)

#==============================================================================#
#               Add total loci data and calculate correction factor
#==============================================================================#

### add loci data frame from analysis functions final.R script
df_loci <- df_loci %>% rename(loci_total = loci)

### Join total loci count data
bottle_data <- left_join(data, df_loci[, 15:16], by = c("runcode" = "run_code"), multiple = "first")

### calculate correction factor
bottle_data$c_factor <- bottle_data$loci_total/(bottle_data$loci * 1000)

### recalculate year and Ne
bottle_data <- bottle_data %>% 
  mutate("c_year" = if_else(method %in% c("epos", "stairway"), year * c_factor, year),
         "c_Ne" = if_else(method %in% c("epos", "stairway"), Ne * c_factor, Ne))

###save bottle dataset with corrected year and Ne values
write_rds(bottle_data, file = "/data/scratch/isobel/bottle_data_corrected.rds")

###add runcode column to data
data$runcode <- paste0(data$model, data$pop_init)
df$runcode <- paste0(df$model, df$pop_init)
data$run_set <- paste0(data$runnum, data$set)

#load bottle data
bottle_data <- read_rds("/data/scratch/isobel/bottle_data_corrected.rds")


#==============================================================================#
#                Interpolate data for analysis efficiency
#==============================================================================#

### split into list by method and run
interp_list <- split(bottle_data, list(bottle_data$run_set, bottle_data$method))

### Filter for short runs
short_bottle <- Filter(function(x) max(x$c_year) < 250, interp_list)
interp_list2 <- Filter(function(x) max(x$c_year) >= 250, interp_list)

### interpolate values at each year/generation
library(parallel)
interp_out <- mclapply(1:length(interp_list2), function(i) {
  x <- interp_list2[[i]]
  maxyr <- ifelse(max(interp_list2[[i]]$c_year) < 350, max(interp_list2[[i]]$c_year), 350)
  df <- data.frame(approx(x$c_year, x$c_Ne, xout = 1:maxyr))
  colnames(df) <- c("year", "Ne")
  df[, 3:9] <- x[2, c(1, 4:7, 9, 11)]
  return(df);
},
mc.cores = 10)

names(interp_out) <- names(interp_list2)

### extract simulation data for comparison
interp_simNe <- mclapply(1:length(interp_out), function(i) {
  num <- grep(paste0(substr(names(interp_out)[i], 1, 18), "sim"), names(interp_out))
  sim_data <- subset(interp_out[[num]], interp_out[[num]]$year <= max(interp_out[[i]]$year))
  interp_out[[i]]$simNe <- sim_data$Ne
  interp_out[[i]]$simrunset <- sim_data$run_set
  return(interp_out[[i]])
}, mc.cores = 10)

names(interp_simNe) <- names(interp_out)

### Remove simulation data frames from list
interp_simNe2 <- Filter(function(x) x$method[1] != 
                          "sim", interp_simNe)


#==============================================================================#
#                Calculate relative absolute error
#==============================================================================#

### calculate relative absolute error at each year
error_abs <- mclapply(interp_simNe2, function(X) {
  df <- X[1, 3:10]
  
  X$error_sq <- abs(X$simNe - X$Ne)#/X$simNe
  return(X)
}, mc.cores = 10)

### convert to data frame
error_boxplot_df <- do.call("rbind", error_abs) #%>% filter(year < 20)

### calculate mean and standard error of relative absolute error
error_meanSE <- error_boxplot_df %>% 
  group_by(method,model, year, ss) %>% 
  summarise(mean = mean(error_sq, na.rm = T),
            count = sum(!is.na(error_sq)),
            SE = sd(error_sq, na.rm = T)/sqrt(sum(!is.na(error_sq))))


#==============================================================================#
#            Restrict to >15ybp and calculate root mean square error
#==============================================================================#

###Calculate root mean square error over desired range for filtered runs
error_rms <- mclapply(interp_simNe2, function(X) {
  df <- X[1, 3:11]
  ###subset for ybp in 15:250
  X <- subset(X, X$year %in% 15:350)
  n <- nrow(X)
  X$error_sq <- ((X$simNe - X$Ne)^2)/n
  rmse <-  X %>% summarise(RMSE = sqrt(sum(error_sq)))
  df$rmse <- unlist(rmse)
  if (n == 0) {
    df$rmse <- NA
  }
  df$n <- n
  return(df)
}, mc.cores = 20)  


##combine output into one dataframe
df_b_rmse <- do.call("rbind", error_rms)

df_b_rmse$run_code <- paste0(df_b_rmse$model, df_b_rmse$pop_init)



#==============================================================================#
#                   Trajectory detection analysis
#==============================================================================#

### load data
interp_simNe2 <- read_rds("/data/scratch/isobel/bottle_interp_simNe2.rds")

### Calculate crash midpoint for each model - use bottledf from simulation script
bottle_df <- bottledf %>% mutate(midpoint = 2301 - (crash_st - (crash_st-recov_st-recov_len)/2))

### perform bottleneck trajectory detection analysis
detect_bottleneck <- lapply(interp_simNe2, function(X) {
  ### subset X to range
  X <- subset(X, X$year >= 15)
  df <- X[1, 3:10]
  
  ### find midpoint value
  mp <- X$model[1]
  p0 <- X$pop_init[1]
  midpoint <- subset(bottle_df, bottle_df$model == mp & bottle_df$pop_init == p0)$midpoint
  
  ### Find minimum value within +/- 50 years of midpoint
  df_range <- subset(X, (X$year <= (midpoint + 50)) & (X$year >= (midpoint - 50)))
  min_Ne_yr <- subset(df_range, df_range$Ne == min(df_range$Ne))$year
  
  ### if multiple minimum points, take earliest one
  if (length(min_Ne_yr > 1) ){
    min_Ne_yr <- min_Ne_yr[length(min_Ne_yr)]
  }
  ###calculate regression statistics
  crash <- subset(X, X$year >= min_Ne_yr)
  recov <- subset(X, X$year <= min_Ne_yr)
  
  lm_c <- lm(Ne~year, data = crash)
  mu_c <- lm_c$coefficients[2]
  p_c <- summary(lm_c)$coefficients[1,4]
  
  lm_r <- lm(Ne ~ year, data = recov)
  mu_r <- lm_r$coefficients[2]
  p_r <- summary(lm_r)$coefficients[1,4]
  
  ### crash detection - positive slope and significant p
  crash_val <- ifelse((mu_c > 1) & (p_c < 0.05), 1, 0)
  
  ### recovery detection - negative slope and significant p
  recov_val <- ifelse((mu_r < -1) & (p_r < 0.05), 1, 0)
  
  ### extract simulated and inferred Ne values at key points
  final <- subset(X, X$year == 15)$Ne
  t_final <- subset(X, X$year == 15)$simNe
  
  historic <- subset(X, X$year == max(X$year))$Ne
  t_historic <- subset(X, X$year == max(X$year))$simNe
  
  ### test bottleneck detection
  overall <- crash_val*recov_val
  ###look for decline/positive slope in crash range
  df <- cbind(df, mu_c, p_c, mu_r, p_r, crash_val, recov_val, midpoint, overall, historic, t_historic, final, t_final)
  return(df)
})

detect_bottle_df <- do.call("rbind", detect_bottleneck)

### summarise bottleneck
detect_bottle_sum <- detect_bottle_df %>% group_by(method, model, loci, ss) %>% 
  summarise(mean_prop = sum(overall, na.rm = T)/n(),
            count_sample = n())

#==============================================================================#
#                Relative error absolute distributions at key points
#==============================================================================#

detect_bottle_df <- read_csv("/data/scratch/isobel/bottle_detection_df (1).csv")

### summarise and calculate relative absolute error
bottle_Ne_err <- detect_bottle_df %>% mutate("abs_final" = abs(t_final - final)/t_final,
                                           "abs_historic" = abs(t_historic - historic)/t_historic) %>% 
  gather(key = "Ne_yr", value = "abs_Ne_yr", abs_final, abs_historic)

### Violin Plot
ggplot(data = bottle_Ne_err %>% filter(Ne_yr == "abs_final", method %in% c("epos", "stairway")), aes(x = as.factor(loci), y = abs_Ne_yr, colour = as.factor(method), fill = as.factor(method))) +
  geom_violin() +
  labs(y = "Relative absolute error - Final Ne",
       x = "Loci (10^3)") +
  facet_grid(ss ~ model, labeller = labeller(model = model.labs)) +
  theme_bw() +
  ylim(c(0,2)) +
  theme(legend.position = "bottom") +
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Method", end = 0.8, labels = method.labs) +
  scale_colour_viridis(discrete = TRUE, option = "D", name = "Method", end = 0.8, labels = method.labs)


