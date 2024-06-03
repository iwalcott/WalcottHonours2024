##############################################################################
#                Decline, expansion, stable output data analysis
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

###set results folder
res_folder1 <- "/data/scratch/isobel/results_gadi/"
res_folder2 <- "/data/scratch/isobel/results2/"

### list files in each folder
fn <- list.files(res_folder1, pattern = "*.rds", full.names = TRUE)
fn2 <- list.files(res_folder2, pattern = "*.rds", full.names = TRUE)

### Load data for each fn list
dat_list = sapply(fn2, function (x) data.table(readRDS(x)))

## Convert list names to run numbers ##
names(dat_list) <- unlist(lapply(names(dat_list), function(x) substr(x, 32, 45)))

#==============================================================================#
#                      Load logged simulated data
#==============================================================================#

###set logfile folder
log_folder <- "/data/scratch/isobel/log_std_set2/"

## Generate list of potential files from folder
sim_list <- list.files(log_folder, pattern = "*.txt", full.names = TRUE)

##Generate list of simulated data
sim_files <- mclapply(names(dat_list), function(x) {
  x <- read.delim2(grep(paste0("Run2_", substr(x, 6, 10)), sim_list, value = T), header = T, sep = ",")}, 
  mc.cores = 10)

sim_files
names(sim_files) <- names(dat_list)

###Add simulated data to main data list
for (i in 1:length(dat_list)) {
  dat_list[i][[1]][4] <- sim_files[i]
  dat_list[i][[1]][[4]]$year <- 2201 - dat_list[i][[1]][[4]]$cycle
}

#==============================================================================#
#                               Tidy data
#==============================================================================#

###subset trimmed dataset to only necessary runs
data_trimmed <- trim(dat_list, ncores = 10)

##convert data to one long dataframe
data_long2 <- par_data_long(data_trimmed, ncores = 10)

###add set number column to dataframe
data_long2$set <- "set2"

###allocate new name to df
data_long1

#combine datasets
data_long <- rbind(data_long2, data_long1)

method_list <- split(data_long, data_long$method)
sim_data <- method_list[3]$sim
method_list <- method_list[c(1,2,4)]

# Assuming combinations is your dataframe

#Get combinations from std_simulation set
combinations <- unique(df[, c(1:4, 6)])

combinations <- combinations %>%
  mutate(model_code = recode(model, "decline" = 1, "expansion" = 2, "stable" = 3))

#==============================================================================#
#        Add extracted total loci information and apply correction factor 
#==============================================================================#

### load loci data
loci_totals <- read_rds("~/R/loci_totals.rds")

### add run_code column for joining data frames
loci_totals$runnum <- rownames(loci_totals)
loci_totals$run_code <- substr(loci_totals$runnum, 5, 9)
colnames(loci_totals)[1] <- "nloc"

### add run_code column
data_long <- data_long %>% mutate(run_code = if_else(set == "set1", substr(runnum, 5, 9), substr(runnum, 6, 10)))

### add loci totals into data frame
new_data <- left_join(data_long, loci_totals[, c(1,3)], by = "run_code")

### calculate correction factor
new_data$c_factor <- new_data$nloc/(new_data$loci * 1000)

### recalculate year and Ne
new_data <- new_data %>% 
  mutate("c_year" = if_else(method %in% c("epos", "stairway"), year * c_factor, year),
         "c_Ne" = if_else(method %in% c("epos", "stairway"), Ne * c_factor, Ne))

### new_data saved in this format
write_rds(new_data, file = "/data/scratch/isobel/new_data_updated.rds")

#create new runnum with set number
new_data$setnum <- paste0(new_data$runnum, new_data$set)

#==============================================================================#
#                Interpolate data for analysis efficiency
#==============================================================================#

### separate into list by method and run
interp_list <- split(new_data, list(new_data$setnum, new_data$method))

### extract short runs
short_runs <- Filter(function(x) max(x$c_year) < 225, interp_list)

### Filter analysis list to runs >= 225
interp_list2 <- Filter(function(x) max(x$c_year) >= 225, interp_list)

### interpolate for each year/generation
library(parallel)
interp_out <-mclapply(1:length(interp_list2), function(i) {
  x <- interp_list2[[i]]
  ####identify max year - interpolate to that value
  maxint <- ifelse(max(x$c_year) > 250, 250, max(x$c_year)) 
  
  df <- data.frame(approx(x$c_year, x$c_Ne, xout = 1:maxint))
  colnames(df) <- c("year", "Ne")
  df[, 3:12] <- x[2, c(1, 4:11, 18)]
  return(df);
  },
  mc.cores = 10)

names(interp_out) <- names(interp_list2)

#==============================================================================#
#                Add simulated data to each run
#==============================================================================#

source("~/R/geohippos/add_sim_data.R")
interp_simNe <- mclapply(interp_out, function(X) {
  popinit <- X$pop_init[1]
  model <- X$model[1]
  tl <- X$tl[1]
  ts <- X$ts[1]
  cp <- as.numeric(X$cp[1])
  for (i in 1:nrow(X)) {
    X$simNe[i] <- model_plot(model, popinit, cp, tl, ts, X$year[i])
  }
  X$medNe <- min(X$simNe) + (max(X$simNe) - min(X$simNe))/2
  X$finalNe <- X$simNe[1]
  return(X)
}, mc.cores = 10)

names(interp_simNe) <- names(interp_out)

#==============================================================================#
#                Calculate relative absolute error at each year
#==============================================================================#

### calculate relative absolute error
error_abs_std <- mclapply(interp_simNe, function(X) {
  df <- X[1, 3:11]
  
  X$error_sq <- abs(X$simNe - X$Ne)/X$simNe
  return(X)
}, mc.cores = 10)

### combine into one dataframe
error_boxplot_df <- do.call("rbind", error_abs_std) #%>% filter(year < 20)

### calculate mean and standard error of relative absolute error
error_meanSE_std <- error_boxplot_df %>% 
  group_by(method,model, year, loci) %>% 
  summarise(mean = mean(error_sq, na.rm = T),
            count = sum(!is.na(error_sq)),
            SE = sd(error_sq, na.rm = T)/sqrt(sum(!is.na(error_sq))))

#==============================================================================#
#            Restrict to >15ybp and calculate root mean square error
#==============================================================================#

### calculate root mean square error across analysis range
error_rms <- mclapply(interp_simNe, function(X) {
  df <- X[1, 3:11]
  ###subset for ybp in 15:250
  X <- subset(X, X$year %in% 15:250)
  n <- nrow(X)
  X$error_sq <- ((X$simNe - X$Ne)^2)/n
  rmse <- X %>% summarise(RMSE = sqrt(sum(error_sq)))
  df$rmse <- unlist(rmse)
  return(df)
}, mc.cores = 20)  

### combine into one data frame
rmse_df <- do.call("rbind", error_rms)

#==============================================================================#
#                   Trajectory detection analysis
#==============================================================================#

### load data
interp_out <- read_rds("/data/scratch/isobel/interp_out_updated.rds")

###create warning when regression is near perfect fit
f <- function(expr) {
  tryCatch(expr, 
           warning = function(w) {
             if (grepl('perfect fit', w))
               return(1) 
             else return(w)
           })  
}

### calculate regression statistics for each run
library(parallel)
model_detect <- mclapply(interp_out, function(X) {
  df <- X[1, 3:ncol(X)]
  df$cp <- as.numeric(df$cp)

  ### calculate regression for entire range
  reg_all <- lm(Ne ~ year, X)
  err <- ifelse(is.numeric(f(summary(reg_all))), "error", "fine")
  
  ### extract slope mean and SE values
  mu <- summary(reg_all)$coefficients[2,1]
  se <- summary(reg_all)$coefficients[2,2]
  low_ci <- mu - 1.96*se
  upp_ci <- mu + 1.96*se
  resid <- summary(reg_all)$sigma
  p <- summary(reg_all)$coefficients[1,4]
  sim_df <- Filter(function(a) (a$runnum[1] == df$runnum) & (a$method[1] == "sim"), interp_out)[[1]]
  
  ###alternative: extract final Ne, compare to historic Ne
  final <- subset(X, X$year == 15)$Ne
  historic <- subset(X, X$year == max(X$year))$Ne
  ts_Ne <- subset(X, X$year == df$ts)$Ne
  
  ###extract target Ne values
  t_final <- subset(sim_df, sim_df$year == 15)$Ne
  t_historic <- subset(sim_df, sim_df$year == max(X$year))$Ne
  
  df <- cbind(df, int_all = reg_all$coefficients[1], slope_all = reg_all$coefficients[2], final, historic, ts_Ne, t_final, t_historic, low_ci, upp_ci, resid, err, p)

  return(df)
  ###declines have positive slope, expansions have negative slope, stable have 0
}, mc.cores = 10)

model_detect_df <- do.call("rbind", model_detect)


### Calculate rates of correct trajectory detection
options(scipen = 900)
detect_out <- mclapply(model_detect, function(X) {
  df <- X
  model <- df$model
  
  ###identify models detected from regression slope
  ###minimum slope value found in simulated cases was ~0.102
  slope_model_all <- ifelse(
    ###stable case
    (abs(df$slope_all) < 0.05) & (df$p < 0.05), 3,
    ifelse(
      #Decline case - slope above threshold and CI doesn't include zero
      (df$slope_all > 0.05) & (df$p < 0.05) & !(0 %in% df$low_ci:df$upp_ci), 1,
      ifelse(
        #Expansion case - slope below threshold and CI doesn't include zero
        (df$slope_all < 0.05) & (df$p < 0.05) & !(0 %in% df$low_ci:df$upp_ci), 2, NA)))
  
  # ####check other stable tests
  # p_test <- ifelse(df$p < 0.05, 1, 0)
  # resids <- ifelse(df$resid < 5, 1, 0)
  # Ne_range <- abs((df$historic - df$final)/df$historic)
  # df$model_detected <- slope_model_all
  
  ### test if detected trajectory matches simulated trajectory
  df$overall <- ifelse(df$model == slope_model_all, 1, 0)
  
  return(df)
}, mc.cores = 20)

options(scipen = 0)


detect_validate_df <- do.call("rbind", detect_out)



####group and summarise to get count
detect_summarised <- detect_validate_df %>% group_by(model, loci, ss, method) %>% 
  summarise(mean_prop = mean(overall, na.rm = T),
            count_sample = n())



#==============================================================================#
#                Relative error absolute distributions at key points
#==============================================================================#

### gather data to plot final and historic on same plot
model_Ne_err <- model_detect_df %>% mutate("abs_final" = abs(t_final - final)/t_final,
                           "abs_historic" = abs(t_historic - historic)/t_historic) %>% 
  gather(key = "Ne_yr", value = "abs_Ne_yr", abs_final, abs_historic)

### Violin Plot
ggplot(data = model_Ne_err %>% filter(Ne_yr == "abs_historic", method %in% c("epos", "stairway")), aes(x = as.factor(loci), y = abs_Ne_yr, colour = as.factor(method), fill = as.factor(method))) +
  geom_violin() +
  labs(y = "Relative absolute error - Historic Ne",
       x = "Loci (10^3)") +
  facet_grid(ss ~ model, labeller = labeller(model = model.labs)) +
  theme_bw() +
  ylim(c(0, 2)) +
  theme(legend.position = "bottom") +
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Method", end = 0.8, labels = method.labs) +
  scale_colour_viridis(discrete = TRUE, option = "D", name = "Method", end = 0.8, labels = method.labs)


