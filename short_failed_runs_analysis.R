###########################################################################
#         Analyse failed and short runs
###########################################################################


##load short std runs
short_runs <- read_rds("/data/scratch/isobel/short_run_list.rds")

short_std_df <- do.call("rbind", short_runs)

##extract unique cases of runs
short_runs_summ <- short_std_df %>% 
  group_by(setnum, method, model, pop_init, tl, ts, ss, cp, loci, set) %>% 
  summarise(max_yr = max(c_year))

##create dataframe of short run counts
short_runs_df <- short_runs_summ %>% group_by(method, model, ss, loci, set) %>% 
  summarise(count = n())

##save
write_csv(short_runs_df, file = "~/R/std_short_runs_table_data.csv")

##summarise short runs by count
ggplot(data = short_runs_summ, aes(x = as.factor(loci), fill = as.factor(ss))) +
  geom_bar(stat = "count", position = position_dodge()) +
  labs(title = "Short runs (<225ybp)",
       y = "Count",
       x = "Loci (10^3)") +
  theme_minimal() +
  facet_grid(method ~ model, labeller = labeller(model = model.labs, method = method.labs)) +
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Sample Size", end = 0.8)

##load failed runs
std_not_run <- read_csv("/data/scratch/isobel/std_notrun.csv")

##count occurrences across parameters
std_nr_sum <- std_not_run %>% 
  group_by(runnumb, model, loci, pop_init, crash_prop, tl, ts, ss) %>% 
  summarise(count = n())

std_nr_sum$run_loci <- paste0(std_nr_sum$runnumb, "_",std_nr_sum$loci)

##create data frame of not run counts
std_nr_df <- std_nr_sum %>% group_by(model, loci, ss) %>% 
  summarise(count = n())

##save
write_csv(std_nr_df, file = "~/R/std_notrun_table_df.csv")


##summarise across loci, sample size and model
df_sum <- df %>% 
  group_by(model, loci, ss) %>% 
  summarise(total_count = n())

##join total and failed summarised data together
std_nr_df$loci <- as.numeric(std_nr_df$loci)

df_failed <- left_join(std_nr_df, std_total_nm, join_by(model, loci, ss)) %>% rename("fail_count" = count.x, "total_count" = count.y)

df_failed$success_count <- df_failed$total_count - df_failed$fail_count



### Create set of all std combinations

#df from sim script
head(df)
df <- crossing(df, loci = c(001, 005, 010, 020, 050, 100))
df <- crossing(df, set = c("set1", "set2"))
df$runset <- paste0(df$runnumb, df$set)

#Save all std cases (before method combinations)
write.csv(df, "/data/scratch/isobel/std_total_df_no_methods.csv")

df_methods <- crossing(df, methods = c("epos", "stairway", "gone"))
write.csv(df_methods, "/data/scratch/isobel/std_total_df_with_methods.csv")

##create total std run summary table (no methods)

std_total_nm <- df %>% group_by(model, loci, ss) %>% 
  summarise(count = n())

write_csv(std_total_nm, file = "~/R/std_total_count_combos_no_methods.csv")



##################Bottleneck Cases######################################
###bottle model/runcode labels
model.labs <- c("Decline", "Expansion", "Stable")
names(model.labs) <- c("croc15000", "croc7500", "frog100", "frog50", "seal2100", "seal4200", "whale6000")


#load data
short_bottle_list <- readRDS("/data/scratch/isobel/short_bottle_runs_250.rds")

#convert to dataframe
short_bottle <- do.call("rbind", short_bottle_list)

##extract unique cases of runs
short_bottle_sum <- short_bottle %>% 
  group_by(run_set, method, model, pop_init, ss, loci, runcode) %>% 
  summarise(max_yr = max(c_year))

##create dataframe of short run counts
short_bottle_runs_df <- short_bottle_sum %>% group_by(method, runcode, ss, loci) %>% 
  summarise(count = n())

##save
write_csv(short_bottle_runs_df, file = "~/R/bottle_short_runs_table_data.csv")

##summarise short runs by count
ggplot(data = short_bottle_sum, aes(x = as.factor(loci), fill = as.factor(ss))) +
  geom_bar(stat = "count", position = position_dodge()) +
  labs(title = "Short bottleneck runs (<250ybp)",
       y = "Count",
       x = "Loci (10^3)") +
  theme_minimal() +
  facet_grid(method ~ runcode, labeller = labeller(method = method.labs)) +
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Sample Size", end = 0.8)

#############################################
#collect all failed bottle data

#reverse calculate from finished runs - start with df pulled from bottleneck script
##create combinations for all loci 
df <- crossing(df, loci = c(001, 005, 010, 020, 050, 100))
df
df$runnum <- paste0(df$runnumb, "_", sprintf("%03d", df$loci))
df <- crossing(df, set = c("set1", "set2", "set3"))

df$runset <- paste0(df$runnum, df$set)

###df now has all possible runset codes represented
notrun_bottle_names <- setdiff(df$runset, bottle_data$run_set)

##extract metadata 
notrun_bottle <- subset(df, df$runset %in% notrun_bottle_names)

##test no 'failed' runs represented in main df
table(unique(bottle_data$run_set) %in% notrun_bottle_names)
################ all false ################################

## Analyse failed bottle runs

bottle_nr_sum <- notrun_bottle %>% 
  group_by(runset, model, loci, pop_init, ss, set) %>% 
  summarise(count = n())

bottle_nr_sum$model_code <- paste0(bottle_nr_sum$model, bottle_nr_sum$pop_init)

ggplot(data = bottle_nr_sum, aes(x = as.factor(loci), fill = as.factor(ss))) +
  geom_bar(stat = "count", position = position_dodge()) +
  labs(title = "Bottleneck Failed Runs",
       y = "Count",
       x = "Loci (10^3)") +
  theme_minimal() +
  facet_wrap(~model_code) +
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Sample Size", end = 0.8)

###create  dataframe for thesis table
bottle_nr_df <- notrun_bottle %>% 
  group_by(model, loci, pop_init, ss, set) %>% 
  summarise(count = n())

###save as csv
write_csv(bottle_nr_df, file = "~/R/bottle_notrun_table_data.csv")


bottle_nr_df
