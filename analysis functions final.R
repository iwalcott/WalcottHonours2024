##############################################################################
#                   Analysis helper functions
##############################################################################


#==============================================================================#
#              Trim dataset to analysis range - single core
#==============================================================================#

trim_single <- function(data_list) {
  trimmed <- list()
  for (i in 1:length(data_list)) {
    print(i)
    curr_run <- data_list[i][[1]]
    ### test more efficient method
    new_run <- curr_run
    
    ### trim dataframes
    new_run[[1]] <- new_run[[1]][new_run[[1]]$X.Time <= 250 #& new_run[[1]]$X.Time >= 1
                                 , ]
    new_run[[2]] <- new_run[[2]][new_run[[2]]$year <= 250 #& new_run[[2]]$year >= 1
                                 , ]
    new_run[[3]] <- new_run[[3]][new_run[[3]]$Generation <= 250  #& new_run[[3]]$Generation >= 1
                                 , ]
    new_run[[4]] <- new_run[[4]][new_run[[4]]$year <= (250 + 1)  #& new_run[[4]]$year >= 1
                                 , ]
    
    ### check dataframes for large gaps - if present, add next available data point

    if (nrow(new_run[[1]]) > 0 & (tail(new_run[[1]]$X.Time, n=1) < 225) & (nrow(new_run[[1]]) < nrow(curr_run[[1]]))) {
      new_run[[1]] <- rbind(new_run[[1]], curr_run[[1]][nrow(new_run[[1]]) + 1, ])
    }
    if (nrow(new_run[[2]]) > 0 & (tail(new_run[[2]]$year, n=1) < 225) & (nrow(new_run[[2]]) < nrow(curr_run[[2]]))) {
      new_run[[2]] <- rbind(new_run[[2]], curr_run[[2]][nrow(new_run[[2]]) + 1, ])
    }
    if (nrow(new_run[[3]]) > 0 & (tail(new_run[[3]]$Generation, n=1) < 225) & (nrow(new_run[[3]]) < nrow(curr_run[[3]]))) {
      new_run[[3]] <- rbind(new_run[[3]], curr_run[[3]][nrow(new_run[[3]]) + 1, ])
    }
    
    ### transfer into output with run name
    trimmed[i] <- data_list[i]
    names(trimmed[i]) <- names(data_list[i])
    trimmed[i][[1]] <- new_run
  }
  names(trimmed) <- names(data_list)
  return(trimmed)
}



#==============================================================================#
#              Trim dataset to analysis range - parallel
#==============================================================================#

### trim_pt = final year to include in range

trim <- function(data_list, ncores = 10, trim_pt = 250) {
  trimmed <- mclapply(1:length(data_list), function(i) {
    print(i)
    curr_run <- data_list[i][[1]]
    new_run <- curr_run
    
    ### subset dataframes to < trim_pt
    new_run[[1]] <- new_run[[1]][new_run[[1]]$X.Time <= trim_pt #& new_run[[1]]$X.Time >= 1
                                 , ]
    new_run[[2]] <- new_run[[2]][new_run[[2]]$year <= trim_pt #& new_run[[2]]$year >= 1
                                 , ]
    new_run[[3]] <- new_run[[3]][new_run[[3]]$Generation <= trim_pt  #& new_run[[3]]$Generation >= 1
                                 , ]
    new_run[[4]] <- new_run[[4]][new_run[[4]]$year <= (trim_pt + 1)  #& new_run[[4]]$year >= 1
                                 , ]
    
    ### check dataframes for large gaps - if present, add next available data point
    if (nrow(new_run[[1]]) > 0 & (tail(new_run[[1]]$X.Time, n=1) < (trim_pt - 25)) & (nrow(new_run[[1]]) < nrow(curr_run[[1]]))) {
      new_run[[1]] <- rbind(new_run[[1]], curr_run[[1]][nrow(new_run[[1]]) + 1, ])
    }
    if (nrow(new_run[[1]][new_run[[1]]$X.Time > trim_pt, ])) {
      new_run[[1]][new_run[[1]]$X.Time > trim_pt, ]$X.Time <- trim_pt
    }
    new_run[[1]][new_run[[1]]$X.Time > trim_pt, ]
    if (nrow(new_run[[2]]) > 0 & (tail(new_run[[2]]$year, n=1) < (trim_pt - 25)) & (nrow(new_run[[2]]) < nrow(curr_run[[2]]))) {
      new_run[[2]] <- rbind(new_run[[2]], curr_run[[2]][nrow(new_run[[2]]) + 1, ])
    }
    if (nrow(new_run[[2]][new_run[[2]]$year > trim_pt, ])) {
      new_run[[2]][new_run[[2]]$year > trim_pt, ]$year <- trim_pt
    }
    if (nrow(new_run[[3]]) > 0 & (tail(new_run[[3]]$Generation, n=1) < (trim_pt - 25)) & (nrow(new_run[[3]]) < nrow(curr_run[[3]]))) {
      new_run[[3]] <- rbind(new_run[[3]], curr_run[[3]][nrow(new_run[[3]]) + 1, ])
    }
    if (nrow(new_run[[3]][new_run[[3]]$Generation > trim_pt, ])) {
      new_run[[3]][new_run[[3]]$Generation > trim_pt, ]$Generation <- trim_pt
    }
    
    return(new_run)
  }, mc.cores = ncores)
  names(trimmed) <- names(data_list)
  return(trimmed)
}



#==============================================================================#
#                              Metadata function  
#==============================================================================#
                               
### add simulation run data to dataframe
metadata <- function(data_list) {
  meta <- data.frame(run = names(data_list));
  meta <- meta %>% mutate(loci = as.numeric(substr(run, 11, 13)));
  for (i in 1:nrow(meta)) {
    model <- sim_files[meta$run[i]][[1]]$model[3]
    pop_init <- sim_files[meta$run[i]][[1]]$pop_init[3]
    ts <- sim_files[meta$run[i]][[1]]$ts[3]
    tl <- sim_files[meta$run[i]][[1]]$tl[3]
    ss <- sim_files[meta$run[i]][[1]]$ss[3]
    cp <- sim_files[meta$run[i]][[1]]$crash_prop[3]
    meta$model[i] <- model
    meta$pop_init[i] <- pop_init
    meta$ts[i] <- ts
    meta$tl[i] <- tl
    meta$ss[i] <- ss
    meta$cp[i] <- cp
  }
  return(meta)
}



#==============================================================================#
#                   Convert to long dataframe (single core)
#==============================================================================#

### convert data to one long dataframe
list_to_long <- function(data_list) {
  dat_long <- data.frame()
  for (i in 1:length(data_list)) {
    epos <- data.frame(method = "epos", year = data_list[[i]][[1]]$X.Time, Ne = data_list[[i]][[1]]$Median)
    stair <- data.frame(method = "stairway", year = data_list[[i]][[2]]$year, Ne = data_list[[i]][[2]]$Ne_median)
    gone <- data.frame(method = "gone", year = data_list[[i]][[3]]$Generation, Ne = data_list[[i]][[3]]$Geometric_mean)
    sim <- data.frame(method = "sim", year = data_list[[i]][[4]]$year, Ne = data_list[[i]][[4]]$nind)

    newrun <- rbind(epos, stair, gone, sim)
    newrun$model <- data_list[[i]][[4]]$model[1]
    newrun$pop_init <- data_list[[i]][[4]]$pop_init[1]
    newrun$tl <- data_list[[i]][[4]]$tl[1]
    newrun$ts <- data_list[[i]][[4]]$ts[1]
    newrun$ss <- data_list[[i]][[4]]$ss[1]
    newrun$cp <- data_list[[i]][[4]]$crash_prop[1]
    newrun$runnum <- names(data_list[i])
    newrun$loci <- as.numeric(substr(names(data_list[i]), 12, 14))
    
    dat_long <- rbind(dat_long, newrun)
  }
  return(dat_long)
}


#==============================================================================#
#           Bottleneck data: Convert to long dataframe (single core)
#==============================================================================#

### convert data to one long dataframe
bottle_long <- function(data_list) {
  dat_long <- data.frame()
  for (i in 1:length(data_list)) {
    epos <- data.frame(method = "epos", year = data_list[[i]][[1]]$X.Time, Ne = data_list[[i]][[1]]$Median)
    stair <- data.frame(method = "stairway", year = data_list[[i]][[2]]$year, Ne = data_list[[i]][[2]]$Ne_median)
    gone <- data.frame(method = "gone", year = data_list[[i]][[3]]$Generation, Ne = data_list[[i]][[3]]$Geometric_mean)
    sim <- data.frame(method = "sim", year = data_list[[i]][[4]]$year, Ne = data_list[[i]][[4]]$nind)
    
    newrun <- rbind(epos, stair, gone, sim)
    newrun$model <- data_list[[i]][[4]]$model[1]
    newrun$pop_init <- data_list[[i]][[4]]$pop_init[1]

    newrun$runnum <- names(data_list[i])
    newrun$loci <- as.numeric(substr(names(data_list[i]), 11, 13))
    
    dat_long <- rbind(dat_long, newrun)
  }
  return(dat_long)
}


#==============================================================================#
#                   Convert to long dataframe - Parallel
#==============================================================================#

### convert data in list form to long dataframe in parallel
par_data_long <- function(data_list, ncores = 10) {
  
  long_list <- mclapply(1:length(data_list), function(i) {
    
    epos <- data.frame(method = "epos", year = data_list[[i]][[1]]$X.Time, Ne = data_list[[i]][[1]]$Median)
    stair <- data.frame(method = "stairway", year = data_list[[i]][[2]]$year, Ne = data_list[[i]][[2]]$Ne_median)
    gone <- data.frame(method = "gone", year = data_list[[i]][[3]]$Generation, Ne = data_list[[i]][[3]]$Geometric_mean)
    sim <- data.frame(method = "sim", year = data_list[[i]][[4]]$year, Ne = data_list[[i]][[4]]$nind)
    
    newrun <- rbind(epos, stair, gone, sim)
    newrun$model <- data_list[[i]][[4]]$model[1]
    newrun$pop_init <- data_list[[i]][[4]]$pop_init[1]
    newrun$tl <- data_list[[i]][[4]]$tl[1]
    newrun$ts <- data_list[[i]][[4]]$ts[1]
    newrun$ss <- data_list[[i]][[4]]$ss[1]
    newrun$cp <- data_list[[i]][[4]]$crash_prop[1]
    newrun$runnum <- names(data_list[i])
    newrun$loci <- as.numeric(substr(names(data_list[i]), 12, 14))
    return(newrun)
  }, mc.cores = ncores)
  names(long_list) <- names(data_list)
  out <- rbindlist(long_list)
  return(data.frame(out))
}


#==============================================================================#
#                   Bottle - Convert to long dataframe - Parallel
#==============================================================================#

par_bottle_long <- function(data_list, ncores = 10) {
  
  long_list <- mclapply(1:length(data_list), function(i) {
    
    epos <- data.frame(method = "epos", year = data_list[[i]][[1]]$X.Time, Ne = data_list[[i]][[1]]$Median)
    stair <- data.frame(method = "stairway", year = data_list[[i]][[2]]$year, Ne = data_list[[i]][[2]]$Ne_median)
    gone <- data.frame(method = "gone", year = data_list[[i]][[3]]$Generation, Ne = data_list[[i]][[3]]$Geometric_mean)
    sim <- data.frame(method = "sim", year = data_list[[i]][[4]]$year, Ne = data_list[[i]][[4]]$nind)
    
    newrun <- rbind(epos, stair, gone, sim)
    newrun$model <- data_list[[i]][[4]]$model[1]
    newrun$pop_init <- data_list[[i]][[4]]$pop_init[1]
    newrun$runnum <- names(data_list[i])
    newrun$loci <- as.numeric(substr(names(data_list[i]), 11, 13))
    return(newrun)
  }, mc.cores = ncores)
  names(long_list) <- names(data_list)
  out <- rbindlist(long_list)
  return(data.frame(out))
}


#==============================================================================#
#   Extract full loci data from vcf files - decline, expansion, stable models
#==============================================================================#

library(geohippos)
library(dartR)
library(parallel)

### get list of simulation vcf output files
fn_loc <- list.files("/data/scratch/isobel/std_vcf_set/", pattern = "*.vcf", full.names = TRUE)

### function to load simulation vcfs as gls files and return total loci counts
glsnloc <- function(i)
{
  gls <- geohippos::gl.read.vcf(fn_loc[i], verbose = 0)
  return(nLoc(gls))
}

### run function over list in parallel
ss <- mclapply(1:length(fn_loc), function(x) glsnloc(x), mc.cores  =14)

### run function over list on single core
ss_loc <- lapply(1:length(fn_loc), function(x) glsnloc(x))

### add run numbers as names
names(ss_loc) <- substr(fn_loc, 28, 36)

data.frame(unlist(ss_loc))


#==============================================================================#
#   Extract full loci data from vcf files - bottleneck models
#==============================================================================#

library(geohippos)
library(dartR)
library(parallel)

### get list of simulation vcf output files
fn_loc <- list.files("/data/scratch/isobel/bottle_vcf_newset/", pattern = "*.vcf", full.names = TRUE)

### function to load simulation vcfs as gls files and return total loci counts
glsnloc <- function(i)
{
  gls <- geohippos::gl.read.vcf(fn_loc[i], verbose = 0)
  return(nLoc(gls))
}

### run function over list in parallel
ss <- mclapply(1:length(fn_loc), function(x) glsnloc(x), mc.cores  =14)

### extract species name
sp_name <- lapply(strsplit(substr(fn_loc, 63, 68), "*_"), function(x) x[[1]])

### extract initial population size number  
p0_name <-lapply(strsplit(substr(fn_loc, 73, 82), "pop"), function(x) strsplit(x[[2]], "_")[[1]][[1]])

### allocate species name and initial pop size for name
names(ss) <- lapply(1:length(sp_name), function(x) {
  sp_name[[x]][[1]]
  p0_name[[x]][[1]]
  return(paste0(sp_name[[x]][[1]], p0_name[[x]][[1]]))
})

df_loci <- data.frame(unlist(ss, use.names = T))

df_loci$loci <- unlist(ss)
df_loci$run_code <- names(ss)


           