##############################################################################
#        Generate decline, expansion, stable model simulations
#############################################################################


#load packages
library(dartR)
library(devtools)
library(slimr)
library(future)
library(geohippos)
library(parallel)
library(plyr)
library(doParallel)
library(tidyr)
library(dplyr)
library(readr)

##If slim_run doesn't work
Sys.setenv(SLIM_HOME='/home/isobel/slim/bin/slim/bin')
#Sys.setenv(SLIM_HOME = "/data/scratch/isobel/win-library/4.1/slimr")
slim_setup()

#==============================================================================#
#                              Simulation script
#==============================================================================#

slim_script(
  slim_block(initialize(),
             {
               
               defineConstant("alpha2", 1e-8);
               defineConstant("L", 5e8);
               defineConstant("model", slimr_template("model"));
               defineConstant("pop_init", slimr_template("pop_init"));
               defineConstant("ts", slimr_template("ts"));
               defineConstant("tl", slimr_template("tl"));
               defineConstant("cp", slimr_template("crash_prop"));
               defineConstant("ss", slimr_template("ss"));
               defineConstant("filename", slimr_template("filename"));
            
               initializeSLiMOptions(nucleotideBased=T);
               initializeAncestralNucleotides(randomNucleotides(L));
               initializeMutationTypeNuc("m1", 0.5, "f",0.0 );
               m1.mutationStackPolicy = "l";
               m1.convertToSubstitution=T;
               defineConstant("mm2", matrix(c(0.0,alpha2,0.0,0.0, alpha2,0.0,0.0,0.0, 0.0,0.0,0.0,alpha2, 0.0,0.0,alpha2,0.0 ),nrow=4, ncol=4));
               initializeGenomicElementType("g1", m1, 1.0, mm2);  
               initializeGenomicElement(g1, 0, L-1);
               initializeRecombinationRate(rates = 1e-8);
               initializeRecombinationRate(rates = c(1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8), ends=c(1e8-1,1e8,2e8-1,2e8,3e8-1, 3e8,4e8-1, 4e8,5e8-1));        
             }),
  
  ### Initialize populations
  slim_block(1, early(),
             {
               ### Define initial population size
               if (model != 2) {
                 sim.addSubpop("p1", pop_init);
               }
               if (model == 2) {
                 exp_pop_init = cp * pop_init;
                 sim.addSubpop("p1", asInteger(round(exp_pop_init)));
               }
               
               ### Log metadata incrementally throughout simulation
               logfilename = paste0("/data/scratch/isobel/log/sim_ind_" + "_rn" + slimr_template("runnumb") + "m" + model + "_p" + pop_init + "_cp" + cp + "_ts" + ts + "_tl" + tl + "_ss" + ss + "f_set3.txt");
               log = community.createLogFile(logfilename, logInterval=10);
               log.addCycle();
               log.addCustomColumn("nind", "p1.individualCount;");
               log.addCustomColumn("model", "model;");
               log.addCustomColumn("pop_init", "pop_init;");
               log.addCustomColumn("crash_prop", "cp;");
               log.addCustomColumn("ts", "ts;");
               log.addCustomColumn("tl", "tl;");
               log.addCustomColumn("ss", "ss;");

             }),
  
  ### Population size changes
  slim_block(2000,2201, late(),
             
             {
               startyr = 2200 - ts;
               endyr = 2200 - ts + tl;
               ###Abort simulation if trajectory length is not compatible with trajectory start time
               if(tl > ts) {
                 sim.simulationFinished();
                 ###create error message###
               }
               ###separate by trajectory start - delay until start time
               if ((sim.cycle + 1) > startyr) {
                 
                 ###Decline model###
                 if (model == 1) {
                   
                   ###run until end of trajectory length, relative to start time
                   x1 = pop_init;
                   xn = cp * pop_init;
                   yr = sim.cycle - startyr;
                   r = log(x1 - xn)/tl;
                   p1.setSubpopulationSize(asInteger(round(
                     (x1 - xn)*exp((-r) * yr) + xn)
                   )); 
                   
                 }
                 
                 ###Expansion model###
                 if (model == 2) {
                   
                   ###run until end of trajectory length, relative to start time
                   if (sim.cycle < endyr) {
                     x1 = cp * pop_init;
                     xn =  pop_init;
                     yr = sim.cycle - startyr;
                     r = log(xn/x1)/tl;
                     
                     p1.setSubpopulationSize(asInteger(round(
                       x1*exp(r*yr))
                     ));
                   }
                   if (sim.cycle > (endyr - 1)) {
                     p1.setSubpopulationSize(asInteger(pop_init));
                   }
                 }
               
             }
               }),
  
  slim_block(2201,late(),
             {
               nn = filename;
               p1.outputVCFSample(sampleSize=ss, replace=F,  outputMultiallelics=F,filePath=nn,  simplifyNucleotides=T);
               sim.simulationFinished();   
             })
) -> script_1

###select folder for vcf files to save to
outdir <- "/data/scratch/isobel/vcf_std_set3/"

###create dataframe of all parameter combinations
df <- expand.grid(rep =1:10,
                  pop_init = c(1000, 500, 100, 50),
                  crash_prop = c(0.5, 0.1, 0.05),
                  tl = c(200, 100, 50, 30), #trajectory length (from ts)
                  ts = c(200, 100, 50, 30), #trajectory start (ybp)
                  ss = c(200, 100, 50, 20),
                  model = c("decline", "expansion", "stable"))

###allocate run numbers
df$runnumb <- paste0("Run_",sprintf("%05d",1:nrow(df)))




#==============================================================================#
# Remove combinations of incompatible parameters and extra stable replicates  
#==============================================================================#

### Remove extra stable replicates
stable <- df %>% filter(model == "stable")
dup_nums <- stable[duplicated(stable[c(1,2,6)]), ]$runnumb
df <- df %>% filter(!(runnumb %in% dup_nums))
rm(stable, dup_nums)

### Filter out dataframes where ss > final pop
invalid_ss <- df %>% filter(as.numeric(ss) > as.numeric(pop_init))
if (nrow(invalid_ss) > 0) {
  df <- df %>% filter(!(runnumb %in% invalid_ss$runnumb))
}

invalid_ss_dec <- df %>% filter(model == "decline") %>% filter(as.numeric(ss) > round((pop_init * crash_prop)))
if (nrow(invalid_ss_dec) > 0) {
  df <- df %>% filter(!(runnumb %in% invalid_ss_dec$runnumb))
}

rm(invalid_ss, invalid_ss_dec)

### Filter out rows where trajectory length is longer than trajectory start
invalid_tstl <- df %>% filter(tl > ts)
if (nrow(invalid_tstl) > 0) {
  df <- df %>% filter(!(runnumb %in% invalid_tstl$runnumb))
}

rm(invalid_tstl)

### Generate final run numbers and filenames
df$filename <- paste0(outdir, "_", df$runnumb, "_slim_rep", df$rep, "_pop", df$pop_init, "_tl", df$tl, "_crash", df$crash_prop, "_ts-ybp", df$ts, "_model", df$model, "_ss", df$ss,  "_final_std_set3.vcf")

#Write excel sheet
write_csv(df, file = "/data/scratch/isobel/combination_set.csv")


#==============================================================================#
#                             Run simulations in SLiM
#==============================================================================#

###run multiple scripts
library(slimr)
library(future)
Sys.setenv(SLIM_HOME='/home/isobel/slim/bin/slim/bin')
plan(multisession, workers = 35)

### Render script
testscript <- slim_script_render(script_1, template = df, parallel = 30)

### Run simulations
sr <- slim_run(testscript, parallel = T)



