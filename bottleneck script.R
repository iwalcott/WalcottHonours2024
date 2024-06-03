################################################################################
#                  Generate bottleneck model simulations
################################################################################

### load packages
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

###If slim_run doesn't work
Sys.setenv(SLIM_HOME = "C:/Users/Isobel/Documents/R/win-library/4.1/slimr")
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
               defineConstant("pop_crash", slimr_template("pop_crash"));
               defineConstant("pop_recov", slimr_template("pop_recov"));
               defineConstant("crash_st", slimr_template("crash_st"));
               defineConstant("recov_st", slimr_template("recov_st"));
               defineConstant("crash_len", slimr_template("crash_len"));
               defineConstant("recov_len", slimr_template("recov_len"));
               defineConstant("k_decay", slimr_template("k_decay"));
               defineConstant("k_growth", slimr_template("k_growth"));
               defineConstant("ss", slimr_template("ss"));
               defineConstant("runnum", slimr_template("runnumb"));
               defineConstant("filename", slimr_template("filename"));
               
               
               #defineConstant("simlen", 5000);
               #how many outputs (simlen=only once at the end)
               initializeSLiMOptions(nucleotideBased=T);
               initializeAncestralNucleotides(randomNucleotides(L));
               initializeMutationTypeNuc("m1", 0.5, "f",0.0 );
               m1.mutationStackPolicy = "l";
               m1.convertToSubstitution=T;
               defineConstant("mm2", matrix(c(0.0,alpha2,0.0,0.0, alpha2,0.0,0.0,0.0, 0.0,0.0,0.0,alpha2, 0.0,0.0,alpha2,0.0 ),nrow=4, ncol=4));
               #defineConstant("mm2", mmJukesCantor((1e-8)/3));
               initializeGenomicElementType("g1", m1, 1.0, mm2);  
               initializeGenomicElement(g1, 0, L-1);
               #initializeRecombinationRate(rates = 1e-8);
               initializeRecombinationRate(rates = c(1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8,0.5,1e-8), ends=c(1e8-1,1e8,2e8-1,2e8,3e8-1, 3e8,4e8-1, 4e8,5e8-1));        
             }),
  
  ### Initialize populations
  slim_block(1, early(),
             {
               ##Define initial population size
               sim.addSubpop("p1", pop_init);
               
               ### Log metadata incrementally throughout simulation
               logfilename = paste0("/data/scratch/isobel/bottle_log/bottle_sim" + runnum + "mod_" + model + pop_init +".txt");
               log = community.createLogFile(logfilename, logInterval=5);
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
  slim_block(2000,2300, late(), {
    
    
    ###Delay until start time (set up finalise each run at the same time)
    if (sim.cycle > (crash_st - 1)) {
      
      ###Start population decay
      if (sim.cycle < (crash_st + crash_len + 1)) {
        k = k_decay;
        n = sim.cycle - crash_st;
        newsize = asInteger(round(pop_init*exp(-k*n)));
        if (newsize < pop_crash) {
          newsize = pop_crash
        };
        p1.setSubpopulationSize(newsize);
      }
      
      ###Start population recovery
      if (sim.cycle > (recov_st - 1)) {
        k = k_growth;
        n = sim.cycle - recov_st;
        newsize = asInteger(round(pop_crash*exp(k*n)));
        if (newsize > pop_recov) {
          newsize = pop_recov;
        };
        p1.setSubpopulationSize(newsize);
      }
    }
               
}),
  
  slim_block(2301,late(),
             {
               nn = filename;
               p1.outputVCFSample(sampleSize=ss, replace=F,  outputMultiallelics=F,filePath=nn,  simplifyNucleotides=T);
               sim.simulationFinished();   
             })
) -> script_bottle

### select folder for vcf files to go
outdir <- "/data/scratch/isobel/bottle_vcf_newset/"

### create dataframe of all test combinations
bottledf <- data.frame(model = c("croc", "croc", "whale", "whale", "seal", "seal", "frog", "frog"), 
                       pop_init = c(7500, 15000, 6000, 8000, 2100, 4200, 50, 100), 
                       pop_crash = c(909, 909, 83, 83, 872, 872, 10, 10),
                       pop_recov = c(7500, 7500, 4000, 4000, 2100, 2100, 50, 50),
                       crash_st = c(2222, 2222, 2067, 2067, 2027, 2027, 2257, 2257),
                       recov_st = c(2248, 2248, 2249, 2249, 2263, 2263, 2278, 2278),
                       crash_len = c(10, 10, 130, 130, 173, 173, 15, 15),
                       recov_len = c(27, 27, 46, 46, 16, 16, 7, 7),
                       k_decay = c(0.211, 0.280, 0.033, 0.035, 0.005, 0.009, 0.107, 0.154),
                       k_growth = c(0.0782, 0.0782, 0.0842, 0.0842, 0.0549, 0.0549, 0.2299, 0.2299))

### create combinations for all replicates and sample sizes
df <- crossing(bottledf, rep = 1:10, ss = c(200, 100, 50, 20))

### allocate run number
df$runnumb <- paste0("Run_",sprintf("%05d",1:nrow(df)))


#==============================================================================#
# Remove combinations of incompatible parameters and extra stable replicates  
#==============================================================================#

### Filter out combinations where ss > final pop
invalid_ss <- df %>% filter(as.numeric(ss) > as.numeric(pop_recov)) %>% select(runnumb)
if (nrow(invalid_ss) > 0) {
  df <- df %>% filter(!(runnumb %in% invalid_ss$runnumb))
}

### reallocate run numbers
df$runnumb <- paste0("Run_",sprintf("%05d",1:nrow(df)))

### generate filenames
df$filename <- paste0(outdir, "bottle_", df$runnumb, "_model", df$model, "_rep", df$rep, "_pop", df$pop_init,  "_ss", df$ss,  "_final.vcf")

#==============================================================================#
#                             Run simulations in SLiM
#==============================================================================#

### run multiple scripts
library(slimr)
library(future)
Sys.setenv(SLIM_HOME='/home/isobel/slim/bin/slim/bin')
plan(multisession, workers = 10)

### Render scripts
testscript <- slim_script_render(script_bottle, template = df, parallel = 10)

### Run simulations
sr <- slim_run(testscript, parallel = T)

