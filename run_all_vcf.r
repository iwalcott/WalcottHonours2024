##############################################################################
#               Run demographic inference methods on Gadi supercomputer
##############################################################################


library(geohippos)

#empty all objects
res_epos <- res_stair <- res_gone <- NULL

#epos

tt <- Sys.time()
res_epos <- gl.epos(dummy, epos.path = paste0("~/binaries/epos/linux"), l = L, u=mu, boot=40, verbose=0, minbinsize = 1, other.options = " -E 2", cleanup = T)
res_epos$runtime <- round(Sys.time()-tt)

#stairway2

tt <- Sys.time()
res_stair <- gl.stairway2(dummy, verbose = T,stairway.path="~/binaries/stairways/", mu = mu, gentime = 1, run=TRUE, nreps = 40, parallel=1, L=L,  minbinsize =1, cleanup = T, nrand=5)
res_stair$runtime <- round(Sys.time()-tt)

#filter gone per hand to minbinsize=1 (gone uses MAF) so MAF in gone to zero and prefilter here

dummy <- gl.filter.maf(dummy, threshold = 1/(2*nInd(dummy)))
dummy@chromosome <- factor(rep("1",nLoc(dummy)))

#gone

tt <- Sys.time()
res_gone <- gl.gone(dummy,gone.path = paste0("~/binaries/gone/linux"), cleanup = T)
res_gone$runtime <- round(Sys.time()-tt)

#save all three in one output

out <- list(epos=res_epos, stair=res_stair, gone=res_gone)
saveRDS(out, file =  file.path("~/res2",paste0(df$id[x],".rds")))



#delete file if it ran successfully

if (is.data.frame(res_epos) & is.data.frame(res_stair) & is.data.frame(res_gone)) ret <- TRUE else ret<-FALSE







