library(data.table)

dat<-fread("t_rw_comb_outg_perform_reg_depth.txt",header=FALSE)
dm<-as.matrix(dat[,-c(1:2)])
mn<-apply(dm,1,mean)
prop<-apply(dm > 0,1,mean)

keep<-mn > 2 & prop > 0.8

pos<-as.numeric(unlist(dat[keep,2]))

save(list=ls(),file="invar_positions.rdat")

## read in file with potentially variable sites, pre-filtering
## grep ^Sc ~/../gompert-group3/data/timema_clines_rw_SV/variants_rw_plus/t_rw_comb_outg_perform.vcf | cut -f 2 > var_pos.txt 

pos_invar<-pos[-which(pos %in% vp)]
length(pos_invar)
## [1] 60612
## left 60,612 confident, invariant sites

write.table(pos_invar,file="invar_positions_perform_og.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
