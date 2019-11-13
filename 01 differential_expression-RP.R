library(Biobase)
library(GEOquery)
library(RankProd)

load("00 meta-analysis/meta_gse.RData")
load("00 meta-analysis/meta_sample_info.RData")
load("01 differential expression-limma/gpl.RData")

####################################### GSE14763 ########################################

# log2 transform
quantile(exprs(gse14763), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)
exprs_gse14763 <- na.omit(log2(exprs(gse14763)+1))

fl <- as.factor(sample_info_gse14763$Stage)

# acute
exprs_gse14763_acute <- exprs_gse14763[,c(which(fl==c("acute")), which(fl==c("control")))]
cl <- rep(c(0,1),c(10,3))

RP_gse14763_acute <- RankProducts(exprs_gse14763_acute, cl, logged=TRUE, 
             na.rm=FALSE,plot=FALSE, rand=123)

plotRP(RP_gse14763_acute, cutoff = 0.05)


RP_gse14763_acute_df <- RPlist2df(RP_gse14763_acute)
RP_gse14763_acute_df <- merge(RP_gse14763_acute_df, gpl2896, by = "ID")

# acute DEGs
RP.DEG_acute_gse14763 <- RP_gse14763_acute_df %>% filter(pfp_up<0.05 | pfp_down<0.05)

# latent
exprs_gse14763_latent <- exprs_gse14763[,c(which(fl==c("latent")), which(fl==c("control")))]
cl <- rep(c(0,1),c(5,3))

RP_gse14763_latent <- RankProducts(exprs_gse14763_latent, cl, logged=TRUE, 
                                  na.rm=FALSE,plot=FALSE, rand=123)
plotRP(RP_gse14763_latent, cutoff = 0.05)

RP_gse14763_latent_df <- RPlist2df(RP_gse14763_latent)
RP_gse14763_latent_df <- merge(RP_gse14763_latent_df, gpl2896, by = "ID")

# latent DEGs
RP.DEG_latent_gse14763 <- RP_gse14763_latent_df %>% filter(pfp_up<0.05 | pfp_down<0.05)


# match to probe annotation
RP.DEG_acute_gse14763 <- merge(RP.DEG_acute_gse14763, gpl2896, by = "ID")
RP.DEG_latent_gse14763 <- merge(RP.DEG_latent_gse14763, gpl2896, by = "ID")


####################################### GSE27166 ########################################

# log2 transform
quantile(exprs(gse27166), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)
exprs_gse27166 <- na.omit(log2(exprs(gse27166)+1))

fl <- as.factor(sample_info_gse27166_right$Stage)
exprs_gse27166_chronic <- exprs_gse27166[,gse27166$`side:ch1`=="right"]
# remove outlier 
fl <- fl[-c(7,12)]
exprs_gse27166_chronic <- exprs_gse27166_chronic[, -c(7,12)]

cl <- rep(c(0,1,0,1),c(3,3,2,2))

RP_gse27166_chronic <- RankProducts(exprs_gse27166_chronic, cl, logged=TRUE, 
                                  na.rm=FALSE,plot=FALSE, rand=123)

plotRP(RP_gse27166_chronic, cutoff = 0.05)

RP_gse27166_chronic_df <- RPlist2df(RP_gse27166_chronic)
RP_gse27166_chronic_df <- merge(RP_gse27166_chronic_df, gpl2896, by = "ID")

# chronic DEGs
RP.DEG_chronic_gse27166 <- RP_gse27166_chronic_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
RP.DEG_chronic_gse27166 <- merge(RP.DEG_chronic_gse27166, gpl2896, by = "ID")

####################################### GSE49849 #######################################

# log2 transform
quantile(exprs(gse49849), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)
exprs_gse49849 <- na.omit(exprs(gse49849))

fl <- as.factor(sample_info_gse49849$Stage)

# acute
exprs_gse49849_acute <- exprs_gse49849[,c(which(fl==c("acute")), which(fl==c("control")))]
cl <- rep(c(0,1),c(5,10))

RP_gse49849_acute <- RankProducts(exprs_gse49849_acute, cl, logged=TRUE, 
                                  na.rm=FALSE,plot=FALSE, rand=123)

plotRP(RP_gse49849_acute, cutoff = 0.05)

RP_gse49849_acute_df <- RPlist2df(RP_gse49849_acute)
RP_gse49849_acute_df <- merge(RP_gse49849_acute_df, gpl6247, by = "ID")

# acute DEGs
RP.DEG_acute_gse49849 <- RP_gse49849_acute_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
RP.DEG_acute_gse49849 <- merge(RP.DEG_acute_gse49849, gpl6247, by = "ID")


# latent
exprs_gse49849_latent <- exprs_gse49849[,c(which(fl==c("latent")), which(fl==c("control")))]
cl <- rep(c(0,1),c(5,10))

RP_gse49849_latent <- RankProducts(exprs_gse49849_latent, cl, logged=TRUE, 
                                   na.rm=FALSE,plot=FALSE, rand=123)

plotRP(RP_gse49849_latent, cutoff = 0.05)

RP_gse49849_latent_df <- RPlist2df(RP_gse49849_latent)
RP_gse49849_latent_df <- merge(RP_gse49849_latent_df, gpl6247, by = "ID")

# latent DEGs
RP.DEG_latent_gse49849 <- RP_gse49849_latent_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
RP.DEG_latent_gse49849 <- merge(RP.DEG_latent_gse49849, gpl6247, by = "ID")

####################################### GSE73878 ########################################

# log2 transform
quantile(exprs(gse73878), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)
exprs_gse73878 <- na.omit(exprs(gse73878))

fl <- as.factor(sample_info_gse73878_IL$Stage)

exprs_gse73878 <- exprs_gse73878[ ,gse73878$`hemisphere:ch1` == "IL"]
exprs_gse73878 <- exprs_gse73878[,-c(16,11,29,30,31,32)]

# acute
exprs_gse73878_acute <- exprs_gse73878[,c(which(fl==c("acute")), which(fl==c("control")))]
cl <- rep(c(0,1),c(7,11))

RP_gse73878_acute <- RankProducts(exprs_gse73878_acute, cl, logged=TRUE, 
                                  na.rm=FALSE,plot=FALSE, rand=123)
plotRP(RP_gse73878_acute, cutoff = 0.05)

RP_gse73878_acute_df <- RPlist2df(RP_gse73878_acute)
RP_gse73878_acute_df <- merge(RP_gse73878_acute_df, gpl6885, by = "ID")

# acute DEGs
RP.DEG_acute_gse73878 <- RP_gse73878_acute_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
RP.DEG_acute_gse73878 <- merge(RP.DEG_acute_gse73878, gpl6885, by = "ID")

# latent
exprs_gse73878_latent <- exprs_gse73878[,c(which(fl==c("latent")), which(fl==c("control")))]
cl <- rep(c(0,1),c(8,11))

RP_gse73878_latent <- RankProducts(exprs_gse73878_latent, cl, logged=TRUE, 
                                   na.rm=FALSE,plot=FALSE, rand=123)
plotRP(RP_gse73878_latent, cutoff = 0.05)

RP_gse73878_latent_df <- RPlist2df(RP_gse73878_latent)
RP_gse73878_latent_df <- merge(RP_gse73878_latent_df, gpl6885, by = "ID")

# latent DEGs
RP.DEG_latent_gse73878 <- RP_gse73878_latent_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
RP.DEG_latent_gse73878 <- merge(RP.DEG_latent_gse73878, gpl6885, by = "ID")

# chronic
exprs_gse73878_chronic <- exprs_gse73878[,c(which(fl==c("chronic")), which(fl==c("control")))]
cl <- rep(c(0,1),c(8,11))

RP_gse73878_chronic <- RankProducts(exprs_gse73878_chronic, cl, logged=TRUE, 
                                   na.rm=FALSE,plot=FALSE, rand=123)
plotRP(RP_gse73878_chronic,cutoff = 0.05)

RP_gse73878_chronic_df <- RPlist2df(RP_gse73878_chronic)
RP_gse73878_chronic_df <- merge(RP_gse73878_chronic_df, gpl6885, by = "ID")

# chronic DEGs
RP.DEG_chronic_gse73878 <- RP_gse73878_chronic_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
RP.DEG_chronic_gse73878 <- merge(RP.DEG_chronic_gse73878, gpl6885, by = "ID")

################################ GSE88992 ########################################

# log2 transform
quantile(exprs(gse88992), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)

exprs_gse88992 <- na.omit(exprs(gse88992))

fl <- as.factor(sample_info_gse88992$Stage)

cl <- rep(c(1,0,1,0,1,0),c(3,3,3,3,3,2))

RP_gse88992_acute <- RankProducts(exprs_gse88992, cl, logged=TRUE, 
                                    na.rm=FALSE,plot=FALSE, rand=123)
plotRP(RP_gse88992_acute, cutoff = 0.05)

RP_gse88992_acute_df <- RPlist2df(RP_gse88992_acute)
RP_gse88992_acute_df <- merge(RP_gse88992_acute_df, gpl1261, by = "ID")

# acute DEGs
RP.DEG_acute_gse88992 <- RP_gse88992_acute_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
RP.DEG_acute_gse88992 = merge(RP.DEG_acute_gse88992, gpl1261, by = "ID")



#########################################################################################
par(mfrow = c(3,3),cex=0.7)
hist(RP.DEG_acute_gse14763$logFC,breaks = 60)
abline(v=c(0.5,-0.5),col="red")
hist(RP.DEG_latent_gse14763$logFC,breaks = 70)
abline(v=c(0.5,-0.5),col="red")

hist(RP.DEG_chronic_gse27166$logFC,breaks = 50)
abline(v=c(0.5,-0.5),col="red")

hist(RP.DEG_acute_gse49849$logFC,breaks = 52)
abline(v=c(0.5,-0.5),col="red")
hist(RP.DEG_latent_gse49849$logFC,breaks = 52)
abline(v=c(0.5,-0.5),col="red")

hist(RP.DEG_acute_gse73878$logFC,breaks = 70)
abline(v=c(0.5,-0.5),col="red")
hist(RP.DEG_latent_gse73878$logFC,breaks = 70)
abline(v=c(0.5,-0.5),col="red")
hist(RP.DEG_chronic_gse73878$logFC,breaks = 50)
abline(v=c(0.5,-0.5),col="red")

hist(RP.DEG_acute_gse88992$logFC,breaks = 52)
abline(v=c(0.5,-0.5),col="red")

# DEG acute
RP.DEG_acute_gse14763_0.5 <- RP.DEG_acute_gse14763 %>% filter(abs(logFC)>0.5)
RP.DEG_acute_gse49849_0.5 <- RP.DEG_acute_gse49849 %>% filter(abs(logFC)>0.5)
RP.DEG_acute_gse73878_0.5 <- RP.DEG_acute_gse73878 %>% filter(abs(logFC)>0.5)
RP.DEG_acute_gse88992_0.5 <- RP.DEG_acute_gse88992 %>% filter(abs(logFC)>0.5)

write.csv(unique(RP.DEG_acute_gse14763_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_acute_gse14763_0.5.csv")
write.csv(unique(RP.DEG_acute_gse49849_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_acute_gse49849_0.5.csv")
write.csv(unique(RP.DEG_acute_gse73878_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_acute_gse73878_0.5.csv")
write.csv(unique(RP.DEG_acute_gse88992_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_acute_gse88992_0.5.csv")

# DEG latent
RP.DEG_latent_gse14763_0.5 <- RP.DEG_latent_gse14763 %>% filter(abs(logFC)>0.5)
RP.DEG_latent_gse49849_0.5 <- RP.DEG_latent_gse49849 %>% filter(abs(logFC)>0.5)
RP.DEG_latent_gse73878_0.5 <- RP.DEG_latent_gse73878 %>% filter(abs(logFC)>0.5)

write.csv(unique(RP.DEG_latent_gse14763_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_latent_gse14763_0.5.csv")
write.csv(unique(RP.DEG_latent_gse49849_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_latent_gse49849_0.5.csv")
write.csv(unique(RP.DEG_latent_gse73878_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_latent_gse73878_0.5.csv")

# DEG chronic
RP.DEG_chronic_gse27166_0.5 <- RP.DEG_chronic_gse27166 %>% filter(abs(logFC)>0.5)
RP.DEG_chronic_gse73878_0.5 <- RP.DEG_chronic_gse73878 %>% filter(abs(logFC)>0.5)

write.csv(unique(RP.DEG_chronic_gse27166_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_chronic_gse27166_0.5.csv")
write.csv(unique(RP.DEG_chronic_gse73878_0.5$`Gene symbol`), file = "01 differential expression-RP/RP.DEG_chronic_gse73878_0.5.csv")

######################################################################################
RPlist2df <- function(RPlist){
  a <- RPlist
  RPdf <- cbind(a$RPs, a$RPrank, a$pfp, a$pval, a$AveFC)
  colnames(RPdf) <- c("RPs_down", "RPs_up", "RPrank_down", "RPrank_up", "pfp_down", "pfp_up", "pval_down", "pval_up", "AveFC")
  RPdf <- as.data.frame(RPdf)
  RPdf$ID <- rownames(RPdf) 
  return(RPdf)
}

##########################################################################################
save(RP_gse14763_acute_df, RP_gse14763_latent_df, RP_gse27166_chronic_df, RP_gse49849_acute_df, RP_gse49849_latent_df, RP_gse73878_acute_df, RP_gse73878_chronic_df, RP_gse73878_latent_df, RP_gse88992_acute_df, file = "01 differential expression-RP/RP_gse_df.RData")

save(RP.DEG_acute_gse14763, RP.DEG_acute_gse49849, RP.DEG_acute_gse73878, RP.DEG_acute_gse88992, RP.DEG_latent_gse14763, RP.DEG_latent_gse49849, RP.DEG_latent_gse73878, RP.DEG_chronic_gse27166, RP.DEG_chronic_gse73878, file = "01 differential expression-RP/DEG_RP.RData")
