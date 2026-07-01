setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(dplyr)
library(minpack.lm)
library(viridis)



##########################
##########################
###
### code to fit models to Lin et al 2015 dataset and generated associated figures
###
##########################
##########################



## load Lin et al 2015 data
x=read.delim("lin2015data.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
x$spp_site <- paste0(x$Species, "_", x$Location)
x$Plantform <- ifelse(x$Plantform=="savanna", "savanna tree", x$Plantform)
x$ci <- x$CO2S - 1.6*x$Photo/x$Cond ## calculate ci and ci/ca ratio
x$cica <- x$ci/x$CO2S

## calculate gamma star and ci at the V/J transition point
x$citx <- NA;
x$gams <- NA
for(r in 1:nrow(x)) {
   if((x$Tleaf[r] > 0) & (!is.na(x$Tleaf[r]))) {
      TleafK <- x$Tleaf[r] + 273.15
   } else {
      TleafK <- 25 + 273.15
   }
   iRT <- 1/(8.3145e-3*TleafK) #Rgas is in kJ/mol K
   Kc <- exp(38.05 - 79.43*iRT)
   iKo <- exp(-20.30 + 36.38*iRT)  #1/Ko
   VJm <- 50/85 #exp(26.35 - 65.33*iRT)/exp(17.57 - 43.54*iRT)
   gams <- exp(13.49 - 24.46*iRT) # gamma star
   kp <- Kc*(1 + 10*21*iKo) #oxygen converted from % to ppt
   citx <- (kp - 8*gams*VJm)/(4*VJm - 1)
   
   x$citx[r] <- citx
   x$gams[r] <- gams
}


## calculate residual photosynthetic capacity for new model
## this calculation of r assumes PS is light-limited
x$rp <- x$Photo*3*x$gams/(x$ci - x$gams)
x$patm <- 101.325*((1 - 2.25577e-5*x$e)^5.25588) ## atmospheric pressure in kPa
x$deltaw <- 1000*x$VPD/x$patm ## deltaw in mmol/mol


### remove impossible points
### ci cannot be < gams in C3 species unless A<0 and ca<gams, and A is never <0 in this dataset
x <- subset(x, !(ci<gams & Pathway=="C3"))

### eqns above don't work for C4, remove these points
x <- subset(x, Pathway!="C4") 

### store pre-processed df
xoriginal <- x

## make df with one row per species-site combination, and initial estimates of Em and k ('k' here is called 'p' in the paper) 
u <- x %>% group_by(spp_site) %>% 
   summarize(nr=n(), 
             Em.est = max(Cond*deltaw), 
             k.est = 20*mean(rp)/3) ## (2/3)*deltaw.ref*r.ref = implication of m/gsref=0.6 in Oren

## bring these paras into main df
x <- left_join(x, u, by="spp_site")

## toss spp x site combos w/ fewer than 3 points bc these will be overfitted by new model
x <- subset(x, nr>3)
u <- subset(u, nr>3)


## fit USO and new model to data
u$g1 <- NA
u$Em <- NA
u$k <- NA
u$fitmode <- NA
for(ss in 1:nrow(u)) {
  ## subset a given spp x site combo
   xx <- subset(x, spp_site==u$spp_site[ss])
   
   ## solve USO model for g1, get initial estimate of g1 as avg of this quantity applied to the real data
   xx$g1est <- (xx$Cond*xx$CO2S/(1.6*xx$Photo) - 1)*sqrt(xx$VPD)
   g1.est = mean(xx$g1est)

   ## fit the model using nls()
   try(mu <- nls(Cond ~ (1.6*Photo/CO2S)*(g1/sqrt(VPD) + 1), 
            data=xx,
            start=list(g1=3)))
   if(!is.null(mu)) {
      u$g1[ss] = coef(mu)[[1]]
   }

   ## try to fit new model, using initial estimates of Em and k calculated earlier
   Em.est = u$Em.est[ss]
   k.est = u$k.est[ss]
   bob <- try(mn <- nls(Cond ~ Em*rp/(k + rp*deltaw),
                        data=xx,
                        start=list(Em=Em.est, k=k.est)),
              silent=T)
   if(class(bob) != "try-error") { 
      ## if that worked, keep results
      u$Em[ss] = coef(mn)[[1]]
      u$k[ss] = coef(mn)[[2]]
      u$fitmode[[ss]] = 1
   } else { 
      ## if it didn't work, try alternative approaches
      ## fit line to 1/g in new model
      ##  1/g = (p/Em)*(1/r) + (1/Em)*deltaw
      xx$ig <- 1/xx$Cond
      xx$ir <- 1/xx$rp
      m <- lm(ig ~ deltaw + ir, data=xx)
      Em.est2 = 1/coef(m)[[2]]
      k.est2 = coef(m)[[3]]/coef(m)[[2]]

      ## try using these coeffs as initial values in nls()
      bob <- try(mn <- nls(Cond ~ Em*rp/(k + rp*deltaw),
                           data=xx,
                           start=list(Em=Em.est2, k=k.est2)),
                 silent=T)
      if(class(bob) != "try-error") {
         ## if that worked, keep the results
         u$Em[ss] = coef(mn)[[1]]
         u$k[ss] = coef(mn)[[2]]
         u$fitmode[[ss]] = 2
      } else {
         ## otherwise, keep the results from the inverse linear fit
         u$Em[ss] = Em.est2
         u$k[ss] = k.est2
         u$fitmode[[ss]] = 3
      }
   }
   
} ## end for ss

## pull fitted paras into original df
x <- left_join(x, u[,c(1,5,6,7,8)], by="spp_site")

## calculate g and residuals for both models
x$gu <- (1.6*x$Photo/x$CO2S)*(x$g1/sqrt(x$VPD) + 1)
x$gn <- x$Em*x$rp/(x$k + x$rp*x$deltaw)
x$resid.u <- x$Cond - x$gu
x$resid.n <- x$Cond - x$gn

## create separate dfs for USO and new model, then rbind into long format
xu <- x[,c(1:40,42)]; colnames(xu)[40:41] <- c("g", "resid"); xu$model <- "USO"
xn <- x[,c(1:39,41,43)]; colnames(xn)[40:41] <- c("g", "resid"); xn$model <- "NEW"
xu$npara <- 1
xn$npara <- 2
x <- rbind(xu, xn)

## fix labels
x$Tregion <- ifelse(x$Tregion=="Arctic", "arctic", x$Tregion)
x$model <- ifelse(x$model=="NEW", "new model", "USO model")

## remove points for which g<0 in new model
x <- subset(x, !(model=="new model" & g<0))

x$modellabel <- ifelse(x$model=="new model", "(a) new model", "(b) USO model")


#########################
#########################
###
### make Figure 4, showing measured vs modeled gs for each model
###
#########################
#########################

ggp <- ggplot() +
  geom_point(data=x, aes(y=Cond, x=g), alpha=0.1, size=1, color="grey70")+
  geom_smooth(data=x, aes(y=Cond, x=g), method="lm", color="black", linewidth=1, linetype="dashed") + 
  xlim(0,1.25) + ylim(0,1.25) +
  geom_abline(slope=1, intercept=0) + theme_bw() + 
  xlab(expression(paste("measured ",italic("g")["sw"]," (mol ",m^-2,s^-1,")")))+
  ylab(expression(paste("predicted ",italic("g")["sw"]," (mol ",m^-2,s^-1,")"))) +
  guides(color=guide_legend(title="region"))+
  facet_wrap(facets="modellabel", nrow=1, ncol=2); ggp
ggsave("Figure 4.png", dpi=600, device="png")


## force Tregion labels to lowercase
x$Tregion <- tolower(x$Tregion)
xu$Tregion <- tolower(xu$Tregion)
xn$Tregion <- tolower(xn$Tregion)



#########################
#########################
###
### make Fig S2, illustrating why fitted Em and p sometimes get very large
###
#########################
#########################

xn.original <- xn
xu.original <- xu

xn <- subset(xn, (Wregion!="arid") | (Wregion=="arid" & Em<1e3))
xu <- subset(xu, (Wregion!="arid") | (Wregion=="arid" & Em<1e3))

xs <- subset(x, Species=="Eucalyptus aparrerinja_Boulia" | Species=="Eucalyptus terminalis_Boulia" | Species=="Senegalia senegal")
us <- unique(xs$Species)
xs$Species <- ifelse(xs$Species==us[1], "Senegalia senegal",
                     ifelse(xs$Species==us[2], "Eucalyptus aparrerinja", "Eucalyptus terminalis"))
pv <- ggplot() +
   geom_point(data=xs, aes(x=VPD, y=Cond, color=Species, shape=Species), size=2) +
   xlim(0,NA) + ylim(0, NA) + theme_bw() + xlab("VPD (kPa)") +
   ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")")))
pc <- ggplot() +
   geom_point(data=xs, aes(x=CO2S, y=Cond, color=Species, shape=Species), size=2) +
   xlim(0,NA) + ylim(0, NA) + theme_bw() + xlab(expression(paste(italic("c")["a"]," (ppm)"))) +
   ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")")))
pr <- ggplot() +
   geom_point(data=xs, aes(x=rp, y=Cond, color=Species, shape=Species), size=2) +
   xlim(0,NA) + ylim(0, NA) + theme_bw() + xlab(expression(paste(italic("A")["m"],"-",italic("A")," (",mu,"mol ",m^-2,s^-1,")"))) +
   ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")")))
library(ggpubr)
ggarrange(pv, pr, pc, nrow=3, common.legend=T, legend="right", labels=c("(a)", "(b)", "(c)"), label.x=c(0.22,0.22,0.22), label.y=c(0.95, 0.95, 0.95))
ggsave("Figures S2.png", device="png", dpi=600)




#########################
#########################
###
### make Figure 5, showing fitted para values
###
#########################
#########################


#### exclude fitmode=3 (where nls couldn't find solution and we kept initial linear estimates)
xu <- subset(xu, fitmode!=3)
xn <- subset(xn, fitmode!=3)

## summarize by spp x site combo
xxu <- xu %>% group_by(spp_site) %>% summarize(g1. = median(g1),
                                               Wregion. = unique(Wregion)[1],
                                               Plantform. = unique(Plantform)[1],
                                               Type. = unique(Type)[1],
                                               Tregion. = unique(Tregion)[1])
xxn <- xn %>% group_by(spp_site) %>% summarize(Em. = median(Em), k. = median(k),
                                               Wregion. = unique(Wregion)[1],
                                               Plantform. = unique(Plantform)[1],
                                               Type. = unique(Type)[1],
                                               Tregion. = unique(Tregion)[1])
xxun <- cbind(xxu, xxn[,2:3])
xxun <- xxun[,c(1,3,4,5,6,2,7,8)]
library(tidyr)
xxun <- pivot_longer(data=xxun, cols=c(6,7,8), names_to="var", values_to="val")

## summarize further by PFT, aridity, latitudinal zone and gymno v angio
xxu.Wregion <- xxu %>% group_by(Wregion.) %>% summarize(g1 = median(g1.))
xxu.Tregion <- xxu %>% group_by(Tregion.) %>% summarize(g1 = median(g1.))
xxu.Plantform <- xxu %>% group_by(Plantform.) %>% summarize(g1 = median(g1.))
xxu.Type <- xxu %>% group_by(Type.) %>% summarize(g1 = median(g1.))

xxn.Wregion <- xxn %>% group_by(Wregion.) %>% summarize(Em = median(Em.), k = median(k.))
xxn.Tregion <- xxn %>% group_by(Tregion.) %>% summarize(Em = median(Em.), k = median(k.))
xxn.Plantform <- xxn %>% group_by(Plantform.) %>% summarize(Em = median(Em.), k = median(k.))
xxn.Type <- xxn %>% group_by(Type.) %>% summarize(Em = median(Em.), k = median(k.))

xx.Wregion <- cbind(xxu.Wregion, xxn.Wregion[,2:3]); xx.Wregion <- pivot_longer(data=xx.Wregion, cols=-1, names_to="var", values_to="val")
xx.Tregion <- cbind(xxu.Tregion, xxn.Tregion[,2:3]); xx.Tregion <- pivot_longer(data=xx.Tregion, cols=-1, names_to="var", values_to="val")
xx.Plantform <- cbind(xxu.Plantform, xxn.Plantform[,2:3]); xx.Plantform <- pivot_longer(data=xx.Plantform, cols=-1, names_to="var", values_to="val")
xx.Type <- cbind(xxu.Type, xxn.Type[,2:3]); xx.Type <- pivot_longer(data=xx.Type, cols=-1, names_to="var", values_to="val")

xxune <- subset(xxun, var=="Em.")
xxunk <- subset(xxun, var=="k.")
xxung <- subset(xxun, var=="g1.")

xxune$Wregion. <- factor(xxune$Wregion., levels=c("arid", "semi-arid", "dry sub-humid", "humid"))
xxunk$Wregion. <- factor(xxunk$Wregion., levels=c("arid", "semi-arid", "dry sub-humid", "humid"))
xxung$Wregion. <- factor(xxung$Wregion., levels=c("arid", "semi-arid", "dry sub-humid", "humid"))
xxn.Wregion$Wregion. <- factor(xxn.Wregion$Wregion., levels=c("arid", "semi-arid", "dry sub-humid", "humid"))
xxu.Wregion$Wregion. <- factor(xxu.Wregion$Wregion., levels=c("arid", "semi-arid", "dry sub-humid", "humid"))

## make plots for each para and each type of classification
pew <- ggplot() +
   geom_jitter(data=xxune, aes(y=Wregion., x=val, fill=Wregion.), size=0.5, alpha=0.6)+
   geom_violin(data=xxune, aes(y=Wregion., x=val, fill=Wregion.), alpha=0.6, adjust=1)+ xlim(0,40) + theme_bw() +
   geom_errorbar(data=xxn.Wregion, aes(y=Wregion., xmin=0.9999*Em, xmax=1.0001*Em, color=Wregion.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab(expression(paste(italic("E")["m"]," (mmol ",m^-2,s^-1,")"))) + 
   guides(fill="none", color="none") + ylab(""); 

pkw <- ggplot() +
   geom_jitter(data=xxunk, aes(y=Wregion., x=val, fill=Wregion.), size=0.5, alpha=0.6)+
   geom_violin(data=xxunk, aes(y=Wregion., x=val, fill=Wregion.), alpha=0.6, adjust=1)+ xlim(0,300) + theme_bw() +
   geom_errorbar(data=xxn.Wregion, aes(y=Wregion., xmin=0.9999*k, xmax=1.0001*k, color=Wregion.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab(expression(paste(italic("p")," (nmol ",m^-2,s^-1,")"))) +
   guides(fill="none", color="none") + ylab(""); 

pgw <- ggplot() +
   geom_jitter(data=xxung, aes(y=Wregion., x=val, fill=Wregion.), size=0.5, alpha=0.6)+
   geom_violin(data=xxung, aes(y=Wregion., x=val, fill=Wregion.), alpha=0.6, adjust=1)+ xlim(0,13) + theme_bw() +
   geom_errorbar(data=xxu.Wregion, aes(y=Wregion., xmin=0.9999*g1, xmax=1.0001*g1, color=Wregion.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab(expression(paste(italic("g")["1"]," (kP",a^0.5,")"))) +
   guides(fill="none", color="none") + ylab(""); 



pet <- ggplot() +
   geom_jitter(data=xxune, aes(y=Tregion., x=val, fill=Tregion.), size=0.5, alpha=0.6)+
   geom_violin(data=xxune, aes(y=Tregion., x=val, fill=Tregion.), alpha=0.6, adjust=1)+ xlim(0,40) + theme_bw() +
   geom_errorbar(data=xxn.Tregion, aes(y=Tregion., xmin=0.9999*Em, xmax=1.0001*Em, color=Tregion.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab("") + #xlab(expression(paste(italic("E")["m"]," (mmol ",m^-2,s^-1,")"))) + 
   guides(fill="none", color="none") + ylab(""); 

pkt <- ggplot() +
   geom_jitter(data=xxunk, aes(y=Tregion., x=val, fill=Tregion.), size=0.5, alpha=0.6)+
   geom_violin(data=xxunk, aes(y=Tregion., x=val, fill=Tregion.), alpha=0.6, adjust=1)+ xlim(0,300) + theme_bw() +
   geom_errorbar(data=xxn.Tregion, aes(y=Tregion., xmin=0.9999*k, xmax=1.0001*k, color=Tregion.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab("") + #xlab(expression(paste(italic("p")," (nmol ",m^-2,s^-1,")"))) +
   guides(fill="none", color="none") + ylab(""); 

pgt <- ggplot() +
   geom_jitter(data=xxung, aes(y=Tregion., x=val, fill=Tregion.), size=0.5, alpha=0.6)+
   geom_violin(data=xxung, aes(y=Tregion., x=val, fill=Tregion.), alpha=0.6, adjust=1)+ xlim(0,13) + theme_bw() +
   geom_errorbar(data=xxu.Tregion, aes(y=Tregion., xmin=0.9999*g1, xmax=1.0001*g1, color=Tregion.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab("") + #xlab(expression(paste(italic("g")["1"]," (kP",a^0.5,")"))) +
   guides(fill="none", color="none") + ylab(""); 





pep <- ggplot() +
   geom_jitter(data=xxune, aes(y=Plantform., x=val, fill=Plantform.), size=0.5, alpha=0.6)+
   geom_violin(data=xxune, aes(y=Plantform., x=val, fill=Plantform.), alpha=0.6, adjust=1)+ xlim(0,40) + theme_bw() +
   geom_errorbar(data=xxn.Plantform, aes(y=Plantform., xmin=0.9999*Em, xmax=1.0001*Em, color=Plantform.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab("")+#xlab(expression(paste(italic("E")["m"]," (mmol ",m^-2,s^-1,")"))) + 
   guides(fill="none", color="none") + ylab(""); 

pkp <- ggplot() +
   geom_jitter(data=xxunk, aes(y=Plantform., x=val, fill=Plantform.), size=0.5, alpha=0.6)+
   geom_violin(data=xxunk, aes(y=Plantform., x=val, fill=Plantform.), alpha=0.6, adjust=1)+ xlim(0,300) + theme_bw() +
   geom_errorbar(data=xxn.Plantform, aes(y=Plantform., xmin=0.9999*k, xmax=1.0001*k, color=Plantform.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab("")+#xlab(expression(paste(italic("p")," (nmol ",m^-2,s^-1,")"))) +
   guides(fill="none", color="none") + ylab(""); 

pgp <- ggplot() +
   geom_jitter(data=xxung, aes(y=Plantform., x=val, fill=Plantform.), size=0.5, alpha=0.6)+
   geom_violin(data=xxung, aes(y=Plantform., x=val, fill=Plantform.), alpha=0.6, adjust=1)+ xlim(0,13) + theme_bw() +
   geom_errorbar(data=xxu.Plantform, aes(y=Plantform., xmin=0.9999*g1, xmax=1.0001*g1, color=Plantform.), linewidth=1)+
   scale_fill_viridis(discrete=T) +scale_color_viridis(discrete=T) + xlab("")+#xlab(expression(paste(italic("g")["1"]," (kP",a^0.5,")"))) +
   guides(fill="none", color="none") + ylab(""); 





pey <- ggplot() +
   geom_jitter(data=xxune, aes(y=Type., x=val, fill=Type.), size=0.5, alpha=0.6)+
   geom_violin(data=xxune, aes(y=Type., x=val, fill=Type.), alpha=0.6, adjust=1)+ xlim(0,40) + theme_bw() +
   geom_errorbar(data=xxn.Type, aes(y=Type., xmin=0.9999*Em, xmax=1.0001*Em, color=Type.), linewidth=1)+
   scale_fill_viridis(discrete=T, begin=0.25, end=0.75) +scale_color_viridis(discrete=T, begin=0.25, end=0.75) + xlab("")+#xlab(expression(paste(italic("E")["m"]," (mmol ",m^-2,s^-1,")"))) + 
   guides(fill="none", color="none") + ylab(""); 

pky <- ggplot() +
   geom_jitter(data=xxunk, aes(y=Type., x=val, fill=Type.), size=0.5, alpha=0.6)+
   geom_violin(data=xxunk, aes(y=Type., x=val, fill=Type.), alpha=0.6, adjust=1)+ xlim(0,300) + theme_bw() +
   geom_errorbar(data=xxn.Type, aes(y=Type., xmin=0.9999*k, xmax=1.0001*k, color=Type.), linewidth=1)+
   scale_fill_viridis(discrete=T, begin=0.25, end=0.75) +scale_color_viridis(discrete=T, begin=0.25, end=0.75) + xlab("")+#xlab(expression(paste(italic("p")," (nmol ",m^-2,s^-1,")"))) +
   guides(fill="none", color="none") + ylab(""); 

pgy <- ggplot() +
   geom_jitter(data=xxung, aes(y=Type., x=val, fill=Type.), size=0.5, alpha=0.6)+
   geom_violin(data=xxung, aes(y=Type., x=val, fill=Type.), alpha=0.6, adjust=1)+ xlim(0,13) + theme_bw() +
   geom_errorbar(data=xxu.Type, aes(y=Type., xmin=0.9999*g1, xmax=1.0001*g1, color=Type.), linewidth=1)+
   scale_fill_viridis(discrete=T, begin=0.25, end=0.75) +scale_color_viridis(discrete=T, begin=0.25, end=0.75) + xlab("")+#xlab(expression(paste(italic("g")["1"]," (kP",a^0.5,")"))) +
   guides(fill="none", color="none") + ylab(""); 



labels <- c("(a) Em", "(b) p", "(c) g1",
            "(d) Em", "(e) p", "(f) g1",
            "(g) Em", "(h) p", "(i) g1",
            "(j) Em", "(k) p", "(l) g1")
library(egg)

ppp <- egg::ggarrange(
          tag_facet(tag_pool="a",x=33,pep),
          tag_facet(tag_pool="b",x=240,pkp + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())),
          tag_facet(tag_pool="c",x=11,pgp + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())), 
          tag_facet(tag_pool="d",x=32.5,pet),
          tag_facet(tag_pool="e",x=240,pkt + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())),
          tag_facet(tag_pool="f",x=11,pgt + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())), 
          tag_facet(tag_pool="g",x=32,pey),
          tag_facet(tag_pool="h",x=240,pky + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())),
          tag_facet(tag_pool="i",x=11,pgy + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())),
          tag_facet(tag_pool="j",x=33,pew),
          tag_facet(tag_pool="k",x=240,pkw + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())), 
          tag_facet(tag_pool="l",x=11,pgw + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())),
          nrow=4, ncol=3, heights=c(1.2,1,0.65,1)); ppp
ggsave("Figure 5.png", ppp, device="png", dpi=1000)

