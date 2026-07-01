setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("initialize.R")

########################
########################
###
### code to generate Figure 8 (ci/ca vs ca for new model and USO under normal and reduced rubisco)
###
########################
########################



### set parameter values
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}

## load in data from von Caemmerer et al 2004
vv=read.delim("antisense.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
vv$trt <- ifelse(vv$trt=="reduced", "reduced Rubisco", "normal Rubisco")

## run simulations at a range of ca for each model, at normal (Vm=50) and -75% rubisco (Vm=12.5)
camin=50; camax=1000;
Vm25=50; y1 <- g.vs.(CA, PPFD, thetaa=thetaa, s=s, ro=ro, g1=g1, mBB=mBB, g0bb=g0bb, g0uso=g0uso, models=c("NEW", "USO"), noplot=TRUE)
Vm25=12.5; y2 <- g.vs.(CA, PPFD, thetaa=thetaa, s=s, ro=ro, g1=g1, mBB=mBB, g0bb=g0bb, g0uso=g0uso, models=c("NEW", "USO"), noplot=TRUE)
y1$vm = "normal"; y2$vm = "reduced"
y <- rbind(y1, y2)
y <- subset(y, ppfd==1000)
unique(y2$model)
unique(y$model)
y$vm <- factor(y$vm, levels=c("normal", "reduced"))
vv$trt <- factor(vv$trt, levels=c("normal Rubisco", "reduced Rubisco"))
ggplot() +
   geom_point(data=vv, aes(x=ca, y=cica, fill=trt), shape=21,color="transparent", size=2)+
   geom_line(data=y, aes(x=ca, y=ci/ca, color=as.factor(vm), linetype=model)) +
   scale_size_manual(values=c(1, 0.5, 0.5))+
   guides(size="none")+
   scale_shape_manual(values=c(19, 21))+
  scale_color_manual(values=c("black", "red"),
                      labels=c("normal Rubisco", "reduced Rubisco"))+
   theme(legend.text = element_text(hjust = 0))+
   scale_fill_manual(values=c("black", "red"))+
   scale_linetype_manual(values=c("solid", "dashed"))+
   ylim(0.6, 1.1) + xlim(0, 1000)+
   guides(fill=guide_legend(title="data", order=1))+
  guides(color=guide_legend(title="simulations", order=2))+
  guides(linetype=guide_legend(title="model", order=3))+
  labs(linetype=NULL)+
  xlab(varnamesb[2])+
   ylab(expression(paste(italic("c")["i"],"/",italic("c")["a"]))) + theme_bw() + theme(legend.key.width=unit(1.5, "cm"))


ggsave("Figure 8.png", device="png", dpi=1000)
