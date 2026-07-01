setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("initialize.R")

######################
######################
###
### code to generate Figure 1 ('standard' response of gsw to various factors)
###
######################
######################


##### deltaw responses at different Kplants
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}
K25=1; y <- g.vs.(DW, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW")); y <- subset(y, TleafC==25); y$K <- "low"; yk1 <- y
K25=3; y <- g.vs.(DW, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW")); y <- subset(y, TleafC==25); y$K <- "med"; yk2 <- y
K25=5; y <- g.vs.(DW, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW")); y <- subset(y, TleafC==25); y$K <- "high"; yk3 <- y
K25=3 
yk1$grel <- yk1$gsw/yk1$gsw[yk1$deltaw==15]
yk2$grel <- yk2$gsw/yk2$gsw[yk2$deltaw==15]
yk3$grel <- yk3$gsw/yk3$gsw[yk3$deltaw==15]
y <- rbind(yk1, yk2, yk3)
y$K <- factor(y$K, levels=c("low", "med", "high"))
p.dw <- ggplot() +  geom_line(data=subset(y, K != "med"), aes(x=deltaw, y=grel, linetype=K), size=1)+
   scale_color_viridis(discrete=T, option="D", begin=0, end=0.7)+
   scale_linetype_manual(values=c("longdash", "solid", "dotdash"))+
   ylim(0,NA) + xlim(0, 35)+
   ylab("relative stomatal conductance")+
   xlab(varnamesb[4])+
   guides(linetype=guide_legend(title="hydraulic\nconductance")) +
   theme_bw() + theme(legend.key.width=unit(1.5, "cm")); p.dw


## psisoil responses at different K
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}
Ki = K25
K25=Ki; y <- g.vs.(PSISOIL, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"),noplot=TRUE); y <- subset(y, TleafC==25); y$K <- "100%"; y$Kp <- 1; yk1 <- y
K25=0.8*Ki; y <- g.vs.(PSISOIL, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"),noplot=TRUE); y <- subset(y, TleafC==25); y$K <- "80%"; y$Kp <- 0.8;  yk2 <- y
K25=0.6*Ki; y <- g.vs.(PSISOIL, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"),noplot=TRUE); y <- subset(y, TleafC==25); y$K <- "60%"; y$Kp <- 0.6;  yk3 <- y
K25=0.4*Ki; y <- g.vs.(PSISOIL, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"),noplot=TRUE); y <- subset(y, TleafC==25); y$K <- "40%"; y$Kp <- 0.4;  yk4 <- y
K25=0.2*Ki; y <- g.vs.(PSISOIL, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"),noplot=TRUE); y <- subset(y, TleafC==25); y$K <- "20%"; y$Kp <- 0.2;  yk5 <- y
K25=Ki
y <- rbind(yk1, yk2, yk3, yk4, yk5)
y$K <- factor(y$K, levels=c("100%", "80%", "60%", "40%", "20%"))
z <- subset(y, (psisoil==0 & K=="100%") | (round(psisoil,2)==-0.4 & K=="80%") | (psisoil==-0.8 & K=="60%") | (psisoil==-1.2 & K=="40%") | (psisoil==-1.6 & K=="20%") | (psisoil==-2 & K=="100%"))
p.psis <- ggplot() +  geom_line(data=y, aes(x=psisoil, y=gsw, linetype=K), size=1)+
   geom_point(data=z, aes(x=psisoil, y=gsw), size=2, color="red")+
  geom_smooth(data=z, aes(x=psisoil, y=gsw), size=0.5, color="red", method="lm", formula=y~poly(x,3), se=F)+
  scale_color_viridis(discrete=T, option="D", begin=0, end=0.7)+
   ylim(0,NA) + 
   ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")")))+
   xlab(varnamesb[5])+
   guides(linetype=guide_legend(title="hydraulic\nconductance")) +
   theme_bw() + theme(legend.key.width=unit(1.5, "cm")); p.psis


## light responses at different photosynthetic capacities
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}
Vm25=10; Jm25=1.7*Vm25; Rd25=0.01*Vm25; y <- g.vs.(PPFD, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"), noplot=TRUE); y <- subset(y, TleafC==25); y$Vm <- "low"; yk1 <- y
Vm25=50; Jm25=1.7*Vm25; Rd25=0.01*Vm25; y <- g.vs.(PPFD, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"), noplot=TRUE); y <- subset(y, TleafC==25); y$Vm <- "med"; yk2 <- y
Vm25=100; Jm25=1.7*Vm25; Rd25=0.01*Vm25; y <- g.vs.(PPFD, TLEAF, thetaa=thetaa, s=s, ro=ro, models=c("NEW"), noplot=TRUE); y <- subset(y, TleafC==25); y$Vm <- "high"; yk3 <- y
Vm25=50; Jm25=1.7*Vm25; Rd25=0.01*Vm25; 
y <- rbind(yk1, yk2, yk3)
y$Vm <- factor(y$Vm, levels=c("high", "med", "low"))
p.ppfd <- ggplot() +  geom_line(data=y, aes(x=ppfd, y=gsw, linetype=Vm), size=1)+
   scale_color_viridis(discrete=T, option="D", begin=0, end=0.7)+
   scale_linetype_manual(values=c("solid", "longdash", "dotdash"))+
   ylim(0,NA) + 
   xlab(varnamesb[1])+
   ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")")))+
   guides(linetype=guide_legend(title="carboxylation\ncapacity")) +
  theme_bw()+theme(legend.key.width=unit(1.5, "cm")); p.ppfd


## CO2 response at different PPFDs
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}
ppfdmin=100; ppfdmax=1000;
y <- g.vs.(CA, PPFD, thetaa=thetaa, s=s, ro=ro, models=c("NEW"), nx2=2, noplot=TRUE)
y <- y %>% group_by(ppfd) %>% mutate(grel = gsw/gsw[ca==425])
p.ca <- ggplot() +  geom_line(data=y, aes(x=ca, y=grel, linetype=as.factor(ppfd)), size=1)+
   scale_color_viridis(discrete=T, option="D", begin=0, end=0.7)+
   scale_linetype_manual(values=c("longdash", "solid", "dotdash"))+
   ylim(0,2) + 
   xlab(varnamesb[2])+
   ylab("relative stomatal conductance")+
   guides(linetype=guide_legend(title=expression(
      atop(paste("     PPFD     "),paste("(",mu,"mol ",m^-2,s^-1,")")))
      )) +
   theme_bw() + theme(legend.key.width=unit(1.5, "cm")); p.ca
ppfdmin=0; ppfdmax=1000;


### compile four general responses into single figure
ggarrange(p.ppfd, p.ca, p.dw, p.psis, ncol=2, nrow=2, widths=c(1,1), common.legend=FALSE,
          labels=c("(a)", "(b)", "(c)", "(d)"), label.x=c(0.15,0.15,0.15,0.15), label.y=c(0.98,0.98,0.98,0.98),
          font.label=list(face="plain"), align="hv")
ggsave("Figure 1.png", device="png", dpi=600)

