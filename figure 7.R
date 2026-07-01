setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("initialize.R")


#########################
#########################
###
### code to generate figure 7 (gsw vs deltaw at night)
###
#########################
#########################

## set para values
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}
ppfd_=0
psisoil_ = 0; psisoilmin = 0; psisoilmax = 0

## run model
y <- g.vs.(DW, PSISOIL, thetaa=thetaa, s=s, ro=ro, g1=g1, mBB=mBB, g0bb=g0bb, g0uso=g0uso, nopsi=TRUE)

## generate figure
ggplot() +  
  geom_line(data=y, aes(x=deltaw, y=gsw, linetype=model), size=1)+
   # scale_color_manual(values=c("black", "red", "blue"))+
   scale_linetype_manual(values=c("solid", "longdash", "dotdash"))+
   scale_size_manual(values=c(1.2, 0.7, 0.7))+
   ylim(0,NA) + xlim(0, 35)+
   ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")")))+
   xlab(varnamesb[4]) + theme_bw() + theme(legend.key.width=unit(1.2, "cm"))

ggsave("Figure 7.png", device="png", dpi=1000)

