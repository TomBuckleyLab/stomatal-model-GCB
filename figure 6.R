setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("initialize.R")

#######################
#######################
###
### code to generate Figure 6, illustrating nonlinearity of the USO factor,
###  but not the residual photosynthetic capacity, wrt light
###
#######################
#######################


#### redefine code for USO model to take input values of thetaj and phij; otherwise same as before
gasxUSO <- function(PPFD, ca, TleafC, D, Vm25, Jm25, Rd25, thetaa, g1, g0, thetaj.input=0, phij.input=0) {
  h = 1 - 10*D/(6.1078*exp(17.27*TleafC/(237.3 + TleafC)))

  if(h>=0) {
    TleafK <- TleafC + 273.15
    iRT <- 1/(8.3145e-3*TleafK) #Rgas is in kJ/mol K
    Vm <- Vm25*exp(26.35 - 65.33*iRT)
    Jm <- Jm25*exp(17.57 - 43.54*iRT)
    Kc <- exp(38.05 - 79.43*iRT)
    iKo <- exp(-20.30 + 36.38*iRT)  #1/Ko
    kp <- Kc*(1 + 10*21*iKo) #oxygen converted from % to ppt
    gams <- exp(13.49 - 24.46*iRT) # gamma star
    Rd <- Rd25*exp(18.72 - 46.39*iRT)
    phiPSII <- 0.352 + 0.022*TleafC - 3.4e-4*TleafC*TleafC #electrons per photon
    phij <- ifelse(phij.input != 0, phij.input, phiPSII*0.9*0.5) # 0.9 is assumed absorptance; 0.5 is photon fraction
    thetaj <- ifelse(thetaj.input != 0, thetaj.input, 0.76 + 0.018*TleafC - 3.7e-4*TleafC*TleafC) #unitless curvature parameter
    i <- PPFD
    J <- ((-(-Jm-phij*i) - sqrt((-Jm-phij*i)*(-Jm-phij*i) - 4*thetaj*Jm*phij*i))/(2*thetaj))

    ## in BB, A/ca is multiplied by product of slope m and RH h
    ## in USO, the equivalent factor is 1.6*(1 + g1/sqrt(D)), where D is VPD in kPa
    ## models are otherwise equivalent
    mh = 1.6*(1 + g1/sqrt(D))

    if(J>0) {
      ## quadratic coeffs for V-lim case
      av = -mh*(Vm - Rd) - g0*ca
      bv = mh*(Vm - Rd)*ca + mh*(Vm*gams + Rd*kp) + g0*ca*ca - g0*ca*kp - 1.6*ca*(Vm - Rd)
      cv = -mh*(Vm*gams + Rd*kp) + g0*ca*ca*kp + 1.6*ca*(Vm*gams + Rd*kp)

      ## quadratic coeffs for J-lim case
      aj = -mh*(0.25*J - Rd) - g0*ca
      bj = mh*(0.25*J - Rd)*ca + mh*(0.25*J*gams + Rd*2*gams) + g0*ca*ca - g0*ca*2*gams - 1.6*ca*(0.25*J - Rd)
      cj = -mh*(0.25*J*gams + Rd*2*gams) + g0*ca*ca*2*gams + 1.6*ca*(0.25*J*gams + Rd*2*gams)

      ## ci solutions
      civ = (0.5/av)*(-bv - sqrt(bv*bv - 4*av*cv))
      cij = (0.5/aj)*(-bj - sqrt(bj*bj - 4*aj*cj))

      ## calculate PS
      Av = Vm*(civ - gams)/(civ + kp) - Rd
      Aj = 0.25*J*(cij - gams)/(cij + 2*gams) - Rd
      A = (0.5/thetaa)*((Av + Aj) - sqrt((Av + Aj)^2 - 4*thetaa*Av*Aj))

      lim = ifelse(Av < Aj, 1, 2)
    } else {
      lim = 2
      A = -Rd
    }
    g = mh*A/ca + g0
    g <- ifelse(g > 0, g, 0)
    ci = ca - 1.6*A/g

    return(c(g, A, ci, lim))
  } else {
    return(c(NA,NA,NA,NA))
  }
} ## end gasxUSO()


### set default parameter values
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; 
   g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 425; camax = 425; 
   TleafC_ = 25; TleafCmin = 25; TleafCmax = 25; deltaw_ = 10; deltawmin = 10; deltawmax = 10;
   psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9;
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin);
   varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)
}

## vm,rd,T,dw from data that gave most-negative USOfac, to get same x range
## run new model for normal Vm and saturating Vm (latter so PPFD always limiting)
Vm25=50; Rd25=3; y <- g.vs.(PPFD, CA, nx1=200, nx2=2, models=c("NEW"), nopsi=TRUE) 
Vm25=1000; yv <- g.vs.(PPFD, CA, nx1=200, nx2=2, models=c("NEW"), nopsi=TRUE)

## calculate USO factor
y$usofac <- 1.6*y$A/(sqrt(y$D)*y$ca)

## generate panel (a)
gb <- ggplot() +
   geom_line(data=yv, aes(x=1.6*A/(sqrt(D)*ca), y=r), size=1, linetype="dashed") + 
   geom_line(data=y, aes(x=1.6*A/(sqrt(D)*ca), y=r), size=1, linetype="solid") + 
   ylim(0,NA) +
   ylab(expression(paste(italic("A")["m"]," - ",italic("A")," (",mu,"mol ",m^-2,s^-1,")")))+
  xlab(expression(paste("USO factor (mol ",m^-2,s^-1,")")))+
   theme_bw()+ theme(text=element_text(size=12)) ; gb


#############################
#############################
###
### load data from Lamour et al 2022
###
#############################
#############################


xg=read.delim("lamour gs.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
xaci=read.delim("lamour aci.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
xaq=read.delim("lamour aq.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
xs=read.delim("lamour samples.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
xs$id <- paste0(xs$Vertical_Profile, "_", xs$Vertical_Elevation, "_", xs$Phenological_Stage, "_", xs$Species_Code)

xgs <- left_join(xg, xs, by="SampleID")
xacis <- left_join(xaci, xs, by="SampleID")
xaqs <- left_join(xaq, xs, by="SampleID")

yacis <- xacis %>% group_by(id) %>% summarize(vm=mean(Vcmax), jm=mean(Jmax), tpu=mean(TPU), rd=mean(Rday))
yaqs <- xaqs %>% group_by(id) %>% summarize(thetaa=mean(Theta), phic=mean(PhiCO2i), ci=mean(Ci))
xgss <- left_join(xgs, yacis, by="id")
xgsss <- left_join(xgss, yaqs, by="id")

x <- subset(xgsss, !is.na(vm) & !is.na(thetaa))
### assume thetaj = thetaa, phij = 0.25*phi_CO2

## generate predictions from new model; will override later
for(r in 1:nrow(x)) x$gswmod[r] <- gasxNEW(x$Qin[r], x$CO2s[r], x$Tleaf[r], x$VPDleaf[r]/x$Patm[r], 0, x$vm[r], x$jm[r], x$rd[r], 
                                           0.99, 10, -3, 0.1, 0, phij.input = x$phic[r]*4, thetaj.input = x$thetaa[r])


## wrapper to drive gasxNEW() w vector of three paras as first arg
paras_to_gasxNEW <- function(pvec, K25, ppfd, ca, tleaf, deltaw, vm, jm, rd, thetaa, phij, thetaj) {
   gasxNEW(ppfd, ca, tleaf, deltaw, 0, vm, jm, rd, thetaa, K25, pvec[1], pvec[2], pvec[3], phij.input=phij, thetaj.input=thetaj)
} 

## function to help fit new model to Lamour data by adjusting psic, s and b0
## (returns sse for given para vals)
fit_lamour <- function(pvec, xx) {
   err2 = 0
   for(r in 1:nrow(xx)) xx$gswmod[r] <- paras_to_gasxNEW(pvec, K25, xx$Qin[r], xx$CO2s[r], xx$Tleaf[r], xx$VPDleaf[r]/xx$Patm[r], xx$vm[r], 
                                                         xx$jm[r], xx$rd[r], 0.99, xx$phic[r]*4, xx$thetaa[r])
   if(!is.na(xx$gswmod[r]) & !is.na(xx$gsw[r])) err2 <- err2 + (xx$gswmod[r] - 0.001*xx$gsw[r])^2
   return(err2)
}

## for each tree in the Lamour data, fit model to the data and generate resulting predictions
y <- data.frame(); oo <- data.frame(); for(i in unique(x$SampleID)) {
   xx <- subset(x, SampleID==i)
   o <- optim(c(psic, s, ro), fit_lamour, xx=xx, method="L-BFGS-B", lower=c(-5, 0.001, 0), upper=c(-0.1, 0.2, 1))
   for(r in 1:nrow(xx)) xx$gswmod[r] <- paras_to_gasxNEW(o$par, K25, xx$Qin[r], xx$CO2s[r], xx$Tleaf[r], xx$VPDleaf[r]/xx$Patm[r], xx$vm[r], 
                                                         xx$jm[r], xx$rd[r], 0.99, xx$phic[r]*4, xx$thetaa[r])
   oo <- rbind(oo,o$par)
   y <- rbind(y, xx)
}


## as above but for USO model
paras_to_gasxUSO <- function(pvec, ppfd, ca, tleaf, VPD, vm, jm, rd, thetaa, phij, thetaj) {
  gasxUSO(ppfd, ca, tleaf, VPD, vm, jm, rd, thetaa, pvec[1], pvec[2], phij.input=phij, thetaj.input=thetaj)
}

fit_lamour_uso <- function(pvec, xx) {
  err2 = 0
  for(r in 1:nrow(xx)) xx$gswmod[r] <- paras_to_gasxUSO(pvec, xx$Qin[r], xx$CO2s[r], xx$Tleaf[r], xx$VPDleaf[r], xx$vm[r], 
                                                        xx$jm[r], xx$rd[r], 0.99, xx$phic[r]*4, xx$thetaa[r])
  if(!is.na(xx$gswmod[r]) & !is.na(xx$gsw[r])) err2 <- err2 + (xx$gswmod[r] - 0.001*xx$gsw[r])^2
  return(err2)
}

g1=3; g0=0.001 ## initial guesses
yu <- data.frame(); oo <- data.frame(); for(i in unique(x$SampleID)) {
  xx <- subset(x, SampleID==i)
  o <- optim(c(g1, g0), fit_lamour_uso, xx=xx, method="L-BFGS-B", lower=c(1, 0), upper=c(10, 0.1))
  for(r in 1:nrow(xx)) xx$gswmod[r] <- paras_to_gasxUSO(o$par, xx$Qin[r], xx$CO2s[r], xx$Tleaf[r], xx$VPDleaf[r], xx$vm[r], 
                                                        xx$jm[r], xx$rd[r], 0.99, xx$phic[r]*4, xx$thetaa[r])
  oo <- rbind(oo,o$par)
  yu <- rbind(yu, xx)
}

## now generate panels b and c and assemble into overall figure

z <- y %>% group_by(Species_Code, SampleID) %>% summarize(n = n())
length(unique(y$SampleID))
length(unique(y$Species_Code))
length(unique(y$Vertical_Elevation))
y$usofac <- 1.6*y$A/(sqrt(y$VPDleaf)*y$CO2s)
p.ga <- ggplot() +
  geom_point(data=y, aes(x=1.6*A/(sqrt(VPDleaf)*CO2s), y=gsw/1000), ) +
  geom_smooth(data=y, aes(x=1.6*A/(sqrt(VPDleaf)*CO2s), y=gsw/1000), method="lm") +
  geom_line(data=y, aes(x=1.6*A/(sqrt(VPDleaf)*CO2s), y=gsw/1000, group=SampleID))+
  ylab(expression(paste("observed ", italic("g")["sw"]," (mol ",m^-2,s^-1,")")))+
  xlab(expression(paste("USO factor (mol ",m^-2,s^-1,")")))+
  theme_bw()+ theme(text=element_text(size=12)); p.ga
p.gg <- ggplot() +
   geom_point(data=y, aes(x=gswmod, y=gsw/1000))+
   geom_smooth(data=y, aes(x=gswmod, y=gsw/1000), method="lm")+
   geom_line(data=y, aes(x=gswmod, y=gsw/1000, group=SampleID))+
   ylab(expression(paste("observed ",italic("g")["sw"]," (mol ",m^-2,s^-1,")")))+
   xlab(expression(paste(italic("g")["sw"]," from Eqn 3 (mol ",m^-2,s^-1,")")))+
   theme_bw()+ theme(text=element_text(size=12)); p.gg


ggarrange(gb, p.ga, p.gg, ncol=3, nrow=1,
          labels=c("(a)", "(b)", "(c)"), label.x=c(0.19,0.22,0.22), label.y=c(0.98,0.98,0.98),
          font.label=list(face="plain"))
ggsave("Figure 6.png", device="png", dpi=1000)
