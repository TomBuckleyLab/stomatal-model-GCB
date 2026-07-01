setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(dplyr)
library(mgcv)
library(ggpubr)
library(pracma)


#########################
#########################
###
### code to generate Figure 3 (test of new model and USO vs heatwave data from Drake/Sicangco)
###
#########################
#########################


### function to fit weibull VC, following sicangco
fitVC <- function(P50, P88) {
  c = -1.118053/(log(P50) - log(P88))
  b = -P50/(0.6931472^(1/c))
  return(c(b, c))
}


### fn to numerically integrate VC to get actual Krel
VC <- function(psileaf, psisoil, b, c) {
  P <- seq(-psisoil, -psileaf, length.out=500) 
  x = exp(-((-P/b)^c)) #local vc
  # Erel = pracma::trapz(P, x) (E relative to value at Kmax)
  # Kplantrel = Erel/(psisoil - psileaf)
  return(pracma::trapz(P, x)/(psisoil - psileaf))
}


## fn to get relative humidity given deltaw and T
getRH <- function(deltaw, temperature) {
  return(1 - deltaw/(6.1078*exp(17.27*temperature/(237.3 + temperature)))) 
}


## fn to get deltaw given relative humidity and T
getDW <- function(h, temperature) {
  return((1 - h)*6.1078*exp(17.27*temperature/(237.3 + temperature))) 
}


## fn to get boundary layer conductance to water vapor given air and leaf T, wind speed, and leaf width
get_gbw <- function(TairC, TleafC, vwind, dleaf) {
  Dh = 2.15e-5 ## diffusivity for heat, m2/s
  TaK = TairC + 273.15 ## air T in kelvins
  gbhw = 0.003*sqrt(vwind/dleaf) ## forced convection; Leuning et al 1995 eqn E1
  Gr = 1.6e8*abs(TleafC - TairC)*(dleaf^3) ## Grashof number; Leuning E4
  gbhf = 0.5*Dh*(Gr^0.25)/dleaf ## free convection, Leuning E3
  gbh = (gbhw + gbhf)*Patm/(8.31446*TaK) ## converted to mol
  gbw = gbh/1.075
  return(gbw)
}


############
### gasxNEW_mod()
###
###    MODIFIED FOR THE DRAKE ET AL 2018 EXPT W/ T RESPONSES GIVEN BY SICANGCO ET AL 2026
###     - due to VC, requires iterative solution, so K25 after VC is an output (Kvc25, 4th element in return vector)
###     - also has numerous differences in PS T responses to mimic plantecophys
###     - also has been updated to include the effect of boundary layer resistance (rbc)
###
gasxNEW_mod <- function(PPFD, ca, TleafC, deltaw, psisoil, rbc, Vm25, Jm25, thetaa, K25, psic, s, ro, Tcrit, T50) {
  h = 1 - deltaw/(6.1078*exp(17.27*TleafC/(237.3 + TleafC)))
  
  if(h>=0) {
    TleafK <- TleafC + 273.15
    iRT <- 1/(8.3145e-3*TleafK) #Rgas is in kJ/mol K
    
    ## T responses used by Sicangco et al; paras from their Table 2; note Hd = "Ed" in table but "Hd" in text
    rgas = 8.314 ## their eqns use rgas in J/mol K; Ea and Hd are in kJ/mol but X 1000 converts them to J/mol
    Eav = 62307; Eaj = 33115
    dSv = 639; dSj = 635
    Hdv = 2e5; Hdj = 2e5
    Vm = Vm25*exp(Eav*(TleafK-298.15)/(298.15*rgas*TleafK))*
      (1 + exp(dSv/rgas - Hdv/(298.15*rgas)))/
      (1 + exp(dSv/rgas - Hdv/(TleafK*rgas)))

    Jm = Jm25*exp(Eaj*(TleafK-298.15)/(298.15*rgas*TleafK))*
      (1 + exp(dSj/rgas - Hdj/(298.15*rgas)))/
      (1 + exp(dSj/rgas - Hdj/(TleafK*rgas)))
    
    r=2/(T50 - Tcrit)
    TC = 1/(1 + exp(-r*(TleafC-T50)))
    Jm = Jm*(1 - TC)

    Rd = 0.92*(1.92^((TleafC-25)/10)) ## formulation used in plantecophys  
    Kc <- 404.9*exp((79430/8.31446)*(1/298.15 - 1/TleafK)) ## formulation used in plantecophys
    Ko <- 278.4*exp((36380/8.31446)*(1/298.15 - 1/TleafK)) ## formulation used in plantecophys
    kp = Kc*(1 + 10*21/Ko)
    gams <- 42.75*exp((37380/8.31446)*(1/298.15 - 1/TleafK)) ## formulation used in plantecophys
    
    phij = 0.24 ##  plantecophys
    thetaj = 0.85 ##  plantecophys
    i <- PPFD
    J <- ((-(-Jm-phij*i) - sqrt((-Jm-phij*i)*(-Jm-phij*i) - 4*thetaj*Jm*phij*i))/(2*thetaj))
    
    B = 0.25*J - Rd + ro
    K = K25/((1.856e-11*exp(4209/TleafK + 0.04527*TleafK - 3.376e-5*TleafK*TleafK))/0.907784) #viscosity correction
    p = s*(psisoil - psic)
    d = s*deltaw/K
    if((J>0) & (p>0)) { ## find cubic solution ## coefficients for V- and J-limited conditions
      x = (Vm - Rd)*ca - (Vm*gams + Rd*kp)
      q3v = 1.6*d + p*rbc ## coefs as in ms but multiplied by s/K, so Em -> p
      q2v = -p*(ca + kp + rbc*(B + Vm - Rd)) - 1.6*(1 + d*(B + Vm - Rd))
      q1v = p*(B*(ca + kp) + x + rbc*B*(Vm - Rd)) + 1.6*(Vm - Rd)*(1 + d*B)
      q0v = -p*B*x
      XV = ifelse(Vm - Rd > B, B, Vm - Rd) ## upper limit of interval of valid A
      x = (0.25*J - Rd)*ca - (0.25*J + 2*Rd)*gams
      q3j = 1.6*d + p*rbc
      q2j = -p*(ca + 2*gams + rbc*(B + 0.25*J - Rd)) - 1.6*(1 + d*(B + 0.25*J - Rd))
      q1j = p*(B*(ca + 2*gams) + x + rbc*B*(0.25*J - Rd)) + 1.6*(0.25*J - Rd)*(1 + d*B)
      q0j = -p*B*x
      XJ = 0.25*J - Rd
      aa. <- c(q3v, q3j)# make vectors with V/J-limited values of inputs and outputs
      bb. <- c(q2v, q2j)
      cc. <- c(q1v, q1j)
      dd. <- c(q0v, q0j)
      XX. <- c(XV, XJ)
      ci. <- c(NA,NA)
      g. <- c(NA,NA)
      A. <- c(NA,NA)
      for(vj in 1:2) {
        aa = aa.[vj]
        bb = bb.[vj]
        cc = cc.[vj]
        dd = dd.[vj]
        XX = XX.[vj]
        if(d==0 & rbc==0) {
          A = (0.5/bb)*(-cc + sqrt(cc*cc - 4*dd*bb))
          g = p*(B - A)/(1 + d*(B - A))
          gtc = g/(1.6 + rbc*g)
          ci = ca - A/gtc
          ci.[vj] = ci; A.[vj] = A; g.[vj] <- g; root = 4
        } else if(ro==0 & vj==2) {
          aaj = 1.6*d + p*rbc
          bbj = -p*(ca + 2*gams + rbc*(0.25*J - Rd)) - 1.6*(1 + d*(0.25*J - Rd))
          ccj = p*((0.25*J - Rd)*ca - (0.25*J + 2*Rd)*gams)
          if(bbj*bbj - 4*ccj*aaj > 0) {
            A = (0.5/aaj)*(-bbj - sqrt(bbj*bbj - 4*ccj*aaj))
          } else {
            A = 0          
          }
          g = p*(B - A)/(1 + d*(B - A))
          gtc = g/(1.6 + rbc*g)
          ci = ca - A/gtc
          if((ci>gams) & (A>-Rd) & (A<XX)) {
            ci.[vj] = ci; A.[vj] = A; g.[vj] <- g; root = 5
          } else {
            ci.[vj] = NA; A.[vj] = NA; g.[vj] <- NA; root = NA
          }
        } else {
          pp = (3*aa*cc - bb*bb)/(3*aa*aa)
          qq = (2*bb*bb*bb - 9*aa*bb*cc + 27*aa*aa*dd)/(27*aa*aa*aa)
          discrim = -pp*pp*pp/27 - qq*qq/4
          bb3aa = bb/(3*aa)
          if(discrim >= 0) {## if discrim is positive, there are three real roots
            tx = 2*sqrt(-pp/3) ## factor multiplied by cos(ty)
            ty = (1/3)*acos((3*qq/(2*pp))*sqrt(-3/pp)) ## argument to cos(), without varible offset (k*2pi/3, k in 0,1,2)
            t2 = tx*cos(ty - 4*pi/3)  ## first calculate root #2
            A2 = t2 - bb3aa
            g2 = p*(B - A2)/(1 + d*(B - A2))
            gtc = g2/(1.6 + rbc*g2)
            ci2 = ifelse(gtc>0, ca - A2/gtc, ca)
            if((ci2 > gams) & (A2 > -Rd) & (A2 < XX)) {
              ci.[vj] = ci2; A.[vj] = A2; g.[vj] <- g2; root = 2
            } else {
              ## if that doesn't work, calculate root #1
              t1 = tx*cos(ty - 2*pi/3)
              A1 = t1 - bb3aa
              g1 = p*(B - A1)/(1 + d*(B - A1))
              gtc = g1/(1.6 + rbc*g1)
              ci1 = ca - A1/gtc
              if((ci1 > gams) & (A1 > -Rd) & (A1 < XX)) {
                ci.[vj] = ci1; A.[vj] = A1; g.[vj] <- g1; root = 1
              } else {
                ## if that doesn't work, calculate root #2
                t0 = tx*cos(ty) ## implicitly tx*cos(ty - 0*2*pi/3)
                A0 = t0 - bb3aa
                g0 = p*(B - A0)/(1 + d*(B - A0))
                gtc = g0/(1.6 + rbc*g0)
                ci0 = ca - A0/gtc
                
                ## if this root is realistic, keep it
                if((ci0 > gams) & (A0 > -Rd) & (A0 < XX)) {
                  ci.[vj] = ci0; A.[vj] = A0; g.[vj] <- g0; root = 0
                }  else {
                  ## if that doesn't work, set everything to NA as a flag
                  ci.[vj] = NA; A.[vj] = NA; g.[vj] = NA; root = NA
                } ## end if root 2
              } ## end if root 1
            } ## end if root 3
          }  else { ## if discriminant is negative, there is only one real root
            sq = sqrt(-discrim)
            u1 = -qq/2 + sq
            u2 = -qq/2 - sq
            t = sign(u1)* (abs(u1)^(1/3)) + sign(u2)* (abs(u2)^(1/3))
            
            bb3aa = bb/(3*aa)
            A = t - bb3aa
            g = p*(B - A)/(1 + d*(B - A))
            gtc = g/(1.6 + rbc*g)
            ci = ca - A/gtc
            root = -1
          }  ## end if discrim
        } ## end if d==0
      } ## end for vj
      
      ## smooth VJ transition
      if(!is.na(A.[1]) & !is.na(A.[2])) {
        A = (0.5/thetaa)*((A.[1] + A.[2]) - sqrt((A.[1] + A.[2])^2 - 4*thetaa*A.[1]*A.[2]))
        lim = ifelse(A.[1] < A.[2], 1, 2)
      } else {
        A = min(A.[1], A.[2], na.rm=T)  
        lim = ifelse(is.na(A.[1]),2,1)
      }
      
    } else {
      ## if J=0 or P<=0, set A to -Rd
      A = -Rd
      lim = 2
      root=-2
    } ## end if J>0, P>0
    
    g = p*(B - A)/(1 + d*(B - A)) ## this is the new model
    g = ifelse(g>=0, g, 0)
    gtc = g/(1.6 + rbc*g)
    ci = ca - A/gtc
    rbw = rbc/1.37
    gtw = g/(1 + rbw*g)
    E = gtw*deltaw
    psileaf = psisoil - gtw*deltaw/K
    return(c(g, A, E, psileaf, ci, root))
  } else {
    ## if h < 0, return NAs
    return(c(NA,NA,NA,NA,NA,NA))
  } ## end if h
} ## end gasxNEW()



#######################
#######################
###
### actual calculations
###
#######################
#######################



## function to compare input and output psileaf, to reconcile VC in new model
test_hydraulics <- function(psileaf, PPFD, ca, TleafC, deltaw, psisoil, rbc,
                            Vm25, Jm25, thetaa, K25max, psic, s, ro, Tcrit, T50) {
  K25 = K25max*VC(psisoil, psileaf, b_VC, c_VC) ## get K implied by input psileaf
  gx <- gasxNEW_mod(PPFD, ca, TleafC, deltaw, psisoil, rbc, Vm25, Jm25, thetaa, K25, psic, s, ro, Tcrit, T50) ## calculate implied gasx
  psileafGX = gx[4] ## get implied psileaf
  return((psileafGX - psileaf)^2)
}


### fn to apply new model iteratively, solving for psileaf that satisfies VC
gswmod_new_withVC <- function(PPFD, ca, TleafC, deltaw, psisoil, rbc, Vm25, Jm25, thetaa,
                              K25max, psic, s, ro, Tcrit, T50) {
  o <- optimize(test_hydraulics, interval=c(-P88, psisoil),
                PPFD=PPFD, ca=ca, TleafC=TleafC, deltaw=deltaw, psisoil=psisoil, rbc=rbc, Vm25=Vm25,
                Jm25=Jm25, thetaa=thetaa, K25max=K25max, psic=psic, s=s, ro=ro, Tcrit=Tcrit, T50=T50)

  K25 = K25max*VC(psisoil, o$minimum, b_VC, c_VC)
  gx <- gasxNEW_mod(PPFD, ca, TleafC, deltaw, psisoil, rbc, Vm25, Jm25, thetaa,
                    K25, psic, s, ro, Tcrit, T50) ## calculate implied gasx
  return(gx)
}


## parameters from Sicangco et al Table 2
Tcrit_C = 43.7; Tcrit_HW = 46.5
T50_C = 48.6; T50_HW = 50.4
g1 = 2.9; g0 = 1e-5
Vm25 = 100.52; Jm25 = 165.53 
P50 = 4.07; P88 = 5.5

Patm = 101325
vwind = 5 ## 5 m/s, not in papers but in their code (analysis_functions.R)
dleaf = 0.025 ## m, in Sicangco

ca=425
thetaa=0.99
ro=0


## get VC parameters
VCpars <- fitVC(P50, P88)
b_VC = VCpars[1]
c_VC = VCpars[2]


## load WTC data
x=read.delim("WTC4_data.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)


## calculate needed internals
x$gbw <- get_gbw(x$Tair, x$Tcan, vwind, dleaf)
x$rbc <- 1.37/x$gbw
x$psisoil = -x$Ps
x$wsata <- 6.1078*exp(17.27*x$Tair/(237.3 + x$Tair))
x$dwa <- 10*x$VPD
x$wa <- x$wsata - x$dwa
x$deltaw <- 6.1078*exp(17.27*x$Tcan/(237.3 + x$Tcan)) - x$wa
x$DateTime_hr <- as.POSIXct(x$DateTime_hr,format="%Y-%m-%d %T",tz="GMT")
x <- subset(x, PPFD >= 500)


## load outputs from sicangco that had USO model predictions
load("control_runs.Rdata")
load("heatwave_runs.Rdata")
colnames(pred.c)
colnames(pred.hw)
pred.c <- pred.c[,-12]
uc <- subset(pred.c, Model=="Medlyn")
uh <- subset(pred.hw, Model=="Medlyn")
uc$HWtrt <- "C"
uh$HWtrt <- "HW"
u <- rbind(uc, uh)
u <- left_join(u, x[,c("DateTime_hr", "E", "A", "gs", "Ps", "Tair")], by="Tair")
colnames(u)[c(5, 8, 9, 14:17)] <- c("Euso", "guso", "Auso", "E", "A", "gs", "Ps")

xc <- subset(x, HWtrt=="C")
xh <- subset(x, HWtrt!="C")


### get gam predictions for control, then heatwave
nlink=6
Egamc <- gam(E ~ s(Tcan, k=nlink), data=xc);# plot.gam(Egamc)
Egamc.pred <- predict(Egamc, xc, type="link", se.fit=T)
xc$Egam <- Egamc.pred$fit; xc$Egam.se <- Egamc.pred$se.fit
Agamc <- gam(A ~ s(Tcan, k=nlink), data=xc);# plot.gam(Agamc)
Agamc.pred <- predict(Agamc, xc, type="link", se.fit=T)
xc$Agam <- Agamc.pred$fit; xc$Agam.se <- Agamc.pred$se.fit

Egamh <- gam(E ~ s(Tcan, k=nlink), data=xh); #plot.gam(Egamh)
Egamh.pred <- predict(Egamh, xh, type="link", se.fit=T)
xh$Egam <- Egamh.pred$fit; xh$Egam.se <- Egamh.pred$se.fit
Agamh <- gam(A ~ s(Tcan, k=nlink), data=xh); #plot.gam(Agamh)
Agamh.pred <- predict(Agamh, xh, type="link", se.fit=T)
xh$Agam <- Agamh.pred$fit; xh$Agam.se <- Agamh.pred$se.fit


## calculate predictions from new model; psic and s have been manually adjusted to give good fit to WTC data
arglistc=list(psic=-P50/2.3, s=0.06, K25max=4., ca=ca, Vm25=Vm25, Jm25=Jm25, thetaa=thetaa, ro=ro, Tcrit=Tcrit_C, T50=T50_C)
arglisth=list(psic=-P50/2.3, s=0.06, K25max=4., ca=ca, Vm25=Vm25, Jm25=Jm25, thetaa=thetaa, ro=ro, Tcrit=Tcrit_HW, T50=T50_HW)

yc <- data.frame(t(mapply(gswmod_new_withVC, xc$PPFD, xc$Tcan, xc$deltaw, xc$psisoil, xc$rbc, MoreArgs=arglistc)))
xc$gmod <- yc[,1]/2.9;xc$Amod <- yc[,2]/2.9;xc$Emod <- yc[,3]/2.9; xc$psileaf <- yc[,4]; xc$ci <- yc[,5]

yh <- data.frame(t(mapply(gswmod_new_withVC, xh$PPFD, xh$Tcan, xh$deltaw, xh$psisoil, xh$rbc, MoreArgs=arglisth)))
xh$gmod <- yh[,1]/2.9;xh$Amod <- yh[,2]/2.9;xh$Emod <- yh[,3]/2.9; xh$psileaf <- yh[,4]; xh$ci <- yh[,5]

xx <- rbind(xc, xh)
xx$Kplant <- NA
for(i in 1:nrow(xx)) {xx$Kplant[i] <- arglistc[3][[1]]*VC(xx$psisoil[i], xx$psileaf[i], b_VC, c_VC)}


## prepare for plotting
Alabel = "assimilation~rate~(mu~mol~m^{-2}~s^{-1})"
Elabel = "transpiration~rate~(mmol~m^{-2}~s^{-1})"

xxe <- xx[, c("Tcan", "Emod", "HWtrt")]; colnames(xxe) <- c("Tleaf", "flux", "HWtrt"); xxe$var <- Elabel
xxa <- xx[, c("Tcan", "Amod", "HWtrt")]; colnames(xxa) <- c("Tleaf", "flux", "HWtrt"); xxa$var <- Alabel
xxea <- rbind(xxe, xxa); xxea$source <- "new model"

ue <- u[,c("Tleaf", "Euso", "HWtrt")]; colnames(ue) <- c("Tleaf", "flux", "HWtrt"); ue$var <- Elabel
ua <- u[,c("Tleaf", "Auso", "HWtrt")]; colnames(ua) <- c("Tleaf", "flux", "HWtrt"); ua$var <- Alabel
uea <- rbind(ue, ua); uea$source <- "USO model"

de <- xx[, c("Tcan", "E", "HWtrt")]; colnames(de) <- c("Tleaf", "flux", "HWtrt"); de$var <- Elabel
da <- xx[, c("Tcan", "A", "HWtrt")]; colnames(da) <- c("Tleaf", "flux", "HWtrt"); da$var <- Alabel
dea <- rbind(de, da); dea$source <- "measured"

xall <- rbind(dea, uea, xxea)
xall$HWtrt <- ifelse(xall$HWtrt=="C", "control", "heatwave")

ge <- xx[, c("Tcan", "Egam", "Egam.se", "HWtrt")]; colnames(ge) <- c("Tleaf", "flux", "se", "HWtrt"); ge$var <- Elabel
ga <- xx[, c("Tcan", "Agam", "Agam.se", "HWtrt")]; colnames(ga) <- c("Tleaf", "flux", "se", "HWtrt"); ga$var <- Alabel
gea <- rbind(ge, ga); gea$data <- "GAM"
gea$HWtrt <- ifelse(gea$HWtrt=="C", "control", "heatwave")


xall$source <- factor(xall$source, levels=c("measured", "USO model", "new model"))
xall$xpos <- ifelse(xall$HWtrt=="control", 14, 17);
xall$ypos <- ifelse(xall$var==Alabel, 12.5, 3.75);
xall$label <- ifelse(xall$HWtrt=="control" & xall$var==Alabel, "(a)",
                     ifelse(xall$HWtrt=="heatwave" & xall$var==Alabel, "(b)",
                            ifelse(xall$HWtrt=="control" & xall$var==Elabel, "(c)", "(d)")))

ggplot() +
  geom_line(data=gea, aes(x=Tleaf, y=flux), color="black", linewidth=1) +
  geom_ribbon(data=gea, aes(x=Tleaf, ymin=flux-2*se, ymax=flux+2*se), alpha=0.2) +
  geom_point(data=xall, aes(x=Tleaf, y=flux, color=source), size=1, shape=1, alpha=0.5)+
  geom_text(data=xall, aes(x=xpos, y=ypos, label=label), family="Calibri")+ ylab("")+
  scale_color_manual(values=c("#777777", "#FF0000", "#0000FF"))+
  facet_grid(rows=vars(var), cols=vars(HWtrt), 
             labeller = label_parsed,
             scales="free", switch="y")+  theme_bw() +
   theme(strip.background=element_blank(), strip.placement="outside", strip.position="left")+
  xlab(expression(paste("leaf temperature /[degrees C]"))) +
  guides(color=guide_legend(override.aes=list(size=2))) 

ggsave("Figure 3.png", dpi=600, device="png")




##############################
##############################
###
###  now make SI figure showing water potentials for same
###
##############################
##############################


library(lubridate)

{
  ## load data
  wp=read.delim("WTC_TEMP-PARRA_CM_WATERPOTENTIAL-HEATWAVE_20161019-20161107_L0.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
  wp <- subset(wp, Timing=="midday" & tissue=="leaf" &
               as.POSIXct(Date) >= as.POSIXct("2016-10-29") &
               as.POSIXct(Date) < as.POSIXct("2016-11-12"))
  wpd <- wp %>% group_by(Date, HW_treatment) %>%
  summarise(LWP. = mean(LWP), LWP.se = SE(LWP))
  
  colnames(wpd)[2] <- "HWtrt"
  wpd$HWtrt <- ifelse(wpd$HWtrt=="control", "C", "HW")
  wpd$measured <- ""
  
  xx$date <- date(xx$DateTime_hr)
  xxmd <- subset(xx, hour(DateTime_hr) > 12 & hour(DateTime_hr) < 15)
  xxmdd <- xxmd %>% group_by(date, HWtrt) %>%
  summarise(LWP. = mean(psileaf), LWP.se = SE(psileaf))
  
  u$psisoil <- -u$Ps
  u$psileafhat <- -u$P
  u$psileaf <- u$psisoil - (u$psisoil - u$psileafhat)*1.5/4
  u$date <- date(u$DateTime_hr)
  umd <- subset(u, hour(DateTime_hr) > 12 & hour(DateTime_hr) < 15)
  umdd <- umd %>% group_by(date, HWtrt) %>%
  summarise(LWP. = mean(psileaf), LWP.se = SE(psileaf))
  
  xxmdd$model <- "NEW"
  umdd$model <- "USO"
  xupsi <- rbind(xxmdd, umdd)
  
  xupsi$HWtrt <- ifelse(xupsi$HWtrt=="C", "control", "heatwave")
  wpd$HWtrt <- ifelse(wpd$HWtrt=="C", "control", "heatwave")
  
  xupsi$xpos <- as.POSIXct("2016-10-30")
  xupsi$ypos <- -0.05
  xupsi$label <- ifelse(xupsi$HWtrt=="control", "(a)", "(b)")
  
}; ggplot() +
  geom_line(data=xupsi, aes(x=date, y=LWP., color=model), size=1) +
  geom_text(data=xupsi, aes(x=xpos, y=ypos, label=label), family="Calibri")+ ylab("")+
  scale_color_manual(values=c("#0000FF", "#FF0000")) +
  geom_point(data=wpd, aes(x=date(Date), y=LWP., shape=measured), color="black", size=3) +
  geom_errorbar(data=wpd, aes(x=date(Date), ymin=LWP.-LWP.se, ymax=LWP.+LWP.se), color="black", width=0.3) +
  theme_bw() + facet_grid(cols=vars(HWtrt)) + ylim(-2, 0) + theme(strip.background=element_blank())+
  ylab("midday leaf water potential /MPa") +
  guides(shape=guide_legend())

ggsave("Figure S1.png", dpi=600, device="png")
