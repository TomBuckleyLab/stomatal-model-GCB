##load necessary libraries
library(doBy)
library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggpubr)

options("scipen"=100, "digits"=8)

#########################
#########################
###
### code to define basic default functions needed to generate figures in
###  Buckley et al (2026), Global Change Biology, "A simple stomatal model that unifies the metabolic and hydraulic
###  control of carbon and water flux"
###
#########################
#########################



## fn to get relative humidity given deltaw and T
getRH <- function(deltaw, temperature) {
   return(1 - deltaw/(6.1078*exp(17.27*temperature/(237.3 + temperature)))) 
}

## fn to get deltaw given relative humidity and T
getDW <- function(h, temperature) {
   return((1 - h)*6.1078*exp(17.27*temperature/(237.3 + temperature))) 
}

############
### gasxNEW()
# calculates gas exchange from the new model
# Vm, Jm, Rd, PPFD ~ umol m-2 s-1
# ca ~ umol/mol
# deltaw ~ mmol/mol
# K25 = plant hydraulic conductance at 25C, mmol m-2 s-1 MPa-1
# psisoil, pi ~ MPa; note pi > 0
# m = fitted parameter, ~ mol air umol-1CO2 MPa-1
# thetaa = photosynthetic colimitation (smoothing) curvature parameter, dimensionless, close to 1 but <=1
## nb "ro" is labelled as "bo" in the actual paper
gasxNEW <- function(PPFD, ca, TleafC, deltaw, psisoil, Vm25, Jm25, Rd25, thetaa, K25, psic, s, ro, thetaj.input=0, phij.input=0) {
      h = 1 - deltaw/(6.1078*exp(17.27*TleafC/(237.3 + TleafC)))
   
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
     B = 0.25*J - Rd + ro
     K = K25/((1.856e-11*exp(4209/TleafK + 0.04527*TleafK - 3.38e-5*TleafK*TleafK))/0.907784) #viscosity correction
     p = s*(psisoil - psic)
     d = s*deltaw/K
     if((J>0) & (p>0)) { ## find cubic solution ## coefficients for V- and J-limited conditions
       x = (Vm - Rd)*ca - (Vm*gams + Rd*kp)
       q3v = 1.6*d 
       q2v = -p*(ca + kp) - 1.6*(1 + d*(B + Vm - Rd))
       q1v = p*(B*(ca + kp) + x) + 1.6*(Vm - Rd)*(1 + d*B)
       q0v = -p*B*x
       XV = ifelse(Vm - Rd > B, B, Vm - Rd) ## upper limit of interval of valid A
       x = (0.25*J - Rd)*ca - (0.25*J + 2*Rd)*gams
       q3j = 1.6*d 
       q2j = -p*(ca + 2*gams) - 1.6*(1 + d*(B + 0.25*J - Rd))
       q1j = p*(B*(ca + 2*gams) + x) + 1.6*(0.25*J - Rd)*(1 + d*B)
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
         if(d==0) {
            A = (0.5/bb)*(-cc + sqrt(cc*cc - 4*dd*bb))
            g = p*(B - A)/(1 + d*(B - A))
            ci = ca - 1.6*A/g
            ci.[vj] = ci; A.[vj] = A; g.[vj] <- g; root = 4
         } else if(ro==0 & vj==2) {
            aaj = 1.6*d 
            bbj = -p*(ca + 2*gams) - 1.6*(1 + d*(0.25*J - Rd))
            ccj = p*((0.25*J - Rd)*ca - (0.25*J + 2*Rd)*gams)
            A = (0.5/aaj)*(-bbj - sqrt(bbj*bbj - 4*ccj*aaj))
            g = p*(B - A)/(1 + d*(B - A))
            ci = ca - 1.6*A/g
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
               ci2 = ca - 1.6*A2/g2
               if((ci2 > gams) & (A2 > -Rd) & (A2 < XX)) {
                  ci.[vj] = ci2; A.[vj] = A2; g.[vj] <- g2; root = 2
               } else {
                  ## if that doesn't work, calculate root #1
                  t1 = tx*cos(ty - 2*pi/3)
                  A1 = t1 - bb3aa
                  g1 = p*(B - A1)/(1 + d*(B - A1))
                  ci1 = ca - 1.6*A1/g1
                  if((ci1 > gams) & (A1 > -Rd) & (A1 < XX)) {
                     ci.[vj] = ci1; A.[vj] = A1; g.[vj] <- g1; root = 1
                  } else {
                     ## if that doesn't work, calculate root #2
                     t0 = tx*cos(ty) ## implicitly tx*cos(ty - 0*2*pi/3)
                     A0 = t0 - bb3aa
                     g0 = p*(B - A0)/(1 + d*(B - A0))
                     ci0 = ca - 1.6*A0/g0
                     
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
               ci = ca - 1.6*A/g
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
     ci = ca - 1.6*A/g 
     return(c(g, A, ci, B-A, lim, root))
   } else {
     ## if h < 0, return NAs
     return(c(NA,NA,NA,NA,NA,NA))
  } ## end if h
} ## end gasxNEW()


### as gasxNEW(), but for BB model (gsw = m*h*A/ca + g0)
gasxBB <- function(PPFD, ca, TleafC, h, Vm25, Jm25, Rd25, thetaa, m, g0) {
   if(h>0) {
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
      phij <- phiPSII*0.9*0.5 #0.9 is assumed absorptance; 0.5 is photon fraction
      thetaj <- 0.76 + 0.018*TleafC - 3.7e-4*TleafC*TleafC #unitless curvature parameter
      i <- PPFD
      J <- ((-(-Jm-phij*i) - sqrt((-Jm-phij*i)*(-Jm-phij*i) - 4*thetaj*Jm*phij*i))/(2*thetaj))
      mh = m*h  ## product of slope m and RH h
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
      g = mh*A/ca + g0 ## this is the NEW model
      g <- ifelse(g > 0, g, 0)
      ci = ca - 1.6*A/g
   
      return(c(g, A, ci, lim))
   } else {
      return(c(NA,NA,NA,NA))
   }
} ## end gasxBB()


### as gasxBB(), but for USO model (gsw = 1.6*(1 + g1/sqrt(D))*A/ca + g0)
##   needs D in kPa
gasxUSO <- function(PPFD, ca, TleafC, D, Vm25, Jm25, Rd25, thetaa, g1, g0) {
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
      phij <- phiPSII*0.9*0.5 #0.9 is assumed absorptance; 0.5 is photon fraction
      thetaj <- 0.76 + 0.018*TleafC - 3.7e-4*TleafC*TleafC #unitless curvature parameter
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



##### function to generate output from models
### var1 = PPFD (1), or CA (2), etc., see a few lines earlier
### var2 = similar: 2nd var to repeat sims by
### defaults to nx1=100 and nx2=5 values of vars 1 and 2, but can override 
###
### can also specify allvars=TRUE instead of specifying var1 and var2; this creates permutations of combos for all vars,
###   with n values per var (default 10). If nopsi=TRUE, psisoil is omitted from vars.
###
### prints figure to screen, return df with results
g.vs. <- function(var1=0, var2=0, allvars=FALSE, nx1=101, nx2=3, n=10, yvar="gsw", models=c("NEW", "BB", "USO"), nopsi=FALSE,
                  mBB = 10, g1 = 5, s = 0.02, g0bb = 0.02, g0uso=0.02, ro = 0.3, thetaa=0.99, noplot=FALSE, simple=FALSE) {

   ## create df to contain inputs and outputs
   if(allvars==TRUE) {
      ppfdrange <- seq(from=ppfdmin, to=ppfdmax, by=(ppfdmax - ppfdmin)/(n-1))
      carange <- seq(from=camin, to=camax, by=(camax - camin)/(n-1))
      TleafCrange <- seq(from=TleafCmin, to=TleafCmax, by=(TleafCmax - TleafCmin)/(n-1))
      deltawrange <- seq(from=deltawmin, to=deltawmax, by=(deltawmax - deltawmin)/(n-1))
      psisoilrange <- seq(from=psisoilmin, to=psisoilmax, by=(psisoilmax - psisoilmin)/(n-1))
      hrange <- seq(from=hmin, to=hmax, by=(hmax - hmin)/(n-1))
   } else {
      ppfdrange <- ppfd_ 
      carange <- ca_
      TleafCrange <- TleafC_
      deltawrange <- deltaw_
      psisoilrange <- psisoil_
      hrange <- h_
      
      if(var1==PPFD) {ppfdrange <- seq(from=ppfdmin, to=ppfdmax, by=(ppfdmax - ppfdmin)/(nx1-1))}
      if(var1==CA) {carange <- seq(from=camin, to=camax, by=(camax - camin)/(nx1-1))}
      if(var1==TLEAF) {TleafCrange <- seq(from=TleafCmin, to=TleafCmax, by=(TleafCmax - TleafCmin)/(nx1-1))}
      if(var1==DW) {deltawrange <- seq(from=deltawmin, to=deltawmax, by=(deltawmax - deltawmin)/(nx1-1))}
      if(var1==PSISOIL) {psisoilrange <- seq(from=psisoilmin, to=psisoilmax, by=(psisoilmax - psisoilmin)/(nx1-1))}
      if(var1==H) {hrange <- seq(from=hmin, to=hmax, by=(hmax - hmin)/(nx1-1))}
      
      if(var2==PPFD) {ppfdrange <- seq(from=ppfdmin, to=ppfdmax, by=(ppfdmax - ppfdmin)/(nx2-1))}
      if(var2==CA) {carange <- seq(from=camin, to=camax, by=(camax - camin)/(nx2-1))}
      if(var2==TLEAF) {TleafCrange <- seq(from=TleafCmin, to=TleafCmax, by=(TleafCmax - TleafCmin)/(nx2-1))}
      if(var2==DW) {deltawrange <- seq(from=deltawmin, to=deltawmax, by=(deltawmax - deltawmin)/(nx2-1))}
      if(var2==PSISOIL) {psisoilrange <- seq(from=psisoilmin, to=psisoilmax, by=(psisoilmax - psisoilmin)/(nx2-1))}
      if(var2==H) {hrange <- seq(from=hmin, to=hmax, by=(hmax - hmin)/(nx2-1))}
   } 

   if(nopsi==TRUE) psisoilrange=psisoil_

   ## calculate RH for BB model unless RH is specified as one of the x vars,
   ##  otherwise calculate deltaw from RH
   if((allvars==TRUE) | (var1!=H & var2!=H)) {
      x <- expand.grid(ppfdrange, carange, TleafCrange, deltawrange, psisoilrange)
      colnames(x) <- c("ppfd", "ca", "TleafC", "deltaw", "psisoil")
      x$h = mapply(getRH, x$deltaw, x$TleafC) 
   } else {
      x <- expand.grid(ppfdrange, carange, TleafCrange, psisoilrange, hrange)
      colnames(x) <- c("ppfd", "ca", "TleafC", "psisoil", "h")
      x$deltaw = mapply(getDW, x$h, x$TleafC)
      x <- x[,c(1,2,3,6,4,5)] ## reorder columns for consistency
   }
   ## calculate VPD for USO model
   x$D <- x$deltaw/10 ## assumes atm P of 100 kPa

   xNEW <- NULL
   xBB <- NULL
   xUSO <- NULL
   xCF <- NULL
   
   if("NEW" %in% models) {
      ## simulate NEW model
      arglist = list(Vm25=Vm25, Jm25=Jm25, Rd25=Rd25, K25=K25, psic=psic, s=s, ro=ro, thetaa=thetaa)
      if(!simple) {
        y <- mapply(gasxNEW, x[,1], x[,2], x[,3], x[,4], x[,5], MoreArgs=arglist)
      } else {
        y <- mapply(gasxNEW.simple, x[,1], x[,2], x[,3], x[,4], x[,5], MoreArgs=arglist)
      }
      y <- data.frame(t(y))
      
      xNEW <- x
      xNEW$model <- "NEW"
      xNEW$modnum <- 1
      xNEW$gsw <- y[,1]
      xNEW$A <- y[,2]
      xNEW$ci <- y[,3]
      xNEW$lim <- y[,5]
      xNEW$r <- y[,4]
   }

   if("BB" %in% models) {
      ## simulate BB model
      arglist = list(Vm25=Vm25, Jm25=Jm25, Rd25=Rd25, m=mBB, g0=g0bb, thetaa=thetaa)
      y <- mapply(gasxBB, x[,1], x[,2], x[,3], x[,6], MoreArgs=arglist) ## x[,4] (deltaw) replaced with x[,6] (h); x[,5] (psisoil) omitted
      y <- data.frame(t(y))
      
      xBB <- x
      xBB$model <- "BB"
      xBB$modnum <- 2
      xBB$gsw <- y[,1]
      xBB$A <- y[,2]
      xBB$ci <- y[,3]
      xBB$lim <- y[,4]
      xBB$r <- 0
   }
   
   if("USO" %in% models) {
      ## simulate USO model
      arglist = list(Vm25=Vm25, Jm25=Jm25, Rd25=Rd25, g1=g1, g0=g0uso, thetaa=thetaa)
      y <- mapply(gasxUSO, x[,1], x[,2], x[,3], x[,7], MoreArgs=arglist) ## x[,4] (deltaw) replaced with x[,7] (VPD/kPa); x[,5] (psisoil) omitted
      y <- data.frame(t(y))
      
      xUSO <- x
      xUSO$model <- "USO"
      xUSO$modnum <- 3
      xUSO$gsw <- y[,1]
      xUSO$A <- y[,2]
      xUSO$ci <- y[,3]
      xUSO$lim <- y[,4]
      xUSO$r <- 0
   }
   
   if(allvars==TRUE) {
      ## combine all results in wide format
      y <- x
      if("NEW" %in% models) {for(i in 1:5) colnames(xNEW)[i+9] <- paste0(colnames(xNEW)[i+9], "_NEW"); y <- cbind(y, xNEW[,c(10:14)])}
      if("BB" %in% models)  {for(i in 1:5) colnames(xBB)[i+9] <-  paste0(colnames(xBB)[i+9],  "_BB");  y <- cbind(y, xBB[, c(10:14)])}
      if("USO" %in% models) {for(i in 1:5) colnames(xUSO)[i+9] <- paste0(colnames(xUSO)[i+9], "_USO"); y <- cbind(y, xUSO[,c(10:14)])}
   } else {
      ## combine all results in long format
      y <- rbind(xNEW, xBB)
      y <- rbind(y, xUSO)

      ## force new model to appear first in plots
      y$model <- factor(y$model, c("NEW", "USO", "BB"))
      # y$size <- as.integer(ifelse(y$model=="NEW", 1, 0.5))
      if(!noplot) {
         p <- ggplot() +  geom_line(data=y, 
                                    aes(x=.data[[varnames[var1]]], 
                                        color=as.factor(.data[[varnames[var2]]]), 
                                        y=.data[[yvar]],
                                        linetype=model, size=model))+
            scale_color_viridis(discrete=T, option="D", begin=0, end=0.7)+
            scale_linetype_manual(values=c("solid", "longdash", "dotdash"))+
            scale_size_manual(values=c(1.2, 0.7, 0.7))+
            ylim(0,NA) + xlim(varlow[var1], varhi[var1])+
            ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")")))+
            xlab(varnamesb[var1])+
            guides(color=guide_legend(title=varnamesc[var2]))
         print(p)
      }
   }
   
   return(y)   
} ## end g.vs.()


