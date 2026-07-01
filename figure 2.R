setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##load libraries
library(dplyr)
library(ggplot2)

options("scipen"=100, "digits"=4)

################
################
###
### code to generate Figure 2 (comparison of models for seasonal drought in woody crops)
###
################
################


##### 
#### override original code for models because we have species-specific photosynthetic parameter T responses
####
gasxNEW <- function(PPFD, ca, TleafC, deltaw, psisoil, Vm25, Jm25, gm25, phij, thetaj, thetaa, K25, psic, s, ro, spp) {
   h = 1 - deltaw/(6.1078*exp(17.27*TleafC/(237.3 + TleafC)))
   
   if(h>=0) {
      TlK = TleafC + 273.15
      TleafK=TlK
      itlk = 1 / TlK
      irt = 0.120273 * itlk * 0.0033557 #0.12... is 1/R, .0033557 is 1/298
      iit = 298.15 * itlk * itlk
      
      KC25 = 272.38
      KO25 = 165.82
      gams25 = KC25 * 210 * 0.21 / (2 * KO25) #210 = 21% * 10 to get mf
      dKc = exp(80990 * (TlK - 298.15) * irt); Kc = KC25 * dKc
      dKo = exp(23720 * (TlK - 298.15) * irt); Ko = KO25 * dKo 
      
      if(spp=="grapevine") {
         TlK = TleafC + 273.15
         TleafK=TlK
         KC25 = 272.38
         KO25 = 165.82
         gams25 = KC25 * 210 * 0.21 / (2 * KO25) #210 = 21% * 10 to get mf
         
         itlk = 1 / TlK
         irt = 0.120273 * itlk * 0.0033557 #0.12... is 1/R, .0033557 is 1/298
         iit = 298.15 * itlk * itlk
         
         #''Vm T response is for grapevine
         dvm = exp(61115 * (TlK - 298.15) * irt)
         Vm = Vm25 * dvm
         drd = exp(44790 * (TlK - 298.15) * irt) #'activation here is from Antonio
         Rd = 0.0089 * Vm25 * drd                 # ' Rd25 vs Vm25 is from dePury
         dKc = exp(80990 * (TlK - 298.15) * irt)
         Kc = KC25 * dKc
         dKo = exp(23720 * (TlK - 298.15) * irt)
         Ko = KO25 * dKo
         kp = Kc*(1 - 210/Ko)
         dgams = exp(24460 * (TlK - 298.15) * irt)
         gams = gams25 * dgams
         
         #''' gm temperature response for grapevine
         dgm = exp(-0.5 * (log(TleafC / 36.7476) / 0.83902) ^ 2)
         gm = gm25 * dgm
         rm = 1/gm
         
         #''' Jm temperature response for grapevine
         ea. = 55791.8767
         s. = 367.1633804
         ed. = 114804.1868
         rgas = 8.3143
         djm = exp(ea. * (TlK - 298.15) / (rgas * 298.15 * TlK))
         djm = djm / (1 + exp((s. * TlK - ed.) / (rgas * TlK)))
         djm = djm * (1 + exp((s. * 298.15 - ed.) / (rgas * 298.15)))
         Jm = Jm25 * djm
      } else if(spp=="almond") {
         TlK = TleafC + 273.15
         TleafK=TlK
         KC25 = 272.38
         KO25 = 165.82
         gams25 = KC25 * 210 * 0.21 / (2 * KO25) # 210 = 21% * 10 to get mf
         
         itlk = 1 / TlK
         irt = 0.120273 * itlk * 0.0033557 # 0.12... is 1/R, .0033557 is 1/298
         iit = 298.15 * itlk * itlk
         
         rgas = 0.0083143
         cvm = 31.58
         eavm = 78.1
         svm = 0.4964
         edvm = 154.21
         
         cjm = 13.775
         eajm = 34.11
         sjm = 0.5824
         edjm = 184.21
         
         cgm = 19.966
         eagm = 49.48
         sgm = 1.3599
         edgm = 426.4
         
         rd25 = 0.75
         q10rd = 2
         
         dvm = exp(cvm - eavm / (rgas * TlK)) / (1 + exp((svm * TlK - edvm) / (rgas * TlK)))
         Vm = Vm25 * dvm
         djm = exp(cjm - eajm / (rgas * TlK)) / (1 + exp((sjm * TlK - edjm) / (rgas * TlK)))
         Jm = Jm25 * djm
         dgm = exp(cgm - eagm / (rgas * TlK)) / (1 + exp((sgm * TlK - edgm) / (rgas * TlK)))
         gm = gm25 * dgm
         drd = q10rd ^ ((TleafC - 25) / 10)
         Rd = rd25 * drd
         dKc = exp(80990 * (TlK - 298.15) * irt)
         Kc = KC25 * dKc
         dKo = exp(23720 * (TlK - 298.15) * irt)
         Ko = KO25 * dKo
         dgams = exp(24460 * (TlK - 298.15) * irt)
         gams = gams25 * dgams
         
         kp = Kc * (1 + 210 / Ko)
         
         rm = 1/gm
      } else if(spp=="olive") {
         TlK = TleafC + 273.15
         TleafK=TlK
         KC25 = 272.38
         KO25 = 165.82
         gams25 = KC25 * 210 * 0.21 / (2 * KO25) # 210 = 21% * 10 to get mf
         
         itlk = 1 / TlK
         irt = 0.120273 * itlk * 0.0033557 # 0.12... is 1/R, .0033557 is 1/298
         iit = 298.15 * itlk * itlk
         
         dvm = exp(46525 * (TlK - 298.15) * irt)
         Vm = Vm25 * dvm
         drd = exp(44790 * (TlK - 298.15) * irt) #'activation here is from Antonio
         Rd = 0.0089 * Vm25 * drd        #          ' Rd25 vs Vm25 is from dePury
         dKc = exp(80990 * (TlK - 298.15) * irt)
         Kc = KC25 * dKc
         dKo = exp(23720 * (TlK - 298.15) * irt)
         Ko = KO25 * dKo
         dgams = exp(24460 * (TlK - 298.15) * irt)
         gams = gams25 * dgams
         
         #''' gm temperature response; parameters fitted from Antonio's lab
         ea = 79999.99606
         s. = 381.3603738
         ed = 113148.9517
         rgas = 8.3143
         dgm = exp(ea * (TlK - 298.15) / (rgas * 298.15 * TlK))
         dgm = dgm / (1 + exp((s. * TlK - ed) / (rgas * TlK)))
         dgm = dgm * (1 + exp((s. * 298.15 - ed) / (rgas * 298.15)))
         gm = gm25 * dgm
         
         #''' Jm temperature response; parameters fitted from Antonio's lab
         ea = 48410.13
         s. = 357.96
         ed = 113149.02
         rgas = 8.3143
         djm = exp(ea * (TlK - 298.15) / (rgas * 298.15 * TlK))
         djm = djm / (1 + exp((s. * TlK - ed) / (rgas * TlK)))
         djm = djm * (1 + exp((s. * 298.15 - ed) / (rgas * 298.15)))
         Jm = Jm25 * djm
         
         kp = Kc * (1 + 210 / Ko)
         
         rm = 1/gm
      }
      
      i <- PPFD
      J <- ((-(-Jm-phij*i) - sqrt((-Jm-phij*i)*(-Jm-phij*i) - 4*thetaj*Jm*phij*i))/(2*thetaj))
      B = 0.25*J - Rd + ro
      K = K25/((1.856e-11*exp(4209/TleafK + 0.04527*TleafK - 3.38e-5*TleafK*TleafK))/0.907784) #viscosity correction
      p = s*(psisoil - psic)
      d = s*deltaw/K + p*rm ## NOTE: this is the value of d giving gtc
      if((J>0) & (p>0)) { ## find cubic solution ## coefficients for V- and J-limited conditions
         x. = (Vm - Rd)*ca - (Vm*gams + Rd*kp)
         q3v = 1.6*d 
         q2v = -p*(ca + kp) - 1.6*(1 + d*(B + Vm - Rd))
         q1v = p*(B*(ca + kp) + x.) + 1.6*(Vm - Rd)*(1 + d*B)
         q0v = -p*B*x.
         XV = ifelse(Vm - Rd > B, B, Vm - Rd) ## upper limit of interval of valid A
         x. = (0.25*J - Rd)*ca - (0.25*J + 2*Rd)*gams
         q3j = 1.6*d 
         q2j = -p*(ca + 2*gams) - 1.6*(1 + d*(B + 0.25*J - Rd))
         q1j = p*(B*(ca + 2*gams) + x.) + 1.6*(0.25*J - Rd)*(1 + d*B)
         q0j = -p*B*x.
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
      
      d = s*deltaw/K ## d was earlier defined as this + p/gm to calculate gtc
      g = p*(B - A)/(1 + d*(B - A)) ## this is the new model
      g = ifelse(g>=0, g, 0)
      gtc=g/(1.6 + g*rm)
      ci = ca - 1.6*A/gtc 
      
      return(c(g, A, ci, lim, root))
   } else {
      ## if h < 0, return NAs
      return(c(NA,NA,NA,NA,NA))
   } ## end if h
} ## end gasxNEW()

gasxUSO <- function(PPFD, ca, TleafC, deltaw, Vm25, Jm25, gm25, phij, thetaj, thetaa, g1, g0, spp) {
  h = 1 - deltaw/(6.1078*exp(17.27*TleafC/(237.3 + TleafC)))
  
  if(h>=0) {
    if(spp=="grapevine") {
      TlK = TleafC + 273.15
      TleafK=TlK
      KC25 = 272.38
      KO25 = 165.82
      gams25 = KC25 * 210 * 0.21 / (2 * KO25) #210 = 21% * 10 to get mf
      
      itlk = 1 / TlK
      irt = 0.120273 * itlk * 0.0033557 #0.12... is 1/R, .0033557 is 1/298
      iit = 298.15 * itlk * itlk
      
      #''Vm T response is for grapevine
      dvm = exp(61115 * (TlK - 298.15) * irt)
      Vm = Vm25 * dvm
      drd = exp(44790 * (TlK - 298.15) * irt) #'activation here is from Antonio
      Rd = 0.0089 * Vm25 * drd                 # ' Rd25 vs Vm25 is from dePury
      dKc = exp(80990 * (TlK - 298.15) * irt)
      Kc = KC25 * dKc
      dKo = exp(23720 * (TlK - 298.15) * irt)
      Ko = KO25 * dKo
      kp = Kc*(1 - 210/Ko)
      dgams = exp(24460 * (TlK - 298.15) * irt)
      gams = gams25 * dgams
      
      #''' gm temperature response for grapevine
      dgm = exp(-0.5 * (log(TleafC / 36.7476) / 0.83902) ^ 2)
      gm = gm25 * dgm
      rm = 1/gm
      
      #''' Jm temperature response for grapevine
      ea. = 55791.8767
      s. = 367.1633804
      ed. = 114804.1868
      rgas = 8.3143
      djm = exp(ea. * (TlK - 298.15) / (rgas * 298.15 * TlK))
      djm = djm / (1 + exp((s. * TlK - ed.) / (rgas * TlK)))
      djm = djm * (1 + exp((s. * 298.15 - ed.) / (rgas * 298.15)))
      Jm = Jm25 * djm
    } else if(spp=="almond") {
      TlK = TleafC + 273.15
      TleafK=TlK
      KC25 = 272.38
      KO25 = 165.82
      gams25 = KC25 * 210 * 0.21 / (2 * KO25) # 210 = 21% * 10 to get mf
      
      itlk = 1 / TlK
      irt = 0.120273 * itlk * 0.0033557 # 0.12... is 1/R, .0033557 is 1/298
      iit = 298.15 * itlk * itlk
      
      rgas = 0.0083143
      cvm = 31.58
      eavm = 78.1
      svm = 0.4964
      edvm = 154.21
      
      cjm = 13.775
      eajm = 34.11
      sjm = 0.5824
      edjm = 184.21
      
      cgm = 19.966
      eagm = 49.48
      sgm = 1.3599
      edgm = 426.4
      
      rd25 = 0.75
      q10rd = 2
      
      dvm = exp(cvm - eavm / (rgas * TlK)) / (1 + exp((svm * TlK - edvm) / (rgas * TlK)))
      Vm = Vm25 * dvm
      djm = exp(cjm - eajm / (rgas * TlK)) / (1 + exp((sjm * TlK - edjm) / (rgas * TlK)))
      Jm = Jm25 * djm
      dgm = exp(cgm - eagm / (rgas * TlK)) / (1 + exp((sgm * TlK - edgm) / (rgas * TlK)))
      gm = gm25 * dgm
      drd = q10rd ^ ((TleafC - 25) / 10)
      Rd = rd25 * drd
      dKc = exp(80990 * (TlK - 298.15) * irt)
      Kc = KC25 * dKc
      dKo = exp(23720 * (TlK - 298.15) * irt)
      Ko = KO25 * dKo
      dgams = exp(24460 * (TlK - 298.15) * irt)
      gams = gams25 * dgams
      
      kp = Kc * (1 + 210 / Ko)
      
      rm = 1/gm
    } else if(spp=="olive") {
      TlK = TleafC + 273.15
      TleafK=TlK
      KC25 = 272.38
      KO25 = 165.82
      gams25 = KC25 * 210 * 0.21 / (2 * KO25) # 210 = 21% * 10 to get mf
      
      itlk = 1 / TlK
      irt = 0.120273 * itlk * 0.0033557 # 0.12... is 1/R, .0033557 is 1/298
      iit = 298.15 * itlk * itlk
      
      dvm = exp(46525 * (TlK - 298.15) * irt)
      Vm = Vm25 * dvm
      drd = exp(44790 * (TlK - 298.15) * irt) #'activation here is from Antonio
      Rd = 0.0089 * Vm25 * drd        #          ' Rd25 vs Vm25 is from dePury
      dKc = exp(80990 * (TlK - 298.15) * irt)
      Kc = KC25 * dKc
      dKo = exp(23720 * (TlK - 298.15) * irt)
      Ko = KO25 * dKo
      dgams = exp(24460 * (TlK - 298.15) * irt)
      gams = gams25 * dgams
      
      #''' gm temperature response; parameters fitted from Antonio's lab
      ea = 79999.99606
      s. = 381.3603738
      ed = 113148.9517
      rgas = 8.3143
      dgm = exp(ea * (TlK - 298.15) / (rgas * 298.15 * TlK))
      dgm = dgm / (1 + exp((s. * TlK - ed) / (rgas * TlK)))
      dgm = dgm * (1 + exp((s. * 298.15 - ed) / (rgas * 298.15)))
      gm = gm25 * dgm
      
      #''' Jm temperature response; parameters fitted from Antonio's lab
      ea = 48410.13
      s. = 357.96
      ed = 113149.02
      rgas = 8.3143
      djm = exp(ea * (TlK - 298.15) / (rgas * 298.15 * TlK))
      djm = djm / (1 + exp((s. * TlK - ed) / (rgas * TlK)))
      djm = djm * (1 + exp((s. * 298.15 - ed) / (rgas * 298.15)))
      Jm = Jm25 * djm
      
      kp = Kc * (1 + 210 / Ko)
      
      rm = 1/gm
    }
    i <- PPFD
    J <- ((-(-Jm-phij*i) - sqrt((-Jm-phij*i)*(-Jm-phij*i) - 4*thetaj*Jm*phij*i))/(2*thetaj))
    
    ## in BB, A/ca is multiplied by product of slope m and RH h
    ## in USO, the equivalent factor is 1.6*(1 + g1/sqrt(D)), where D is VPD in kPa
    ## models are otherwise equivalent
    
    ## in soln used in this version of the code, the 1.6 is not part of mh
    mh = (1 + g1/sqrt(deltaw/10))
    
    if(J>0) {
      ##### nb.: quadratic is solved in terms of A, not ci in this case
      ## quadratic coeffs for V-lim case
      av = rm
      bv = ca/mh - ca - kp - rm*(Vm - Rd)
      cv = (Vm - Rd)*ca*(1 - 1/mh) - (Vm*gams + Rd*kp)
      
      ## quadratic coeffs for J-lim case
      aj = rm
      bj = ca/mh - ca - 2*gams - rm*(0.25*J - Rd)
      cj = (0.25*J - Rd)*ca*(1 - 1/mh) - (0.25*J*gams + Rd*2*gams)
      
      ## ci solutions
      Av = 0; Aj = 0
      Av = ifelse(bv*bv>4*av*cv, (0.5/av)*(-bv - sqrt(bv*bv - 4*av*cv)), 0)
      Aj = ifelse(bj*bj>4*aj*cj, (0.5/aj)*(-bj - sqrt(bj*bj - 4*aj*cj)), 0)
      
      ## calculate PS
      A = (0.5/thetaa)*((Av + Aj) - sqrt((Av + Aj)^2 - 4*thetaa*Av*Aj))

      lim = ifelse(Av < Aj, 1, 2)
    } else {
      lim = 2
      A = -Rd
    }
    g = 1.6*mh*A/ca + g0 ## stomatal conductance to water vapor
    g <- ifelse(g > 0, g, 0)
    G = g/(rm*g + 1.6) ## total conductance to CO2
    cc = ca - A/G
    
    return(c(g, J, Av, Aj))
  } else {
    return(c(NA,NA,NA,NA))
  }
} ## end gasxUSO()

gasxBB <- function(PPFD, ca, TleafC, deltaw, Vm25, Jm25, gm25, phij, thetaj, thetaa, g1, g0, spp) {
   h = 1 - deltaw/(6.1078*exp(17.27*TleafC/(237.3 + TleafC)))
   
   if(h>=0) {
      if(spp=="grapevine") {
         TlK = TleafC + 273.15
         TleafK=TlK
         KC25 = 272.38
         KO25 = 165.82
         gams25 = KC25 * 210 * 0.21 / (2 * KO25) #210 = 21% * 10 to get mf
         
         itlk = 1 / TlK
         irt = 0.120273 * itlk * 0.0033557 #0.12... is 1/R, .0033557 is 1/298
         iit = 298.15 * itlk * itlk
         
         #''Vm T response is for grapevine
         dvm = exp(61115 * (TlK - 298.15) * irt)
         Vm = Vm25 * dvm
         drd = exp(44790 * (TlK - 298.15) * irt) #'activation here is from Antonio
         Rd = 0.0089 * Vm25 * drd                 # ' Rd25 vs Vm25 is from dePury
         dKc = exp(80990 * (TlK - 298.15) * irt)
         Kc = KC25 * dKc
         dKo = exp(23720 * (TlK - 298.15) * irt)
         Ko = KO25 * dKo
         kp = Kc*(1 - 210/Ko)
         dgams = exp(24460 * (TlK - 298.15) * irt)
         gams = gams25 * dgams
         
         #''' gm temperature response for grapevine
         dgm = exp(-0.5 * (log(TleafC / 36.7476) / 0.83902) ^ 2)
         gm = gm25 * dgm
         rm = 1/gm
         
         #''' Jm temperature response for grapevine
         ea. = 55791.8767
         s. = 367.1633804
         ed. = 114804.1868
         rgas = 8.3143
         djm = exp(ea. * (TlK - 298.15) / (rgas * 298.15 * TlK))
         djm = djm / (1 + exp((s. * TlK - ed.) / (rgas * TlK)))
         djm = djm * (1 + exp((s. * 298.15 - ed.) / (rgas * 298.15)))
         Jm = Jm25 * djm
      } else if(spp=="almond") {
         TlK = TleafC + 273.15
         TleafK=TlK
         KC25 = 272.38
         KO25 = 165.82
         gams25 = KC25 * 210 * 0.21 / (2 * KO25) # 210 = 21% * 10 to get mf
         
         itlk = 1 / TlK
         irt = 0.120273 * itlk * 0.0033557 # 0.12... is 1/R, .0033557 is 1/298
         iit = 298.15 * itlk * itlk
         
         rgas = 0.0083143
         cvm = 31.58
         eavm = 78.1
         svm = 0.4964
         edvm = 154.21
         
         cjm = 13.775
         eajm = 34.11
         sjm = 0.5824
         edjm = 184.21
         
         cgm = 19.966
         eagm = 49.48
         sgm = 1.3599
         edgm = 426.4
         
         rd25 = 0.75
         q10rd = 2
         
         dvm = exp(cvm - eavm / (rgas * TlK)) / (1 + exp((svm * TlK - edvm) / (rgas * TlK)))
         Vm = Vm25 * dvm
         djm = exp(cjm - eajm / (rgas * TlK)) / (1 + exp((sjm * TlK - edjm) / (rgas * TlK)))
         Jm = Jm25 * djm
         dgm = exp(cgm - eagm / (rgas * TlK)) / (1 + exp((sgm * TlK - edgm) / (rgas * TlK)))
         gm = gm25 * dgm
         drd = q10rd ^ ((TleafC - 25) / 10)
         Rd = rd25 * drd
         dKc = exp(80990 * (TlK - 298.15) * irt)
         Kc = KC25 * dKc
         dKo = exp(23720 * (TlK - 298.15) * irt)
         Ko = KO25 * dKo
         dgams = exp(24460 * (TlK - 298.15) * irt)
         gams = gams25 * dgams
         
         kp = Kc * (1 + 210 / Ko)
         
         rm = 1/gm
      } else if(spp=="olive") {
         TlK = TleafC + 273.15
         TleafK=TlK
         KC25 = 272.38
         KO25 = 165.82
         gams25 = KC25 * 210 * 0.21 / (2 * KO25) # 210 = 21% * 10 to get mf
         
         itlk = 1 / TlK
         irt = 0.120273 * itlk * 0.0033557 # 0.12... is 1/R, .0033557 is 1/298
         iit = 298.15 * itlk * itlk
         
         dvm = exp(46525 * (TlK - 298.15) * irt)
         Vm = Vm25 * dvm
         drd = exp(44790 * (TlK - 298.15) * irt) #'activation here is from Antonio
         Rd = 0.0089 * Vm25 * drd        #          ' Rd25 vs Vm25 is from dePury
         dKc = exp(80990 * (TlK - 298.15) * irt)
         Kc = KC25 * dKc
         dKo = exp(23720 * (TlK - 298.15) * irt)
         Ko = KO25 * dKo
         dgams = exp(24460 * (TlK - 298.15) * irt)
         gams = gams25 * dgams
         
         #''' gm temperature response; parameters fitted from Antonio's lab
         ea = 79999.99606
         s. = 381.3603738
         ed = 113148.9517
         rgas = 8.3143
         dgm = exp(ea * (TlK - 298.15) / (rgas * 298.15 * TlK))
         dgm = dgm / (1 + exp((s. * TlK - ed) / (rgas * TlK)))
         dgm = dgm * (1 + exp((s. * 298.15 - ed) / (rgas * 298.15)))
         gm = gm25 * dgm
         
         #''' Jm temperature response; parameters fitted from Antonio's lab
         ea = 48410.13
         s. = 357.96
         ed = 113149.02
         rgas = 8.3143
         djm = exp(ea * (TlK - 298.15) / (rgas * 298.15 * TlK))
         djm = djm / (1 + exp((s. * TlK - ed) / (rgas * TlK)))
         djm = djm * (1 + exp((s. * 298.15 - ed) / (rgas * 298.15)))
         Jm = Jm25 * djm
         
         kp = Kc * (1 + 210 / Ko)
         
         rm = 1/gm
      }
      i <- PPFD
      J <- ((-(-Jm-phij*i) - sqrt((-Jm-phij*i)*(-Jm-phij*i) - 4*thetaj*Jm*phij*i))/(2*thetaj))
      
      ## in soln used in this version of the code, the 1.6 is not part of mh
      mh = g1*h
      
      if(J>0) {
         ##### nb.: quadratic is solved in terms of A, not ci in this case
         ## quadratic coeffs for V-lim case
         av = rm
         bv = ca/mh - ca - kp - rm*(Vm - Rd)
         cv = (Vm - Rd)*ca*(1 - 1/mh) - (Vm*gams + Rd*kp)
         
         ## quadratic coeffs for J-lim case
         aj = rm
         bj = ca/mh - ca - 2*gams - rm*(0.25*J - Rd)
         cj = (0.25*J - Rd)*ca*(1 - 1/mh) - (0.25*J*gams + Rd*2*gams)
         
         ## ci solutions
         Av = 0; Aj = 0
         Av = ifelse(bv*bv>4*av*cv, (0.5/av)*(-bv - sqrt(bv*bv - 4*av*cv)), 0)
         Aj = ifelse(bj*bj>4*aj*cj, (0.5/aj)*(-bj - sqrt(bj*bj - 4*aj*cj)), 0)
         
         ## calculate PS
         A = (0.5/thetaa)*((Av + Aj) - sqrt((Av + Aj)^2 - 4*thetaa*Av*Aj))
         
         lim = ifelse(Av < Aj, 1, 2)
      } else {
         lim = 2
         A = -Rd
      }
      g = 1.6*mh*A/ca + g0 ## stomatal conductance to water vapor
      g <- ifelse(g > 0, g, 0)
      G = g/(rm*g + 1.6) ## total conductance to CO2
      cc = ca - A/G
      
      return(c(g, A, cc, lim))
   } else {
      return(c(NA,NA,NA,NA))
   }
} ## end gasxBB()


SE <- function(x, na.rm=FALSE) {
   if (na.rm) x <- na.omit(x)
   sqrt(var(x)/length(x))
}



### load data
species <- c("grapevine", "olive", "almond")
xall <- data.frame(matrix(nrow=0, ncol=22))
for(ss in species) {
   filename=paste0("celia ", ss, ".csv") ## data from (Celia) Dominguez-Rodriguez et al 2016 PCE
   x=read.delim(filename, header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
   x$spp <- ss
   x$TleafK <- x$tl + 273.15
   x$K25 <- x$k*((1.856e-11*exp(4209/x$TleafK + 0.04527*x$TleafK - 3.38e-5*x$TleafK*x$TleafK))/0.907784) #viscosity correction
   ## get midday K25
   if(ss!="olive") {
      y <- subset(x, time>0.5 & time <0.6)
   } else {
      y <- subset(x, time==0.5)
   }
   colnames(y)[22] <- "K25."
   x <- left_join(x, y[,c(1,22)], by="tree")
   x$K25 <- NULL
   colnames(x)[22] <- "K25"
   xall <- rbind(xall, x)
}
colnames(xall)[4] <- "deltaw"

x <- xall
x$trt <- substr(x$tree, 2, 4)
x$gswmod <- NA
x$A <- NA
x$ci <- NA
x$lim <- NA
ncol(x)


### generate default predictions, will override later
arglist=list(psic=-2.5, s=0.012, ro=3)
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
xall <- x

### code to fit new model, returns sse
fitmod <- function(pvec, xx) {
   arglist=list(psic=pvec[1], s=pvec[2], ro=pvec[3])
   y <- data.frame(t(mapply(gasxNEW, xx$ppfd, xx$ca, xx$tl, xx$deltaw, xx$psis, xx$vm, xx$jm, xx$gm, xx$phi, xx$thetaj, xx$thetaa, xx$K25, x$spp, MoreArgs=arglist)))
   y$gsw <- xx$gsw
   return(sum((y[,1]-y$gsw)^2))
}


### fit new model for each treatment x species combination (reg=WW, seq=WS)

x <- subset(xall, trt=="reg" & spp=="grapevine")
o <- optim(c(-3.5, 0.04, 3), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.rg <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yrg <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yrg$trt <- "reg"; yrg$spp <- "grapevine"

x <- subset(xall, trt=="seq" & spp=="grapevine")
o <- optim(c(-2.5, 0.04, 3), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.sg <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ysg <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ysg$trt <- "seq"; ysg$spp <- "grapevine"

x <- subset(xall, trt=="reg" & spp=="almond")
o <- optim(c(-3.3, 0.005, 10), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.ra <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yra <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yra$trt <- "reg"; yra$spp <- "almond"

x <- subset(xall, trt=="seq" & spp=="almond")
o <- optim(c(-2.7, 0.01, 3), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.sa <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ysa <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ysa$trt <- "seq"; ysa$spp <- "almond"

x <- subset(xall, trt=="reg" & spp=="olive")
o <- optim(c(-3.3, 0.005, 10), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.ro <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yro <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yro$trt <- "reg"; yro$spp <- "olive"

x <- subset(xall, trt=="seq" & spp=="olive")
o <- optim(c(-2.7, 0.01, 3), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.so <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yso <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yso$trt <- "seq"; yso$spp <- "olive"


ys <- rbind(yrg, ysg, yra, ysa, yro, yso); ys$mode="parameters fitted separately\nfor each treatment"


### fit new model for each species, for both treatments combined
x <- subset(xall, spp=="grapevine")
o <- optim(c(-3.5, 0.04, 3), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.rg <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yg <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yg$spp <- "grapevine"

x <- subset(xall, spp=="almond")
o <- optim(c(-3.3, 0.005, 10), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.ra <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ya <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ya$spp <- "almond"

x <- subset(xall, spp=="olive")
o <- optim(c(-3.3, 0.005, 10), fitmod, xx=x, method="Nelder-Mead")
arglist=list(psic=o$par[1], s=o$par[2], ro=o$par[3]); par.ro <- o$par
y <- data.frame(t(mapply(gasxNEW, x$ppfd, x$ca, x$tl, x$deltaw, x$psis, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$K25, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yo <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yo$spp <- "olive"

yt <- rbind(yg, ya, yo); yt$mode="both treatments have\nsame parameters"

y <- rbind(ys, yt)

y$mode <- factor(y$mode, levels=c("parameters fitted separately\nfor each treatment", "both treatments have\nsame parameters"))
y$trt <- ifelse(y$trt=="reg", "control", "drought")
y$trt <- factor(y$trt, levels=c("control", "drought"))
y$spp <- factor(y$spp, levels=c("grapevine", "olive", "almond"))
y$time <- 24*y$time
ynew <- y
ynew$model <- "NEW"



#################
###############
### repeat above for USO model

### load data, calculate Delta, Kplant and deltapsi
species <- c("grapevine", "olive", "almond")
xall <- data.frame(matrix(nrow=0, ncol=22))
for(ss in species) {
   filename=paste0("celia ", ss, ".csv")
   x=read.delim(filename, header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
   x$spp <- ss
   x$TleafK <- x$tl + 273.15
   x$K25 <- x$k*((1.856e-11*exp(4209/x$TleafK + 0.04527*x$TleafK - 3.38e-5*x$TleafK*x$TleafK))/0.907784) #viscosity correction
   ## get midday K25
   if(ss!="olive") {
      y <- subset(x, time>0.5 & time <0.6)
   } else {
      y <- subset(x, time==0.5)
   }
   colnames(y)[22] <- "K25."
   x <- left_join(x, y[,c(1,22)], by="tree")
   x$K25 <- NULL
   colnames(x)[22] <- "K25"
   xall <- rbind(xall, x)
}

colnames(xall)[4] <- "deltaw"

x <- xall
x$trt <- substr(x$tree, 2, 4)
x$gswmod <- NA
x$A <- NA
x$ci <- NA
x$lim <- NA
ncol(x)

arglist=list(g1=3, g0=0.02)
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
xall <- x


fitmod <- function(pvec, xx) {
   arglist=list(g1=pvec[1], g0=pvec[2])
   y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
   y$gsw <- xx$gsw
   return(sum((y[,1]-y$gsw)^2, na.rm=T))
}


x <- subset(xall, trt=="reg" & spp=="grapevine")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.rg <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yrg <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yrg$trt <- "reg"; yrg$spp <- "grapevine"

x <- subset(xall, trt=="seq" & spp=="grapevine")
o <- optim(c(1, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.sg <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ysg <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ysg$trt <- "seq"; ysg$spp <- "grapevine"

x <- subset(xall, trt=="reg" & spp=="almond")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ra <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yra <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yra$trt <- "reg"; yra$spp <- "almond"

x <- subset(xall, trt=="seq" & spp=="almond")
o <- optim(c(1, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.sa <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ysa <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ysa$trt <- "seq"; ysa$spp <- "almond"

x <- subset(xall, trt=="reg" & spp=="olive")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ro <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yro <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yro$trt <- "reg"; yro$spp <- "olive"

x <- subset(xall, trt=="seq" & spp=="olive")
o <- optim(c(1, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.so <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yso <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yso$trt <- "seq"; yso$spp <- "olive"


ys <- rbind(yrg, ysg, yra, ysa, yro, yso); ys$mode="parameters fitted separately\nfor each treatment"


x <- subset(xall, spp=="grapevine")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.g <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yg <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yg$spp <- "grapevine"

x <- subset(xall, spp=="almond")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ra <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ya <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ya$spp <- "almond"

x <- subset(xall, spp=="olive")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ro <- o$par
y <- data.frame(t(mapply(gasxUSO, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yo <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yo$spp <- "olive"

yt <- rbind(yg, ya, yo); yt$mode="both treatments have\nsame parameters"

y <- rbind(ys, yt)
y$mode <- factor(y$mode, levels=c("parameters fitted separately\nfor each treatment", "both treatments have\nsame parameters"))
y$trt <- ifelse(y$trt=="reg", "control", "drought")
y$trt <- factor(y$trt, levels=c("control", "drought"))
y$spp <- factor(y$spp, levels=c("grapevine", "olive", "almond"))
y$time <- 24*y$time
y$model <- "USO"
yuso <- y




#################
###############
### repeat for BB

### load data, calculate Delta, Kplant and deltapsi
species <- c("grapevine", "olive", "almond")
xall <- data.frame(matrix(nrow=0, ncol=22))
for(ss in species) {
   filename=paste0("celia ", ss, ".csv")
   x=read.delim(filename, header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
   x$spp <- ss
   x$TleafK <- x$tl + 273.15
   x$K25 <- x$k*((1.856e-11*exp(4209/x$TleafK + 0.04527*x$TleafK - 3.38e-5*x$TleafK*x$TleafK))/0.907784) #viscosity correction
   ## get midday K25
   if(ss!="olive") {
      y <- subset(x, time>0.5 & time <0.6)
   } else {
      y <- subset(x, time==0.5)
   }
   colnames(y)[22] <- "K25."
   x <- left_join(x, y[,c(1,22)], by="tree")
   x$K25 <- NULL
   colnames(x)[22] <- "K25"
   xall <- rbind(xall, x)
}

colnames(xall)[4] <- "deltaw"

x <- xall
x$trt <- substr(x$tree, 2, 4)
x$gswmod <- NA
x$A <- NA
x$ci <- NA
x$lim <- NA
ncol(x)

arglist=list(g1=3, g0=0.02)
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
xall <- x


fitmod <- function(pvec, xx) {
   arglist=list(g1=pvec[1], g0=pvec[2])
   y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
   y$gsw <- xx$gsw
   return(sum((y[,1]-y$gsw)^2, na.rm=T))
}


x <- subset(xall, trt=="reg" & spp=="grapevine")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.rg <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yrg <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yrg$trt <- "reg"; yrg$spp <- "grapevine"

x <- subset(xall, trt=="seq" & spp=="grapevine")
o <- optim(c(1, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.sg <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ysg <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ysg$trt <- "seq"; ysg$spp <- "grapevine"

x <- subset(xall, trt=="reg" & spp=="almond")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ra <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yra <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yra$trt <- "reg"; yra$spp <- "almond"

x <- subset(xall, trt=="seq" & spp=="almond")
o <- optim(c(1, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.sa <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ysa <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ysa$trt <- "seq"; ysa$spp <- "almond"

x <- subset(xall, trt=="reg" & spp=="olive")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ro <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yro <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yro$trt <- "reg"; yro$spp <- "olive"

x <- subset(xall, trt=="seq" & spp=="olive")
o <- optim(c(1, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.so <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yso <- x %>% group_by(time) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yso$trt <- "seq"; yso$spp <- "olive"


ys <- rbind(yrg, ysg, yra, ysa, yro, yso); ys$mode="parameters fitted separately\nfor each treatment"


x <- subset(xall, spp=="grapevine")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.g <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yg <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yg$spp <- "grapevine"

x <- subset(xall, spp=="almond")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ra <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
ya <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
ya$spp <- "almond"

x <- subset(xall, spp=="olive")
o <- optim(c(3, 0.02), fitmod, xx=x, method="Nelder-Mead")
arglist=list(g1=o$par[1], g0=o$par[2]); par.ro <- o$par
y <- data.frame(t(mapply(gasxBB, x$ppfd, x$ca, x$tl, x$deltaw, x$vm, x$jm, x$gm, x$phi, x$thetaj, x$thetaa, x$spp, MoreArgs=arglist)))
x[,c(24:27)] <- y[,c(1:4)]
yo <- x %>% group_by(time, trt) %>%  summarise(meangsw=mean(gsw), SEgsw = SE(gsw), meangswmod = mean(gswmod), SEgswmod=SE(gswmod))
yo$spp <- "olive"

yt <- rbind(yg, ya, yo); yt$mode="both treatments have\nsame parameters"

y <- rbind(ys, yt)
y$mode <- factor(y$mode, levels=c("parameters fitted separately\nfor each treatment", "both treatments have\nsame parameters"))
y$trt <- ifelse(y$trt=="reg", "control", "drought")
y$trt <- factor(y$trt, levels=c("control", "drought"))
y$spp <- factor(y$spp, levels=c("grapevine", "olive", "almond"))
y$time <- 24*y$time
y$model <- "BB"
ybb <- y



#############
### combine all results together and generate figure


yall <- rbind(ynew, yuso)
yall <- rbind(yall, ybb)
yall$model <- factor(yall$model, levels=c("NEW", "USO", "BB"))
yall$linesize <- ifelse(yall$model=="NEW", "A", "B")


pp <- ggplot() +
   geom_line(data=yall, aes(x=time, y=meangswmod, color=trt, linetype=model))+
   geom_line(data=ynew, aes(x=time, y=meangswmod, color=trt), linewidth=1)+
   geom_point(data=yall, aes(x=time, y=meangsw, color=trt), shape=15)+
   geom_errorbar(data=yall, aes(x=time, ymin=meangsw-SEgsw, ymax=meangsw+SEgsw, color=trt), width=0.01)+
   facet_grid(rows=vars(spp), cols=vars(mode), scales="free_y") +
   guides(color=guide_legend(title="treatment",
                             override.aes = list(size=0.1))) +
   scale_color_manual(values=c("blue", "red"))+
   scale_linetype_manual(values=c("solid", "longdash", "dotdash"))+
   ylab(expression(paste(italic("g")["sw"]," (mol ",m^-2,s^-1,")"))) +
   xlab("time (H GMT)") + theme_bw() + xlim(6.5,20.5); pp

library(egg)
pp <- tag_facet(pp)
pp <- pp + theme(strip.text = element_text(size=9, margin=margin(t=3, b=3, r=3, l=3)), 
                   strip.background = element_rect(fill="grey90"))
pp

ggsave("Figure 2.png", device="png", dpi=1000)

