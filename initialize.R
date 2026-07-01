source("main functions.R")

## general paras for sims
PPFD=1; CA=2; TLEAF=3; DW=4; PSISOIL=5; H=6
varnames <- c("ppfd", "ca", "TleafC", "deltaw", "psisoil", "h")
varnamesb <- c(expression(paste("PPFD (",mu,"mol ",m^-2,s^-1,")")),
               expression(paste(italic("c")["a"]," (",mu,"mol ",mol^-1,")")),
               expression(paste(italic("T")["leaf"]," (deg C)")),
               expression(paste(Delta,italic("w")," (mmol/mol)")),
               expression(paste(psi["soil"], " (MPa)")),
               "RH")

varnamesc <- c(expression(atop(paste("PPFD"), paste("(",mu,"mol ",m^-2,s^-1,")"))),
               expression(atop(paste(italic("c")["a"]), paste("\n(",mu,"mol ",mol^-1,")"))),
               expression(atop(paste(italic("T")["leaf"]), paste("\n(deg C)"))),
               expression(atop(paste(Delta,italic("w")), paste("\n(mmol/mol)"))),
               expression(atop(paste(psi["soil"]), paste("\n(MPa)"))),
               "RH")


### set default parameter values and ranges
Vm25 = 50
Jm25 = 85
Rd25 = 0.01*Vm25
thetaa=0.99

## key paras for this model
s = 0.02 ## a for new model
ro = 0.3 ## b for new model
K25 = 3
psic = -2

ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000 ### default, min, max
ca_ = 425; camin = 150; camax = 700
TleafC_ = 25; TleafCmin = 5; TleafCmax = 45
deltaw_ = 15; deltawmin = 5; deltawmax = 45
psisoil_ = 0; psisoilmin = psic; psisoilmax = 0
h_ = 0.5; hmin = 0.1; hmax = 0.9
varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin)
varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)


### fit BB and USO models to new model at two points:
## darkness and 1000 ppfd
## ca=425, T=25, DW=15, psisoil=0

#### first BB
### g0 = glo + m*h*Rd/ca
### g0 = ghi - m*h*Ahi/ca
### -> m = (ghi - glo)/((h/ca)*(Ahi + Rd))
### solve iteratively by adjusting m to satisfy the last eqn

gnew.lo <- gasxNEW(0, 425, 25, 15, 0, Vm25, Jm25, Rd25, thetaa, K25, psic, s, ro)[1]
gnew.hi <- gasxNEW(1000, 425, 25, 15, 0, Vm25, Jm25, Rd25, thetaa, K25, psic, s, ro)[1]

fitBB <- function(m.in) {
   h <- getRH(15, 25)
   go.in <- gnew.lo + m.in*h*Rd25/ca_
   ABB.hi <- gasxBB(1000, 425, 25, h, Vm25, Jm25, Rd25, thetaa, m.in, go.in)[2]
   m.out = (gnew.hi - gnew.lo)/((h/425)*(ABB.hi + Rd25))
   return((m.out - m.in)^2)
}
o <- optimize(fitBB, interval=c(5,15))
mBB.o <- o$minimum
g0bb.o <- gnew.lo + mBB.o*getRH(15, 25)*Rd25/425


#### now USO
### g0 = glo + (1.6*Rd/ca)*(1 + g1/sqrt(D))
### g0 = ghi - (1.6*Ahi/ca)*(1 + g1/sqrt(D))
### -> g1 = [(ghi - glo)/((1.6/ca)*(Ahi + Rd)) - 1]*sqrt(D)
### solve iteratively by adjusting m to satisfy the last eqn

fitUSO <- function(g1.in) {
   D=1.5
   go.in <- gnew.lo + (1.6*Rd25/425)*(1 + g1.in/sqrt(D))
   AUSO.hi <- gasxUSO(1000, 425, 25, D, Vm25, Jm25, Rd25, thetaa, g1.in, go.in)[2]
   g1.out <- ((gnew.hi - gnew.lo)/((1.6/425)*(AUSO.hi + Rd25)) - 1)*sqrt(D)
   return((g1.out - g1.in)^2)
}
o <- optimize(fitUSO, interval=c(1,5))
g1.o <- o$minimum
g0uso.o <- gnew.lo + (1.6*Rd25/425)*(1 + g1.o/sqrt(1.5))

g1 <- g1.o
g0uso <- g0uso.o
mBB <- mBB.o
g0bb <- g0bb.o

### code to reset all parameters to initial default values
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}

