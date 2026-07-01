setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("initialize.R")


################
################
###
### code to generate Figure 9 (response of stomata to temperature at constant deltaw)
###
################
################

## set para values
{  Vm25 = 50; Jm25 = 85; Rd25 = 0.01*Vm25; thetaa = 0.99; s = 0.02; ro = 0.3; K25 = 3; psic = -2; g1 <- g1.o; g0uso <- g0uso.o; mBB <- mBB.o; g0bb <- g0bb.o
   ppfd_ = 1000; ppfdmin = 0; ppfdmax = 1000; ca_ = 425; camin = 150; camax = 700; TleafC_ = 25; TleafCmin = 5; TleafCmax = 45; deltaw_ = 15; deltawmin = 5; deltawmax = 45; psisoil_ = 0; psisoilmin = psic; psisoilmax = 0; h_ = 0.5; hmin = 0.1; hmax = 0.9
   varlow <- c(ppfdmin, camin, TleafCmin, deltawmin, psisoilmin, hmin); varhi <- c(ppfdmax, camax, TleafCmax, deltawmax, psisoilmax, hmax)}
ppfdmin=100; ppfdmax=700

## generate model predictions for low and high PPFD
y <- g.vs.(TLEAF, PPFD, g1=g1, mBB=mBB, thetaa=thetaa, s=s, ro=ro, g0bb=g0bb, g0uso=g0uso, nx2=2)

## group by model and PPFD level and express gs relative to its value at 25C
y <- y %>% group_by(ppfd, model) %>% mutate(grel=gsw/gsw[TleafC==25])


#### load in previously published data for the response, from Mills et al 2024 PCE
x=read.delim("mills 2024 data.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
x$response <- as.factor(x$response)
x$rnum <- as.factor(x$rnum)
x <- x %>% group_by(response) %>% mutate(npoints = n())
z <- data.table(x)

## fit line to each response in that dataset and use it to express gs wrt 25C in each response
znew <- z[,	list(intercept = coef(lm(gsw ~ t))[1], slope = coef(lm(gsw ~ t))[2]), by = response]
znew$y25 <- znew$slope*25 + znew$intercept
z$y25 <- znew$y25[match(z$response, znew$response)]
z <- z %>% group_by(response) %>% mutate(grel=gsw/y25)

## plot simulations and data at high and low light
p1 <- ggplot() +
   geom_line(data=subset(y, ppfd==700), aes(x=TleafC, y=grel, linetype=model), size=1)+
  scale_linetype_manual(values=c("solid", "dashed", "dotdash"))+
  scale_size_manual(values=c(1.2, 0.7, 0.7))+
   xlab(varnamesb[3])+
   ylab("relative stomatal conductance") +
  geom_smooth(method="lm", formula = y ~ x + I(x^2), data=z, aes(x=t, y=grel, color="red"),
              size=0, alpha=0.1, fill="red", se=TRUE) +
  geom_line(stat="smooth", method="lm", formula = y ~ x + I(x^2), data=z, aes(x=t, y=grel, color="red"),
            size=0.5, alpha=1, se=FALSE) +
  ylim(0,2.5) + theme_bw() + theme(legend.key.width=unit(1.5, "cm")) +
  scale_color_manual(values=c("red"), labels=c("")) + 
  guides(color=guide_legend(title="observations")); p1

p2 <- ggplot() +
   geom_line(data=subset(y, ppfd==100), aes(x=TleafC, y=grel, linetype=model), size=1)+
   scale_linetype_manual(values=c("solid", "dashed", "dotdash"))+
   scale_size_manual(values=c(1.2, 0.7, 0.7))+
   xlab(varnamesb[3])+
  ylab("relative stomatal conductance") +
  geom_smooth(method="lm", formula = y ~ x + I(x^2), data=z, aes(x=t, y=grel, color="red"),
             size=0, alpha=0.1, fill="red", se=TRUE) +
   geom_line(stat="smooth", method="lm", formula = y ~ x + I(x^2), data=z, aes(x=t, y=grel, color="red"),
               size=0.5, alpha=1, se=FALSE) +
   ylim(0,2.5) + theme_bw() + theme(legend.key.width=unit(1.5, "cm")) +
  scale_color_manual(values=c("red"), labels=c("")) + 
  guides(color=guide_legend(title="observations")); p2

ggarrange(p1, p2, ncol=1, widths=c(1,1), common.legend=TRUE, legend="right",
          labels=c("(a)", "(b)"), label.x=c(0.15,0.15), label.y=c(0.98,0.98),
          font.label=list(face="plain"))
ggsave("Figure 9.png", device="png", dpi=1000)
