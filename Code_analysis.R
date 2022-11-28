### Wild selection line recruits - Data visualization and analysis
### ---------------------------------------------------------------------------
###
###
### Author R script: M Lindner (M.Lindner@nioo.knaw.nl), v2022


### Set-up 
### ---------------------------------------------------------------------------

### environment and R packages
options(width=200)

library(lubridate)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

library(rethinking)
library(cmdstanr)

setwd("/Users/melanielindner/Documents/NIOO/Projects/ERC_phenotypic_data/WIld_selection_lines/EggSwapping_AllYears")

## INFO
# Sexe:       0 - unknown; 1 - female; 2 - male
# BroodType:  0 - first clutch; 1 - replacement clutch; 2 - second clutch after successful first clutch; 6 - probably replacement clutch of first clutch
# Include:    0 - fitness NOT ok; 1 - fitness ok


### load data
# selection line recruits
Recruits <- read.table("Tbl_SelectLineRecruitsAtHV_new.csv", sep=",", header=TRUE)
Recruits[duplicated(Recruits),]

# remove duplicated row
Recruits <- Recruits[!duplicated(Recruits),]

# selection line fledglings
Fledglings <- read.table("Tbl_SelectLineChicksAtHV_new.csv", sep=",", header=TRUE)
Fledglings$DeviceNumber <- toupper(Fledglings$DeviceNumber)
# remove duplicate entry from data: transponder was changed for 522982 -> double entries
Fledglings[Fledglings$IndivID==522982,]
Fledglings <- Fledglings[Fledglings$DeviceNumber!="01101721A8",]

# selection line eggs:
Eggs <- read.table("Tbl_SelectLineEggsAtHV_new.csv", sep=",", header=TRUE)

# fitness data of (non-selection line) local recruits
Fitness <- read.table("Tbl_Fitness_GT_HV_AllClutches.csv", sep=",", header=TRUE)

# local (non-selection line) fledglings
FledglingsHV_new <- read.table("Tbl_AllChicksAtHV_new.csv", sep=",", header=TRUE)
FledglingsHV <- FledglingsHV_new[FledglingsHV_new$BroodYear!=2021,]


### Introduction of selection line eggs from the aviaries to the wild
### ---------------------------------------------------------------------------

### Explorative analysis

Eggs.temp <- Eggs[Eggs$BroodYear>=2017,]
Eggs.temp$SelectionLineNumeric <- ifelse(grepl("E", Eggs.temp$SelectionLine), 1, 2)
Eggs.temp$count <- 1

# get number of eggs per line and year:
Eggs.temp.noF5 <- Eggs.temp[!grepl("F5", Eggs.temp$SelectionLine),]
N_Eggs <- aggregate(count ~ BroodYear + SelectionLineNumeric, data=Eggs.temp.noF5, sum)
# Table S34
write.table(N_Eggs, "out/N_Eggs.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t") 

# get number of social broods per (line and) year: mixed broods and broods composed of eggs from one selection line only
Eggs.temp.HvBroodID <- Eggs.temp.noF5[!duplicated(Eggs.temp.noF5[,c(9,8,14,15)]), c(9,8,14,15)]
mixed <- Eggs.temp.HvBroodID[duplicated(Eggs.temp.HvBroodID$BrIDatHV),]$BrIDatHV
Eggs.temp.HvBroodID_mixed <- Eggs.temp.HvBroodID[Eggs.temp.HvBroodID$BrIDatHV %in% mixed,]
Eggs.temp.HvBroodID_one.line <- Eggs.temp.HvBroodID[!Eggs.temp.HvBroodID$BrIDatHV %in% mixed,]

# mixed broods
N_HVBroodID_mixed <- aggregate(count ~ BroodYear, data=Eggs.temp.HvBroodID_mixed, sum)

# broods composed of eggs from one selection line only
N_HVBroodID_one.line <- aggregate(count ~ BroodYear + SelectionLineNumeric, data=Eggs.temp.HvBroodID_one.line, sum)
# Table S33
write.table(N_HVBroodID_one.line, "out/N_HVBroodID_one.line.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t") 



###  Local recruitment probability of fledglings at the study site
### ---------------------------------------------------------------------------

### Prepare data

# 1. selection line fledglings (here often 'ERC' birds..)

# which years included
unique(Fledglings$BroodYear)
# remove birds born earlier than 2017 and chicks that did not fledge
ERC_Fledglings <- Fledglings[Fledglings$BroodYear>=2017 & Fledglings$Fledgling==1,]

# all locations HV locations? Location of brood, ringing capture and d15 capture in data frame
ERC_Fledglings[!grepl("HV", ERC_Fledglings$BroodLocationatHV),] # Yes
ERC_Fledglings[!grepl("HV", ERC_Fledglings$RingLocationAtHV),] # Yes
ERC_Fledglings[!grepl("HV", ERC_Fledglings$D15LocationAtHV),] # Yes

# sex info ok?
table(ERC_Fledglings$Sexe) ## one bird with no sex info
ERC_Fledglings$Sexe[ERC_Fledglings$Sexe==0] <- NA # set to NA rather than 0

# missing data?
help_ERC <- ERC_Fledglings[!complete.cases(ERC_Fledglings), c(2,6,12:17)] # used below for comparison with all data from study site
help_ERC$type <- rep("SL", nrow(help_ERC))

# 2. local (non-selection line) fledglings from HV 

# remove selection line birds
ERCs_in_wild <- intersect(FledglingsHV$IndivID, Fledglings$IndivID) # get all selection line birds
# note: number of birds (724) matches number of birds in data.frame 'Fledglings' (726) when two birds from 2016 are removed (and 2016 is not considered here)
FledglingsHV.0 <- FledglingsHV[!FledglingsHV$IndivID %in% ERCs_in_wild,] # remove all selection line birds from data

# which years included
unique(FledglingsHV.0$BroodYear)

# keep birds born 2017-2019 and remove chicks that did not fledge
HV_Fledglings.0 <- FledglingsHV.0[FledglingsHV.0$BroodYear>=2017 & FledglingsHV.0$BroodYear<=2019 & FledglingsHV.0$Fledgling==1,]

# all birds from HV? Location of brood, ringing capture and d15 capture in data frame
HV_Fledglings.0[!grepl("HV", HV_Fledglings.0$BroodLocationatHV),] # Yes
HV_Fledglings.0[!grepl("HV", HV_Fledglings.0$RingLocationAtHV),] # Yes
HV_Fledglings.0[!grepl("HV", HV_Fledglings.0$D15LocationAtHV),] # Yes

# sex info ok?
table(HV_Fledglings.0$Sexe) 
# sex info not available for many birds; expected as chicks with local genetic parents were not sexed in the lab and hence sex info is only available for birds that recruited into the study area

# missing data problems:
# in BAS (our data base) some capture data were missing or incomplete, fixed in BAS now but here the old file version was used and hence these data issues are fixed manually here
nrow(HV_Fledglings.0[!complete.cases(HV_Fledglings.0),]) # 27 birds with missing data now
nrow(HV_Fledglings.0[!complete.cases(HV_Fledglings.0[,-(13:16)]),-(13:16)]) # all related to day 15 capture

# fix problems
# 1. remove two chicks that were found dead
HV_Fledglings.1 <- HV_Fledglings.0[HV_Fledglings.0$Condition!="5 - nestling found dead in nestbox", ]
# 2. correct weight data for chicks (also for chicks that do not have missing data):
HV_Fledglings.1$Weight[HV_Fledglings.1$RingNumber=="BD...59071"] <- 1721
HV_Fledglings.1$Weight[HV_Fledglings.1$RingNumber=="BD...88439"] <- 1690
HV_Fledglings.1$Weight[HV_Fledglings.1$RingNumber=="BK...03197"] <- 1614
HV_Fledglings.1$Weight[HV_Fledglings.1$RingNumber=="BK...05405"] <- 1966
HV_Fledglings.1$Weight[HV_Fledglings.1$RingNumber=="BK...05865"] <- 1981
HV_Fledglings.1$Weight[HV_Fledglings.1$RingNumber=="BD...86983"] <- 2000

# combine data from selection line and local (non-selection line) birds
unique(ERC_Fledglings$SelectionLine)
ERC_Fledglings$SelectionLineNumeric <- ifelse(grepl("E", ERC_Fledglings$SelectionLine),1,2)
HV_Fledglings.1$SelectionLineNumeric <- rep(0, nrow(HV_Fledglings.1))

FINAL_Dat_Recruitment.0 <- rbind(ERC_Fledglings[,c(2,4,20,6,10,14:17,19,12)], HV_Fledglings.1[,c(2,3,18,4,8,12:15,17,10)])

# some more checks..
# are fledglings marked as 'recruited' in data.frame 'Fledglings' present present in data.frame 'Fitness'?
Recruited_in_fledgling_table <- FINAL_Dat_Recruitment.0[FINAL_Dat_Recruitment.0$Recruit==1,]
Not_in_fitness <- Recruited_in_fledgling_table[!Recruited_in_fledgling_table$RingNumber %in% c(Fitness$RingNumberFemale, Fitness$Father),]; nrow(Not_in_fitness)
# 17 birds not in data.frame 'Fitness'
# of these four are selection line recruits that recruited into another study site but the Hoge Veluwe
# The 'recruited' tag in data.frame 'Fledglings' is considering any recruitment irrespective of the study area. The data.frame 'Fitness' only includes local recruits, i.e. only includes birds recruiting as breeding bird into the Hoge Veluwe study site

# here we only consider local recruits (so probability of *local* recruitment)
# set 'recruit' of birds breeding in some other location to '0'
FINAL_Dat_Recruitment.good <- FINAL_Dat_Recruitment.0[!FINAL_Dat_Recruitment.0$RingNumber %in% Not_in_fitness$RingNumber,]
FINAL_Dat_Recruitment.bad <- FINAL_Dat_Recruitment.0[FINAL_Dat_Recruitment.0$RingNumber %in% Not_in_fitness$RingNumber,]
FINAL_Dat_Recruitment.bad$Recruit <- rep(0, nrow(FINAL_Dat_Recruitment.bad))
FINAL_Dat_Recruitment <- rbind(FINAL_Dat_Recruitment.good , FINAL_Dat_Recruitment.bad)

# save data:
save(FINAL_Dat_Recruitment, file="data/FINAL_Dat_Recruitment.RData")


### Exploratory analysis

# get number of selection line fledglings that recruited per line, sex and year:
FINAL_Dat_Recruitment.temp <- FINAL_Dat_Recruitment[FINAL_Dat_Recruitment$SelectionLineNumeric!=0,]
FINAL_Dat_Recruitment.temp$count <- 1
N_Recruit <- aggregate(Recruit ~ BroodYear + SelectionLineNumeric + Sexe, data=FINAL_Dat_Recruitment.temp, sum)
N_Total <- aggregate(count ~ BroodYear + SelectionLineNumeric + Sexe, data=FINAL_Dat_Recruitment.temp, sum)
if(unique(as.vector(N_Recruit[,1:3]==N_Total[,1:3]))) N_Recruit$Total <- N_Total$count
N_Recruit$Propportion_Recruit <- N_Recruit$Recruit/N_Recruit$Total
N_Recruit.ordered <- N_Recruit[order(N_Recruit$BroodYear, N_Recruit$SelectionLineNumeric, N_Recruit$Sexe),]
# Table S35
write.table(N_Recruit.ordered, "out/N_Recruit.Selection.lines.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# get April date for d15 measurement (as proxy for fledgling date):
FINAL_Dat_Recruitment$D15Catch_date <- dmy(FINAL_Dat_Recruitment$D15Catch, tz="EST")
FINAL_Dat_Recruitment$D15CatchApril <- as.numeric(yday(FINAL_Dat_Recruitment$D15Catch_date)-yday(as.POSIXct(paste(year(FINAL_Dat_Recruitment$D15Catch_date), "-03-31", sep=""), format='%Y-%m-%d')))


### Statistical analysis

# fix brood location ID format; use intergers from 1 to n broods rather than character strings
FINAL_Dat_Recruitment_all <- FINAL_Dat_Recruitment
temp <- FINAL_Dat_Recruitment_all[,c(4,11)]
temp.0 <- temp[!duplicated(temp),]
temp.0$Brood <- 1:nrow(temp.0)
For_Stats_Recruitment_all <- merge(FINAL_Dat_Recruitment_all, temp.0, by=c("BroodYear", "D15LocationAtHV"))

# use 3 (instead of 0 for 'wild' fledglings)
For_Stats_Recruitment_all$help <- ifelse(For_Stats_Recruitment_all$SelectionLineNumeric==0,3,0)
For_Stats_Recruitment_all$SelectionLineNumeric2 <- For_Stats_Recruitment_all$SelectionLineNumeric+For_Stats_Recruitment_all$help

# MODEL for fledgling weight (Eq. S8) of selection line and local (non-selection line) fledglings:

# prepare variables for model
d <- For_Stats_Recruitment_all[complete.cases(For_Stats_Recruitment_all$Weight),] # remove fledglings with no weight info
d$Date_std <- standardize(d$D15CatchApril) # standardize April date of d15 measurements (relative timing)
d$Weight_std <- standardize(d$Weight) # standardize weight
d$Line <- as.factor(d$SelectionLineNumeric2) # make sure is factor/integer starting at 1
d$Year <- as.factor(d$BroodYear-2016) # make sure is factor/integer starting at 1
d$Brood_factor <- as.factor(d$Brood) # make sure is factor/integer starting at 1
table(d$Brood_factor)[table(d$Brood_factor)==1] # 7 broods with one chick only. Might cause trouble..

# input data:
dat <- list(
  R=d$Weight_std,
  A=d$Date_std,
  L=d$Line,
  B=d$Brood_factor,
  Y=d$Year
)

# fit model
mFW.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + x[B]*sigma_g + bL[L] + bA*A,
    z[Y] ~ dnorm(0, 1),
    x[B] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5),
    sigma ~ dhalfnorm(0,1),
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bA ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[B]:g <<- x*sigma_g
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))

# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
# Warnings about hitting the maximum treedepth are not as serious as warnings about divergent transitions. While divergent transitions are a validity concern, hitting the maximum treedepth is an efficiency concern.
# fix: increase max_treedepth

mFW.1b <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + x[B]*sigma_g + bL[L] + bA*A,
    z[Y] ~ dnorm(0, 1),
    x[B] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5),
    sigma ~ dhalfnorm(0,1),
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bA ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[B]:g <<- x*sigma_g
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(max_treedepth=15, adapt_delta=0.99))
save(mFW.1b, file="temp/mFW.1b.multilevel.RData") # save model output!

## no warnings about maximum tree-depth !

# Information Criteria and Pareto-Smoothed Importance Sampling Cross-Validation (PSIS)
PSIS.FW <- PSIS(mFW.1b, pointwise=TRUE)
d[PSIS.FW$k>0.5,]
nrow(PSIS.FW[PSIS.FW$k>0.5,]) # 46 values with k>0.5, but no obvious mistake in data, just really extreme values for fledgling weight whih is not unusual

# visualize PSIS outcome
par(mfrow = c(1,1))
par(mar = c(3.5,3.5,1.2,1.2), mgp=c(2.0, 0.5, 0))
plot(NULL, xlim=c(-0.1,1.15), ylim=c(0,8.1),
     xlab="Pareto k", ylab="PSIS penalty", cex.main=1)
points(PSIS.FW$k, PSIS.FW$penalty, pch=1)
abline(v=0.5, lty=2)

# get model summary:
(out <- precis(mFW.1b, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S15
write.table(out.0, "out/Bayes_FledglingWeight_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mFW.1b)), row.names(precis(mFW.1b, depth=2))[grepl("z", row.names(precis(mFW.1b, depth=2)))], row.names(precis(mFW.1b, depth=2))[grepl("bL", row.names(precis(mFW.1b, depth=2)))])
par(mfrow = c(1,1))
# Fig. S29
pdf(file = "Plots/Bayes_FledglingWeight_traceplots.pdf", width = 6, height = 8)
traceplot(mFW.1b, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# Fig. S30
pdf(file = "Plots/Bayes_FledglingWeight_trankplot.pdf", width = 6, height = 8)
trankplot(mFW.1b, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior mean and ci for wegiht of selection line and local (non-selection line) fledglings
# prepare data for prediction; use unique combinations of year, selection line and brood location ID as in fitted data
dat_help.0 <- d[!duplicated(d[,19:21]), 19:21]
dat_help <- dat_help.0[order(dat_help.0$Line, dat_help.0$Year, dat_help.0$Brood_factor),]
d_pred <- list(L=dat_help$Line, Y=dat_help$Year, B=dat_help$Brood_factor, A=rep(0, times=nrow(dat_help)))

# get posterior
p_post <- link(mFW.1b, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:3) {
  dat.temp <- p_post[,dat_help$Line==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
diff_e_w <- new.temp[,3]-new.temp[,1]
d_e_w_mu <- mean(diff_e_w); d_e_w_ci <- PI(diff_e_w)
diff_l_w <- new.temp[,3]-new.temp[,2]
d_l_w_mu <- mean(diff_l_w); d_l_w_ci <- PI(diff_l_w)

post.out <- data.frame(line=c(1,2,3,"diff_e_l", "diff_e_w", "diff_l_w"), mu=c(p_mu, d_mu, d_e_w_mu, d_l_w_mu), ci.5.5=c(p_ci[1,],d_ci[1],d_e_w_ci[1],d_l_w_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2],d_e_w_ci[2],d_l_w_ci[2]))
# Table S16
write.table(post.out, "out/Bayes_FledglingWeight_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction:
# prepare plot
n_line <- data.frame(table(d$L)[-3]); names(n_line)[2] <- "n"; n_line$lab <- ifelse(n_line$Var1==1, "early", "late")
d <- For_Stats_Recruitment_all[complete.cases(For_Stats_Recruitment_all$Weight) & For_Stats_Recruitment_all$SelectionLineNumeric!=0,]
d$R <- standardize(d$Weight) # standardize weight
d$L <- as.factor(d$SelectionLineNumeric)

# Fig. S5
pdf(file = "Plots/Main_Bayes_FledgingWeight.pdf", width = 7, height = 3)
layout(t(matrix(c(1,2,3))), width=c(1.25,0.9, 1.5))

# A
par(mar = c(4,4,1.2,1), mgp=c(2, 0.5, 0))
dat_out <- post.out[1:2,]; dat_out$line <- as.numeric(dat_out$line)
plot(NULL, xlim=c(0.7,2.3), ylim=c(-5,2.15),
     xlab="", ylab="Standardized fledgling weight",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(n_line$lab, " (n=", n_line$n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-5,2,1), labels=seq(-5,2,1), lwd = 0, lwd.ticks = 1)
points(jitter(as.numeric(d$L),0.5), d$R, pch=1, cex=1.5, col=ifelse(d$L==1, "darkorange", "darkslategrey"))
arrows(x0=dat_out$line, y0=dat_out$ci.5.5, x1=dat_out$line, y1=dat_out$ci.94.5, code=3, angle=90, length=0.2, col="black", lwd=1)
points(dat_out$line, dat_out$mu, pch=21, cex=3, col="black", bg=ifelse(dat_out$line==1, "darkorange", "darkslategrey"))
text(0.7+0.03,2.15, "A", cex=1.2, pos=1, offset=0, font=2)

# B
par(mar = c(4,4,1.2,1), mgp=c(2, 0.5, 0))
dat_out_diff <- post.out[4,]; dat_out$line <- as.numeric(dat_out$line)
plot(NULL, xlim=c(0.7,1.3), ylim=c(-0.01,0.2),
     xlab="", ylab="Difference in standardized\nfledgling weight (late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(0,0.2,0.05), labels=seq(0,0.2,0.05), lwd = 0, lwd.ticks = 1)
arrows(x0=1, y0=dat_out_diff$ci.5.5, x1=1, y1=dat_out_diff$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(1, dat_out_diff$mu, pch=19, cex=3, col="black")
abline(h=0, lty=2)
text(0.7+0.02,0.2, "B", cex=1.2, pos=1, offset=0, font=2)

# C
par(mar = c(4.2,1,1.2,1), mgp=c(2, 0.5, 0))
plot(NULL, xlim=range(new.temp), ylim=c(0,4500),
     xlab="Standardized fledgling weight", ylab="",
     xaxt="n", yaxt="n")
hist(new.temp[,1], col=col.alpha("darkorange", 0.75), border="white", add=TRUE)
hist(new.temp[,2], col=col.alpha("darkslategrey", 0.75), border="white", add=TRUE)
hist(new.temp[,3], col=col.alpha("white", 0.5), border="black", add=TRUE)
axis(1, at=seq(-0.3,0.1, 0.1), labels=seq(-0.3,0.1, 0.1), lwd = 0, lwd.ticks = 1)
abline(v=post.out[3,2], lty=2)
text(range(new.temp)[1],4500, "C", cex=1.2, pos=1, offset=0, font=2)
dev.off()


# MODEL for fledging date (Eq. S7) of selection line and local (non-selection line) fledglings:

# run only for selection lines:
d <- For_Stats_Recruitment_all[complete.cases(For_Stats_Recruitment_all$Weight),]
d$Date_std <- standardize(d$D15CatchApril) # standardize April date of d15 measurements (relative timing)
d$Line <- as.factor(d$SelectionLineNumeric2) # make sure is factor/integer starting at 1
d$Year <- as.factor(d$BroodYear-2016) # make sure is factor/integer starting at 1
d <- d[d$Line!=3,]

# make input data:
dat <- list(
  R=d$Date_std,
  L=d$Line,
  Y=d$Year
)

# fit model
mFD.1c <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + bL[L],
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5),
    sigma ~ dhalfnorm(0,1),
    sigma_a ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(max_treedepth=15, adapt_delta=0.99))
save(mFD.1c, file="temp/mFD.1c.multilevel.RData") # save model output!

# get model summary:
(out <- precis(mFD.1c, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S13
write.table(out.0, "out/Bayes_FledgingDate_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

PSIS.FD <- PSIS(mFD.1c, pointwise=TRUE) # good!

# prepare data for prediction
d_pred <- list(L=rep(1:2,times=3), Y=rep(1:3, each=2))

# get posterior
p_post <- link(mFD.1c, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S14
write.table(post.out, "out/Bayes_FledgingDate_posterior_prediction.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mFD.1c)), row.names(precis(mFD.1c, depth=2))[grepl("z", row.names(precis(mFD.1c, depth=2)))], row.names(precis(mFD.1c, depth=2))[grepl("bL", row.names(precis(mFD.1c, depth=2)))])
par(mfrow = c(1,1))
# Fig. S27
pdf(file = "Plots/Bayes_FledglingDate_traceplots.pdf", width = 6, height = 8)
traceplot(mFD.1c, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# Fig. S28
pdf(file = "Plots/Bayes_FledglingDate_trankplot.pdf", width = 6, height = 8)
trankplot(mFD.1c, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# plot posterior prediction:
# prepare plot
post.out$x.lab <- c(1:3)
n <- as.numeric(table(dat$L))[1:2]

# Fig. S4
pdf(file = "Plots/Bayes_FledgingDate_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[-3,]
plot(NULL, xlim=c(0.7,2.3), ylim=c(-0.15,0.1), xlab="",
     ylab="standardized fledging date",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-0.12,0.08, 0.04), labels=seq(-0.12,0.08, 0.04), lwd = 0, lwd.ticks = 1)
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
text(0.7, 0.1, "A", cex=1.2, pos=1, offset=0, font=2)

plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-0.05,0.2),
     xlab="", ylab="difference in standardized\nfledging date (late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(-0.04,0.16, 0.04), labels=seq(-0.04,0.16, 0.04), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025, 0.2, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()


## MODEL for local recruitment probability (Eq. 4) of selection line and local (non-selection line) fledglings

# prepare variables for model
d <- For_Stats_Recruitment_all[complete.cases(For_Stats_Recruitment_all$Weight),] # remove fledglings without weight info
d$Date_std <- standardize(d$D15CatchApril) # standardize April date of d15 measurements (relative timing)
d$Weight_std <- standardize(d$Weight) # standardize weight
d$Line <- as.factor(d$SelectionLineNumeric2) # make sure is factor/integer starting at 1
d$Year <- as.factor(d$BroodYear-2016) # make sure is factor/integer starting at 1
d$Brood_factor <- as.factor(d$Brood) # make sure is factor/integer starting at 1

table(d$Brood_factor)[table(d$Brood_factor)==1] # 7 broods with one chick only. Might cause trouble..

# make input data:
dat <- list(
  R=d$Recruit,
  A=d$Date_std,
  W=d$Weight_std,
  L=d$Line,
  B=d$Brood_factor,
  Y=d$Year
)

# fit model
mR.1a <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + x[B]*sigma_g + bL[L] + bA*A + bW*W,
    z[Y] ~ dnorm(0, 1),
    x[B] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bA ~ dnorm(0, 1.5),
    bW ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[B]:g <<- x*sigma_g
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mR.1a, file="temp/mR.1a.multilevel.RData")

# get model summary:
(out <- precis(mR.1a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S11 
write.table(out.0, "out/Bayes_Recruit_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mR.1a)), row.names(precis(mR.1a, depth=2))[grepl("z", row.names(precis(mR.1a, depth=2)))], row.names(precis(mR.1a, depth=2))[grepl("bL", row.names(precis(mR.1a, depth=2)))])
par(mfrow = c(1,1))
# Fig. S23
pdf(file = "Plots/Bayes_Recruit_traceplots.pdf", width = 6, height = 8)
traceplot(mR.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# Fig. S24
pdf(file = "Plots/Bayes_Recruit_trankplot.pdf", width = 6, height = 8)
trankplot(mR.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior mean and ci for recruitment probability of selection line and local (non-selection line) fledglings
# prepare data for prediction; use unique combinations of year, selection line and brood location ID as in fitted data
dat_help.0 <- d[!duplicated(d[,19:21]), 19:21]
dat_help <- dat_help.0[order(dat_help.0$Line, dat_help.0$Year, dat_help.0$Brood_factor),]
d_pred <- list(L=dat_help$Line, Y=dat_help$Year, B=dat_help$Brood_factor, A=rep(0, times=nrow(dat_help)), W=rep(0, times=nrow(dat_help)))

# get posterior
p_post <- link(mR.1a, data=d_pred)

new.temp <- NULL
for(i in 1:3) {
  dat.temp <- p_post[,dat_help$Line==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
diff_e_w <- new.temp[,3]-new.temp[,1]
d_e_w_mu <- mean(diff_e_w); d_e_w_ci <- PI(diff_e_w)
diff_l_w <- new.temp[,3]-new.temp[,2]
d_l_w_mu <- mean(diff_l_w); d_l_w_ci <- PI(diff_l_w)

post.out <- data.frame(line=c(1,2,3,"diff_e_l", "diff_e_w", "diff_l_w"), mu=c(p_mu, d_mu, d_e_w_mu, d_l_w_mu), ci.5.5=c(p_ci[1,],d_ci[1],d_e_w_ci[1],d_l_w_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2],d_e_w_ci[2],d_l_w_ci[2]))
# Table S12
write.table(post.out, "out/Bayes_Recruit_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction:

# prepare plot
# 1. selection line fledglings (A)
post.out$x.lab <- c(1:2,0,3,0,0)
n <- as.numeric(table(d$SelectionLineNumeric[d$SelectionLineNumeric!=0]))

# get mean and sd of April day 15 days post hatching (B)
d2 <- d[d$SelectionLineNumeric!=0,]
April_mean <- tapply(d2$Date_std, d2$SelectionLineNumeric, mean)
April_ci <- tapply(d2$Date_std, d2$SelectionLineNumeric, PI)
April_day_lines <- data.frame(mu_April=April_mean, ci.5.5_April=c(April_ci[[1]][1],April_ci[[2]][1]), ci.94.5_April=c(April_ci[[1]][2],April_ci[[2]][2]))
# combine in data frame with recruitment probability
for_common_plot <- cbind(post.out[-(3:6),-5], April_day_lines)

# 1. local (non-selection line) fledglings (B)
dat_help.0 <- d[!duplicated(d[,19:21]), 19:21]; dat_help.0 <- dat_help.0[dat_help.0$Line==3,]
dat_help <- dat_help.0[order(dat_help.0$Year, dat_help.0$Brood_factor),]
range <- range(d$Date_std)
xseq_base <- seq(range[1]-0.05, range[2]+0.05, length.out=50)

d_pred <- list(Y=rep(dat_help$Year, times=50), B=rep(dat_help$Brood_factor, times=50), A=rep(xseq_base, each=nrow(dat_help)), W=rep(0, nrow(dat_help)*50), L=rep(3, nrow(dat_help)*50))
p_post <- link(mR.1a, data=d_pred) # Quite memory intensive. Depending on your local machine, better run on HPC.
save(p_post, file="temp/mR.1a.post.RData")

d_plot <- d[d$Line==3,]
range <- range(d_plot$Date_std)
xseq <- seq(range[1], range[2], length.out=50)

new.temp <- NULL
for(i in 1:50){
  dat.temp <- p_post[,(nrow(dat_help)*i-nrow(dat_help)+1):(nrow(dat_help)*i)]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)

# make plot
pdf(file = "Plots/Bayes_Recruit.pdf", width = 3.5, height = 3.5)
par(mfrow = c(1,1))
par(mar = c(4.2,3.5,1.2,2.5), mgp=c(2.5, 0.5, 0))

plot(NULL, xaxt="n", yaxt="n", xlim=range(xseq_base), ylim=c(0,0.2),
     xlab="standardized April date\nat 15 days post hatching", ylab="", cex.main=1)
axis(1, at=seq(-1,4,1), labels=seq(-1,4,1), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,0.2,0.1), labels=seq(0,0.2,0.1), lwd = 0, lwd.ticks = 1)
axis(4, at=c(0,0.2), labels=c(0,1), lwd = 0, lwd.ticks = 1, col.axis="grey")
mtext("local recruitment probability", side=2, line=1.75, cex=1, col="black")
mtext("recruited", side=4, line=0.5, cex=1, col="grey")
points(d_plot$Date_std, d_plot$Recruit/5, pch=1, col="grey")
lines(xseq, p_mu, lwd=2)
shade(p_ci, xseq, col=col.alpha("black", 0.1))

# add selection line data
arrows(x0=for_common_plot$mu_April, y0=for_common_plot$ci.5.5, x1=for_common_plot$mu_April, y1=for_common_plot$ci.94.5, code=3, angle=90, length=0.1, col=ifelse(for_common_plot$line==1, "darkorange", "darkslategrey"), lwd=2)
arrows(x0=for_common_plot$ci.5.5_April, y0=for_common_plot$mu, x1=for_common_plot$ci.94.5_April, y1=for_common_plot$mu, code=3, angle=90, length=0.1, col=ifelse(for_common_plot$line==1, "darkorange", "darkslategrey"), lwd=2)
points(for_common_plot$mu_April, for_common_plot$mu, pch=19, cex=1.5, col=ifelse(for_common_plot$line==1, "darkorange", "darkslategrey"))
dev.off()

# Fig. 3
pdf(file = "Plots/Main_Bayes_Recruit_v2.pdf", width = 6.5, height = 3.5)
layout(t(matrix(c(1,2))), width=c(1.25,1.5))

# A
par(mar = c(4,4,1.2,1.2), mgp=c(2, 0.5, 0))
plot(NULL, xlim=c(0.7,2.3), ylim=c(0.05,0.115), xlab="",
     ylab="Local recruitment probability",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0.05,0.11, 0.02), labels=seq(0.05,0.11, 0.02), lwd = 0, lwd.ticks = 1)
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
text(0.7+0.04, 0.115, "A", cex=1.2, pos=1, offset=0, font=2)

# B
par(mar = c(4,3,1.2,1.2), mgp=c(2, 0.5, 0))
plot(NULL, xaxt="n", yaxt="n", xlim=range(xseq_base), ylim=c(0,0.2),
     xlab="Standardized April date (d15)", ylab="Local recruitment probability", cex.main=1)
axis(1, at=seq(-1,4,1), labels=seq(-1,4,1), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,0.2,0.1), labels=seq(0,0.2,0.1), lwd = 0, lwd.ticks = 1)
shade(p_ci, xseq, col=col.alpha("black", 0.1))
lines(xseq, p_mu, lwd=1)
text(range(xseq_base)[1]+0.05,0.2, "B", cex=1.2, pos=1, offset=0, font=2)

# add selection line data
arrows(x0=for_common_plot$mu_April, y0=for_common_plot$ci.5.5, x1=for_common_plot$mu_April, y1=for_common_plot$ci.94.5, code=3, angle=90, length=0.1, col=ifelse(for_common_plot$line==1, "darkorange", "darkslategrey"), lwd=1)
arrows(x0=for_common_plot$ci.5.5_April, y0=for_common_plot$mu, x1=for_common_plot$ci.94.5_April, y1=for_common_plot$mu, code=3, angle=90, length=0.1, col=ifelse(for_common_plot$line==1, "darkorange", "darkslategrey"), lwd=1)
points(for_common_plot$mu_April, for_common_plot$mu, pch=19, cex=1.5, col=ifelse(for_common_plot$line==1, "darkorange", "darkslategrey"))
dev.off()


## MODEL for local recruitment probability (Eq. 4 but with an additional sex-effect) of selection line fledglings

# fix brood location IDs; make it range from 1 to n broods rather than character as used now
FINAL_Dat_Recruitment_all <- FINAL_Dat_Recruitment[FINAL_Dat_Recruitment$SelectionLineNumeric!=0 & !is.na(FINAL_Dat_Recruitment$Sexe),]
temp <- FINAL_Dat_Recruitment_all[,c(4,11)]
temp.0 <- temp[!duplicated(temp),]
temp.0$Brood <- 1:nrow(temp.0)
For_Stats_Recruitment_SL <- merge(FINAL_Dat_Recruitment_all, temp.0, by=c("BroodYear", "D15LocationAtHV"))

# prepare variables for model
d <- For_Stats_Recruitment_SL[complete.cases(For_Stats_Recruitment_SL$Weight),] # remove fledglings without weight information
d$Date_std <- standardize(d$D15CatchApril) # standardize April date of d15 measurements (relative timing)
d$Weight_std <- standardize(d$Weight) # standardize weight
d$Line <- as.factor(d$SelectionLineNumeric) # make sure is factor/integer starting at 1
d$Year <- as.factor(d$BroodYear-2016) # make sure is factor/integer starting at 1
d$Brood_factor <- as.factor(d$Brood) # make sure is factor/integer starting at 1
d$Sex <- as.factor(d$Sexe)

table(d$Brood_factor)[table(d$Brood_factor)==1] # two broods with one chick only. Might cause trouble..

# with weight and date
# make input data:
dat <- list(
  R=d$Recruit,
  A=d$Date_std,
  W=d$Weight_std,
  L=d$Line,
  B=d$Brood_factor,
  Y=d$Year,
  S=d$Sex
)

# fit model
mR.2 <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + x[B]*sigma_g + bL[L] + bS[S] + bA*A + bW*W,
    z[Y] ~ dnorm(0, 1),
    x[B] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bS[S] ~ dnorm(0,1.5),
    bA ~ dnorm(0, 1.5),
    bW ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[B]:g <<- x*sigma_g
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mR.2, file="temp/mR.2.multilevel.RData")

# get model summary:
(out <- precis(mR.2, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# table S17
write.table(out.0, "out/Bayes_Recruit_SL_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mR.2)), row.names(precis(mR.2, depth=2))[grepl("z", row.names(precis(mR.2, depth=2)))], row.names(precis(mR.2, depth=2))[grepl("bL", row.names(precis(mR.2, depth=2)))], row.names(precis(mR.2, depth=2))[grepl("bS", row.names(precis(mR.2, depth=2)))])
par(mfrow = c(1,1))
# Fig. S25
pdf(file = "Plots/Bayes_Recruit_SL_traceplots.pdf", width = 6, height = 8)
traceplot(mR.2, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# Fig. S26
pdf(file = "Plots/Bayes_Recruit_SL_trankplot.pdf", width = 6, height = 8)
trankplot(mR.2, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior mean and ci for recruitment probability of female and male selection line fledglings
dat_help.0 <- d[!duplicated(d[,17:20]), 17:20]
dat_help <- dat_help.0[order(dat_help.0$Line, dat_help.0$Sex, dat_help.0$Year, dat_help.0$Brood_factor),]
d_pred <- list(L=dat_help$Line,S=dat_help$Sex, Y=dat_help$Year, B=dat_help$Brood_factor, A=rep(0, times=nrow(dat_help)), W=rep(0, times=nrow(dat_help)))

# get posterior
p_post <- link(mR.2, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,dat_help$Line==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S18 (1)
write.table(post.out, "out/Bayes_Recruit_SL_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# get mean for female and male , respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,dat_help$Sex==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
post.out <- data.frame(line=c("f","m","diff_f_m"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S18 (2)
write.table(post.out, "out/Bayes_Recruit_SL_sex_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# test for differenc ebetween sexes for early and late selection line birds:

# early
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,dat_help$Line==1 & dat_help$Sex==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
post.out.e <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))


# late
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,dat_help$Line==2 & dat_help$Sex==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out.l <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))

temp <- rbind(post.out.e[3,], post.out.l[3,])
temp$line <- c("early", "late")

# Table S18 (3)
write.table(temp, "out/Bayes_Recruit_SL_sex_int_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# get posterior distribution 
# early: 
dat.temp <- p_post[,dat_help$Line==1 & dat_help$Sex==1]
e_f <- apply(dat.temp, 1, mean)

# late:
dat.temp <- p_post[,dat_help$Line==2 & dat_help$Sex==1]
l_f <- apply(dat.temp, 1, mean)

Keep.recruitment.prob <- list(
  early=e_f, late=l_f
)


### MODEL comparison (in line with model presented in Eq. 4): local recruitment probability (selection line and local (non-selection) fledglings); weight and/or date included?

# Models: with weight and date (a); with weight and without date (b); with date and without weight (c)

# prepare data for model:
For_Stats_Recruitment_all$help <- ifelse(For_Stats_Recruitment_all$SelectionLineNumeric==0,3,0)
For_Stats_Recruitment_all$SelectionLineNumeric2 <- For_Stats_Recruitment_all$SelectionLineNumeric+For_Stats_Recruitment_all$help
# get all variables as needed in model input
d <- For_Stats_Recruitment_all[complete.cases(For_Stats_Recruitment_all$Weight),]
d$Date_std <- standardize(d$D15CatchApril) # standardize April date of d15 measurements (relative timing)
d$Weight_std <- standardize(d$Weight) # standardize weight
d$Line <- as.factor(d$SelectionLineNumeric2) # make sure is factor/integer starting at 1
d$Year <- as.factor(d$BroodYear-2016) # make sure is factor/integer starting at 1
d$Brood_factor <- as.factor(d$Brood) # make sure is factor/integer starting at 1

table(d$Brood_factor)[table(d$Brood_factor)==1] # 7 broods with one chick only. Might cause trouble..

# (a) with weight and date
# make input data:
dat <- list(
  R=d$Recruit,
  A=d$Date_std,
  W=d$Weight_std,
  L=d$Line,
  B=d$Brood_factor,
  Y=d$Year
)

# fit model
mR.1a <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + x[B]*sigma_g + bL[L] + bA*A + bW*W,
    z[Y] ~ dnorm(0, 1),
    x[B] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bA ~ dnorm(0, 1.5),
    bW ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[B]:g <<- x*sigma_g
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mR.1a, file="temp/mR.1a.multilevel.RData")

# (b) with weight only
# make input data:
dat <- list(
  R=d$Recruit,
  W=d$Weight_std,
  L=d$Line,
  B=d$Brood_factor,
  Y=d$Year
)

# fit model
mR.1b <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + x[B]*sigma_g + bL[L] + bW*W,
    z[Y] ~ dnorm(0, 1),
    x[B] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bW ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[B]:g <<- x*sigma_g
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mR.1b, file="temp/mR.1b.multilevel.RData")

# (c) with date only
dat <- list(
  R=d$Recruit,
  A=d$Date_std,
  L=d$Line,
  B=d$Brood_factor,
  Y=d$Year
)

# fit model
mR.1c <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + x[B]*sigma_g + bL[L] + bA*A,
    z[Y] ~ dnorm(0, 1),
    x[B] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bA ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[B]:g <<- x*sigma_g
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mR.1c, file="temp/mR.1c.multilevel.RData")

(model_comparison <- compare(mR.1a, mR.1b, mR.1c, func=PSIS))
model_comparison$Model <- row.names(model_comparison)
# Table S36
write.table(model_comparison, "temp/Bayes_Recruit_model_comparison.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")



###   Lay dates of female selection line recruits and local (non-selection line) female recruits
### ---------------------------------------------------------------------------

### Prepare data

# some checks 
(keep <- unique(Recruits$SelectionLine)[-c(3,6)])
(Years <- sort(unique(Recruits[Recruits$SelectionLine %in% keep,]$YearOfBreeding)))

# only great tits at HV?
unique(Fitness$Species)
unique(Fitness$Area)

# all ERC recruits in fitness table?
ERC <- unique(Recruits$RingNumber)
Recruits[!Recruits$RingNumber %in% ERC,]
FitTbl <- c(unique(Fitness$Mother), unique(Fitness$Father))

FitTbl[FitTbl %in% ERC]
ERC_birds_missing <- ERC[!ERC %in% FitTbl]

unique(Recruits[!Recruits$RingNumber %in% FitTbl,]$RingNumber)==ERC[!ERC %in% FitTbl]
Recruits[!Recruits$RingNumber %in% FitTbl,] ## missing birds are HK recruits !!!
Recruits_noHK <- Recruits[!grepl("HK", Recruits$UserPlaceName),]

# 1. get females selection line recruits
temp_sel <- Recruits_noHK[Recruits_noHK$SelectionLine %in% keep,]
# remove HK birds
temp_sel$SelectionLineNumeric <- ifelse(grepl("E", temp_sel$SelectionLine), 1, 2)
# select females and 'first clutch'-broods
dat_ERC_females <- temp_sel[temp_sel$Sexe==1 & temp_sel$BroodType==0,]
ERC_Females <- unique(dat_ERC_females[,c(2,19)])

# 2. get mean and standard deviation for standardization (from local (non-selection line) female recruits)
  # only keep first clutches from 2018-2020; use all non-selection line females (local and immigrants) 
ForLayDates <- Fitness[Fitness$BroodType==0 & Fitness$YearOfBreeding>=2018 & Fitness$YearOfBreeding<=2020,]
ForLayDates_females <- unique(ForLayDates$RingNumberFemale); ForLayDates_females <- ForLayDates_females[ForLayDates_females!="0000000000"]

# later when number of total fledglings is related to lay date, the first record of lay date is used for all females that were breeding in 2018-2020 
# for this, we need to get all years in which females with recorded lay date in 2018-2020 have a lay date in an earlier year and also calculate mean and sd for all non-selection line females in those years
# for the analysis of lay dates (i.e. presented here) onle mean and sd from 2018-2020 are used.
Years_with_lay_date_records <- unique(Fitness[Fitness$RingNumberFemale %in% ForLayDates_females,]$YearOfBreeding) # include 2013-2020
ForLayDates.1.0 <- Fitness[Fitness$BroodType==0 & Fitness$YearOfBreeding>=2013 & Fitness$YearOfBreeding<=2020,]

# remove selection line females (must not be included in calculating mean and standard deviation for standardization)
ForLayDates.1 <- ForLayDates.1.0[!ForLayDates.1.0$RingNumberFemale %in% ERC_Females$RingNumber,]
mean <- tapply(ForLayDates.1$LayDateApril, ForLayDates.1$YearOfBreeding, mean)
sd <- tapply(ForLayDates.1$LayDateApril, ForLayDates.1$YearOfBreeding, sd)
n <- tapply(ForLayDates.1$LayDateApril, ForLayDates.1$YearOfBreeding, length)
Summary_lay <- data.frame(YearOfBreeding=names(mean), Mean=mean, Sd=sd, n_BreedingFemales=n); row.names(Summary_lay) <- 1:nrow(Summary_lay)
keep_LayDates_for_std <- Summary_lay

# standardization for local (non-selection line) females
wild_ForLayDates <- ForLayDates.1[ForLayDates.1$YearOfBreeding>=2018, c(4,7,12)]
wild_ForLayDates$SelectionLineNumeric <- rep(0, nrow(wild_ForLayDates))
wild_ForLayDates.0 <- merge(wild_ForLayDates, keep_LayDates_for_std, by="YearOfBreeding")
wild_ForLayDates.0$LayDateApril_std_zscore <- (wild_ForLayDates.0$LayDateApril-wild_ForLayDates.0$Mean)/wild_ForLayDates.0$Sd
wild_ForLayDates.0$LayDateApril_std_meancentered <- wild_ForLayDates.0$LayDateApril-wild_ForLayDates.0$Mean

# standardization for selection line females
# get lay date of first year of breeding for selection line recruits; sort selection line females by ring number and year of breeding
temp <- dat_ERC_females[order(dat_ERC_females$RingNumber, dat_ERC_females$YearOfBreeding),]
temp.0 <- temp[!duplicated(temp$RingNumber), c(6, 2, 11, 19)]; names(temp.0)[2] <- "RingNumberFemale"
ERC_wild_ForLayDates.0 <- merge(temp.0, keep_LayDates_for_std, by="YearOfBreeding")
ERC_wild_ForLayDates.0$LayDateApril_std_zscore <- (ERC_wild_ForLayDates.0$LayDateApril-ERC_wild_ForLayDates.0$Mean)/ERC_wild_ForLayDates.0$Sd
ERC_wild_ForLayDates.0$LayDateApril_std_meancentered <- ERC_wild_ForLayDates.0$LayDateApril-ERC_wild_ForLayDates.0$Mean

# combine data
FINAL_Dat_LayDates <- rbind(ERC_wild_ForLayDates.0, wild_ForLayDates.0)
save(FINAL_Dat_LayDates, file="temp/Dat_lay_dates.RData")


### Statistical analysis

## MODEL for lay dates (Eq. 1) of female selection line recruits and local (non-selection line) female recruits

# prepare data for model:
d.0 <- FINAL_Dat_LayDates
d.0$help <- ifelse(d.0$SelectionLineNumeric==0,3,0)
d.0$SelectionLineNumeric2 <- d.0$SelectionLineNumeric+d.0$help
d <- d.0
d$R <- d$LayDateApril_std_zscore
d$L <- as.factor(d$SelectionLineNumeric2)

# input data:
dat <- list(
  R=d$R,
  L=d$L
)

# run model:
mL.1c <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a+b[L],
    a ~ dnorm(0,1.5),
    b[L] ~ dnorm(0,1.5),
    sigma ~ dhalfnorm(0,1)
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mL.1c, file="temp/mL.1c.RData")

dat <- list(L=1:3)
p_post <- link(mL.1c, data=dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)
diff <- p_post[,2]-p_post[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
diff_e_w <- p_post[,3]-p_post[,1]
d_e_w_mu <- mean(diff_e_w); d_e_w_ci <- PI(diff_e_w)
diff_l_w <- p_post[,3]-p_post[,2]
d_l_w_mu <- mean(diff_l_w); d_l_w_ci <- PI(diff_l_w)
post.out <- data.frame(line=c(1,2,3,"diff_e_l", "diff_e_w", "diff_l_w"), mu=c(p_mu, d_mu, d_e_w_mu, d_l_w_mu), ci.5.5=c(p_ci[1,],d_ci[1],d_e_w_ci[1],d_l_w_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2],d_e_w_ci[2],d_l_w_ci[2]))

# plot posterior prediction
n_line <- data.frame(table(d$L)[-3]); names(n_line)[2] <- "n"; n_line$lab <- ifelse(n_line$Var1==1, "early", "late")
d <- FINAL_Dat_LayDates[FINAL_Dat_LayDates$SelectionLineNumeric!=0,]
d$R <- d$LayDateApril_std_zscore
d$L <- as.factor(d$SelectionLineNumeric)

# Fig. S12
pdf(file = "Plots/Bayes_Lay_date.c.pdf", width = 8, height = 5)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0))

dat_out <- post.out[1:2,]; dat_out$line <- as.numeric(dat_out$line)
plot(NULL, xlim=c(0.7,2.3), ylim=c(-1.2,2.7),
     xlab="", ylab="standardized lay date",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(n_line$lab, " (n=", n_line$n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-1,2.6, 0.6), labels=seq(-1,2.6, 0.6), lwd = 0, lwd.ticks = 1)
arrows(x0=dat_out$line, y0=dat_out$ci.5.5, x1=dat_out$line, y1=dat_out$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(dat_out$line, dat_out$mu, pch=19, cex=3, col=ifelse(dat_out$line==1, "darkorange", "darkslategrey"))
points(jitter(as.numeric(d$L),0.5), d$LayDateApril_std_zscore, pch=1, cex=1.5, col=ifelse(d$L==1, "darkorange", "darkslategrey"))
text(0.7,2.7, "A", cex=1.2, pos=1, offset=0, font=2)

dat_out_diff <- post.out[4,]; dat_out$line <- as.numeric(dat_out$line)
plot(NULL, xlim=c(0.7,1.3), ylim=c(-0.1,1.7),
     xlab="", ylab="difference in standardized lay date (late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(0,1.6, 0.4), labels=seq(0,1.6, 0.4), lwd = 0, lwd.ticks = 1)
arrows(x0=1, y0=dat_out_diff$ci.5.5, x1=1, y1=dat_out_diff$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(1, dat_out_diff$mu, pch=19, cex=3, col="black")
abline(h=0, lty=2)
text(0.7+0.02,1.7, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()

## check whether female is a potential outlier
# prepare data for model:
d <- FINAL_Dat_LayDates[FINAL_Dat_LayDates$SelectionLineNumeric!=0,]
d$R <- d$LayDateApril_std_zscore
d$L <- as.factor(d$SelectionLineNumeric)

dat <- list(
  R=d$R,
  L=d$L
)

# run model:
mL.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a+b[L],
    a ~ dnorm(0,1.5),
    b[L] ~ dnorm(0,1.5),
    sigma ~ dhalfnorm(0,1)
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mL.1a, file="temp/mL.1a.RData")

(PSIS.a <- PSIS(mL.1a, pointwise=TRUE))
PSIS.a$RingNumberFemale <- d$RingNumberFemale
(WAIC.a <- WAIC(mL.1a, pointwise=TRUE)) # same picture as with PSIS

# plot PSIS outcome
# Fig. S11
pdf(file = "Plots/Bayes_Lay_date_PSIS.a.pdf", width = 3.5, height = 3.5)
par(mfrow = c(1,1))
par(mar = c(3.5,3.5,1.2,1.2), mgp=c(2.0, 0.5, 0))
plot(NULL, xlim=c(0,0.7), ylim=c(0,1.8),
     xlab="Pareto k", ylab="PSIS penalty", cex.main=1)
points(PSIS.a$k, PSIS.a$penalty, pch=ifelse(PSIS.a$RingNumberFemale!="AU...91626",1,19))
abline(v=0.5, lty=2)
dev.off()


## repeat analysis but remove outlier:
d.0 <- FINAL_Dat_LayDates[FINAL_Dat_LayDates$RingNumberFemale!="AU...91626",]
d.0$help <- ifelse(d.0$SelectionLineNumeric==0,3,0)
d.0$SelectionLineNumeric2 <- d.0$SelectionLineNumeric+d.0$help
d <- d.0
d$R <- d$LayDateApril_std_zscore
d$L <- as.factor(d$SelectionLineNumeric2)

dat <- list(
  R=d$R,
  L=d$L
)

# run model:
mL.1d <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a+b[L],
    a ~ dnorm(0,1.5),
    b[L] ~ dnorm(0,1.5),
    sigma ~ dhalfnorm(0,1)
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mL.1d, file="temp/mL.1d.RData")

(out <- precis(mL.1d, depth=2))
pars <- row.names(out)
out.0 <- data.frame(parameter=pars, out)

# Table S1
write.table(out.0, "out/Bayes_Lay_date_precis.d.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

pars <- row.names(precis(mL.1d, depth=2))
par(mfrow = c(1,1))
# Fig. S13
pdf(file = "Plots/Bayes_Lay_date_Line_traceplots.pdf", width = 6, height = 6)
traceplot(mL.1d, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# Fig. S14
pdf(file = "Plots/Bayes_Lay_date_Line_trankplot.pdf", width = 6, height = 6)
trankplot(mL.1d, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# get posterior predictions
dat <- list(L=1:3)
p_post <- link(mL.1d, data=dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)
diff <- p_post[,2]-p_post[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
diff_e_w <- p_post[,3]-p_post[,1]
d_e_w_mu <- mean(diff_e_w); d_e_w_ci <- PI(diff_e_w)
diff_l_w <- p_post[,3]-p_post[,2]
d_l_w_mu <- mean(diff_l_w); d_l_w_ci <- PI(diff_l_w)
post.out <- data.frame(line=c(dat$L,"diff_e_l", "diff_e_w", "diff_l_w"), mu=c(p_mu, d_mu, d_e_w_mu, d_l_w_mu), ci.5.5=c(p_ci[1,],d_ci[1],d_e_w_ci[1],d_l_w_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2],d_e_w_ci[2],d_l_w_ci[2]))
# Table S2
write.table(post.out, "out/Bayes_Lay_date_posterior_prediction.d.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior predictions:
# prepare plot
n_line <- data.frame(table(d$L)[-3]); names(n_line)[2] <- "n"; n_line$lab <- ifelse(n_line$Var1==1, "early", "late")
d <- FINAL_Dat_LayDates[FINAL_Dat_LayDates$SelectionLineNumeric!=0,]
d$R <- d$LayDateApril_std_zscore
d$L <- as.factor(d$SelectionLineNumeric)
d$Fem <- d$RingNumberFemale

# Fig. 1
pdf(file = "Plots/Main_Bayes_Lay_date_v2.pdf", width = 6, height = 3)
layout(t(matrix(c(1,2))), width=c(1.25, 1.5))

# A
par(mar = c(4,4,1.2,1), mgp=c(2, 0.5, 0))
dat_out <- post.out[1:2,]; dat_out$line <- as.numeric(dat_out$line)
plot(NULL, xlim=c(0.7,2.3), ylim=c(-1.2,2.7),
     xlab="", ylab="Standardized lay date",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(n_line$lab, " (n=", n_line$n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-1,2.6, 0.6), labels=seq(-1,2.6, 0.6), lwd = 0, lwd.ticks = 1)
arrows(x0=dat_out$line, y0=dat_out$ci.5.5, x1=dat_out$line, y1=dat_out$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(dat_out$line, dat_out$mu, pch=19, cex=3, col=ifelse(dat_out$line==1, "darkorange", "darkslategrey"))
points(jitter(as.numeric(d$L),0.5), d$LayDateApril_std_zscore, pch=ifelse(d$Fem=="AU...91626",13,1), cex=1.5, col=ifelse(d$L==1, "darkorange", "darkslategrey"))
text(0.7+0.05,2.7, "A", cex=1.2, pos=1, offset=0, font=2)

# B
par(mar = c(4.2,1,1.2,1), mgp=c(2, 0.5, 0))
plot(NULL, xlim=range(p_post)+0.05, ylim=c(0,4500),
     xlab="Standardized lay date", ylab="",
     xaxt="n", yaxt="n")
hist(p_post[,1], col=col.alpha("darkorange", 0.75), border="white", add=TRUE)
hist(p_post[,2], col=col.alpha("darkslategrey", 0.75), border="white", add=TRUE)
hist(p_post[,3], col="gray30", border="gray30", add=TRUE)
axis(1, at=seq(-1.5,2, 0.5), labels=c("",-1,"",0,"",1,"",2), lwd = 0, lwd.ticks = 1)
abline(v=post.out[3,2], lty=2)
text(range(p_post)[1]+0.06,4500, "B", cex=1.2, pos=1, offset=0, font=2)

dev.off()


### Life-time number of fledglings produced by female selection line recruits and local (non-selection line) females
### ---------------------------------------------------------------------------

### Prepare data

# 1. get mean and standard deviation needed for standardization 
# only count fledglings from broods between 2018-2020 ; use all females (local and immigrants)
# remove females without ring number 
ForFledglings <- Fitness[Fitness$RingNumberFemale!="0000000000" & Fitness$YearOfBreeding>=2018 & Fitness$Include==1, c(4,7,11,12,19,23)]
# remove selection line females
ForFledglings.1 <- ForFledglings[!ForFledglings$RingNumberFemale %in% ERC_Females$RingNumber,]

# for the standardization we first calculate the total number of fledglings a female produced within a year and use the sum to do a within year standardization
# later; total number of fledglings within a year is standardized within year and the total number of within-year standradized fledglings across years is calculated

# so, (1) get the sum within year
ForFledglings.2 <- split(ForFledglings.1[,c(1:3,5)], ForFledglings.1$YearOfBreeding)

# check whether females have more than one brood within a year: YES!
lapply(ForFledglings.2, dim)
lapply(ForFledglings.2, function(x) length(unique(x$RingNumberFemale)))

# get sum of fledglings for each female in within a year (there mus be a MUUUCH easier way to do this..)
sum <- lapply(ForFledglings.2, function(x) as.data.frame(tapply(x$NoFledged, x$RingNumberFemale, sum)))
sum.0 <- lapply(sum, function(x) data.frame(RingNumberFemale=row.names(x), SumFledged=x[,1]))
lapply(sum.0, dim)
sum.1 <- data.frame(YearOfBreeding=c(rep(names(sum.0)[1], as.numeric(lapply(sum.0, nrow)[1])),
                                     rep(names(sum.0)[2], as.numeric(lapply(sum.0, nrow)[2])),
                                     rep(names(sum.0)[3], as.numeric(lapply(sum.0, nrow)[3]))),
                    RingNumberFemale=unlist(sapply(sum.0, "[", 1)),
                    SumFledged=unlist(sapply(sum.0, "[", 2))); row.names(sum.1) <- 1:nrow(sum.1)
nrow(sum.1)==sum(as.numeric(lapply(sum.0, nrow)))

mean <- tapply(sum.1$SumFledged, sum.1$YearOfBreeding, mean)
sd <- tapply(sum.1$SumFledged, sum.1$YearOfBreeding, sd)
n <- tapply(sum.1$SumFledged, sum.1$YearOfBreeding, length)
Summary_fledge <- data.frame(YearOfBreeding=names(mean), Mean=mean, Sd=sd, n_BreedingFemales=n); row.names(Summary_fledge) <- 1:nrow(Summary_fledge)
keep_NoFledged_for_std <- Summary_fledge

# (2) standardize data

# 1. selection line recruits:
table(temp_sel$YearBorn) # data from lay date analysis
table(temp_sel$YearOfBreeding)
table(temp_sel$BroodType)

# exclude females with biased brood success (make sure to exclude all broods of females that have any biased brood)
out_ERC <- temp_sel[temp_sel$Include==0,]$RingNumber 
temp_sel[temp_sel$RingNumber %in% out_ERC, c(2,3,5,6,8,14,19)]

# the two males only bred once, exclude from data (mean of wild population cannot be used to impute number of fledglings as we expect the selection lines to differ from the wild population)
temp_sel.0 <- temp_sel[!temp_sel$RingNumber %in% out_ERC & temp_sel$YearOfBreeding>2017, c(2,3,5,6,8,14,19)]

# calculate total number of fledglings produced per year and bird
# split data by year
temp_sel.2 <- split(temp_sel.0[,c(1,4,6)], temp_sel.0$YearOfBreeding)
# are there birds with >1 brood per year? Yes, in 2019.
lapply(temp_sel.2, dim)
lapply(temp_sel.2, function(x) length(unique(x$RingNumber))) 

# get sum of fledglings for each bird in each year
sum <- lapply(temp_sel.2, function(x) as.data.frame(tapply(x$NoFledged, x$RingNumber, sum)))
sum.0 <- lapply(sum, function(x) data.frame(RingNumber=row.names(x), SumFledged=x[,1]))
lapply(sum.0, dim)
sum.1 <- data.frame(YearOfBreeding=c(rep(names(sum.0)[1], as.numeric(lapply(sum.0, nrow)[1])),
                                     rep(names(sum.0)[2], as.numeric(lapply(sum.0, nrow)[2])),
                                     rep(names(sum.0)[3], as.numeric(lapply(sum.0, nrow)[3]))),
                    RingNumber=unlist(sapply(sum.0, "[", 1)),
                    NoFledged_SumPerYear=unlist(sapply(sum.0, "[", 2))); row.names(sum.1) <- 1:nrow(sum.1)
nrow(sum.1)==sum(as.numeric(lapply(sum.0, nrow)))

# standardize total number of fledglings within a year
ERC_NoFledged.0 <- merge(sum.1, keep_NoFledged_for_std[,-4], by="YearOfBreeding")
ERC_NoFledged.0$NoFledged_SumPerYear_zscore <- (ERC_NoFledged.0$NoFledged_SumPerYear-ERC_NoFledged.0$Mean)/ERC_NoFledged.0$Sd
ERC_NoFledged.0$NoFledged_SumPerYear_meancentered <- ERC_NoFledged.0$NoFledged_SumPerYear-ERC_NoFledged.0$Mean

# get total number of fledglings summed across years
sum_zscore <- tapply(ERC_NoFledged.0$NoFledged_SumPerYear_zscore, ERC_NoFledged.0$RingNumber, sum)
sum_meancentered <- tapply(ERC_NoFledged.0$NoFledged_SumPerYear_meancentered, ERC_NoFledged.0$RingNumber, sum)

ERC_NoFledged.1 <- data.frame(SumNoFledged_zscore=sum_zscore, SumNoFledged_meancentered=sum_meancentered); ERC_NoFledged.1$RingNumber <- row.names(ERC_NoFledged.1); row.names(ERC_NoFledged.1) <- 1:nrow(ERC_NoFledged.1)
length(unique(ERC_NoFledged.1$RingNumber))==nrow(ERC_NoFledged.1)

# combine with bird-specific data
temp_sel.1 <- temp_sel.0[!duplicated(temp_sel.0$RingNumber), c(1:3,7)]
ERC_NoFledged <- merge(temp_sel.1, ERC_NoFledged.1, by="RingNumber")

# 2. local (non-selection line recruits)
table(Fitness$YearOfBreeding)

# combine males and females in one row and remove unknown individuals
temp_fitness <- gather(Fitness, names(Fitness)[c(7,9)], key = "Sex", value = "RingNumber"); temp_fitness$Sexe <- ifelse(temp_fitness$Sex=="RingNumberFemale",1,2)
temp_fitness$YearBorn <- rep(NA, nrow(temp_fitness))
temp_fitness$SelectionLineNumeric <- rep(0, nrow(temp_fitness))
temp_fitness.0 <- temp_fitness[temp_fitness$RingNumber!="0000000000" & temp_fitness$YearOfBreeding>=2018, c(23:25,4,9,17,26,21)]
temp_fitness.0$NoFledged[temp_fitness.0$Include==0] <- NA

# remove selection line birds; difference in two birds are the two selection line males that are excluded from analysis (see above)
# note that for local (non-selection line) individuals, the within-year total number of fledglings is imputed when data is biased (see below)
ERCs_recruits_in_wild <- intersect(temp_fitness.0$RingNumber, temp_sel$RingNumber)
temp_fitness.1 <- temp_fitness.0[!temp_fitness.0$RingNumber %in% ERCs_recruits_in_wild,]

# calculate within-year total number of fledglings produced 
# split data by year
temp_fitness.2 <- split(temp_fitness.1[,c(1,4,6)], temp_fitness.1$YearOfBreeding)
# are there birds with >1 brood per year? YES!
lapply(temp_fitness.2, dim)
lapply(temp_fitness.2, function(x) length(unique(x$RingNumber))) 
# ---> all years have birds breeding >1

# get within-year sum of fledglings for each bird
sum <- lapply(temp_fitness.2, function(x) as.data.frame(tapply(x$NoFledged, x$RingNumber, sum)))
sum.0 <- lapply(sum, function(x) data.frame(RingNumber=row.names(x), SumFledged=x[,1]))
lapply(sum.0, dim)
sum.1 <- data.frame(YearOfBreeding=c(rep(names(sum.0)[1], as.numeric(lapply(sum.0, nrow)[1])),
                                     rep(names(sum.0)[2], as.numeric(lapply(sum.0, nrow)[2])),
                                     rep(names(sum.0)[3], as.numeric(lapply(sum.0, nrow)[3]))),
                    RingNumber=unlist(sapply(sum.0, "[", 1)),
                    NoFledged_SumPerYear=unlist(sapply(sum.0, "[", 2))); row.names(sum.1) <- 1:nrow(sum.1)
nrow(sum.1)==sum(as.numeric(lapply(sum.0, nrow)))

nrow(sum.1)
nrow(sum.1[complete.cases(sum.1),])
head(sum.1[!complete.cases(sum.1),], 20)

wild_NoFledged.0 <- merge(sum.1, keep_NoFledged_for_std[,-4], by="YearOfBreeding")
wild_NoFledged.0$NoFledged_SumPerYear_zscore <- (wild_NoFledged.0$NoFledged_SumPerYear-wild_NoFledged.0$Mean)/wild_NoFledged.0$Sd
wild_NoFledged.0$NoFledged_SumPerYear_meancentered <- wild_NoFledged.0$NoFledged_SumPerYear-wild_NoFledged.0$Mean

# impute within-year total number of fledglings for birds with missing data
missing_fitness_data <- unique(sum.1[!complete.cases(sum.1),]$RingNumber)
dat <- wild_NoFledged.0[wild_NoFledged.0$RingNumber %in% missing_fitness_data,]
dat_summary <- data.frame(table(dat$RingNumber))

dat_out <- NULL
for(i in 1:length(missing_fitness_data)) {
  bird_x <- missing_fitness_data[i]
  dat_bird_x <- dat[dat$RingNumber==bird_x,]
  
  if(nrow(dat_bird_x)==1) {
    dat_bird_x$NoFledged_SumPerYear_zscore <- 0
    dat_bird_x$NoFledged_SumPerYear_meancentered <- 0
  } else {
    if(nrow(dat_bird_x)==2) {
      if(nrow(dat_bird_x[!complete.cases(dat_bird_x),])==2) {
        dat_bird_x$NoFledged_SumPerYear_zscore <- c(0,0)
        dat_bird_x$NoFledged_SumPerYear_meancentered <- c(0,0)
      } else {
        zscore <- dat_bird_x$NoFledged_SumPerYear_zscore[!is.na(dat_bird_x$NoFledged_SumPerYear_zscore)]
        meancentered <- dat_bird_x$NoFledged_SumPerYear_meancentered[!is.na(dat_bird_x$NoFledged_SumPerYear_meancentered)]
        dat_bird_x$NoFledged_SumPerYear_zscore[is.na(dat_bird_x$NoFledged_SumPerYear_zscore)] <- zscore
        dat_bird_x$NoFledged_SumPerYear_meancentered[is.na(dat_bird_x$NoFledged_SumPerYear_meancentered)] <- meancentered
      }} else {
        if(nrow(dat_bird_x)==3) {
          if(nrow(dat_bird_x[!complete.cases(dat_bird_x),])==3) {
            dat_bird_x$NoFledged_SumPerYear_zscore <- c(0,0,0)
            dat_bird_x$NoFledged_SumPerYear_meancentered <- c(0,0,0)
          } else {
            if(nrow(dat_bird_x[!complete.cases(dat_bird_x),])==2) {
              zscore <- dat_bird_x$NoFledged_SumPerYear_zscore[!is.na(dat_bird_x$NoFledged_SumPerYear_zscore)]
              meancentered <- dat_bird_x$NoFledged_SumPerYear_meancentered[!is.na(dat_bird_x$NoFledged_SumPerYear_meancentered)]
              dat_bird_x$NoFledged_SumPerYear_zscore[is.na(dat_bird_x$NoFledged_SumPerYear_zscore)] <- c(zscore, zscore)
              dat_bird_x$NoFledged_SumPerYear_meancentered[is.na(dat_bird_x$NoFledged_SumPerYear_meancentered)] <- c(meancentered, meancentered)
            } else {
              zscore <- mean(dat_bird_x$NoFledged_SumPerYear_zscore[!is.na(dat_bird_x$NoFledged_SumPerYear_zscore)])
              meancentered <- mean(dat_bird_x$NoFledged_SumPerYear_meancentered[!is.na(dat_bird_x$NoFledged_SumPerYear_meancentered)])
              dat_bird_x$NoFledged_SumPerYear_zscore[is.na(dat_bird_x$NoFledged_SumPerYear_zscore)] <- zscore
              dat_bird_x$NoFledged_SumPerYear_meancentered[is.na(dat_bird_x$NoFledged_SumPerYear_meancentered)] <- meancentered
            }}}}}
  dat_out <- rbind(dat_out, dat_bird_x)
}

dat_ok <- wild_NoFledged.0[!wild_NoFledged.0$RingNumber %in% missing_fitness_data,]
dat_imputed_fitness <- rbind(dat_ok, dat_out)
nrow(dat_imputed_fitness[complete.cases(dat_imputed_fitness$NoFledged_SumPerYear_zscore),])==nrow(dat_imputed_fitness)

# get total number of fledgling across years
sum_zscore <- tapply(dat_imputed_fitness$NoFledged_SumPerYear_zscore, dat_imputed_fitness$RingNumber, sum)
sum_meancentered <- tapply(dat_imputed_fitness$NoFledged_SumPerYear_meancentered, dat_imputed_fitness$RingNumber, sum)

Wild_NoFledged.1 <- data.frame(SumNoFledged_zscore=sum_zscore, SumNoFledged_meancentered=sum_meancentered); Wild_NoFledged.1$RingNumber <- row.names(Wild_NoFledged.1); row.names(Wild_NoFledged.1) <- 1:nrow(Wild_NoFledged.1)
length(unique(dat_imputed_fitness$RingNumber))==nrow(Wild_NoFledged.1)

# combine with bird-specific data
temp_fitness.3 <- temp_fitness.1[!duplicated(temp_fitness.1$RingNumber), c(1:3,7)]
Wild_NoFledged <- merge(temp_fitness.3, Wild_NoFledged.1, by="RingNumber")

# combine ERC and wild data
FINAL_Dat_NoFledged <- rbind(ERC_NoFledged, Wild_NoFledged)

# 4. combine data on number of fledglings with lay date data: for females only of course..
# get lay dates of all birds that laid in 2018-2020 (recorded in any year) and only use the first recorded lay date (to make better comparable with lay dates from ERC birds)
lay_dates_temp.0 <- ForLayDates.1.0[ForLayDates.1.0$RingNumberFemale %in% ForLayDates_females, c(4,7,12)]
# get z-scored and mean-centered lay dates
lay_dates_temp.0_std <- merge(lay_dates_temp.0, keep_LayDates_for_std[,-4], by="YearOfBreeding")
lay_dates_temp.0_std$LayDateApril_std_zscore <- (lay_dates_temp.0_std$LayDateApril-lay_dates_temp.0_std$Mean)/lay_dates_temp.0_std$Sd
lay_dates_temp.0_std$LayDateApril_std_meancentered <- lay_dates_temp.0_std$LayDateApril-lay_dates_temp.0_std$Mean

# sort data to keep only first lay date record of each female
lay_dates_temp.1 <- lay_dates_temp.0_std[order(lay_dates_temp.0_std$RingNumberFemale, lay_dates_temp.0_std$YearOfBreeding, lay_dates_temp.0_std$LayDateApril),]
lay_dates_temp.2 <- lay_dates_temp.1[!duplicated(lay_dates_temp.1$RingNumberFemale),]

# all females in data frame with number of fledglings?
females_with_lay_date_and_fitness <- intersect(FINAL_Dat_NoFledged$RingNumber, lay_dates_temp.2$RingNumber)
# all females that have a lay date also have number of fledglings
# 10 females with number of fledglings have no lay date

# check which females to be sure all ok here:
check <- FINAL_Dat_NoFledged[FINAL_Dat_NoFledged$Sexe==1 & !(FINAL_Dat_NoFledged$RingNumber %in% females_with_lay_date_and_fitness),]
# two of them are selection line females; makes sense as 18 females in lay date analysis and 20 females have data on fledglings
nrow(FINAL_Dat_NoFledged[FINAL_Dat_NoFledged$Sexe==1 & FINAL_Dat_NoFledged$SelectionLineNumeric!=0,])
Fitness[Fitness$RingNumberFemale %in% check$RingNumber & Fitness$YearOfBreeding>=2018,] # all have Brood Type 6 or 1 (replacement or second clutch)

# all ok, so merge data:
dot_merge <- lay_dates_temp.2[,c(2,6,7)]; names(dot_merge)[1] <- "RingNumber"
FINAL_Dat_NoFledged_LayDate <- merge(FINAL_Dat_NoFledged, dot_merge, by="RingNumber", all=FALSE)

# save formatted data
save(FINAL_Dat_NoFledged, FINAL_Dat_NoFledged_LayDate,
     file="data/FINAL_Dat_Number_fledglings.RData")


### Statistical analysis

## MODEL for life-time number of fledglings produced (Eq. 6 and Eq. 5) for female selection line fledglings and female local (non-selection line) fledglings

# prepare data for model
d.0 <- FINAL_Dat_NoFledged_LayDate[FINAL_Dat_NoFledged_LayDate$RingNumber!="AU...91626",]
d.0$help <- ifelse(d.0$SelectionLineNumeric==0,3,0)
d.0$SelectionLineNumeric2 <- d.0$SelectionLineNumeric+d.0$help
d <- d.0
d$R <- d$SumNoFledged_zscore
d$L <- as.factor(d$SelectionLineNumeric2)
d$A <- d$LayDateApril_std_zscore

table(d$L)

# run analysis for female local (non-selection line) recruits (Eq. 6):
d.local <- d[d$L==3,]
dat <- list(
  R=d.local$R,
  A=d.local$A
)

mF.local <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a+bA*A,
    a ~ dnorm(0,1.5),
    bA ~ dnorm(0,1.5),
    sigma ~ dhalfnorm(0,1)
  ), data=dat, chains=4, cores=4, iter=10000)
save(mF.local, file="temp/mF.local.RData")

(out <- precis(mF.local, depth=2))
pars <- row.names(out)
out.0 <- data.frame(parameter=pars, out)
# Table S21
write.table(out.0, "out/Bayes_NoFLedged_Local_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

pars <- row.names(precis(mF.local, depth=2))
par(mfrow = c(1,1))
# Fig. S33
pdf(file = "Plots/Bayes_NoFledged_Local_traceplots.pdf", width = 6, height = 6)
traceplot(mF.local, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()
# Fig. S34
pdf(file = "Plots/Bayes_NoFledged_Local_trankplot.pdf", width = 6, height = 6)
trankplot(mF.local, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# run analysis for selection line birds only (Eq. 5):
d.SL <- d[d$L!=3,]
dat <- list(
  R=d.SL$R,
  L=d.SL$L
)

mF.SL <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a+bL[L],
    a ~ dnorm(0,1.5),
    bL[L] ~ dnorm(0,1.5),
    sigma ~ dhalfnorm(0,1)
  ), data=dat, chains=4, cores=4, iter=10000)
save(mF.SL, file="temp/mF.SL.RData")

(out <- precis(mF.SL, depth=2))
pars <- row.names(out)
out.0 <- data.frame(parameter=pars, out)
# Table S19
write.table(out.0, "out/Bayes_NoFLedged_SL_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

pars <- row.names(precis(mF.SL, depth=2))
par(mfrow = c(1,1))
# Fig. S31
pdf(file = "Plots/Bayes_NoFledged_SL_traceplots.pdf", width = 6, height = 6)
traceplot(mF.SL, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()
# Fig. S32
pdf(file = "Plots/Bayes_NoFledged_SL_trankplot.pdf", width = 6, height = 6)
trankplot(mF.SL, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# get posterior mean and ci
dat <- list(L=1:2)
p_post <- link(mF.SL, data=dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)
diff <- p_post[,2]-p_post[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1:2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
keep_post.out <- post.out
# Table S20
write.table(post.out, "out/Bayes_number_fledglings_SL_posterior_prediction.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

############### small sidetrack start..
# get total fitness: recruitment probability of selection line fledglings * life-time number of offspring produced by selection line fledglings
Keep.no.fledglings <- list(early=p_post[,1], late=p_post[,2])
Combined_fitness <- matrix(c(Keep.recruitment.prob$early*Keep.no.fledglings$early,
                         Keep.recruitment.prob$late*Keep.no.fledglings$late),
                         ncol=2)
p_mu <- apply(Combined_fitness, 2, mean)
p_ci <- apply(Combined_fitness, 2, PI)
diff <- Combined_fitness[,2]-Combined_fitness[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
post.out <- data.frame(line=c(1:2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S32
write.table(post.out, "out/Bayes_Combined_fitness_posterior_prediction.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
############### .. small sidetrack end

# plot posterior prediction 
# prepare plot
post.out <- keep_post.out[-3,]
post.out$x.lab <- c(1:2)
post.out$n <- as.numeric(table(d$SelectionLineNumeric[d$SelectionLineNumeric!=0]))
plot.dat_SL <- post.out

# get lay dates for plot below (Table S2):
post.out.lay <- read.table("out/Bayes_Lay_date_posterior_prediction.d.txt", header=TRUE, sep="\t") 
names(post.out.lay) <- paste(names(post.out.lay), "_April", sep="")
# combine in data frame with recruitment probability
for_common_plot <- cbind(post.out[1:2,-(5:6)], post.out.lay[1:2,-c(1)])

# plot posterior prediction: total number of fledglings relative to standardized April date of lay date
# prepare plot
d_plot <- d.local
range <- range(d_plot$A)
xseq <- seq(range[1], range[2], length.out=50)
d_pred <- list(A=xseq)
p_post <- link(mF.local, data=d_pred)

p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)

# Fig. 4
pdf(file = "Plots/Main_Bayes_number_fledglings_v2_revised.pdf", width = 6.5, height = 3.5)
layout(t(matrix(c(1,2))), width=c(1.25,1.5))

# A
par(mar = c(4,4,1.2,1), mgp=c(2, 0.5, 0))
plot(NULL, xlim=c(0.7,2.3), ylim=c(-1.8,1.7), xlab="",
     ylab="Standardized fledglings produced",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", post.out$n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-1.5,1.5, 0.5), labels=seq(-1.5,1.5, 0.5), lwd = 0, lwd.ticks = 1)
arrows(x0=post.out$x.lab, y0=post.out$ci.5.5, x1=post.out$x.lab, y1=post.out$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(post.out$x.lab, post.out$mu, pch=19, cex=3, col=ifelse(post.out$x.lab==1, "darkorange", "darkslategrey"))
points(jitter(as.numeric(d$L),0.5), d$R, pch=1, cex=1.5, col=ifelse(d$L==1, "darkorange", "darkslategrey"))
text(0.7+0.04, 1.7, "A", cex=1.2, pos=1, offset=0, font=2)

# B
par(mar = c(4,4,1.2,1), mgp=c(2, 0.5, 0))
plot(NULL, xaxt="n", yaxt="n", xlim=range(xseq)+0.05, ylim=c(-4.3,5.7),
     xlab="Standardized lay date", ylab="Standardized fledglings produced", cex.main=1)
axis(1, at=seq(-1,4,1), labels=seq(-1,4,1), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-4,6,2), labels=seq(-4,6,2), lwd = 0, lwd.ticks = 1)
# points(d_plot$A, d_plot$R, pch=1, col="grey")
shade(p_ci, xseq, col=col.alpha("black", 0.1))
lines(xseq, p_mu, lwd=1, lty=2)
text(range(xseq)[1]+0.09,5.7, "B", cex=1.2, pos=1, offset=0, font=2)

# add selection line data
arrows(x0=for_common_plot$mu_April, y0=for_common_plot$ci.5.5, x1=for_common_plot$mu_April, y1=for_common_plot$ci.94.5, code=3, angle=90, length=0.1, col=ifelse(for_common_plot$line==1, "black", "black"), lwd=1)
arrows(x0=for_common_plot$ci.5.5_April, y0=for_common_plot$mu, x1=for_common_plot$ci.94.5_April, y1=for_common_plot$mu, code=3, angle=90, length=0.1, col=ifelse(for_common_plot$line==1, "black", "black"), lwd=1)
points(for_common_plot$mu_April, for_common_plot$mu, pch=19, cex=1.5, col=ifelse(for_common_plot$line==1, "darkorange", "darkslategrey"))
dev.off()


###  Selective survival during nest building
### ---------------------------------------------------------------------------

### load and prepare data

Nest_building_records.lay <- read.table("out/Nestbuilding_Records_females.txt", header=TRUE, row.names=NULL, sep="\t")

# exclude birds with <=5 visits
temp.0 <- Nest_building_records.lay[Nest_building_records.lay$Visits>5,]

# get nest building state from books
Nb.dates <- read.table("temp/Nest.building.locations.csv", sep=",", header=TRUE)

if(unique(temp.0$YearDetected==Nb.dates$YearDetected & temp.0$HVLocation==Nb.dates$HVLocation)) temp.0$StartNest <- Nb.dates$StartNest
if(unique(temp.0$YearDetected==Nb.dates$YearDetected & temp.0$HVLocation==Nb.dates$HVLocation)) temp.0$StatusNest <- Nb.dates$StatusNest

help <- data.frame(StatusNest=sort(unique(temp.0$StatusNest)), CorrectionNest=c(4,2,5,0))
temp.0 <- merge(temp.0, help, by="StatusNest"); temp.0$CorrectedStartNest <- temp.0$StartNest-temp.0$CorrectionNest

temp.0$SelectionLine <- sub("_F[45]", "", temp.0$SelectionLine); temp.0$Number.Females <- 1

# check for females with multiple occurrences
temp.0[duplicated(temp.0[,c(1:3)]),] # one female found in same nest box in same year
temp.0[temp.0$RingNumber=="BD...70561",] # most likely same female, but first breeding attempt failed.

temp.0[duplicated(temp.0[,c(1,3)]),] # another female found in two nest boxes in same year
temp.0[temp.0$RingNumber=="BD...59338" & temp.0$Year==2018,]

# get females that did end up laying in the same nest box
temp.1 <- temp.0[!is.na(temp.0$RingNumberFemale),]
same.box <- temp.1[temp.1$RingNumber==temp.1$RingNumberFemale,]; same.box$lay <- 1
n.same.box <- aggregate(Number.Females ~ YearDetected + SelectionLine, data=same.box, FUN="sum"); n.same.box <- n.same.box[order(n.same.box$YearDetected, n.same.box$SelectionLine),]
n.same.box$lay <- 1

# get females that did not end up in the nest box and where no female was laying in that same nest box
no.lay <- temp.0[is.na(temp.0$RingNumberFemale),]; no.lay$lay <- 0
n.no.lay <- aggregate(Number.Females ~ YearDetected + SelectionLine, data=no.lay, FUN="sum")
n.no.lay$lay <- 0

# get females that did not end up in the nest box and where another female was laying in that same nest box
temp.1 <- temp.0[!is.na(temp.0$RingNumberFemale),]
other.box <- temp.1[temp.1$RingNumber!=temp.1$RingNumberFemale,]; other.box$lay <- 0
n.other.box <- aggregate(Number.Females ~ YearDetected + SelectionLine, data=other.box, FUN="sum"); n.other.box <- n.other.box[order(n.other.box$YearDetected, n.other.box$SelectionLine),]
n.other.box$lay <- 1

nest.building.laid <- rbind(same.box, other.box, no.lay)
Summary_nest.building <- rbind(n.same.box, n.other.box, n.no.lay)


### Statistical analysis

## MODEL for local recruitment probability (Eq. S2) when female selection line recruits awere identified during nest building

# make input data:
d <- nest.building.laid
dat <- list(
  R=d$lay,
  A=standardize(d$CorrectedStartNest),
  L=ifelse(d$SelectionLine=="E",1,2),
  Y=d$YearDetected-2017
)

# fit model
mNB.1a <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + bL[L] + bA*A,
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bA ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mNB.1a, file="temp/mNB.1a.multilevel.RData")

# fit model again but with line-specific slopes
mNB.1b <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + bL[L] + bAL[L]*A,
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    bAL[L] ~ dnorm(0, 1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mNB.1b, file="temp/mNB.1b.multilevel.RData")

# get model summary:
(out <- precis(mNB.1b, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S9
write.table(out.0, "out/Bayes_NestBuilding_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

pars <- row.names(precis(mNB.1b, depth=2)); pars <- pars[-c((length(pars)-2):length(pars))]
par(mfrow = c(1,1))
# Fig. S21
pdf(file = "Plots/Bayes_NestBuilding_traceplots.pdf", width = 6, height = 6)
traceplot(mNB.1b, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()
# Fig. S22
pdf(file = "Plots/Bayes_NestBuilding_trankplot.pdf", width = 6, height = 6)
trankplot(mNB.1b, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# posterior prediction
# prepare data for prediction
d_pred <- list(L=rep(1:2, each=3), Y=rep(1:3, times=2), A=rep(0, times=6))

# get posterior
p_post <- link(mNB.1b, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S10
write.table(post.out, "out/Bayes_NestBuilding_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction for fledglings genetically derived from selection line breeding pairs
# prepare plot
post.out$x.lab <- c(1:3)
n <- as.numeric(table(dat$L))

# Fig. S3
pdf(file = "Plots/Bayes_NestBuilding_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(3,5,2,0.5), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[-3,]
plot(NULL, xlim=c(0.7,2.3), ylim=c(0.4,0.95), xlab="",
     ylab="Local recruitment probability\n when nest building",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0.5,0.9, 0.2), labels=seq(0.5,0.9, 0.2), lwd = 0, lwd.ticks = 1)
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
text(0.7, 0.95, "A", cex=1.2, pos=1, offset=0, font=2)

par(mar = c(3,5,2,2), mgp=c(2, 0.5, 0))
plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-0.4,0.2),
     xlab="", ylab="Difference in ocal recruitment probability\nwhen nest building (late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(-0.4,0.2, 0.2), labels=seq(-0.4,0.2, 0.2), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025,0.2, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()



###  Selective survival based on roost inspections
### ---------------------------------------------------------------------------

### load and prepare data
Roost.raw <- read.table("data/Qry_FindSelectionLineBirdsIn AVC_2.csv", sep=",", header=TRUE) 
Roost.raw <- Roost.raw[Roost.raw$YearOfAVC<2021,]

### Statistical analysis

## MODEL for local recruitment probability (Eq. S1) when female selection line recruits awere identified during roosting inspections in january

dat <- list(
  R=ifelse(!is.na(Roost.raw$YearOfBreeding), 1, 0),
  L=ifelse(grepl("E", Roost.raw$SelectionLine),1,2),
  Y=Roost.raw$YearOfAVC-2017
)

# fit model
mRC.1a <- ulam(
  alist(
    R ~ dbinom(1,p),
    logit(p) <- a_bar + z[Y]*sigma_a + bL[L],
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mRC.1a, file="temp/mRC.1a.multilevel.RData")

# get model summary:
(out <- precis(mRC.1a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S7
write.table(out.0, "out/Bayes_RoostCheck_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

pars <- row.names(precis(mRC.1a, depth=2)); pars <- pars[-c((length(pars)-2):length(pars))]
par(mfrow = c(1,1))
# Fig. S19
pdf(file = "Plots/Bayes_RoostCheck_traceplots.pdf", width = 6, height = 6)
traceplot(mRC.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()
# Fig. S20
pdf(file = "Plots/Bayes_RoostCheck_trankplot.pdf", width = 6, height = 6)
trankplot(mRC.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# get posterior prediction
# prepare 'new' data for prediction; use unique combinations of year, selection line and brood location ID as in fitted data
d_pred <- list(L=rep(1:2, each=3), Y=rep(1:3, times=2))

# get posterior
p_post <- link(mRC.1a, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S8
write.table(post.out, "out/Bayes_RoostCheck_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction 
# prepare plot
post.out$x.lab <- c(1:3)
n <- as.numeric(table(dat$L))

# Fig. S2
pdf(file = "Plots/Bayes_RoostCheck_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(3,5,2,0.5), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[-3,]
plot(NULL, xlim=c(0.7,2.3), ylim=c(0.27,0.9), xlab="",
     ylab="Local recruitment probability\nwhen roosting",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0.3,0.9, 0.3), labels=seq(0.3,0.9, 0.3), lwd = 0, lwd.ticks = 1)
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
text(0.7, 0.9, "A", cex=1.2, pos=1, offset=0, font=2)

par(mar = c(3,5,2,2), mgp=c(2, 0.5, 0))
plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-0.2,0.47),
     xlab="", ylab="Difference in local recruitment\nprobability when roosting (late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(-0.2,0.4, 0.2), labels=seq(-0.2,0.4, 0.2), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025,0.47, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()


###  Phenological mismatch
### ---------------------------------------------------------------------------

# get lay date and food peak data to calculate optimal lay date
lay.date.dat <- FINAL_Dat_LayDates[,1:4]
food <- read.csv("data/tbl_PeakDate_Biomass_HVLim.csv",header=T,sep=';')

Optimal.dal.date <- merge(lay.date.dat, food[,2:3], by.x="YearOfBreeding", by.y="Year")
Optimal.dal.date$OptimalLayDateApril <- Optimal.dal.date$MidDate-33
Optimal.dal.date$Diff_LD_OptLD <- Optimal.dal.date$LayDateApril-Optimal.dal.date$OptimalLayDateApril

### Statistical analysis

## MODEL for phenological mismatch (Eq. S6) of local (non-selection line) female recruits

# prepare data for model:
d <- Optimal.dal.date[Optimal.dal.date$SelectionLineNumeric==0,] 
mean_mm <- mean(d$Diff_LD_OptLD)
sd_mm <- sd(d$Diff_LD_OptLD)

# run model:
dat <- list(
  R=d$Diff_LD_OptLD,
  Y=as.factor(d$YearOfBreeding-2017)
)

# run model:
mMM.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a,
    sigma ~ dhalfnorm(0,1),
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 2.5), 
    sigma_a ~ dhalfnorm(0,1),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mMM.1a, file="temp/mMM.1a.RData")

# get model summary:
(out <- precis(mMM.1a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S5
write.table(out.0, "out/Bayes_Mismatch.1_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

pars <- row.names(precis(mMM.1a, depth=2)); pars <- pars[-c((length(pars)-2):length(pars))]
par(mfrow = c(1,1))
# Fig. S17
pdf(file = "Plots/Bayes_Mismatch.1_traceplots.pdf", width = 6, height = 6)
traceplot(mMM.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()
# Fig. S18
pdf(file = "Plots/Bayes_Mismatch.1_trankplot.pdf", width = 6, height = 6)
trankplot(mMM.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# get posterior prediction
# prepare data for prediction
d_pred <- list(Y=1:3)

# get posterior for each year
p_post <- link(mMM.1a, data=d_pred)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)
diff <- p_post[,2]-p_post[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
diff_e_w <- p_post[,3]-p_post[,1]
d_e_w_mu <- mean(diff_e_w); d_e_w_ci <- PI(diff_e_w)
diff_l_w <- p_post[,3]-p_post[,2]
d_l_w_mu <- mean(diff_l_w); d_l_w_ci <- PI(diff_l_w)
post.out <- data.frame(year=c(2018:2020,"diff_2018_2019", "diff_2018_2020", "diff_2019_2020"), mu=c(p_mu, d_mu, d_e_w_mu, d_l_w_mu), ci.5.5=c(p_ci[1,],d_ci[1],d_e_w_ci[1],d_l_w_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2],d_e_w_ci[2],d_l_w_ci[2]))
# Table S6
write.table(post.out, "out/Bayes_Mismatch.1_per.year_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Fig. S1
pdf(file = "Plots/Main_Bayes_Mismatch_per.year.pdf", width = 5, height = 4.5)
layout(matrix(c(1,2,3)), heights=c(1,1,1.3))

cap <- c("2018", "2019", "2020")
d.sl <- Optimal.dal.date[Optimal.dal.date$SelectionLineNumeric!=0,]
lim <- range(Optimal.dal.date$Diff_LD_OptLD)

for(i in 1:3) {
  d.sub <- d.sl[d.sl$YearOfBreeding==(i+2017),]
  if(i!=3) {
    par(mar = c(1,1,1.2,1), mgp=c(2, 0.5, 0))
    plot(NULL, xlim=lim+0.05, ylim=c(0,5000),
         xlab="", ylab="",
         xaxt="n", yaxt="n")
    hist(p_post[,i], col="gray30", border="gray30", add=TRUE)
    axis(1, at=seq(-15,25, 5), labels=FALSE, lwd = 0, lwd.ticks = 1)
    abline(v=0, lty=2)
    text(lim[1]+0.03,5000, cap[i], cex=1.2, pos=1, offset=0, font=2)
    # add post distributions of lines
    points(d.sub$Diff_LD_OptLD, jitter(rep(4500, nrow(d.sub)), 3), pch=ifelse(d.sub$RingNumberFemale=="AU...91626",13,19), cex=1.5, col=ifelse(d.sub$SelectionLineNumeric==1, "darkorange", "darkslategrey"))
    
  } else {
    par(mar = c(4.2,1,1.2,1), mgp=c(2, 0.5, 0))
    plot(NULL, xlim=lim+0.05, ylim=c(0,5000),
         xlab="mismatch (observed - optimal lay date)", ylab="",
         xaxt="n", yaxt="n")
    hist(p_post[,i], col="gray30", border="gray30", add=TRUE)
    axis(1, at=seq(-15,25, 5), labels=seq(-15,25, 5), lwd = 0, lwd.ticks = 1)
    abline(v=0, lty=2)
    text(lim[1]+0.03,5000, cap[i], cex=1.2, pos=1, offset=0, font=2)
    # add post distributions of lines
    points(d.sub$Diff_LD_OptLD, jitter(rep(4500, nrow(d.sub)), 3), pch=ifelse(d.sub$RingNumberFemale=="AU...91626",13,19), cex=1.5, col=ifelse(d.sub$SelectionLineNumeric==1, "darkorange", "darkslategrey"))
    
    
    
  }}
dev.off()


## MODEL for phenological mismatch (Eq. 2) of female selection line recruits

# get within-year standardized mismatch (in line with lay date analysis)
# get mean and standard deviation for standardization: 
ForMM.1 <- Optimal.dal.date[Optimal.dal.date$SelectionLineNumeric==0,,]
# remove selection line females (must not be included in calculating mean and standard deviation for standardization)
mean <- tapply(ForMM.1$Diff_LD_OptLD, ForMM.1$YearOfBreeding, mean)
sd <- tapply(ForMM.1$Diff_LD_OptLD, ForMM.1$YearOfBreeding, sd)
n <- tapply(ForMM.1$Diff_LD_OptLD, ForMM.1$YearOfBreeding, length)
Summary_MM <- data.frame(YearOfBreeding=names(mean), Mean=mean, Sd=sd, n_BreedingFemales=n); row.names(Summary_lay) <- 1:nrow(Summary_lay)
keep_MM_for_std <- Summary_MM

# standardization
Optimal.dal.date_std <- merge(Optimal.dal.date, keep_MM_for_std, by="YearOfBreeding")
Optimal.dal.date_std$Mismatch_std <- Optimal.dal.date_std$Diff_LD_OptLD/Optimal.dal.date_std$Sd

# run model:
d <- Optimal.dal.date_std[Optimal.dal.date_std$RingNumberFemale!="AU...91626",]
d$help <- ifelse(d$SelectionLineNumeric==0,3,0)
d$SelectionLineNumeric2 <- d$SelectionLineNumeric+d$help

dat <- list(
  R=d$Mismatch_std,
  L=as.factor(d$SelectionLineNumeric2)
)

# run model:
mMM.2a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a + b[L],
    sigma ~ dhalfnorm(0,1),
    a ~ dnorm(0, 1.5), 
    b[L] ~ dnorm(0,1.5)
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mMM.2a, file="temp/mMM.2a.RData")

# get model summary:
(out <- precis(mMM.2a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S3
write.table(out.0, "out/Bayes_Mismatch_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

pars <- row.names(precis(mMM.2a, depth=2))
par(mfrow = c(1,1))
# Fig. S15
pdf(file = "Plots/Bayes_Mismatch_traceplots.pdf", width = 6, height = 6)
traceplot(mMM.2a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()
# Fig. S16
pdf(file = "Plots/Bayes_Mismatch_trankplot.pdf", width = 6, height = 6)
trankplot(mMM.2a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=2)
dev.off()

# get posterior prediction
# prepare data for prediction
d_pred <- list(L=1:3)

# get posterior
p_post <- link(mMM.2a, data=d_pred)

# get mean for early and late selection line, respectively
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)
diff <- p_post[,2]-p_post[,1]
d_mu <- mean(diff); d_ci <- PI(diff)
diff_e_w <- p_post[,3]-p_post[,1]
d_e_w_mu <- mean(diff_e_w); d_e_w_ci <- PI(diff_e_w)
diff_l_w <- p_post[,3]-p_post[,2]
d_l_w_mu <- mean(diff_l_w); d_l_w_ci <- PI(diff_l_w)
post.out <- data.frame(line=c(1:3,"diff_e_l", "diff_e_w", "diff_l_w"), mu=c(p_mu, d_mu, d_e_w_mu, d_l_w_mu), ci.5.5=c(p_ci[1,],d_ci[1],d_e_w_ci[1],d_l_w_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2],d_e_w_ci[2],d_l_w_ci[2]))
# Table S4
write.table(post.out, "out/Bayes_Mismatch_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction
n_line <- data.frame(table(dat$L)[-3]); names(n_line)[2] <- "n"; n_line$lab <- ifelse(n_line$Var1==1, "early", "late")
d <- Optimal.dal.date_std[Optimal.dal.date_std$SelectionLineNumeric!=0,]
d$R <- d$Mismatch_std
d$L <- as.factor(d$SelectionLineNumeric)
d$Fem <- d$RingNumberFemale

# Fig. 2
pdf(file = "Plots/Bayes_Mismatch.distributions_v2.pdf", width = 3.3, height = 3)
par(mfrow=c(1,1), mar = c(4,1,1.2,1), mgp=c(2, 0.5, 0))

plot(NULL, xlim=range(p_post), ylim=c(0,6000),
     xlab="Standardized mismatch", ylab="",
     xaxt="n", yaxt="n")
hist(p_post[,1], col=col.alpha("darkorange", 0.75), border="white", add=TRUE)
hist(p_post[,2], col=col.alpha("darkslategrey", 0.75), border="white", add=TRUE)
hist(p_post[,3], col="gray30", border="gray30", add=TRUE)
abline(v=0, lty=2)
axis(1, at=seq(-2,5, 1), labels=seq(-2,5, 1), lwd = 0, lwd.ticks = 1)
dev.off()


###  Fledgling weight, tarsus length and P3 length of fledglings produced by female selection line recruits
### ---------------------------------------------------------------------------

Chicks.sl.raw <- read.table("Qry_FindBiometryFledglingsOfRecruitedSelectionLineBirdsAtHV_4.csv", sep=",", header=TRUE)
Chicks.sl.weight <- Chicks.sl.raw[!is.na(Chicks.sl.raw$Weight),]
Chicks.sl.weight <- Chicks.sl.weight[!duplicated(Chicks.sl.weight),]
Dup.chicks <- Chicks.sl.weight[duplicated(Chicks.sl.weight$IndivID),]$IndivID

# quite a few duplicates; check why and remove
Dup.chicks.dat <- Chicks.sl.weight[Chicks.sl.weight$IndivID %in% Dup.chicks,]
NoDup.chicks.dat <- Chicks.sl.weight[!Chicks.sl.weight$IndivID %in% Dup.chicks,] # keep data for uniwue ring numbers

nrow(Chicks.sl.weight)==nrow(NoDup.chicks.dat)+nrow(Dup.chicks.dat)
unique(Dup.chicks.dat$YearBorn)
# only 2018; weight taken d7 and d15!
# remove and keep d15 data for duplicated ring numbers
Dup.chicks.sorted <- Dup.chicks.dat[order(Dup.chicks.dat$RingNumber, Dup.chicks.dat$CaptureDate, decreasing=TRUE),]
Dup.chicks.no_dupl <- Dup.chicks.sorted[!duplicated(Dup.chicks.sorted$IndivID),]

# combine
Chicks.noDup <- rbind(NoDup.chicks.dat, Dup.chicks.no_dupl)
nrow(Chicks.noDup)==nrow(Chicks.sl.weight[!duplicated(Chicks.sl.weight$IndivID),]) #ok

# get females and year for which we need chick wieght:
females.with.young <- FINAL_Dat_NoFledged_LayDate[FINAL_Dat_NoFledged_LayDate$SelectionLineNumeric!=0, c(1,3,4)]
females.with.lay.date <- FINAL_Dat_LayDates[FINAL_Dat_LayDates$SelectionLineNumeric!=0, c(1,2,4)]
sort(females.with.young$RingNumber)==sort(females.with.lay.date$RingNumberFemale)

# only use chick weight of first year of breeding:
laying.females <- FINAL_Dat_LayDates[FINAL_Dat_LayDates$SelectionLineNumeric!=0, 1:4]

chicks.laying.females <- merge(laying.females, Chicks.noDup[,c(2,4,7,10,12,14,15)], by.x=c("YearOfBreeding", "RingNumberFemale"), by.y=c("YearBorn", "Mother"))
range(chicks.laying.females$Weight)
range(chicks.laying.females$P3_Length)
range(chicks.laying.females$Tarsus)

fledglings.laying.females <- merge(chicks.laying.females, FledglingsHV[,c(2,3,12)], by="RingNumber")
fledglings.laying.females[fledglings.laying.females$Fledgling!=1,] # 2 chicks did not fledge; remove
fledglings.laying.females <- fledglings.laying.females[fledglings.laying.females$Fledgling==1,]

## cross check this with total number of chick given above !!!
# get number of fledglings per female (add year info) from capture data
fledglings.per.female <- as.data.frame(table(fledglings.laying.females$RingNumberFemale)); names(fledglings.per.female) <- c("RingNumber", "NumberFledged")
brood.year <- fledglings.laying.females[!duplicated(fledglings.laying.females[,2:3]),2:3]
brood.year.0 <- brood.year[match(fledglings.per.female$RingNumber, brood.year$RingNumberFemale),]
if(unique(brood.year.0$RingNumberFemale==fledglings.per.female$RingNumber)) fledglings.per.female$YearOfBreeding  <- brood.year.0$YearOfBreeding 

# add total number of fledglings from brood data:
fledglings.per.female.0 <- merge(fledglings.per.female, ERC_NoFledged.0[,1:3], by=c("YearOfBreeding", "RingNumber"))
fledglings.per.female.0$difference <- fledglings.per.female.0$NumberFledged-fledglings.per.female.0$NoFledged_SumPerYear

# 3 females have different number..
problems <- fledglings.per.female.0[fledglings.per.female.0$difference!=0,]
Recruit.problems <- merge(Recruits_noHK[,c(6,9,15,18)], problems, by.x=c("YearOfBreeding", "Mother"), by.y=c("YearOfBreeding", "RingNumber"))

### GET BACK TO BOOKS AND CHECK WHAT TABLE TELLS THE TRUTH  

# for now keep as is and prepare data for analysis:
# get d15 capture as April date
chicks.laying.females <- chicks.laying.females[chicks.laying.females$RingNumber %in% fledglings.laying.females$RingNumber,]
chicks.laying.females$CaptureDate_formatted <- as.POSIXct(chicks.laying.females$CaptureDate, format='%Y-%m-%d', tz="EST") 
chicks.laying.females$CaptureDateApril <- as.numeric(difftime(chicks.laying.females$CaptureDate_formatted, as.POSIXct(paste(chicks.laying.females$YearOfBreeding ,"03-31", sep="-"), format='%Y-%m-%d', tz="EST"), units="days"))

# how many females have more than one brood?
broods.per.female <- chicks.laying.females[!duplicated(chicks.laying.females[,c(2,11)]), c(2,11)]
trouble.chicks <- chicks.laying.females[chicks.laying.females$RingNumberFemale=="BD...70661",]
trouble.chicks$CaptureDateApril <- 84
good.chicks <- chicks.laying.females[chicks.laying.females$RingNumberFemale!="BD...70661",]

chicks.laying.females.0 <- rbind(good.chicks, trouble.chicks)
broods.per.female <- chicks.laying.females.0[!duplicated(chicks.laying.females.0[,c(2,11)]), c(2,11)]

FINAL_chick.data <- chicks.laying.females.0[,c(1:5,7:9,11)]
save(FINAL_chick.data, file="data/FINAL_chick.data.RData")

### Statistical analysis

## MODEL for fledgling weight, tarsus length and P3 length (Eq. S9) of fledglings produced by female selection line recruits 

# prepare data for model:
d <- FINAL_chick.data
# use female ring number as integer:
females.as.interger <- NULL
r.no <- d$RingNumberFemale
r.no.unique <- unique(d$RingNumberFemale)
for(i in 1:length(unique(d$RingNumberFemale))) {
  females.as.interger <- append(females.as.interger, rep(i, length(r.no[r.no==r.no.unique[i]])), length(females.as.interger))
}

dat <- list(
  R=standardize(d$Weight),
  L=as.factor(d$SelectionLineNumeric),
  Y=as.factor(d$YearOfBreeding-2017),
  F=as.factor(females.as.interger),
  A=standardize(d$CaptureDateApril)
)

# run model for fledgling weight:
mR.FW.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + x[F]*sigma_g + bL[L] + bA*A,
    sigma ~ dhalfnorm(0,1),
    z[Y] ~ dnorm(0, 1),
    x[F] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1),
    bA ~ dnorm(0, 1),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[F]:g <<- x*sigma_g
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mR.FW.1a, file="temp/mR.FW.1a.RData")

# get model summary:
(out <- precis(mR.FW.1a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S22
write.table(out.0, "out/Bayes_Recruit.Fledgling.weight_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mR.FW.1a)), row.names(precis(mR.FW.1a, depth=2))[grepl("z", row.names(precis(mR.FW.1a, depth=2)))], row.names(precis(mR.FW.1a, depth=2))[grepl("bL", row.names(precis(mR.FW.1a, depth=2)))])
par(mfrow = c(1,1))
# Fig. S35
pdf(file = "Plots/Bayes_Recruit.Fledgling.weight_traceplots.pdf", width = 6, height = 8)
traceplot(mR.FW.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()
# Fig. S36
pdf(file = "Plots/Bayes_Recruit.Fledgling.weight_trankplot.pdf", width = 6, height = 8)
trankplot(mR.FW.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior prediction
# prepare data for prediction
dat.0 <- data.frame(dat$F, dat$Y)
dat_help.0 <- dat.0[!duplicated(dat.0),]
dat_help <- dat_help.0[order(dat_help.0$dat.F, dat_help.0$dat.Y),]
d_pred <- list(L=rep(1:2, each=nrow(dat_help)), Y=rep(dat_help$dat.Y, times=2), F=rep(dat_help$dat.F, times=2), A=rep(0, times=nrow(dat_help)*2))

# get posterior
p_post <- link(mR.FW.1a, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S23
write.table(post.out, "out/Bayes_Recruit.Fledgling.weight_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction for fledglings genetically derived from selection line breeding pairs
# prepare plot
post.out$x.lab <- c(1:2,3)
n <- as.numeric(table(d$SelectionLineNumeric[d$SelectionLineNumeric!=0]))

# Fig. S6
pdf(file = "Plots/Bayes_Recruit.Fledgling.weight_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[1:2,]
plot(NULL, xlim=c(0.7,2.3), ylim=c(-4.2,2.75), xlab="",
     ylab="standardized fledgling weight",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-4,2, 1), labels=seq(-4,2, 1), lwd = 0, lwd.ticks = 1)
points(jitter(as.numeric(dat$L),0.5), dat$R, pch=1, cex=1.5, col=ifelse(dat$L==1, "darkorange", "darkslategrey"))
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.2, col="black", lwd=1.5)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
text(0.7, 2.75, "A", cex=1.2, pos=1, offset=0, font=2)

plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-0.8,0.8),
     xlab="", ylab="difference in standardized fledgling weight\n(late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(-0.6,0.6, 0.3), labels=seq(-0.6,0.6, 0.3), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025,0.8, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()

# run model for P3:
d <- FINAL_chick.data
dat <- list(
  R=standardize(d$P3_Length),
  L=as.factor(d$SelectionLineNumeric),
  Y=as.factor(d$YearOfBreeding-2017),
  F=as.factor(females.as.interger),
  A=standardize(d$CaptureDateApril)
)

# run model:
mR.FP3.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + x[F]*sigma_g + bL[L] + bA*A,
    sigma ~ dhalfnorm(0,1),
    z[Y] ~ dnorm(0, 1),
    x[F] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1),
    bA ~ dnorm(0, 1),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[F]:g <<- x*sigma_g
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mR.FP3.1a, file="temp/mR.FP3.1a.RData")


# get model summary:
(out <- precis(mR.FP3.1a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S26
write.table(out.0, "out/Bayes_Recruit.Fledgling.P3_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mR.FP3.1a)), row.names(precis(mR.FP3.1a, depth=2))[grepl("z", row.names(precis(mR.FP3.1a, depth=2)))], row.names(precis(mR.FP3.1a, depth=2))[grepl("bL", row.names(precis(mR.FP3.1a, depth=2)))])
par(mfrow = c(1,1))
# Fig. S39
pdf(file = "Plots/Bayes_Recruit.Fledgling.P3_traceplots.pdf", width = 6, height = 8)
traceplot(mR.FP3.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()
# Fig. S40
pdf(file = "Plots/Bayes_Recruit.Fledgling.P3_trankplot.pdf", width = 6, height = 8)
trankplot(mR.FP3.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior mean and ci for recruitment probability of chicks genetically derived from the early and late selection line and chicks genetically derived from local breeding birds
# prepare 'new' data for prediction; use unique combinations of year, selection line and brood location ID as in fitted data
dat.0 <- data.frame(dat$F, dat$Y)
dat_help.0 <- dat.0[!duplicated(dat.0),]
dat_help <- dat_help.0[order(dat_help.0$dat.F, dat_help.0$dat.Y),]
d_pred <- list(L=rep(1:2, each=nrow(dat_help)), Y=rep(dat_help$dat.Y, times=2), F=rep(dat_help$dat.F, times=2), A=rep(0, times=nrow(dat_help)*2))

# get posterior
p_post <- link(mR.FP3.1a, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S27
write.table(post.out, "out/Bayes_Recruit.Fledgling.P3_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction for fledglings genetically derived from selection line breeding pairs
# prepare plot
post.out$x.lab <- c(1:2,3)
n <- as.numeric(table(d$SelectionLineNumeric[d$SelectionLineNumeric!=0]))

# Fig. S8
pdf(file = "Plots/Bayes_Recruit.Fledgling.P3_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[1:2,]
plot(NULL, xlim=c(0.7,2.3), ylim=c(-3.4,2.1), xlab="",
     ylab="standardized fledgling P3",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-3,2, 1), labels=seq(-3,2, 1), lwd = 0, lwd.ticks = 1)
points(jitter(as.numeric(dat$L),0.5), dat$R, pch=1, cex=1.5, col=ifelse(dat$L==1, "darkorange", "darkslategrey"))
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.2, col="black", lwd=1.5)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
text(0.7, 2.1, "A", cex=1.2, pos=1, offset=0, font=2)

plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-0.4,1.1),
     xlab="", ylab="difference in standardized fledgling P3\n(late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(-0.3,0.9, 0.3), labels=seq(-0.3,0.9, 0.3), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025,1.1, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()


# run model for tarsus length:
d <- FINAL_chick.data
dat <- list(
  R=standardize(d$Tarsus),
  L=as.factor(d$SelectionLineNumeric),
  Y=as.factor(d$YearOfBreeding-2017),
  F=as.factor(females.as.interger),
  A=standardize(d$CaptureDateApril)
)

# run model:
mR.FT.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + x[F]*sigma_g + bL[L] + bA*A,
    sigma ~ dhalfnorm(0,1),
    z[Y] ~ dnorm(0, 1),
    x[F] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1), 
    sigma_a ~ dhalfnorm(0,1),
    sigma_g ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1),
    bA ~ dnorm(0, 1),
    gq> vector[Y]:a <<- a_bar+z*sigma_a,
    gq> vector[F]:g <<- x*sigma_g
  ), data=dat, chains=4, cores=4, iter=10000, log_lik=TRUE)
save(mR.FT.1a, file="temp/mR.FT.1a.RData")

# get model summary:
(out <- precis(mR.FT.1a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S24
write.table(out.0, "out/Bayes_Recruit.Fledgling.tarsus_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mR.FT.1a)), row.names(precis(mR.FT.1a, depth=2))[grepl("z", row.names(precis(mR.FT.1a, depth=2)))], row.names(precis(mR.FT.1a, depth=2))[grepl("bL", row.names(precis(mR.FT.1a, depth=2)))])
par(mfrow = c(1,1))
# Fig. S37
pdf(file = "Plots/Bayes_Recruit.Fledgling.tarsus_traceplots.pdf", width = 6, height = 8)
traceplot(mR.FT.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()
# Fig. S38
pdf(file = "Plots/Bayes_Recruit.Fledgling.tarsus_trankplot.pdf", width = 6, height = 8)
trankplot(mR.FT.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior prediction
# prepare data for prediction
dat.0 <- data.frame(dat$F, dat$Y)
dat_help.0 <- dat.0[!duplicated(dat.0),]
dat_help <- dat_help.0[order(dat_help.0$dat.F, dat_help.0$dat.Y),]
d_pred <- list(L=rep(1:2, each=nrow(dat_help)), Y=rep(dat_help$dat.Y, times=2), F=rep(dat_help$dat.F, times=2), A=rep(0, times=nrow(dat_help)*2))

# get posterior
p_post <- link(mR.FT.1a, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S25
write.table(post.out, "out/Bayes_Recruit.Fledgling.tarsus_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction for fledglings genetically derived from selection line breeding pairs
# prepare plot
post.out$x.lab <- c(1:2,3)
n <- as.numeric(table(d$SelectionLineNumeric[d$SelectionLineNumeric!=0]))

# Fig. S7
pdf(file = "Plots/Bayes_Recruit.Fledgling.tarsus_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[1:2,]
plot(NULL, xlim=c(0.7,2.3), ylim=c(-4.5,2.1), xlab="",
     ylab="standardized fledgling tarsus length",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-4,2, 1), labels=seq(-4,2, 1), lwd = 0, lwd.ticks = 1)
points(jitter(as.numeric(dat$L),0.5), dat$R, pch=1, cex=1.5, col=ifelse(dat$L==1, "darkorange", "darkslategrey"))
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.2, col="black", lwd=1.5)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
text(0.7, 2.1, "A", cex=1.2, pos=1, offset=0, font=2)

plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-0.2,1.1),
     xlab="", ylab="difference in standardized fledgling tarsus length\n(late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(0,0.9, 0.3), labels=seq(0,0.9, 0.3), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025,1.1, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()


### Feeding frequency
### ---------------------------------------------------------------------------

### load and prepare data
keep_FeedingfrequencyData_All <- read.table("/Users/melanielindner/Documents/NIOO/Projects/ERC_phenotypic_data/WIld_selection_lines/Feeding.frequency.txt", header=TRUE)
Breeding.pairs.transponder <- read.table("/Users/melanielindner/Documents/NIOO/Projects/ERC_phenotypic_data/WIld_selection_lines/Breeding.pairs.transponder.txt", header=TRUE)

# only females:
FF.females <- keep_FeedingfrequencyData_All[keep_FeedingfrequencyData_All$Sex==1,]
Bird.dat <- Breeding.pairs.transponder[Breeding.pairs.transponder$Parent=="Mother",]

FF.females.0 <- merge(FF.females, Bird.dat[,c(1:6,8)], by.x=c("Year", "UserPlaceName"), by.y=c("YearOfBreeding", "UserPlaceName"))

### Statistical analysis

## MODEL for daily feeding frequency (Eq. S4) of  female selection line recruits 

# prepare data for model:
dat <- list(
  R=standardize(FF.females.0$FeedingCounts),
  L=as.factor(FF.females.0$SelectionLineNumericMother),
  Y=as.factor(FF.females.0$Year-2017)
)

# fit model

mFF.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + bL[L],
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))

precis(mFF.1a, depth=2)
PSIS.out <- PSIS(mFF.1a, pointwise=TRUE)
FF.females.0[PSIS.out$k>0.5,] # clear outlier..

# remove outlier and rerun
FF.females.1 <- FF.females.0[!PSIS.out$k>0.5,]

# make input data:
dat <- list(
  R=standardize(FF.females.1$FeedingCounts),
  L=as.factor(FF.females.1$SelectionLineNumericMother),
  Y=as.factor(FF.females.1$Year-2017)
)

mFF.1b <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + bL[L],
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mFF.1b, file="temp/mFF.1b.multilevel.RData")

# get model summary:
(out <- precis(mFF.1b, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S30
write.table(out.0, "out/Bayes_FF_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mFF.1b)), row.names(precis(mFF.1b, depth=2))[grepl("z", row.names(precis(mFF.1b, depth=2)))], row.names(precis(mFF.1b, depth=2))[grepl("bL", row.names(precis(mFF.1b, depth=2)))])
par(mfrow = c(1,1))
# Fig. S43
pdf(file = "Plots/Bayes_FF_traceplots.pdf", width = 6, height = 8)
traceplot(mFF.1b, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()
# Fig. S44
pdf(file = "Plots/Bayes_FF_trankplot.pdf", width = 6, height = 8)
trankplot(mFF.1b, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior mean and ci
d_pred <- list(L=rep(1:2, each=3), Y=rep(1:3, times=2))
# get posterior
p_post <- link(mFF.1b, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S31
write.table(post.out, "out/Bayes_FF_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction for fledglings genetically derived from selection line breeding pairs
# prepare plot
post.out$x.lab <- c(1:3)
n <- as.numeric(table(FF.females.0[-10,]$SelectionLineNumericMother))

# Fig. S10
pdf(file = "Plots/Bayes_FF_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[1:2,]
dat1 <- FF.females.0[-10,]
dat1$FeedingCounts_std <- standardize(dat1$FeedingCounts)

plot(NULL, xlim=c(0.7,2.3), ylim=c(-1.6,2.1), xlab="",
     ylab="standardized feeding frequency",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-1.5,2.0, 0.5), labels=seq(-1.5,2.0, 0.5), lwd = 0, lwd.ticks = 1)
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
points(jitter(dat1$SelectionLineNumericMother,0.5), dat1$FeedingCounts_std, pch=1, cex=1.5, col=ifelse(dat1$SelectionLineNumericMother==1, "darkorange", "darkslategrey"))
text(0.7, 2.1, "A", cex=1.2, pos=1, offset=0, font=2)

plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-1,1),
     xlab="", ylab="difference in standardized\nfeeding frequency (late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(-1,1,0.4), labels=seq(-1,1,0.4), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025,1, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()


### Daily energy expenditure
### ---------------------------------------------------------------------------

### load and prepare data
DEE.all.years <- read.table("data/DEE.txt", header=TRUE)
Breeding.pairs.transponder <- read.table("/Users/melanielindner/Documents/NIOO/Projects/ERC_phenotypic_data/WIld_selection_lines/Breeding.pairs.transponder.txt", header=TRUE)

# only females:
DEE.females <- DEE.all.years[DEE.all.years$Sex=="Female",]
Bird.dat <- Breeding.pairs.transponder[Breeding.pairs.transponder$Parent=="Mother",]

DEE.females.0 <- merge(DEE.females[,c(3,14)], Bird.dat[,c(1:6,8)], by.x=c("RingNo"), by.y=c("RingNumber"))

### Statistical analysis

## MODEL for daily energy expenditure (Eq. S3) of  female selection line recruits 

# prepare data for model:
dat <- list(
  R=standardize(DEE.females.0$DEE_Speakman),
  L=as.factor(DEE.females.0$SelectionLineNumericMother),
  Y=as.factor(DEE.females.0$YearOfBreeding-2017)
)

mDEE.1a <- ulam(
  alist(
    R ~ dnorm(mu,sigma),
    mu <- a_bar + z[Y]*sigma_a + bL[L],
    z[Y] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5), 
    sigma_a ~ dhalfnorm(0,1),
    sigma ~ dhalfnorm(0,1),
    bL[L] ~ dnorm(0,1.5),
    gq> vector[Y]:a <<- a_bar+z*sigma_a
  ), data=dat, chains=4 , cores=4, iter=10000, log_lik=TRUE, control=list(adapt_delta=0.99))
save(mDEE.1a, file="temp/mDEE.1a.multilevel.RData")

# get model summary:
(out <- precis(mDEE.1a, depth=2))
pars <- row.names(out)#; pars <- sub("a\\[1\\]", "a\\[2017\\]", pars); pars <- sub("a\\[2\\]", "a\\[2018\\]", pars); pars <- sub("a\\[3\\]", "a\\[2019\\]", pars)
out.0 <- data.frame(parameter=pars, out)
# Table S28
write.table(out.0, "out/Bayes_DEE_precis.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# check whether convergence looks ok: trace- and trankplots for posterior distribution of parameters estimated
pars <- c(row.names(precis(mDEE.1a)), row.names(precis(mDEE.1a, depth=2))[grepl("z", row.names(precis(mDEE.1a, depth=2)))], row.names(precis(mDEE.1a, depth=2))[grepl("bL", row.names(precis(mDEE.1a, depth=2)))])
par(mfrow = c(1,1))
# Fig. S41
pdf(file = "Plots/Bayes_DEE_traceplots.pdf", width = 6, height = 8)
traceplot(mDEE.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()
# Fig. S42
pdf(file = "Plots/Bayes_DEE_trankplot.pdf", width = 6, height = 8)
trankplot(mDEE.1a, pars=pars,
          col=c(rangi2, "black", "gray60", "gray30"), n_cols=3)
dev.off()

# get posterior mean and ci
d_pred <- list(L=rep(1:2, each=3), Y=rep(1:3, times=2))
# get posterior
p_post <- link(mDEE.1a, data=d_pred)

# get mean for early and late selection line, respectively
new.temp <- NULL
for(i in 1:2) {
  dat.temp <- p_post[,d_pred$L==i]
  mu_over_broods <- apply(dat.temp, 1, mean)
  new.temp <- cbind(new.temp, mu_over_broods)
}
p_mu <- apply(new.temp, 2, mean)
p_ci <- apply(new.temp, 2, PI)
diff <- new.temp[,2]-new.temp[,1]
d_mu <- mean(diff); d_ci <- PI(diff)

post.out <- data.frame(line=c(1,2,"diff_e_l"), mu=c(p_mu, d_mu), ci.5.5=c(p_ci[1,],d_ci[1]), ci.94.5=c(p_ci[2,],d_ci[2]))
# Table S29
write.table(post.out, "out/Bayes_DEE_posterior_predictions.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# plot posterior prediction
# prepare plot
post.out$x.lab <- c(1:3)
n <- as.numeric(table(DEE.females.0$SelectionLineNumericMother))

# Fig. S9
pdf(file = "Plots/Bayes_DEE_Selection.lines.pdf", width = 7, height = 4)
layout(t(matrix(c(1,2))), width=c(1.5,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0))

plot.dat1 <- post.out[1:2,]
dat1 <- DEE.females.0
dat1$DEE_Speakman_std <- standardize(dat1$DEE_Speakman)

plot(NULL, xlim=c(0.7,2.3), ylim=c(-1.5,2), xlab="",
     ylab="standardized daily energy expenditure",
     xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=paste(c("early", "late"), " (n=", n, ")", sep=""), lwd = 0, lwd.ticks = 1)
axis(2, at=seq(-1.5,2.0, 0.5), labels=seq(-1.5,2.0, 0.5), lwd = 0, lwd.ticks = 1)
arrows(x0=plot.dat1$x.lab, y0=plot.dat1$ci.5.5, x1=plot.dat1$x.lab, y1=plot.dat1$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat1$x.lab, plot.dat1$mu, pch=19, cex=3, col=ifelse(plot.dat1$x.lab==1, "darkorange", "darkslategrey"))
points(jitter(dat1$SelectionLineNumericMother,0.5), dat1$DEE_Speakman_std, pch=1, cex=1.5, col=ifelse(dat1$SelectionLineNumericMother==1, "darkorange", "darkslategrey"))
text(0.7, 2, "A", cex=1.2, pos=1, offset=0, font=2)

plot.dat2 <- post.out[3,]
plot(NULL, xlim=c(2.7,3.3), ylim=c(-1.4,0.5),
     xlab="", ylab="difference in standardized\ndaily energy expenditure (late - early)",
     xaxt="n", yaxt="n")
axis(2, at=seq(-1.3,0.3,0.4), labels=seq(-1.3,0.3,0.4), lwd = 0, lwd.ticks = 1)
abline(h=0, lty=2)
arrows(x0=plot.dat2$x.lab, y0=plot.dat2$ci.5.5, x1=plot.dat2$x.lab, y1=plot.dat2$ci.94.5, code=3, angle=90, length=0.1, col="black", lwd=1)
points(plot.dat2$x.lab, plot.dat2$mu, pch=19, cex=3, col="black")
text(2.7+0.025,0.5, "B", cex=1.2, pos=1, offset=0, font=2)
dev.off()






