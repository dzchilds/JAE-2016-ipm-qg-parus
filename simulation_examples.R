
library(ggplot2)
library(tidyr)
library(dplyr)

source("demog_functions.R")

load(file="model_param.rda")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Example simulations + some summary figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---------------------------------------------------#
# Example Part 1 -- Set up various control variables #
# ---------------------------------------------------#

# set up control parameters for all simulations
iPar <- 
  list(minx=-30, maxx=+30, msizex=50, # l-state mesh (residual)
       ming=-30, maxg=+30, msizeg=50, # q-state mesh (just 'g')
       maxA=8)                        # max age allowed

# !!!! must use same upper/lower bound for 'g' (i.e. zero centred) in the above !!!!

# construct everything we need to numericall implement the model
cPar <- calcControlParams(iPar)
# copy the stored model parameter into a new matrix
mPar <- mParStore
# number of years in initial run in period
nStepInt <- 30
# number of simulation years
nStepSim <- 200
# set the secular trend effect (assume it is fixed for simplicity)
yTrendSim <- 2010
# number of replicate simulations
nSim <- 250
# lists to store the state distributions of interest -- 'g'enotype and 'p'henotype
gSim <- pSim <- list()
# progress bar
pb <- txtProgressBar(width = 25)

# ---------------------------------------------#
# Example Part 2 -- Carry out the simulations  #
# ---------------------------------------------#

# loop over simulations
for (i in 1:nSim) {

  #
  # Simulation Part 1 -- Run in phase
  #
  
  # uncouple parent-offpsring genotype 
  getOffgDist <- getOffgDistNoE
  # random sample of years (iid assumption)
  ySet <- sample(seq.int(1, ncol(mPar)), size = nStepInt, replace = TRUE)
  # construct the environment (year type, secular trend effect, period)
  eSet <- list(ySet = ySet, yTrend = 2010 , period = 0)
  # initial state distribution
  nt0 <- makeInitState(cPar, mPar, Nt = 240, sigmaG = 0.5)
  # simulate dynamics without evolution
  simSys <- iterateEcoEvoDynamics(mPar, cPar, eSet, nt0, keepAll = FALSE)
  # keep the final state for use in the next step
  nt0 <- list(F = simSys$F[[nStepInt]], M = simSys$M[[nStepInt]])
  
  #
  # Simulation Part 2 -- Simulate 'eco-evo' dynamics 
  #

  # couple parent-offpsring genotype (i.e. let the system evolve) 
  getOffgDist <- getOffgDistEvo
  # random sample of years (iid assumption)
  ySet <- sample(seq.int(1, ncol(mPar)), size = nStepSim, replace = TRUE)
  # construct the environment (year type, secular trend effect, period)
  eSet <- list(ySet = ySet, yTrend = yTrendSim, period = 0)
  # simulate dynamics with evolution
  simSys <- iterateEcoEvoDynamics(mPar, cPar, eSet, nt0, keepAll = FALSE)
  
  #
  # Simulation Part 3a -- Calculate ASPE components for genotype 
  #
  
  # single sex terms
  ASPE <- calcPriceDecompTerms(cPar, simSys, vrType = "g")
  # two sex terms (just a reweighting)
  ts.ASPE <- lapply(ASPE, function(x) calcPrice2Sx(x)) 
  # housekeeping to get it all into one data frame
  var <- expand.grid(Sim.No = i, Age = seq.int(cPar$maxA), Sex = c("F", "M"), Time = seq.int(nStepSim))
  dat <- lapply(ts.ASPE, function(lst) bind_rows(lst$F, lst$M))
  dat <- bind_rows(dat)
  gSim[[i]] <- bind_cols(var, dat)
  
  #
  # Simulation Part 3b -- Calculate ASPE components for phenotype 
  #
  
  # single sex terms
  ASPE <- calcPriceDecompTerms(cPar, simSys, vrType = "p")
  # two sex terms (just a reweighting)
  ts.ASPE <- lapply(ASPE, function(x) calcPrice2Sx(x)) 
  # housekeeping to get it all into one data frame
  var <- expand.grid(Sim.No = i, Age = seq.int(cPar$maxA), Sex = c("F", "M"), Time = seq.int(nStepSim))
  dat <- lapply(ts.ASPE, function(lst) bind_rows(lst$F, lst$M))
  dat <- bind_rows(dat)
  pSim[[i]] <- bind_cols(var, dat)
  #
  setTxtProgressBar(pb, i/nSim)
}
# bind everything together into one big data frame 
gSim <- bind_rows(gSim)
pSim <- bind_rows(pSim)

# optionally save the results (simulations take ~ 1-2 hours)
# save(gSim, pSim, file = "sim_out.rda")

# -----------------------------------------#
# Example Part 3 -- Visualise the results  #
# -----------------------------------------#

# vector containing nicer labels
name.map <- c(SD.S = "Survival Selection", SD.F = "Recruitment Selection",
              ON = "Ontogeny", PO = "Inheritance", AF = "Demography")

#
# Trends, summing over ages classes
# -> density, breading value, age structure 

sim.dens <- 
  gSim %>% filter(Sex == "F", Time <= 500) %>% 
  select(Sim.No, Time, N) %>% distinct

plt.mean.dens <- 
  sim.dens %>% 
  group_by(Time) %>% 
  summarise(N = mean(N))

plt.dens <- 
  sim.dens %>% 
  filter(Sim.No == 1)

ggplot(plt.dens, aes(x = Time, y = N)) +
  geom_line(alpha = 0.2, aes(group = Sim.No)) + 
  geom_smooth(data = plt.mean.dens,
              method = "gam", formula = y ~ s(x), se = FALSE, gamma = 10) + 
  theme_bw() + xlab("Time (years)") + ylab("Breeding Pair Density")

#
# Trends, summing over ages classes
# -> mean breading value

sim.mu.g <- 
  gSim %>% 
  filter(Sex == "F", Time <= 500) %>% 
  group_by(Sim.No, Time) %>% 
  summarise(mu = sum(c.A * Z.1))

plt.mean.mu.g <- sim.mu.g %>% group_by(Time) %>% summarise(mu = mean(mu))
plt.mu.g <- sim.mu.g %>% filter(Sim.No <= 10)

ggplot(plt.mu.g, aes(x = Time, y = mu)) +
  geom_line(alpha = 0.5, aes(group = Sim.No)) + 
  geom_smooth(data = plt.mean.mu.g,
              method = "gam", formula = y ~ s(x), se = FALSE, gamma = 10) + 
  theme_bw() + xlab("Time (years)") + ylab("Mean Breeding Value")

#
# Trends, summing over ages classes
# -> age structure 

sim.ages <- 
  gSim %>% 
  filter(Sex == "F", Time <= 500)

plt.mean.age <- 
  sim.ages %>% 
  group_by(Time, Age) %>% 
  summarise(pA = mean(c.A))

ggplot(plt.mean.age, aes(x = Time, y = pA, group = Age,)) +
  geom_smooth(aes(colour = Age), method = "gam", 
              formula = y ~ s(x), se = FALSE, gamma = 10)
  
#
# Trends, summing over ages classes
# -> breeding values and phenotype

# Extract trends in ASPE terms -> G

ASPE.terms <- 
  gSim %>% 
  group_by(Sim.No, Sex, Time) %>% 
  summarise(SD.S = sum(SD.S), SD.F = sum(SD.F), 
            ON = sum(ON), PO = sum(PO), AF = sum(AF.F + AF.1 + AF.0)) %>% 
  gather(Term, Contrib, SD.S:AF) %>% 
  filter(Sex == "F", Term != "ON")
# mean trend
plt.means.g <- ASPE.terms %>% group_by(Time, Term) %>% summarise(Contrib= mean(Contrib))
# individual simulation
plt.1sim.g <- ASPE.terms %>% filter(Sim.No == 1, Time %% 10 == 0, 
                                    Contrib > -0.1, Contrib < +0.1)
rng <- max(abs(range(plt.1sim.g$Contrib)))
plt.1sim.g <- mutate(plt.1sim.g, 
                     Contrib = ifelse(Contrib == min(Contrib), -rng, Contrib),
                     Contrib = ifelse(Contrib == max(Contrib), +rng, Contrib))

# Extract trends in ASPE terms -> P

ASPE.terms <- 
  pSim %>% 
  group_by(Sim.No, Sex, Time) %>% 
  summarise(SD.S = sum(SD.S), SD.F = sum(SD.F), 
            ON = sum(ON), PO = sum(PO), AF = sum(AF.F + AF.1 + AF.0)) %>% 
  gather(Term, Contrib, SD.S:AF) %>% 
  filter(Sex == "F")
# mean trend
plt.means.p <- ASPE.terms %>% group_by(Time, Term) %>% summarise(Contrib= mean(Contrib))
# individual simulation
plt.1sim.p <- ASPE.terms %>%  filter(Sim.No == 1, Time %% 10 == 0,
                                     Contrib > -4, Contrib < +4) # remove extreme values for plots
#
rng <- max(abs(range(plt.1sim.p$Contrib)))
plt.1sim.p <- mutate(plt.1sim.p, 
                     Contrib = ifelse(Contrib == min(Contrib), -rng, Contrib),
                     Contrib = ifelse(Contrib == max(Contrib), +rng, Contrib))
  
# Plot trends in ASPE terms -> P

#
plt.means <- bind_rows(mutate(Trait = "Phenotype",       plt.means.p),
                       mutate(Trait = "Breeding Value", plt.means.g))
#
plt.1sim  <- bind_rows(mutate(Trait = "Phenotype",       plt.1sim.p),
                       mutate(Trait = "Breeding Value", plt.1sim.g))
#
plt.1sim  <- mutate(plt.1sim,  Term = name.map[Term])
plt.means <- mutate(plt.means, Term = name.map[Term])
# 
ggplot(plt.1sim, aes(x = Time, y = Contrib, colour = Term, shape = Term)) +
  geom_point(alpha = 0.5, size = 1.8) + 
  geom_smooth(data = plt.means, method = "gam", formula = y ~ s(x), se = FALSE) + 
  theme_bw() + xlab("Time (years)") + ylab("Contribution") +
  facet_wrap(~Trait, ncol = 2, scales = "free_y")

