# Script to make predictions from the boral models
# These are mainly for reproducing Figure 2 and 3 in the main article

nsim <- 500

# generate counterfactual age gradient
Xv <- seq(min(R.nb.boral$X[,1]), max(R.nb.boral$X[,1]), length.out = nsim)
Xv.org <- Xv * X_attr$scale + X_attr$center
age.new <- cbind(age = Xv)

## Prediction with observed species' traits
R.nb.pred <-
    predict(
        R.nb.boral,
        newX = age.new,
        predict.type = "marginal",
        noroweff = TRUE,
        lv.mc = 1000,
        return.alllinpred = TRUE
    )

## Prediction with newtraits ("imaginary species")
# define "low" and "high" trait values
trait.select <- apply(scale(R.nb.boral$traits), 2, quantile, probs = c(0.1, 0.9))
# extract traits that strongly moderate species responses to age
sig.trait.R <- c("SDM", "H95", "LPC2")
## Single trait
# Counterfactual prediction with hypothetical species that vary by only one trait
# create a list of newtrait matrices, 
# one for each traits while holding other traits at zero
imaginary.spp.unitrait <- list()
for (i in seq_len(ncol(trait.select))) {
    imaginary.spp.unitrait[[i]] <- trait.select
    imaginary.spp.unitrait[[i]][, -i] <- 0   # set other traits to zero
}
names(imaginary.spp.unitrait) <- colnames(trait.select)
imaginary.spp.unitrait.R <- imaginary.spp.unitrait[sig.trait.R]
imaginary.spp.unitrait.R <- do.call(rbind, imaginary.spp.unitrait.R)
rownames(imaginary.spp.unitrait.R) <- NULL  # reset just in case
# now make single newtrait predictions
R.nb.unitrait.pred <-
    predict(
        R.nb.boral,
        newX = age.new,
        predict.type = "marginal",
        newtraits = imaginary.spp.unitrait.R,
        lvcoefs.type = "median",
        noroweff = TRUE,
        lv.mc = 1000,
        return.alllinpred = TRUE
    )
## Multiple traits
# similar to above, but each hypothetical species is a unique, factorial 
# combination of the high and low values of each traits
tmp <- trait.select; tmp[, !colnames(tmp)%in%sig.trait.R] <- 0
imaginary.spp.multtrait.R <- 
    expand.grid(WD = tmp[, "WD"],
                SDM = tmp[, "SDM"],
                H95 = tmp[, "H95"],
                LPC1 = tmp[, "LPC1"],
                LPC2 = tmp[, "LPC2"]) %>% 
    distinct(.keep_all = TRUE)
imaginary.spp.multtrait.R <- as.matrix(imaginary.spp.multtrait.R)
# add a ninth species with all traits at their mean
imaginary.spp.multtrait.R <- 
    rbind(imaginary.spp.multtrait.R, rep(0, ncol(imaginary.spp.multtrait.R)))
# now make multiple-newtrait predictions
R.nb.multtrait.pred <-
    predict(
        R.nb.boral,
        newX = age.new,
        predict.type = "marginal",
        newtraits = imaginary.spp.multtrait.R,
        lvcoefs.type = "median",
        noroweff = TRUE,
        lv.mc = 1000,
        return.alllinpred = TRUE
    )

## Ditto the above for mortality
S.pred <-
    predict.boral(
        S.boral,
        newX = age.new,
        predict.type = "marginal",
        noroweff = TRUE,
        lv.mc = 1000,
        return.alllinpred = TRUE
    )
sig.trait.S <- "LPC2"
imaginary.spp.unitrait.S <- imaginary.spp.unitrait[sig.trait.S]
imaginary.spp.unitrait.S <- do.call(rbind, imaginary.spp.unitrait.S)
rownames(imaginary.spp.unitrait.S) <- NULL
S.unitrait.pred <-
  predict(
    S.boral,
    newX = age.new,
    predict.type = "marginal",
    newtraits = imaginary.spp.unitrait.S,
    lvcoefs.type = "median",
    noroweff = TRUE,
    lv.mc = 1000,
    return.alllinpred = TRUE
  )
S.multtrait.pred <-
  predict(
    S.boral,
    newX = age.new,
    predict.type = "marginal",
    # for multitrait prediction we use the same trait mixture as recruitment
    # because we need the same set of traits for Figure 3
    newtraits = imaginary.spp.multtrait.R, 
    lvcoefs.type = "median",
    noroweff = TRUE,
    lv.mc = 1000,
    return.alllinpred = TRUE
  )


## Simulate population dynamics
# this is to reproduce Figure 3

set.seed(101)

R.mcmc <- get.mcmcsamples(R.nb.boral)
S.mcmc <- get.mcmcsamples(S.boral)

# some index for loops
age.per.step <- max(Xv.org) / n.step

# some index for loops
n.mcmc <- dim(R.nb.multtrait.pred$all.linpred)[3]
n.spp  <- dim(R.nb.multtrait.pred$all.linpred)[2]
n.step <- dim(R.nb.multtrait.pred$all.linpred)[1]

## Get the overdispersion parameters for neg. binom Recruitment model
get_phis <- matrix(R.mcmc[1, grep("lv.coefs", colnames(R.mcmc))], nrow = R.nb.boral$p)
get_phis <- get_phis[, ncol(get_phis)] 
get_phis_hypo <- median(get_phis)

## initial number of individuals per species
mean.log.N.init <- rep(0, n.spp)
N.init <- mean.log.N.init

# simulate Recruitment across forest age
# because Recruitment is independent of N (population size) in previous step
# we can simulate Recruitment across age at once
# (see below: mortality has to be simulated for each time step separately,
# because mortality depends on N in previous step)
all.simR <- R.nb.multtrait.pred$all.linpred
for (k in 1:n.mcmc) {
    for (j in 1:n.spp) {
        all.simR[, j, k] <-
            rnbinom(
                n.step,
                mu = exp(R.nb.multtrait.pred$all.linpred[, j, k]) * age.per.step,
                size = 1 / get_phis_hypo
            )
    }
}

# empty count array to fill in
N.pred.hypo <- array(0, dim(R.nb.multtrait.pred$all.linpred))

# begin simulation
for (k in seq_len(n.mcmc)) {
    N.pred.hypo[1, , k] <- N.init  # initial abundance
    for (i in 2:n.step) {
        # Recruitment had been simulated (see above)
        R.sim <- all.simR[i-1, , k]
        # Simulate mortality now, using previous-step N as trial size
        M.sim <- 
            rbinom(
                n.spp,
                size = N.pred.hypo[i-1, , k],
                prob = inv.logit.func(S.multtrait.pred$all.linpred[i-1, , k]) * age.per.step
            )
        N.pred.hypo[i, , k] <- N.pred.hypo[i-1, , k] - M.sim + R.sim
    }
}

# relative abundance
N.rel.pred <- N.pred.hypo
for (k in seq_len(dim(N.pred.hypo)[3])) {
    N.rel.pred[,,k] <- N.rel.pred[,,k] / rowSums(N.rel.pred[,,k])
}

# summarise to median and CIs
N.pred.summary.hypo <-
    list(
        median = apply(N.pred.hypo, 1:2, median),
        mean   = apply(N.pred.hypo, 1:2, mean),
        mean.log = apply(N.pred.hypo, 1:2, function(x) mean(log(x + 1))),
        upper  = apply(N.pred.hypo, 1:2, quantile, probs = 0.975),
        lower  = apply(N.pred.hypo, 1:2, quantile, probs = 0.025),
        median.rel = apply(N.rel.pred, 1:2, median, na.rm = TRUE),
        mean.rel   = apply(N.rel.pred, 1:2, mean, na.rm = TRUE),
        upper.rel  = apply(N.rel.pred, 1:2, quantile, probs = 0.975, na.rm = TRUE),
        lower.rel  = apply(N.rel.pred, 1:2, quantile, probs = 0.025, na.rm = TRUE)
    )
# paste species names to columns
# and convert to dataframes
N.pred.summary.hypo <- 
    lapply(N.pred.summary.hypo, function(x) {
        colnames(x) <- as.character(1:n.spp)
        return(x)
    }) %>% 
    reshape2::melt() %>% 
    spread(L1, value) %>% 
    rename(step = Var1,
           specid = Var2) %>% 
    left_join(
        reshape2::melt(age.new) %>% 
            select(-Var2) %>% 
            rename(step = Var1,
                   age = value)
    ) %>% 
    # backscale age
    mutate(age = age * X_attr$scale + X_attr$center,
           specid = rep(as.character(1:n.spp), n.step))  %>% 
    left_join(imaginary.spp.multtrait.R %>% 
                  as.data.frame %>% 
                  rownames_to_column("specid") %>% 
                  mutate(specid = as.character(specid)))
