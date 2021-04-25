# Script to reproduce figures in Lai et al. (2021) Successional syndromes 
# of saplings in tropical secondary forests emerge from environment-dependent 
# trait−demography relationships. Ecology Letters.

# Please contact Hao Ran Lai [hrlai.ecology@gmail.com] for issues
# Some codes are clunky and inconsistent in style as the project spanned a few
# years, apologies in advance!

# A special note for trait data:
# The leaf trait data provided here are the two principal components that went
# into the modelling process. For raw leaf traits, please visit another repository 
# (https://figshare.com/articles/Functional_Leaf_Traits_55_spp_in_central_Panama_/1402253)
# If you opt to use the raw leaf traits, please instead cite the other repository
# and references therein. 

# We believe that the sharing of datasets is important for advancing ecology. 
# At the same time, for data sharing to be successful and sustainable, it is 
# important that those individuals whose time and efforts generated the data 
# are acknowledged. Therefore, when you use the data or model output in your 
# original research or meta-analysis, we appreciate if the our work is 
# acknowledge :-) All the best.



# Load packages -----------------------------------------------------------

rm(list = ls())
library(boral)
library(mvtnorm)
library(tidyverse)
library(ggpubr)
library(viridis)
library(scales)

# function to backtransform logit to prob
inv.logit.func <- function(x) {1 / (1 + 1 / (exp(x)))}  

# modified prediction function from the boral package
source("R/modifiedpredictboral.R")



# Read JSDM models and housekeeping ----------------------------------------

# recruitment model
R.nb.boral <- readRDS(file = "data/recruitment_model.rds")

# sapling mortality model
S.boral    <- readRDS(file = "data/mortality_model.rds")

# also import dimension names for the input matrices of JSDM
dim_names  <- readRDS(file = "data/dim_names.rds")

# also import predictions made from the boral models above
# these files are mostly to reproduce the figures
# load(file = "R.nb.pred.Rdata")
# load(file = "S.pred.Rdata")
# load(file = "R.nb.unitrait.pred.Rdata")
# load(file = "S.unitrait.pred.Rdata")
# load(file = "R.nb.multtrait.pred.Rdata")
# load(file = "S.multtrait.pred.Rdata")
# load(file = "imaginary.spp.unitrait.R.Rdata")
# load(file = "imaginary.spp.unitrait.S.Rdata")
# load(file = "imaginary.spp.multtrait.R.Rdata")

# also load population simulation results
# load(file = "N.pred.summary.hypo.RData")

# also import mean and scale of age covariate and traits for backscaling
X_attr <- readRDS(file = "data/X_attr.rds")
Z_attr <- readRDS(file = "data/Z_attr.rds")

# put dimension names back to the boral model
rownames(R.nb.boral$traits.coefs.median)[-1] <- dim_names$covariate
colnames(R.nb.boral$traits.coefs.median)[-c(1,ncol(R.nb.boral$traits.coefs.median))] <- dim_names$trait
rownames(R.nb.boral$hpdintervals$traits.coefs)[-1] <- dim_names$covariate
colnames(R.nb.boral$hpdintervals$traits.coefs)[-c(1,ncol(R.nb.boral$hpdintervals$traits.coefs[,,1]))] <- dim_names$trait
dimnames(R.nb.boral$traits) <- list(dim_names$species, dim_names$trait)
dimnames(R.nb.boral$X) <- list(dim_names$sampling_unit, dim_names$covariate)
rownames(S.boral$traits.coefs.median)[-1] <- dim_names$covariate
colnames(S.boral$traits.coefs.median)[-c(1,ncol(S.boral$traits.coefs.median))] <- dim_names$trait
rownames(S.boral$hpdintervals$traits.coefs)[-1] <- dim_names$covariate
colnames(S.boral$hpdintervals$traits.coefs)[-c(1,ncol(S.boral$hpdintervals$traits.coefs[,,1]))] <- dim_names$trait
dimnames(S.boral$traits) <- list(dim_names$species, dim_names$trait)
dimnames(S.boral$X) <- list(dim_names$sampling_unit, dim_names$covariate)





# Example use of our model and data ---------------------------------------

# Here we show a few codes that extracts data from the boral objects
# Feel free to contact Hao Ran if you are unable to extract the desired data
# or parameters
# We are going to use the recruitment model as an example, but the following
# code should work for the sapling mortality model too by swapping
# the object names

# 1. Extract forest age data for plot
# showing first six rows
# rownames are concatenated by site ID, plot ID, and census ID
head(R.nb.boral$X)
# age column is scaled forest age, to unscale use our mean and standard deviation
R.nb.boral$X * X_attr$scale + X_attr$center

# 2. Extract trait data for species
head(R.nb.boral$traits)
# traits were scaled, to unscale use our mean and standard deviation
sweep(sweep(R.nb.boral$traits, 2, Z_attr$scale, "*"), 2, Z_attr$center, "+")

# 3. Extract the median estimate of fourth-corner coefficients
# the helpfile of boral provides lots of great tips, check out ?boral
R.nb.boral$traits.coefs.median





# Reproduce Figure 1 ------------------------------------------------------

# prepare a fourth-corner coefficient dataframe for plotting
R.nb.trait.coefs <- 
    as.data.frame(R.nb.boral$traits.coefs.median) %>% 
    rownames_to_column("Env") %>% 
    gather(Trait, Mean, -Env) %>% 
    left_join(
        as.data.frame(R.nb.boral$hpdintervals$traits.coefs[,,1]) %>% 
            rownames_to_column("Env") %>% 
            gather(Trait, LCI, -Env)
    ) %>% 
    left_join(
        as.data.frame(R.nb.boral$hpdintervals$traits.coefs[,,2]) %>% 
            rownames_to_column("Env") %>% 
            gather(Trait, UCI, -Env)
    ) %>% 
    filter(!Trait %in% c("sigma")) %>% 
    mutate(Trait = factor(Trait, 
                          levels = rev(c("kappa0","LPC1","LPC2","WD","SDM","H95"))),
           Env = factor(Env, 
                        levels = c("beta0","age"), 
                        labels = c("italic(beta)[0]","italic(beta)[Age]")),
           Sig = sign(LCI)==sign(UCI))
S.trait.coefs <- 
    as.data.frame(S.boral$traits.coefs.median) %>% 
    rownames_to_column("Env") %>% 
    gather(Trait, Mean, -Env) %>% 
    left_join(
        as.data.frame(S.boral$hpdintervals$traits.coefs[,,1]) %>% 
            rownames_to_column("Env") %>% 
            gather(Trait, LCI, -Env)
    ) %>% 
    left_join(
        as.data.frame(S.boral$hpdintervals$traits.coefs[,,2]) %>% 
            rownames_to_column("Env") %>% 
            gather(Trait, UCI, -Env)
    ) %>% 
    filter(!Trait %in% c("sigma")) %>% 
    mutate(Trait = factor(Trait, 
                          levels = rev(c("kappa0","LPC1","LPC2","WD","SDM","H95"))),
           Env = factor(Env, 
                        levels = c("beta0","age"), 
                        labels = c("italic(beta)[0]","italic(beta)[Age]")),
           Sig = sign(LCI)==sign(UCI))

# make figure 1!
Fig1theme <- 
    theme(plot.title = element_text(size=12),
          legend.position="none",
          axis.text.y = element_text(colour = "black"), 
          plot.margin = margin(2,2,2,-12))
y.labels <- c('kappa0' = expression(italic(T)[0]),
              'LPC1'   = "LPC1",
              'LPC2'   = "LPC2",
              'WD'     = "WD",
              'SDM'    = "SDM",
              'H95'    = expression(H[max]))
Fig1a <-
    ggplot(R.nb.trait.coefs, aes(Mean, Trait)) +
    facet_wrap(~ Env, nrow = 1, scale = "free_x", labeller = label_parsed) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = LCI, xmax = UCI, colour = Sig), height=0) +
    geom_point(aes(colour = Sig)) +
    labs(y = "", x = "Effect sizes", title = "Recruitment rate") +
    scale_colour_grey(start=0.6, end=0) +
    scale_y_discrete(labels = y.labels) +
    theme_classic() + Fig1theme
Fig1b <-
    ggplot(S.trait.coefs, aes(Mean, Trait)) +
    facet_wrap(~ Env, nrow = 1, scale = "free_x", labeller = label_parsed) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = LCI, xmax = UCI, colour = Sig), height=0) +
    geom_point(aes(colour = Sig)) +
    labs(y = "", x = "Effect sizes", title = "Mortality rate") +
    scale_colour_grey(start=0.6, end=0) +
    scale_y_discrete(labels = y.labels) +
    theme_classic() + Fig1theme
ggarrange(Fig1a, Fig1b,
          labels = c("(a)", "(b)"),
          hjust = -.5, vjust = 1,
          font.label = list(size = 12),
          ncol = 2, nrow = 1)




# Reproduce Figure 2 ------------------------------------------------------

# run script that makes predictions from the boral models
# this will take a while, go for a walk and relax!
source("R/prepare_predictions.R")

# some indices
n.mcmc <- dim(R.nb.unitrait.pred$all.linpred)[3]
n.spp  <- dim(R.nb.unitrait.pred$all.linpred)[2]
n.spp.obs  <- dim(R.nb.pred$all.linpred)[2]
n.step <- dim(R.nb.unitrait.pred$all.linpred)[1]

# generate counterfactual age gradient
Xv <- seq(min(R.nb.boral$X[,1]), max(R.nb.boral$X[,1]), length.out = n.step)
Xv.org <- Xv * X_attr$scale + X_attr$center

# traits
trait <- R.nb.boral$traits %>% as.data.frame %>% rownames_to_column("specid")

# tidy prediction of hypothetical species
# first convert predictions from link to response scale
R.nb.hypo.respred <- exp(R.nb.unitrait.pred$all.linpred)
R.nb.hypo.respred <- 
    apply(R.nb.hypo.respred, 1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
    reshape::melt(varnames = c("probs", "age.step", "specid")) %>% 
    spread(probs, value) %>% 
    mutate(specid = as.character(specid),
           age = rep(Xv.org, each = nrow(imaginary.spp.unitrait.R))) %>% 
    rename(median = `50%`,
           lower  = `2.5%`,
           upper  = `97.5%`) %>% 
    left_join(imaginary.spp.unitrait.R %>% 
                  as.data.frame %>% 
                  rownames_to_column("specid") %>% 
                  mutate(specid = as.character(specid))) %>% 
    # trick to label trait
    gather(Trait, TraitVal, WD:LPC2) %>% 
    filter(TraitVal != 0)
# do the same for mortality rate
S.hypo.respred <- inv.logit.func(S.unitrait.pred$all.linpred)
S.hypo.respred <- 
    apply(S.hypo.respred, 1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
    reshape::melt(varnames = c("probs", "age.step", "specid")) %>% 
    spread(probs, value) %>% 
    mutate(specid = as.character(specid),
           age = rep(Xv.org, each = nrow(imaginary.spp.unitrait.S))) %>% 
    rename(median = `50%`,
           lower  = `2.5%`,
           upper  = `97.5%`) %>% 
    left_join(imaginary.spp.unitrait.S %>% 
                  as.data.frame %>% 
                  rownames_to_column("specid") %>% 
                  mutate(specid = as.character(specid))) %>% 
    # trick to label trait
    gather(Trait, TraitVal, WD:LPC2) %>% 
    filter(TraitVal != 0)
# tidy predictions for multitrait hypothetical species
R.nb.hypo.mult.respred <- exp(R.nb.multtrait.pred$all.linpred)
R.nb.hypo.mult.respred <- 
    apply(R.nb.hypo.mult.respred, 1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
    reshape::melt(varnames = c("probs", "age.step", "specid")) %>% 
    spread(probs, value) %>% 
    mutate(specid = as.character(specid),
           age = rep(Xv.org, each = nrow(imaginary.spp.multtrait.R))) %>% 
    rename(median = `50%`,
           lower  = `2.5%`,
           upper  = `97.5%`) %>% 
    left_join(imaginary.spp.multtrait.R %>% 
                  as.data.frame %>% 
                  rownames_to_column("specid") %>% 
                  mutate(specid = as.character(specid)))
S.hypo.mult.respred <- inv.logit.func(S.multtrait.pred$all.linpred)
S.hypo.mult.respred <- 
    apply(S.hypo.mult.respred, 1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
    reshape::melt(varnames = c("probs", "age.step", "specid")) %>% 
    spread(probs, value) %>% 
    mutate(specid = as.character(specid),
           age = rep(Xv.org, each = nrow(imaginary.spp.multtrait.R))) %>% 
    rename(median = `50%`,
           lower  = `2.5%`,
           upper  = `97.5%`) %>% 
    left_join(imaginary.spp.multtrait.R %>% 
                  as.data.frame %>% 
                  rownames_to_column("specid") %>% 
                  mutate(specid = as.character(specid)))

# Plot figure 2!
Fig2theme <- 
    theme(legend.position = "none",
          axis.ticks.length=unit(-1, "mm"),
          axis.text.x = element_text(size = 8, margin = margin(t = unit(5, "mm"))),
          axis.text.y = element_text(size = 8, margin = margin(r = unit(5, "mm"))),
          axis.title = element_text(size = 9),
          plot.margin = margin(2,2,2,2),
          plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))
Fig2ylim <- range(as.matrix(R.nb.hypo.respred[, c("lower","upper")]))

Fig2b <- 
    ggplot(R.nb.hypo.respred %>% filter(Trait == "LPC2"), 
           aes(age, median, group = specid)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = TraitVal), 
                alpha = 0.2) +
    geom_line(aes(colour = TraitVal), size = 1) +
    scale_colour_viridis(option = "inferno", begin = 0.3, end = 0.7, 
                         aesthetics = c("colour", "fill")) +
    coord_cartesian(ylim = Fig2ylim) +
    scale_y_log10(label = trans_format("log10", math_format(10^.x))) +
    labs(x = "", y = "", title = "Leaf PC2") +
    theme_classic() + Fig2theme
Fig2c <- 
    ggplot(R.nb.hypo.respred %>% filter(Trait == "H95"), 
           aes(age, median, group = specid)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                alpha = 0.15) +
    geom_line(aes(alpha = as.factor(TraitVal)), 
              size = 1) +
    scale_alpha_manual(values = c(0.4, 1)) +
    coord_cartesian(ylim = Fig2ylim) +
    scale_y_log10(label = trans_format("log10", math_format(10^.x))) +
    labs(x = "Forest age (yr)", y = "", title = expression(bold(H[max]))) +
    theme_classic() + Fig2theme
Fig2d <- 
    ggplot(R.nb.hypo.respred %>% filter(Trait == "SDM"), 
           aes(age, median, group = specid)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                alpha = 0.15) +
    geom_line(aes(linetype = as.factor(TraitVal)), size = 1, lineend = "round") +
    scale_linetype_manual(values = c("12", "solid")) +
    coord_cartesian(ylim = Fig2ylim) +
    scale_y_log10(label = trans_format("log10", math_format(10^.x))) +
    labs(x = "Forest age (yr)", y = "", title = "log SDM") +
    theme_classic() + Fig2theme
Fig2_R_mult <- 
    ggplot() +
    geom_line(data = R.nb.hypo.mult.respred %>% filter(LPC2 == 0), 
              aes(age, median),
              size = 0.3, colour = "grey", linetype = 4) +
    geom_line(data = R.nb.hypo.mult.respred %>% filter(LPC2 != 0), 
              aes(age, median, 
                  group = specid, 
                  colour = LPC2,
                  linetype = as.factor(SDM),
                  alpha = as.factor(H95)),
              size = 1, lineend = "round") +
    scale_colour_viridis(option = "inferno", begin = 0.3, end = 0.7) +
    scale_linetype_manual(values = c("12", "solid")) +
    scale_alpha_manual(values = c(0.3, 1)) +
    coord_cartesian(ylim = Fig2ylim) +
    scale_y_log10(label = trans_format("log10", math_format(10^.x))) +
    labs(x = "", 
         y = expression(paste("Recruitment rate (count ", yr^-1, ")")), 
         title = "Trait mixtures") +
    theme_classic() + Fig2theme

Fig2f <- 
    ggplot(S.hypo.respred %>% filter(Trait == "LPC2"), 
           aes(age, median, group = specid)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = TraitVal), 
                alpha = 0.2) +
    geom_line(aes(colour = TraitVal), size = 1) +
    scale_colour_viridis(option = "inferno", begin = 0.3, end = 0.7, 
                         aesthetics = c("colour", "fill")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Forest age (yr)", y = "", title = "Leaf PC2") +
    theme_classic() + Fig2theme
Fig2_S_mult <- 
    ggplot() +
    geom_line(data = S.hypo.mult.respred %>% filter(LPC2 == 0), 
              aes(age, median),
              size = 0.3, colour = "grey", linetype = 4) +
    geom_line(data = S.hypo.mult.respred %>% filter(LPC2 != 0), 
              aes(age, median, 
                  group = specid, 
                  colour = LPC2,
                  linetype = as.factor(SDM),
                  alpha = as.factor(H95)),
              size = 1, lineend = "round") +
    scale_colour_viridis(option = "inferno", begin = 0.3, end = 0.7) +
    scale_linetype_manual(values = c("12", "solid")) +
    scale_alpha_manual(values = c(0.3, 1)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Forest age (yr)",
         y = expression(paste("Mortality rate (% ", yr^-1, ")")),
         title = "Trait mixtures") +
    theme_classic() + Fig2theme

Fig2_legend <-
    ggplot() +
    coord_fixed(xlim = c(-0.5, 2), ylim = c(0.2, 2)) +
    geom_text(aes(x = rep(0, 3), y = c(2, 1.3, 0.7), 
                  label = c("Colour",
                            "Colour intensity",
                            "Line type")),
              size = 3, hjust = 0, fontface = "bold") +
    geom_ribbon(aes(x = c(0, 0.2), ymin = c(1.7, 1.7), ymax = c(1.9, 1.9)),
                fill = inferno(1, 1, 0.3)) +
    geom_ribbon(aes(x = c(0, 0.2), ymin = c(1.4, 1.4), ymax = c(1.6, 1.6)),
                fill = inferno(1, 1, 0.7)) +
    geom_segment(aes(x = 0, y = 1.1, xend = 0.2, yend = 1.1), size = 1) +
    geom_segment(aes(x = 0, y = 0.9, xend = 0.2, yend = 0.9), size = 1, colour = "grey") +
    geom_segment(aes(x = 0, y = 0.5, xend = 0.2, yend = 0.5), size = 1) +
    geom_segment(aes(x = 0, y = 0.3, xend = 0.2, yend = 0.3), size = 1, linetype = "12") +
    geom_text(aes(x = rep(0.25, 6), y = c(0.3,0.5,0.9,1.1,1.5,1.8),
                  label = c("Light SDM",
                            "Heavy SDM",
                            "Short Hmax",
                            "Tall Hmax",
                            "High leaf PC2",
                            "Low leaf PC2")),
              size = 3, hjust = 0) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
ggarrange(Fig2_R_mult, Fig2b, Fig2c, Fig2d, 
          Fig2_S_mult, Fig2f, Fig2_legend,
          labels = paste0("(", letters[1:6], ")"),
          hjust = -1, vjust = 1,
          font.label = list(size = 12),
          ncol = 4, nrow = 2)




# Reproduce Figure 3 ------------------------------------------------------

fig3a <- 
    ggplot(N.pred.summary.hypo %>% filter(LPC2 == min(LPC2)), 
           aes(age, exp(mean.log)-1, group = specid)) +
    facet_wrap(~ factor(LPC2, labels = "bold('(a)')~Low~Leaf~PC2~(low~A[area])"),
               labeller = label_parsed) +
    geom_line(data = N.pred.summary.hypo %>% filter(LPC2 == 0), 
              aes(age, exp(mean.log)-1),
              colour = "grey", linetype = 4) +
    geom_line(aes(alpha = factor(H95, labels = c("Short", "Tall")),
                  linetype = factor(SDM, labels = c("Light", "Heavy"))),
              colour = inferno(1, begin = 0.3, end = 0.3),
              size = 1, lineend = "round") +
    annotate("text", x = 31, y = 10, label = "'Late-\nsuccessional'", 
             angle = 30, vjust = 0, colour = inferno(1, begin = 0.3, end = 0.3)) +
    scale_alpha_manual(values = c(0.3, 1), name = expression(H[max])) +
    scale_linetype_manual(values = c("12", "solid")) +
    scale_y_log10() +
    coord_cartesian(ylim = c(0.015, 30)) +
    labs(x = "Forest age (yr)", 
         y = expression(paste("Simulated mean plot abundance, ", italic(hat(N)))),
         colour = expression(H[max]),
         linetype = "log SDM") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.margin = margin(0,10,0,10),
          strip.text = element_text(size = 10, hjust = 0),
          strip.background = element_blank(),
          legend.text = element_text(size = 10))
fig3b <- 
    ggplot(N.pred.summary.hypo %>% filter(LPC2 == max(LPC2)), 
           aes(age, exp(mean.log)-1, group = specid)) +
    facet_wrap(~ factor(LPC2, labels = "bold('(b)')~High~Leaf~PC2~(high~A[area])"),
               labeller = label_parsed) +
    geom_line(data = N.pred.summary.hypo %>% filter(LPC2 == 0), 
              aes(age, exp(mean.log)-1),
              colour = "grey", linetype = 4) +
    geom_line(aes(alpha = factor(H95, labels = c("Short", "Tall")),
                  linetype = factor(SDM, labels = c("Light", "Heavy"))),
              colour = inferno(1, begin = 0.7, end = 0.7),
              size = 1, lineend = "round") +
    annotate("text", x = 30, y = 0.1, label = "'Pioneers'", 
             colour = inferno(1, begin = 0.7, end = 0.7, alpha = 0.4), angle = -30) +
    annotate("text", x = 30, y = 1, label = "'Long-lived pioneers'",
             angle = -15, colour = inferno(1, begin = 0.7, end = 0.7)) +
    scale_alpha_manual(values = c(0.3, 1), name = expression(H[max])) +
    scale_linetype_manual(values = c("12", "solid")) +
    scale_y_log10() +
    coord_cartesian(ylim = c(0.015, 30)) +
    labs(x = "Forest age (yr)", y = "",
         colour = expression(H[max]),
         linetype = "log SDM") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.margin = margin(0,10,0,10),
          strip.text = element_text(size = 10, hjust = 0),
          strip.background = element_blank(),
          legend.text = element_text(size = 10))
ggarrange(fig3a, fig3b,
          ncol = 2,
          common.legend = TRUE, legend = "bottom")
