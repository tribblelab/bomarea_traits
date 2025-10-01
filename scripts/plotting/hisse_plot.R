library(coda)
library(phytools)
library(ggplot2)
library(dplyr)
library(RevGadgets)

# Make sure you are where your files are supposed to be
setwd("~/Desktop/bomarea_traits/output/")

# With coda
# mcmc_run1 <- readTrace(path = "hisse_binary_type_run_1.log", burnin = 0.1)
# mcmc_trace <- as.mcmc(mcmc_run1[[1]])
# traceplot(mcmc_trace)
# effectiveSize(mcmc_trace)

# If you run into issues use it as below
# mcmc_trace[[1]] <- as.mcmc(mcmc_trace[[1]]

# ggplot+revgadgets
mcmc_run1 <- readTrace(path = "hisse_run_1.log", burnin = 0.1)
mcmc_run1 <- data.frame(mcmc_run1)
mcmc_run1 <- cbind(mcmc_run1, run = rep("run 1", length(mcmc_run1$Iteration)))

mcmc_run2 <- readTrace(path = "hisse.log", burnin = 0.1)
mcmc_run2 <- data.frame(mcmc_run2)
mcmc_run2 <- cbind(mcmc_run2, run = rep("run 2", length(mcmc_run2$Iteration)))

mcmc_table <- rbind(mcmc_run1, mcmc_run2)

trace_plot <- ggplot(mcmc_table, aes(x = Iteration,y = Posterior, group = run)) +
              geom_line(aes(color = run)) +
              theme_classic()
trace_plot

# Make sure you have cut sufficiently the burn-in

# ggplot2

traitcols <- c("#3D348B", "#7678ED", "#F18701", "#F35B04")

hisse <- read.table("hisse_run_1.log", header = TRUE)
hisse <- hisse[-seq(1, 400000, 1), ] # make sure you are cutting the burn in!

transition_rates <- data.frame(dens = c(hisse$q_01A, hisse$q_01B,
                                        hisse$q_10A, hisse$q_10B),
                               rate = rep(c("q_01A","q_01B", "q_10A", "q_10B"),
                               each = length(hisse$q_01A)))

violin_transitions <- ggplot(transition_rates,
                             aes(x = rate, y = dens, fill = rate)) +
  geom_violin(trim = FALSE) +
  labs(title = "Transition Rates") +
  scale_fill_manual(values = traitcols) +
  xlab("Transition type") +
  ylab("Rate") +
  theme_classic()
violin_transitions
ggsave("hisse_results/hisse_violin_transitions.png")

## Hidden state transitions

hidden_rates <- data.frame(dens = c(hisse$hidden_rate1,
                                    hisse$hidden_rate2),
                           rate = rep(c("alpha", "beta"),
                           each = length(hisse$hidden_rate1)))

violin_hidden <- ggplot(hidden_rates, aes(x = rate, y = dens, fill = rate)) +
  geom_violin(trim = FALSE) +
  labs(title = "Hidden state transitions ") +
  scale_fill_manual(values = traitcols[1:2]) +
  xlab("Transition type") +
  ylab("Rate") +
  theme_classic()
violin_hidden
ggsave("hisse_results/hisse_violin_hidden.png")

divcols <- c("#E63946", "#F3A5AB", "#1D3557", "#457B9D")

# In RevBayes 1=0A, 2=1A, 3=0B, and 4=1B
netdiversification_rates<- data.frame(dens=c(hisse$speciation.1.-hisse$extinction.1.,
                                             hisse$speciation.2.-hisse$extinction.2., 
                                             hisse$speciation.3.-hisse$extinction.3.,
                                             hisse$speciation.4.-hisse$extinction.4.),
                                      rate=rep(c("net_div_0A","net_div_1A","net_div_0B","net_div_1B"),
                                      each=length(hisse$speciation.1.)))

violin_diversification <- ggplot(netdiversification_rates, aes(x = rate,y = dens, fill = rate))+
  geom_violin(trim = FALSE) +
  labs(title = "Rates of Transition") +
  scale_fill_manual(values = divcols) +
  xlab("") +
  ylab("Rate") +
  theme_classic()
violin_diversification
ggsave("hisse_results/hisse_violin_diversification.png")

difcols <- c("#FF006E","#FFC2DC")
T_diff <- data.frame(dens=c((hisse$speciation.1.-hisse$extinction.1.)-(hisse$speciation.2.-hisse$extinction.2.),
                            (hisse$speciation.3.-hisse$extinction.3.)-(hisse$speciation.4.-hisse$extinction.4.)),
                     difference=rep(c("T_A","T_B"),
                     each=length(hisse$speciation.1.)))

violin_difference <- ggplot(T_diff, aes(x = difference,y = dens, fill = difference)) +
  geom_violin(trim = FALSE) +
  labs(title = "Test statistics") +
  scale_fill_manual(values = difcols) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 1) +
  xlab("Test statistic") +
  ylab("Difference") +
  theme_classic()

violin_difference
ggsave("hisse_results/hisse_violin_difference.png")
# 0 crosses these posterior distributions, but we have to check the probability.



# dplyr

quantile_diff <- T_diff %>%
  group_by(difference) %>%
  reframe(res = quantile(dens, probs = c(0.025, 0.975)))
quantile_diff

#ggplot2
difcols <- c("#1DB32C", "#BFEEC3")
T_diff <- data.frame(dens = c((hisse$speciation.1.-hisse$extinction.1.) -
                              (hisse$speciation.3.-hisse$extinction.3.),
                              (hisse$speciation.2.-hisse$extinction.2.) -
                              (hisse$speciation.4.-hisse$extinction.4.)),
                     difference=rep(c("T_0","T_1"),
                     each = length(hisse$speciation.1.)))

violin_difference <- ggplot(T_diff, aes(x = difference, y = dens, fill = difference)) +
  geom_violin(trim = FALSE) +
  labs(title = "Test statistics") +
  scale_fill_manual(values = difcols) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 1) +
  xlab("Test statistic") +
  ylab("Differences") +
  theme_classic()

violin_difference
ggsave("hisse_results/hisse_violin_differences")

# dplyr

quantile_diff <- T_diff %>%
  group_by(difference) %>%
  reframe(res = quantile(dens, probs = c(0.025, 0.975)))
quantile_diff

anc_states <- processAncStates(path ="asr_hisse_binary_type.tree",
                               state_labels=c("0"="simple A","1"="compound A","2"="simple B", "3"="compound B"))
plotAncStatesMAP(t = anc_states, tree_layout = "rectangular",
                 state_transparency = 0.5,
                 node_size = c(0.1, 5),
                 tip_labels_size = 2,
                 tip_states_size=2,
                 node_color = c("#D9081D", "#F8D3D8", "#072AC8", "#A2D6F9"))
ggsave("hisse_results/plotAncStatesMAP.png")
