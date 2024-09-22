# Load catch data (updated to 2024)

library(tidyverse)
library(MARSS)
theme_set(theme_bw())

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

catch <- read.csv("./Data/GOA_salmon_catch.csv") 

# Process catch data
catch %>%
  pivot_longer(cols = c(2:5), names_to = "species") %>% # pivot longer
  rename(catch = value) %>%
  mutate(lag = case_when((species == "sockeye") ~ 2, # add in lag values for ocean entry
                         (species %in% c("pink", "coho")) ~ 1,
                         (species == "chum") ~ 3),
         lagged_year = year - lag) %>% # lag the years
  filter(lagged_year > 1964) %>% # limit to 1965 onward
  group_by(lagged_year, species) %>%
  reframe(total_catch = sum(catch), # summarize
          log_catch = log10(total_catch)) %>%
  dplyr::select(!total_catch) %>%
  pivot_wider(., names_from = c(species), values_from = "log_catch") -> catch.wide 
  #filter(lagged_year < 2022) add in if running PCA so no NAs
  
 
# scale the data
pc.scale <- scale(catch.wide[1:57,-1])

# run pca
pc.catch <- prcomp(pc.scale)

# Get pc1
catch.pc1 <- pc.catch$x[,1]


# Combine pc1 with year
data.frame(year = 1965:2021, catch = catch.pc1) -> pc1.dat


# Read in old pc1 data from Litzow et al. 2018
old.pc1.dat <- read.csv("./Data/GOA_salmon_catch_PC1_SST.csv") %>%
  rename(catch = catch_pc1) %>%
  mutate(type = "Old PC1") 

# Modify new pc1.dat
new.pc1.dat <- pc1.dat %>%
  mutate(type = "New PC1")

# Bind new and old for comparison
rbind(old.pc1.dat %>% select(year, catch, type), new.pc1.dat) -> plot.dat

# plot
ggplot(plot.dat, aes(year, catch, color = type))+
    geom_line(size = 1)+
      geom_point(size = 2) + 
  theme_bw()

## run DFA ------------------------------------

# set up data
dfa.dat <- catch.wide %>%
  dplyr::select(-lagged_year) %>%
  t()

colnames(dfa.dat) <- 1965:2023

# and look at correlations
cors <- cor(t(dfa.dat), use = "p")
diag(cors) <- 0

cors

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# changing convergence criterion to ensure convergence
# cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# fit models & store results
m = 1 # single-trend models

for(R in levels.R) {
    
    dfa.model = list(A="zero", R=R, m=m)
    
    kemz = MARSS(dfa.dat, model=dfa.model,
                 form="dfa", z.score=TRUE) #, control=cntl.list)
    
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    
    assign(paste("kemz", m, R, sep="."), kemz)
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data # unconstrained is best

# save model selection table--note that unconstrained models did not converge!
write.csv(model.data, "./Output/dfa_model_selection_table.csv",
          row.names = F)

# fit best model
model.list = list(A="zero", m=1, R="unconstrained") # second-best model - this is the borealization index

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa")


# process loadings and trend

CI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names=rownames(dfa.dat),
                      mean=CI$par$Z[1:4],
                      upCI=CI$par.upCI$Z[1:4],
                      lowCI=CI$par.lowCI$Z[1:4])

dodge <- position_dodge(width=0.9)


plot.CI$names <- reorder(plot.CI$names, CI$par$Z[1:4])


# plot loadings
ggplot(plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[6]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

# looks as you'd expect

ggsave("./Figures/salmon_DFA_loadings.png", width = 4, height = 4, units = 'in')


# plot trend
trend <- data.frame(t=1965:2023,
                    estimate=as.vector(mod$states),
                    conf.low=as.vector(mod$states)-1.96*as.vector(mod$states.se),
                    conf.high=as.vector(mod$states)+1.96*as.vector(mod$states.se))

ggplot(trend, aes(t, estimate)) +
  theme_bw() +
  geom_line() +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_errorbar(aes(x=t, ymin=conf.low, ymax=conf.high)) +
  xlab("") +
  ylab("DFA trend") +
  scale_x_continuous(breaks = seq(1970, 2020, 10))

ggsave("./Figures/salmon_DFA_trend.png", width = 6, height = 4, units = 'in')

# lots of uncertainty in 2022 and 2023 b/c of NAs

# and save loadings and trend
write.csv(plot.CI, "./Output/dfa_loadings.csv", row.names = F)
write.csv(trend, "./Output/dfa_trend.csv", row.names = F)

# compare with PC1 scores
trend <- trend %>%
  rename(year = t, 
         catch = estimate) %>%
  mutate(type = "DFA") %>%
  select(year, catch, type)

plot.dat <- rbind(plot.dat, trend)

ggplot(plot.dat, aes(year, catch, color = type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = cb[c(2,4,6)])

# interesting - similar but different!

ggsave("./Figures/PC1_DFA_comparison.png", width = 6, height = 4, units = 'in')
