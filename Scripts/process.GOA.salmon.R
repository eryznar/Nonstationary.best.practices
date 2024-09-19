# Load catch data (updated to 2024)
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
pc.scale <- scale(catch.wide[,-1])

# run pca
pc.catch <- prcomp(pc.scale)

# Get pc1
catch.pc1 <- pc.catch$x[,1]


# Combine pc1 with year
data.frame(year = catch.wide$lagged_year, catch = catch.pc1) -> pc1.dat


# Read in old pc1 data from Litzow et al. 2018
old.pc1.dat <- read.csv("./Data/GOA_salmon_catch_PC1_SST.csv") %>%
  rename(catch = catch_pc1) %>%
  mutate(type = "Old")


# Modify new pc1.dat
new.pc1.dat <- pc1.dat %>%
  mutate(type = "New")

# Bind new and old for comparison
rbind(old.pc1.dat %>% select(!sst_3yr_running_mean), new.pc1.dat) -> plot.dat

# plot
ggplot(plot.dat, aes(year, catch, color = type))+
    geom_line(size = 1)+
      geom_point(size = 2) + 
  theme_bw()


