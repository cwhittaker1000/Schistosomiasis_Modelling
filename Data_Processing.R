# Data Processing
Schisto_Data <- Schisto_Data
Minimal_Schisto_Data <- Schisto_Data[, c("age", "epg_base", "TAL1_IgE_base")]
na_subsetter <- !is.na(Minimal_Schisto_Data$TAL1_IgE_base)
Ready_Schisto <- Minimal_Schisto_Data[na_subsetter, ]

# Adds Age as a Categorical Variable to the Data
Ready_Schisto <- within(Ready_Schisto,  {
  agecat <- cut(age, breaks = c(seq(0, 45, by=5), max(age)))
})

# Calculating the mean epg and IgE for different age categories
# age_categories <- c("(5,10]","(10,15]","(15,20]",
#                     "(20,25]","(25,30]","(30,35]", "(35, 40]", "(40, 45]", "(45, 50]")

mean_epg_age_group <- aggregate(epg_base ~ agecat, data=Ready_Schisto, FUN="mean")
mean_IgE_age_group <- aggregate(TAL1_IgE_base ~ agecat, data=Ready_Schisto, FUN="mean")
