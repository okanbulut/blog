
# Preparing the original data for measurement invariance analysis
library("dplyr")
library("magrittr")

# Load the data in to a dataframe for later use, please enter the appropriate path for the CSV.
df <- read.csv('NFWBS_PUF_2016_data.csv') %>%
  filter(sample == 1, agecat < 6) %>%
  select(starts_with("FWB1_"), starts_with("FWB2_"), PPGENDER)

# Recode gender
df$PPGENDER = plyr::revalue(factor(df$PPGENDER), c(
  `1` = "male",
  `2` = "female"
))

# Rename variables
df <- rename(df, gender = PPGENDER) %>%
  as.data.frame()

# Update item names
colnames(df)[1:10] <- paste0("item", 1:10)

# Recode missing values for items into NA
df <- mutate_if(df, is.integer, list(~na_if(., -1)))
df <- mutate_if(df, is.integer, list(~na_if(., -4)))

# Remove missing cases (only 11 participants)
df <- na.omit(df)

# Export the dataset
write.csv(df, "finance.csv", quote = FALSE, row.names = FALSE)

