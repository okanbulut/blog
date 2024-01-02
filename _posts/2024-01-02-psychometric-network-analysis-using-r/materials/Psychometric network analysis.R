
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                      EDPY 607 - MEASUREMENT THEORY II                    ----
##                        PSYCHOMETRIC NETWORK MODELS                         ~~
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Install the new packages
install.packages(c("devtools", "shiny", "qgraph", "bootnet", "EGAnet", 
                   "psychonetrics", "psychTools"))

# Activate all required packages

# for data wrangling
library("dplyr")
library("ggcorrplot")  
library("DataExplorer") 

# for network modeling 
library("qgraph")
library("psychonetrics")
library("bootnet")  
library("EGAnet")  

# Ancillary packages
library("mirt") # for simulating binary data
library("shiny") # for running shiny apps
library("psychTools") # for additional tools & data from the psych package

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Example 1                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read the data in
sapa <- read.csv("sapa_data.csv", header = TRUE)

# View the first 6 rows
head(sapa)

# Save the correlation matrix
cormat <- psych::tetrachoric(x = sapa)$rho

# Correlation matrix plot
ggcorrplot::ggcorrplot(corr = cormat, # correlation matrix
                       type = "lower", # print only the lower part of the correlation matrix
                       hc.order = TRUE, # hierarchical clustering
                       show.diag = TRUE, # show the diagonal values of 1
                       lab = TRUE, # add correlation values as labels
                       lab_size = 3) # Size of the labels


#....................Graphical Gaussian Models...................

# Now, let's start with regular network analysis
network_sapa_1 <- bootnet::estimateNetwork(
  data = sapa,
  # Alternatively, "cov" for covariances, "cor" for correlations 
  corMethod = "cor_auto", # for polychoric and polyserial correlations
  # Alternatively, "ggmModSelect" for an unregularized GGM using glasso
  default = "EBICglasso", # for estimating GGM with gLASSO and EBIC
  tuning = 0.5 # EBIC tuning parameter; set to zero for BIC model selection
)

# Print the estimated network
print(network_sapa_1)

# View the estimated network
plot(network_sapa_1, layout = "spring") 
plot(network_sapa_1, layout = "circle") 

# See the weighted partial correlations as a heatmap
cor.net <- network_sapa_1$graph

# Plot the correlation matrix as an interactive heatmap
plotly::plot_ly(x = colnames(cor.net), 
                y = colnames(cor.net),
                z = cor.net, 
                type = "heatmap")

# Centrality indices
qgraph::centralityTable(network_sapa_1)

qgraph::centralityPlot(network_sapa_1,
                       include = c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
                       scale = "raw", # or use "z-scores"
                       orderBy = "Strength")

# Let's also try the same model using psychonetrics
obsvars <- colnames(sapa)

network_sapa_2 <- psychonetrics::ggm(sapa, 
                                     # vars argument gets handy when it comes to selecting
                                     # some variables but not all in the dataset
                                     vars = obsvars) %>%
  psychonetrics::runmodel() # Run the model

# View the model parameters
network_sapa_2 %>% psychonetrics::parameters()

# Look at the model fit
network_sapa_2 %>% psychonetrics::fit()

# Plot the confidence intervals
psychonetrics::CIplot(network_sapa_2,  
                      matrices = "omega",
                      alpha_ci = 0.05)

# We can prune this model to remove insignificant edges
network_sapa_3 <- psychonetrics::ggm(sapa, 
                                     vars = obsvars) %>%
  psychonetrics::runmodel() %>%
  psychonetrics::prune(adjust = "fdr", alpha = 0.05)

# View the model parameters
network_sapa_3 %>% psychonetrics::parameters()

# Look at the model fit
network_sapa_3 %>% psychonetrics::fit()

# Check modification indices
network_sapa_3 %>% psychonetrics::MIs()

# Update the model with new edges
network_sapa_4 <- psychonetrics::ggm(sapa, 
                                     vars = obsvars) %>%
  psychonetrics::runmodel() %>%
  psychonetrics::prune(adjust = "fdr", alpha = 0.05) %>%
  # To automatically add edges at  alpha=0.05 until BIC is no longer be improved
  psychonetrics::stepup(criterion = "bic", alpha = 0.05) %>%
  psychonetrics::modelsearch()

# Look at the model fit
network_sapa_4 %>% psychonetrics::fit()

# Obtain the network plot
net4 <- getmatrix(network_sapa_4, "omega")

qgraph::qgraph(net4, 
               layout = "spring", 
               theme = "colorblind",
               labels = obsvars)

# Compare the models
comparison <- psychonetrics::compare(
  `1. Original model`  = network_sapa_2,
  `2. Sparse Model: Only Pruning` = network_sapa_3,
  `3. Sparse Model: Pruning + Step-up` = network_sapa_4)

print(comparison)

#...............................................................................
#                                                                              .
#                                  Exercises                                   .
#                                                                              .
#...............................................................................


# 1. Run the network model again by changing the tuning parameter to 0 and 0.25.
# How did the model change?

# 2. Change the adjustment method for p-values. See ?p.adjust to see alternative 
# methods.

#....................Confirmatory Fit of GGMs....................

# Split the data into two parts
set.seed(2023)
index <- sample.int(n = nrow(sapa), size = 800)

trainData <- sapa[index,]
testData <- sapa[-index,]

# Form a saturated model
GGM_train <-  psychonetrics::ggm(trainData)

# Perform model selection steps:
GGM_train <- GGM_train %>% 
  psychonetrics::runmodel() %>% 
  psychonetrics::prune(adjust = "fdr", alpha = 0.05) %>%
  psychonetrics::stepup(criterion = "bic", alpha = 0.05) %>%
  psychonetrics::modelsearch()

# Obtain the network for the best model
net_train <- psychonetrics::getmatrix(GGM_train, "omega")
qgraph::qgraph(net_train, theme = "colorblind", layout = "spring")

# Obtain the structure to apply to the test data
structure <- 1*(net_train != 0)

# Form the test model
GGM_test <- psychonetrics::ggm(testData, # test dataset
                omega = structure) %>% 
  psychonetrics::runmodel()

# Obtain the network for the estimated model
net_test <- psychonetrics::getmatrix(GGM_test, "omega")
qgraph::qgraph(net_test, theme = "colorblind", layout = "spring")

# Evaluate the fit
GGM_test %>% psychonetrics::fit()

# Lastly, we can run some simulations to identify the required sample size
# etc. to estimate a stable network
simRes <- bootnet::netSimulator(network_sapa_1$graph,
                                dataGenerator = ggmGenerator(
                                  ordinal = TRUE, 
                                  nLevels = 2),
                                default = "EBICglasso",
                                nCases = c(250,500,1000),
                                tuning = 0.5,
                                nReps = 50,
                                nCores = 8
)

plot(simRes)

plot(simRes, yvar = c("strength","closeness","betweenness"))

#...................Exploratory Graph Analysis...................

# Dimension stability analysis via EGAnet
bootEGA_sapa1 <- EGAnet::bootEGA(
  # we could also provide the cor matrix but then
  # n (i.e., number of rows) must also be specified
  data = sapa, 
  cor = "cor_auto",
  uni.method = "louvain",
  iter = 500, # Number of replica samples to generate
  # resampling" for n random subsamples of the original data
  # parametric" for n synthetic samples from multivariate normal dist.
  type = "parametric", 
  # EGA Uses standard exploratory graph analysis
  # EGA.fit Uses total entropy fit index (tefi) to determine best fit of EGA
  # hierEGA Uses hierarchical exploratory graph analysis
  EGA.type = "EGA", 
  model = "glasso", 
  algorithm = "walktrap", # or "louvain" (better for unidimensional structures)
  # use "highest_modularity", "most_common", or "lowest_tefi"
  consensus.method = "highest_modularity", 
  typicalStructure = TRUE, # typical network of partial correlations
  plot.typicalStructure = TRUE, # returns a plot of the typical network
  ncores = 8 # Number of cores to use in computing results
)

# View the number of communities
bootEGA_sapa1$EGA
bootEGA_sapa1$typicalGraph$typical.dim.variables

# Dimension (i.e., structural) stability results
dim_sapa <- EGAnet::dimensionStability(bootEGA_sapa1)
dim_sapa$dimension.stability

# Item stability results
dim_sapa$item.stability
dim_sapa$item.stability$plot # to see only the plot


# Now, let's try hierarchical EGA
bootEGA_sapa2 <- EGAnet::bootEGA(
  # we could also provide the cor matrix but then
  # n (i.e., number of rows) must also be specified
  data = sapa, 
  cor = "cor_auto",
  uni.method = "louvain",
  # We will reduce it down to 50 to run it fast but
  # 500 or more iterations are recommended
  iter = 50, 
  type = "parametric", 
  # hierEGA Uses hierarchical exploratory graph analysis
  EGA.type = "hierEGA",
  model = "glasso", 
  algorithm = "walktrap", 
  consensus.method = "highest_modularity", 
  typicalStructure = TRUE,
  plot.typicalStructure = TRUE, 
  ncores = 4 
)

#...............................................................................
#                                                                              .
#                                  Exercises                                   .
#                                                                              .
#...............................................................................


# 1. Run the network model again by changing the type to "resampling". Did the 
# results change? Do you think parametric or non-parametric is more suitable for
# our analysis?

# 2. Use EGA.type = "EGA.fit" to run the model again. What is the best fitting model
# based on Total Entropy Fit Index?

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Example 2                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read the data in
finance <- read.csv("finance.csv", header = TRUE)

# View the first six rows
head(finance)

# Select the items only
response <- finance %>%
  dplyr::select(dplyr::starts_with("fwb"))

# Recode missing items
response <- apply(response, # data to apply the function
                  2, # 1 to apply to each row; 2 to apply to each column
                  function(x) ifelse(x %in% c(-1, -4), NA, x)) %>%
  as.data.frame()

# Check response utilization
apply(response, 2, table)

# Check the correlation matrix
cormat_fwb <- psych::polychoric(x = response)$rho

# Create a correlation matrix plot with the new data
ggcorrplot::ggcorrplot(corr = cormat_fwb, # correlation matrix
                       type = "lower", # print only the lower part of the correlation matrix
                       hc.order = TRUE, # hierarchical clustering
                       show.diag = TRUE, # show the diagonal values of 1
                       lab = TRUE, # add correlation values as labels
                       lab_size = 3) # Size of the labels


# For the following analysis, we can either generate our correlation matrix or
# use the original dataset. qgraph automatically detects the best matrix
# for a given data type
qgraph::cor_auto(response)

# Alternatively, we could use the following, especially if the amount of
# missing data was high, because it uses full information maximum likelihood
psych::corFiml(response)

#...................Exploratory Graph Analysis...................

# EGA with random-intercept model for detecting item wording effects
# This estimates the number of dimensions after controlling for wording effects
# EGA is applied to a residual correlation matrix after subtracting a random intercept 
# factor model with equal unstandardized loadings from all the regular and unrecoded 
# reversed items in the data
riEGA_1a <- EGAnet::bootEGA(
  response, 
  uni.method = "LE",
  iter = 500,
  type = "parametric",
  cor = "cor_auto",
  model = "glasso",
  EGA.type = "riEGA", # select random-intercept EGA here
  consensus.method	= "highest_modularity", 
  consensus.iter = 100,
  algorithm="walktrap"
)


# See the results
print(riEGA_1a)

# The same model using riEGA without bootstrapping
?EGAnet::riEGA

riEGA_1b <- EGAnet::riEGA(response, # the data must be unrecoded.
                          cor = "cor_auto",
                          model = "glasso",
                          # use "highest_modularity", "most_common", or "lowest_tefi"
                          consensus.method	= "highest_modularity", 
                          consensus.iter = 100,
                          algorithm="walktrap",
                          plot.EGA = TRUE)

# See the EGA results
print(riEGA_1b)

# Print the random-intercept results
riEGA_1b$RI


# Change consensus method to "highest_modularity"
riEGA_2 <- EGAnet::riEGA(response, 
                         cor = "cor_auto",
                         model = "glasso",
                         # use "highest_modularity", "most_common", or "lowest_tefi"
                         consensus.method	= "most_common", 
                         consensus.iter = 100,
                         algorithm="walktrap",
                         plot.EGA = TRUE)


# See the EGA results
print(riEGA_2)


# Change the algorithm to "louvain"
riEGA_3 <- EGAnet::riEGA(response, 
                         cor = "cor_auto",
                         model = "glasso",
                         # use "highest_modularity", "most_common", or "lowest_tefi"
                         consensus.method	= "highest_modularity", 
                         consensus.iter = 100,
                         algorithm="louvain",
                         plot.EGA = TRUE)

# See the results
print(riEGA_3)

#...............................................................................
#                                                                              .
#                                  Exercises                                   .
#                                                                              .
#...............................................................................


# Use the correlation matrix (from the riEGA model) to create a network plot. You
# can fit a GGM model, or alternatively, use qgraph::EBICglasso() to create a plot
# with qgraph::qgraph afterwards. 

cormat <- riEGA_1b$RI$correlation

# You can change gamma to a different value
parcor <- EBICglasso(S = cormat, n = nrow(response), gamma = 0.5, threshold = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Example 3                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#...............................................................................
#                                                                              .
#  In this example, we will use the Experiences in Close Relationships (ECR)   .
#  scale (Brennan et al., 1998). The ECR scale consists of 36 items measuring  .
#  two higher-order attachment dimensions for adults (18 items per             .
#  dimension): avoidance and anxiety (see Figure 2) 2 . The items are based    .
#  on a 5-point Likert scale (i.e., 1 = strongly disagree to 5 = strongly      .
#  agree). For each subscale (i.e., dimension), higher scores indicate higher  .
#  levels of avoidance (or anxiety). Individuals who score high on either or   .
#  both of these dimensions are assumed to have an insecure adult attachment   .
#  orientation (Wei et al., 2007).                                             .
#                                                                              .
#...............................................................................

# Read the data in
ecr <- read.csv("ecr_data.csv", header = TRUE)

head(ecr)

# Read the item content
ecr_items <- read.csv("ecr_items.csv", header = TRUE)

head(ecr_items)

# Avoidance subscale
avoidance <- ecr[, seq(1, 35, by = 2)]
avoidance_items <- as.character(ecr_items$Label[1:18])

# Anxiety subscale
anxiety <- ecr[, seq(2, 36, by = 2)]
anxiety_items <- as.character(ecr_items$Label[19:36])

# Unique variable analysis to remove redundant items
avoidance_uva <- EGAnet::UVA(
  data = avoidance,
  key = avoidance_items, # item labels
  model = "glasso",
  corr = "cor_auto",
  method = "wTO",
  type = "threshold",
  reduce = TRUE, # reduce the data based on redundant variables
  auto = TRUE, # automatic reduction? If FALSE, manual selection is done
  reduce.method = "latent", # or "remove"
  plot.redundancy = TRUE
)

# View redundant (locally-dependent) items
avoidance_uva$redundancy$redundant
avoidance_uva$redundancy$plot

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Example 4                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Let's simulate binary response data based on a multidimensional model
# Parameters from Reckase (2009), p. 153

set.seed(2023)

a <- matrix(c(
  .7471, .0250, 
  .4595, .0097, 
  .8613, .0067, 
  1.0141, .0080, 
  .5521, .0204, 
  1.3547, .0064, 
  1.3761, .0861, 
  .8525, .0383, 
  1.0113, .0055, 
  .9212, .0119, 
  .3073, .9704, 
  .1819, .4980, 
  .4115,1.1136, 
  .1536,1.7251, 
  .1530, .6688, 
  .2890,1.2419, 
  .1341,1.4882, 
  .0524, .4754, 
  .2139, .4612, 
  .1761,1.1200), 20, 2, byrow=TRUE)*1.702

d <- matrix(
  c(.1826,-.1924,-.4656,-.4336,-.4428,-.5845,-1.0403,
    .6431,.0122,.0912,.6172,-.1955,-.3668,
    -1.7590,-.2434,.4925,-.3410,.2896,.006,.0329), ncol=1)*1.702


data_sim <- simdata(a, d, 2000, itemtype = '2PL')

# Now, let's fit the Ising model using eLASSO
# Check out this article: https://www.nature.com/articles/srep05918
network_simdata <- bootnet::estimateNetwork(
  data = data_sim,
  default = "IsingFit",
  tuning = 0.25 
)

# Print the estimated network
print(network_simdata)

# View the estimated network
plot(network_simdata, layout = "spring") 

# Manually create the network plot
qgraph::qgraph(network_simdata$results$weiadj,
               groups = list(
                 Factor1 = 1:10,
                 Factor2 = 11:20),
               layout = "spring",
               labels = colnames(data_sim),
               color=c("#fc8d59", "#4575b4"))

# Print the parameters
network_simdata$results

# weiadj: The weighted adjacency matrix
# thresholds: Thresholds of the variables
# asymm.weights: Asymmetrical weighted adjacency matrix

# Centrality indices
qgraph::centralityPlot(network_simdata,
                       include = c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
                       scale = "raw", # or use "z-scores"
                       orderBy = "Strength")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                           Additional Information                         ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# The original source: 
# https://github.com/AdelaIsvoranu/SEA
# https://github.com/SachaEpskamp/simulation_exploration_app

shiny::runGitHub("SEA", "okanbulut")


