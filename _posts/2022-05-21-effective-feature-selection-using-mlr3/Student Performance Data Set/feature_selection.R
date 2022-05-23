
library("dplyr")
library("mlr3verse")
library("mlr3fselect")
library("mlr3learners")

data <- read.table("student-por.csv", header = TRUE, sep = ";")

# Binary value based on G3
data$class <- as.factor(ifelse(data$G3 >= 10, "pass", "fail"))

# Remove variables that won't be used in the prediction
data <- select(data, -G1, -G2, -G3)

# Construct a binary classification task
task <- as_task_classif(data, target = "class", positive = "pass")

# Select a classifier
learner <- lrn("classif.ranger", importance = "impurity")

# 5-fold cross-validation
resampling <- rsmp("cv", folds = 5)

# Classification accuracy as the evaluation metric
measure <- msr("classif.ce")
resampling$instantiate(task)

# Termination criteria
terminator = trm("none")

instance = FSelectInstanceSingleCrit$new(
  task = task,
  learner = learner,
  resampling = resampling,
  measure = measure,
  terminator = terminator,
  store_models = TRUE)

fselector = fs("rfe", recursive = FALSE)
fselector$optimize(instance)

as.data.table(instance$archive, 
              exclude_columns = c("runtime_learners", "timestamp", "batch_nr", "resample_result", "uhash"))
