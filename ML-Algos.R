## Data Pre-processing
library(tidyverse)
library(readr)
library(mlr3verse)
library(mlr3tuning)
library(paradox)
library(xgboost)
library(glmnet)
library(mltools)

new_feature_matrix_csv <- read_csv("C:/Users/Frank/Desktop/R Projects/mPCR/new_feature_matrix.csv.txt")

data <- new_feature_matrix_csv %>% drop_na()

# OH Encoded Data
data.oh <- data %>% mutate(expoutput = factor(expoutput)) %>% as.data.table
data.oh <- one_hot(data.oh) %>% select(where(is.numeric)) %>% mutate(
  log.n_reads = log(n_reads + 1)
) %>% select(-n_reads)

# Numeric data
data.numeric <- data %>% select(is.numeric) %>% drop_na()

# Stacked symmetric data set
data.symmetric <- data.numeric 
colnames(data.symmetric) <- gsub("_L$", "_R", colnames(data.symmetric))
data.symmetric <- rbind(data.numeric, data.numeric) 
data.symmetric <- data.symmetric %>% mutate(
  log.n_reads = log(n_reads +1)
) %>% select(-n_reads)

# Stacked non-symmetric
data.stacked.nonsymm <- rbind(data.numeric, data.numeric) 
data.stacked.nonsymm <- data.stacked.nonsymm %>% mutate(
  log.n_reads = log(n_reads +1)
) %>% select(-n_reads)
symmetric.data <- data.numeric["n_reads"]
fun.list <- list(mean = function(L,R){(L + R)/2}, min = pmin, max = pmax)
for(feature in sub("_L$", "", grep("_L$", names(data.numeric), value = TRUE))){
  for(funname in names(fun.list)){
    fun <- fun.list[[funname]]
    symmetric.data[[paste0(feature, "_", funname)]] <- fun(data.numeric[[paste0(feature, "_L")]], 
                                  data.numeric[[paste0(feature, "_R")]])
    
  }
}

symmetric.data <- symmetric.data %>% mutate(
  log.n_reads = log(n_reads +1)
) %>% select(-n_reads)

# X <- select(data.numeric, -n_reads, -strain, -Ampliconoutput, -Amplicon)
# y <- data.numeric$n_reads

# Setting colnames for task creation
colnames(data.numeric) <- make.names(colnames(data.numeric),unique = T)
colnames(data.oh) <- make.names(colnames(data.oh),unique = T)
colnames(symmetric.data)<- make.names(colnames(symmetric.data),unique = T)
colnames(data.stacked.nonsymm)<- make.names(colnames(data.stacked.nonsymm),unique = T)

log.label.data <- data.numeric %>% mutate(
  log.n_reads = log(n_reads +1)
) %>% select(-n_reads)

# mean.log.reads <- mean(log.label.data$log.n_reads)
# mean((mean.log.reads - log.label.data$log.n_reads)^2)

task_mpcr <- TaskRegr$new(id = "mpcr", backend = data.numeric, target = "n_reads")
task_log.label <- TaskRegr$new(id = "mpcr.asymmetric", backend = log.label.data, target = "log.n_reads")
task_cat <- TaskRegr$new(id = "mpcr.cat", backend = data.oh, target = "log.n_reads")
task_symmetric <- TaskRegr$new(id = "mpcr.symmetric.features", backend = data.symmetric, target = "log.n_reads")
task_stacked.nonsymm <- TaskRegr$new(id = "mpcr.asymmetric.stacked", backend = data.symmetric, target = "log.n_reads")
measure = msr("regr.mse")

xgb_learn <- lrn("regr.xgboost")

set.seed(103)
fivefold.cv = rsmp("cv", folds = 5)

param.list <- list(#alpha = p_dbl(lower = 1e-20, upper = 1e10, logscale = TRUE), 
                   #lambda = p_dbl(lower = 1e-30, upper = 100, logscale = TRUE),
                   eta = p_dbl(lower = 0, upper = 1)
                  # max_depth = p_int(2, 10)
)

model.list <- list()
for(model.i in 1:length(param.list)){
  
  param.list.subset <- param.list[model.i]
  search_space <- do.call(ps, param.list.subset)
  
  at <- AutoTuner$new(
    learner = xgb_learn,
    resampling = rsmp("cv", folds = 5),
    measure = measure,
    search_space = search_space,
    terminator = trm("none"),
    tuner = tnr("grid_search", resolution = 15),
    store_tuning_instance = TRUE
  )
  
  at$id = paste0(at$id, ".", names(param.list[model.i]))
  
  model.list[[model.i]] <- at
}

model.list <- c(model.list, list(xgb_learn, lrn("regr.featureless"), lrn("regr.cv_glmnet")))

set.seed(103)
grid <- benchmark_grid(
  task = list(task_log.label, task_symmetric),
  learner = model.list,
  resampling = rsmp("cv", folds =3)
)

bmr <- benchmark(grid, store_models = TRUE)
bmr.data <- bmr$data$as_data_table()
llist<- lapply(1:nrow(bmr.data), function(i)as.data.table(bmr.data$learner[[i]]$learner$param_set$values))
llist$fill=TRUE
do.call(rbind, llist)

bmr.df <- as.data.frame(bmr$score(measure)) %>% select(iteration, regr.mse, learner_id,
                                                       uhash, nr, task_id)
bmr.df <- bmr.df %>% mutate(
  model = str_replace(learner_id, "regr.", ""),
  fold_id = iteration
)

bmr.df <- bmr.df %>% select(model, regr.mse, fold_id, task_id)

ggplot(bmr.df, aes(x = regr.mse, y = task_id)) +
  geom_point(aes(colour = factor(fold_id))) +
  facet_grid(model ~ .) +
  labs(colour ="Fold ID")
#  geom_point(aes(colour = factor(task_id))) +
#  labs(colour = "Training data with entries")













## K fold CV Stuff
n.folds <- 5
set.seed(1)
fold.vec <- rep(sample(1:n.folds), l = nrow(data.numeric))
# fold.vec <- rep(sample(1:n.folds), l = nrow(primer_sets))
for(test.fold in 1:n.folds) {
  
  is.test <- fold.vec == test.fold
  is.train <- !is.test
  
  X.train <- X.mat[is.train, ]
  y.train <- y[is.train]
  X.test <- X.mat[is.test, ]
  y.test <- y[is.test]
  baseline <- mean(y.train)
  
  xgb_train <- xgb.DMatrix(data = X.train, label = y.train)
  xgb_model <- xgboost(data = xgb_train, eval_metric = "poisson-nloglik",
                       nrounds = 149, objective = "count:poisson", 
                       max_depth = 9, mind_child_weight = 33.334,
                       reg_alpha = 44.445, reg_lambda = 0.001)
  
  xgb.pred.reads <- predict(xgb_model, X.test)
  
  # Glmnet
  lasso_reg <- cv.glmnet(X.train, y.train, family = 'poisson')
  lasso.pred.reads <- exp(predict(lasso_reg, X.test))
  
  # Figure for mean poisson loss per model
  model.pred.list <- list("baseline" = baseline, 
                          "glmnet" = lasso.pred.reads,
                          "xgboost" = xgb.pred.reads)
  test.error.list <- NULL
  
  for(model.name in  names(model.pred.list)){
    test.error.list[[test.fold]] <- data.frame(
      test.fold, 
      model.name, 
      poissloss = poisson.loss(y.test, model.pred.list[[model.name]])
      
    )
    
  }
  
  test.error <- do.call(rbind, test.error.list)
}

poisson.loss <- function(x, m) {
  mean(m - x*log(m))
}
