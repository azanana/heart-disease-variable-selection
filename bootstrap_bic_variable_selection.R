library(stats)
library(parallel)

set.seed(123)

data_file <- "data/processed.cleveland.data"
col_names <- c("age","sex","cp","trestbps","chol","fbs","restecg",
               "thalach","exang","oldpeak","slope","ca","thal","num")

data <- read.csv(data_file, header=FALSE, col.names=col_names, 
                 na.strings="?", stringsAsFactors=FALSE)

data$heart_disease <- factor(ifelse(data$num > 0, 1, 0), 
                              levels=c(0,1), labels=c("No","Yes"))
data$num <- NULL

data$sex <- factor(data$sex, levels=c(0,1), labels=c("Female","Male"))
data$cp <- factor(data$cp, levels=1:4, 
                  labels=c("Typical Angina","Atypical Angina",
                          "Non-anginal Pain","Asymptomatic"))
data$fbs <- factor(data$fbs, levels=c(0,1), 
                   labels=c("<=120 mg/dl",">120 mg/dl"))
data$restecg <- factor(data$restecg, levels=0:2, 
                       labels=c("Normal","ST-T Abnormality","LV Hypertrophy"))
data$exang <- factor(data$exang, levels=c(0,1), labels=c("No","Yes"))
data$slope <- factor(data$slope, levels=1:3, 
                     labels=c("Upsloping","Flat","Downsloping"))
data$thal[is.na(data$thal)] <- 3
data$thal <- factor(data$thal, levels=c(3,6,7), 
                    labels=c("Normal","Fixed Defect","Reversible Defect"))
data$ca[is.na(data$ca)] <- 0
data$ca <- as.numeric(data$ca)
data <- data[complete.cases(data), ]

predictors <- c("age","sex","cp","trestbps","chol","fbs","restecg",
                "thalach","exang","oldpeak","slope","ca","thal")

bootstrap_step <- function(seed_val, dat, vars) {
  set.seed(seed_val)
  boot_data <- dat[sample(nrow(dat), replace=TRUE), ]
  boot_data <- droplevels(boot_data)
  
  formula_str <- paste("heart_disease ~", paste(vars, collapse=" + "))
  fit <- try(glm(as.formula(formula_str), data=boot_data, 
                 family=binomial()), silent=TRUE)
  
  if (inherits(fit, "try-error")) return(character(0))
  
  step_fit <- try(step(fit, direction="backward", trace=0, 
                       k=log(nrow(boot_data))), silent=TRUE)
  
  if (inherits(step_fit, "try-error")) return(character(0))
  
  attr(terms(step_fit), "term.labels")
}

n_bootstrap <- 300
seeds <- 123 + seq_len(n_bootstrap)

n_cores <- min(detectCores(), 4)
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, 123)
clusterExport(cl, c("data","predictors","bootstrap_step"), 
              envir=environment())
clusterEvalQ(cl, library(stats))

results <- parLapply(cl, seeds, function(s) {
  bootstrap_step(s, data, predictors)
})
stopCluster(cl)

results <- Filter(function(x) length(x) > 0, results)

all_vars <- unlist(results)
inclusion_freq <- table(all_vars) / length(results)
inclusion_freq <- inclusion_freq[match(predictors, names(inclusion_freq))]
names(inclusion_freq) <- predictors
inclusion_freq[is.na(inclusion_freq)] <- 0

models <- sapply(results, function(x) paste(sort(x), collapse=" + "))
model_counts <- sort(table(models), decreasing=TRUE)

best_bic <- Inf
best_model <- NULL
for (m in names(model_counts)) {
  form <- as.formula(paste("heart_disease ~", m))
  fit <- try(glm(form, data=data, family=binomial()), silent=TRUE)
  if (!inherits(fit, "try-error")) {
    bic_val <- AIC(fit, k=log(nrow(data)))
    if (bic_val < best_bic) {
      best_bic <- bic_val
      best_model <- fit
    }
  }
}

print("Variable inclusion frequencies:")
print(sort(inclusion_freq, decreasing=TRUE))
print(paste("\nBest model BIC:", round(best_bic, 3)))
print("\nBest model summary:")
print(summary(best_model))

