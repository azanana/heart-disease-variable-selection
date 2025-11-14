set.seed(123)

col_names <- c("age","sex","cp","trestbps","chol","fbs","restecg",
               "thalach","exang","oldpeak","slope","ca","thal","num")

data <- read.csv("data/processed.cleveland.data", header=FALSE,
                 col.names=col_names, na.strings="?")

data$heart_disease <- ifelse(data$num > 0, 1, 0)
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

formula <- heart_disease ~ age + sex + cp + trestbps + chol + fbs + 
           restecg + thalach + exang + oldpeak + slope + ca + thal

X <- model.matrix(formula, data)
y <- as.integer(data$heart_disease)

X[, -1] <- scale(X[, -1])

assign_vec <- attr(X, "assign")
term_labels <- attr(terms(formula), "term.labels")
groups <- lapply(seq_along(term_labels), function(i) which(assign_vec == i))
names(groups) <- term_labels

log_likelihood <- function(y, X, beta) {
  z <- as.numeric(X %*% beta)
  log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
  sum(y * z - log1pexp(z))
}

log_prior <- function(beta, mu, sigma2) {
  sum(dnorm(beta[beta != 0], mu, sqrt(sigma2), log=TRUE))
}

rjmcmc <- function(y, X, n_iter, step_size, mu, sigma2, p_within, groups) {
  
  n_pred <- ncol(X) - 1
  gamma <- c(1, rep(0, n_pred))
  beta <- rep(0, n_pred + 1)
  
  ll <- log_likelihood(y, X, beta)
  lp <- log_prior(beta, mu, sigma2)
  
  gamma_trace <- matrix(0, n_iter, n_pred + 1)
  beta_trace <- matrix(0, n_iter, n_pred + 1)
  log_post <- numeric(n_iter)
  
  for (iter in 1:n_iter) {
    
    if (runif(1) < p_within) {
      for (j in which(gamma == 1)) {
        beta_prop <- beta
        beta_prop[j] <- beta[j] + runif(1, -step_size, step_size)
        
        ll_prop <- log_likelihood(y, X, beta_prop)
        lp_prop <- log_prior(beta_prop, mu, sigma2)
        
        log_ratio <- (ll_prop + lp_prop) - (ll + lp)
        
        if (log(runif(1)) < log_ratio) {
          beta <- beta_prop
          ll <- ll_prop
          lp <- lp_prop
        }
      }
    } else {
      group_idx <- sample(length(groups), 1)
      cols <- groups[[group_idx]]
      
      gamma_prop <- gamma
      beta_prop <- beta
      
      for (j in cols) {
        gamma_prop[j] <- 1 - gamma[j]
        if (gamma_prop[j] == 0) {
          beta_prop[j] <- 0
        } else {
          beta_prop[j] <- rnorm(1, mu, sqrt(sigma2))
        }
      }
      
      ll_prop <- log_likelihood(y, X, beta_prop)
      lp_prop <- log_prior(beta_prop, mu, sigma2)
      
      log_ratio <- (ll_prop + lp_prop) - (ll + lp)
      
      if (log(runif(1)) < log_ratio) {
        gamma <- gamma_prop
        beta <- beta_prop
        ll <- ll_prop
        lp <- lp_prop
      }
    }
    
    gamma_trace[iter, ] <- gamma
    beta_trace[iter, ] <- beta
    log_post[iter] <- ll + lp
  }
  
  list(gamma=gamma_trace, beta=beta_trace, log_posterior=log_post)
}

n_iterations <- 20000
burnin <- 10000

output <- rjmcmc(y, X, n_iterations, 
                 step_size=0.1, mu=0, sigma2=1, 
                 p_within=0.8, groups=groups)

gamma_post_burnin <- output$gamma[(burnin+1):n_iterations, ]

inclusion_prob <- numeric(length(groups))
names(inclusion_prob) <- names(groups)

for (i in seq_along(groups)) {
  cols <- groups[[i]]
  inclusion_prob[i] <- mean(rowSums(gamma_post_burnin[, cols, drop=FALSE]) > 0)
}

print("Posterior inclusion probabilities:")
print(sort(inclusion_prob, decreasing=TRUE))

