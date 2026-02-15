IGaGa.MLE <- function(alpha, Beta, Var.tau, x, z, ep = 1e-4) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- length(x)
  
  # 确定τ分布的参数λ
  lambda <- 1 / Var.tau + 2
  
  log_f <- function(x, z, alpha, Beta, lambda) {
    f <- function(data) {
      x <- data[1]
      z <- data[-1]
      mu <- exp(z %*% Beta)
      term1 <- gamma(alpha + lambda) / gamma(alpha) / gamma(lambda)
      term2 <- (alpha / mu)^alpha
      term3 <- (lambda - 1)^lambda * x^(alpha - 1)
      term4 <- (alpha * x / mu + lambda - 1)^(-alpha)
      term5 <- (alpha * x / mu + lambda - 1)^(-lambda)
      term4 * term1 * term5 * term2 * term3
    }
    pdf <- apply(cbind(x, z), 1, f)
    return(sum(log(pdf)))
  }
  
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  # US算法求c+1/(s-1)+log(s)-ψ(s)=0的零点
  US.IGa <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (k <= 1000) {
      s <- s_new
      q <- c + log(s - 1) + 1 / (s - 1) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # GIG的矩
  GIG.m <- function(x, z, alpha, Beta, lambda) {
    compute.b <- function(data) {
      x <- data[1]
      z <- data[-1]
      mu <- exp(t(z) %*% Beta)
      alpha * x / mu + lambda - 1
    }
    a <- alpha + lambda
    b <- apply(cbind(x, z), 1, compute.b)
    
    result1 <- a / b
    result2 <- log(b) - digamma(a)
    
    
    return(cbind(result1, result2))
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  loglikeli_new <- log_f(x, z, alpha, Beta, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    loglikeli <- loglikeli_new
    mu <- exp(z %*% Beta)
    
    gig.vals <- GIG.m(x, z, alpha, Beta, lambda)
    b1 <- gig.vals[, 1]
    d1 <- gig.vals[, 2]
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    Wb1 <- diag(b1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - b1 * x / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    r <- log(mu) + x / mu - 1 / b1
    Beta_new <- solve(t(z) %*% Wb1 %*% z) %*% t(z) %*% Wb1 %*% r
    cat("Beta:", Beta_new, "\n")
    # 更新λ
    c <- 1 - Bb1 - Bd1
    lambda_new <- US.IGa(lambda, c)
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, z, alpha_new, Beta_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

IGauGa.MLE <- function(alpha, Beta, Var.tau, x, z, ep = 1e-4) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- length(x)
  
  # 确定τ分布的参数λ
  lambda <- 1 / Var.tau
  
  # 定义对数似然函数
  log_f <- function(x, z, alpha, Beta, lambda) {
    f <- function(data) {
      x <- data[1]
      z <- data[-1]
      mu <- exp(z %*% Beta)
      term1 <- exp(lambda) * sqrt(2 * lambda / pi)
      term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
      term3 <- (lambda / (lambda + 2 * alpha * x / mu))^(alpha / 2 + 1 / 4)
      term4 <- sqrt(lambda * (lambda + 2 * alpha * x / mu))
      term5 <- besselK(term4, alpha + 1 / 2, expon.scaled = FALSE)
      term1 * term2 * term3 * term5
    }
    pdf <- apply(cbind(x, z), 1, f)
    return(sum(log(pdf)))
  }
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # GIG的矩
  GIG.m <- function(x, z, alpha, Beta, lambda) {
    compute.b <- function(data) {
      x <- data[1]
      z <- data[-1]
      mu <- exp(t(z) %*% Beta)
      2 * alpha * x / mu + lambda
    }
    a <- lambda
    b <- apply(cbind(x, z), 1, compute.b)
    p <- -alpha - 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  loglikeli_new <- log_f(x, z, alpha, Beta, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    loglikeli <- loglikeli_new
    mu <- exp(z %*% Beta)
    
    gig.vals <- GIG.m(x, z, alpha, Beta, lambda)
    a1 <- gig.vals[, 1]
    b1 <- gig.vals[, 2]
    d1 <- gig.vals[, 3]
    Ba1 <- mean(a1)
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    Wb1 <- diag(b1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - b1 * x / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    r <- log(mu) + x / mu - 1 / b1
    Beta_new <- solve(t(z) %*% Wb1 %*% z) %*% t(z) %*% Wb1 %*% r
    cat("Beta:", Beta_new, "\n")
    # 更新λ
    lambda_new <- 1 / (Ba1 + Bb1 - 2)
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, z, alpha_new, Beta_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

RIGGa.MLE <- function(alpha, Beta, Var.tau, x, z, ep = 1e-4) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- length(x)
  
  # 确定τ分布的参数λ
  lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
  
  log_f <- function(x, z, alpha, Beta, lambda) {
    f <- function(data) {
      x <- data[1]
      z <- data[-1]
      mu <- exp(z %*% Beta)
      term1 <- exp(1 / (lambda - 1)) * sqrt(2 * lambda / (pi * (lambda - 1)))
      term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
      term3 <- (((lambda * (lambda - 1))^(-1) + 2 * alpha * x / mu) / (lambda / (lambda - 1)))^(-alpha / 2 + 1 / 4)
      term4 <- sqrt(lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * alpha * x / mu))
      term5 <- besselK(term4, -alpha + 1 / 2, expon.scaled = FALSE)
      term1 * term2 * term3 * term5
    }
    pdf <- apply(cbind(x, z), 1, f)
    return(sum(log(pdf)))
  }
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  # GIG的矩
  GIG.m <- function(x, z, alpha, Beta, lambda) {
    compute.b <- function(data) {
      x <- data[1]
      z <- data[-1]
      mu <- exp(t(z) %*% Beta)
      2 * alpha * x / mu + 1 / (lambda * (lambda - 1))
    }
    a <- lambda / (lambda - 1)
    b <- apply(cbind(x, z), 1, compute.b)
    p <- -alpha + 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  alpha_new <- alpha
  Beta_new <- Beta
  lambda_new <- lambda
  loglikeli_new <- log_f(x, z, alpha, Beta, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    Beta <- Beta_new
    lambda <- lambda_new
    loglikeli <- loglikeli_new
    mu <- exp(z %*% Beta)
    
    gig.vals <- GIG.m(x, z, alpha, Beta, lambda)
    a1 <- gig.vals[, 1]
    b1 <- gig.vals[, 2]
    d1 <- gig.vals[, 3]
    Ba1 <- mean(a1)
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    Wb1 <- diag(b1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - b1 * x / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    r <- log(mu) + x / mu - 1 / b1
    Beta_new <- solve(t(z) %*% Wb1 %*% z) %*% t(z) %*% Wb1 %*% r
    cat("Beta:", Beta_new, "\n")
    # 更新λ
    lambda_new <- (1 + 2 * Bb1 + sqrt((1 + 2 * Bb1)^2 - 4 * Bb1 * (3 - Ba1))) / (6 - 2 * Ba1)
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, z, alpha_new, Beta_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Beta = Beta_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

#===============================================================================
set.seed(20260104)
x <- insurance$charges / 1000
n <- length(x)
ind <- sample(1:n, 1000)
x.train <- x[ind]
x.test <- x[-ind]

z <- cbind(1, insurance$age, insurance$bmi, insurance$children, as.numeric(insurance$sex == "male"), as.numeric(insurance$smoker == "yes"))
z.train <- z[ind, ]
z.test <- z[-ind, ]

alpha <- 1
Beta <- rep(0, 6)
Var.tau <- 0.1
Res1 <- IGaGa.MLE(alpha, Beta, Var.tau, x.train, z.train)


alpha <- 3
Beta <- rep(0, 6)
Var.tau <- 10
Res2 <- IGauGa.MLE(alpha, Beta, Var.tau, x.train, z.train)

alpha <- 1
Beta <- rep(0, 6)
Var.tau <- 1
Res3 <- RIGGa.MLE(alpha, Beta, Var.tau, x.train, z.train)

# Gamma 均值回归
data <- data.frame(x = x.train, z = z.train[, -1])
fit <- glm(x ~ ., data = data, family = Gamma(link = "log"))
(beta_ga <- coef(fit))
(alpha_ga <- 1 / summary(fit)$dispersion)
(loglik_ga <- as.numeric(logLik(fit)))

# 逆高斯均值回归
fit_ig <- glm(x ~ ., data = data, family = inverse.gaussian(link = "log"))
(beta_ig <- coef(fit_ig))
(lambda_ig <- 1 / summary(fit_ig)$dispersion)
(loglik_ig <- as.numeric(logLik(fit_ig)))

Loglikelihood <- c(Res1$Loglikelihood, Res2$Loglikelihood, Res3$Loglikelihood, loglik_ga, loglik_ig)
AIC.value <- -2 * Loglikelihood + 2 * c(8, 8, 8, 7, 7)
BIC.value <- -2 * Loglikelihood + log(length(ind)) * c(8, 8, 8, 7, 7)
Final <- data.frame(Model = c("IGaGa", "IGauGa", "RIGGa", "Gamma", "IGaussian"),
                    Loglikelihood = Loglikelihood %>% round(digits = 4), 
                    AIC = AIC.value %>% round(digits = 4), 
                    BIC = BIC.value %>% round(digits = 4))
Final %>% kable()

# 预测===================================================================================
# IGaGa 预测
IGaGa.loglikeli <- function(x, z, alpha, Beta, lambda) {
  f <- function(data) {
    x <- data[1]
    z <- data[-1]
    mu <- exp(z %*% Beta)
    term1 <- gamma(alpha + lambda) / gamma(alpha) / gamma(lambda)
    term2 <- (alpha / mu)^alpha
    term3 <- (lambda - 1)^lambda * x^(alpha - 1)
    term4 <- (alpha * x / mu + lambda - 1)^(-alpha)
    term5 <- (alpha * x / mu + lambda - 1)^(-lambda)
    term4 * term1 * term5 * term2 * term3
  }
  pdf <- apply(cbind(x, z), 1, f)
  return(sum(log(pdf)))
}
IGaGa.test <- IGaGa.loglikeli(x.test, z.test, Res1$Alpha, Res1$Beta, Res1$Lambda)
IGaGa.MAE <- mean(abs(exp(z.test %*% Res1$Beta) - x.test))

# IGau-Ga 预测
IGauGa.loglikeli <- function(x, z, alpha, Beta, lambda) {
  f <- function(data) {
    x <- data[1]
    z <- data[-1]
    mu <- exp(z %*% Beta)
    term1 <- exp(lambda) * sqrt(2 * lambda / pi)
    term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
    term3 <- (lambda / (lambda + 2 * alpha * x / mu))^(alpha / 2 + 1 / 4)
    term4 <- sqrt(lambda * (lambda + 2 * alpha * x / mu))
    term5 <- besselK(term4, alpha + 1 / 2, expon.scaled = FALSE)
    term1 * term2 * term3 * term5
  }
  pdf <- apply(cbind(x, z), 1, f)
  return(sum(log(pdf)))
}
IGauGa.test <- IGauGa.loglikeli(x.test, z.test, Res2$Alpha, Res2$Beta, Res2$Lambda)
IGauGa.MAE <- mean(abs(exp(z.test %*% Res2$Beta) - x.test))
# RIGau-Ga 预测
RIGGa.loglikeli <- function(x, z, alpha, Beta, lambda) {
  f <- function(data) {
    x <- data[1]
    z <- data[-1]
    mu <- exp(z %*% Beta)
    term1 <- exp(1 / (lambda - 1)) * sqrt(2 * lambda / (pi * (lambda - 1)))
    term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
    term3 <- (((lambda * (lambda - 1))^(-1) + 2 * alpha * x / mu) / (lambda / (lambda - 1)))^(-alpha / 2 + 1 / 4)
    term4 <- sqrt(lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * alpha * x / mu))
    term5 <- besselK(term4, -alpha + 1 / 2, expon.scaled = FALSE)
    term1 * term2 * term3 * term5
  }
  pdf <- apply(cbind(x, z), 1, f)
  return(sum(log(pdf)))
}
RIGGa.test <- RIGGa.loglikeli(x.test, z.test, Res3$Alpha, Res3$Beta, Res3$Lambda)
RIGGa.MAE <- mean(abs(exp(z.test %*% Res3$Beta) - x.test))
# Gamma 预测
Gamma.loglikeli <- function(x, z, alpha, Beta) {
  f <- function(data) {
    x <- data[1]
    z <- data[-1]
    mu <- exp(z %*% Beta)
    alpha^alpha / mu^alpha / gamma(alpha) * x^(alpha - 1) * exp(-alpha * x / mu)
  }
  pdf <- apply(cbind(x, z), 1, f)
  return(sum(log(pdf)))
}
Gamma.test <- Gamma.loglikeli(x.test, z.test, alpha_ga, beta_ga)
Gamma.MAE <- mean(abs(exp(z.test %*% beta_ga) - x.test))
# IGaussian 预测
IGaussian.loglikeli <- function(x, z, lambda, Beta) {
  f <- function(data) {
    x <- data[1]
    z <- data[-1]
    mu <- exp(z %*% Beta)
    sqrt(lambda / 2 / pi) * x^(-1.5) * exp(-lambda * (x - mu)^2 / 2 / mu^2 / x)
  }
  pdf <- apply(cbind(x, z), 1, f)
  return(sum(log(pdf)))
}
IGaussian.test <- IGaussian.loglikeli(x.test, z.test, lambda_ig, beta_ig)
IGaussian.MAE <- mean(abs(exp(z.test %*% beta_ig) - x.test))

Loglikelihood <- c(IGaGa.test, IGauGa.test, RIGGa.test, Gamma.test, IGaussian.test)
AIC.value <- -2 * Loglikelihood + 2 * c(8, 8, 8, 7, 7)
BIC.value <- -2 * Loglikelihood + log(n - length(ind)) * c(8, 8, 8, 7, 7)
MAE.value <- c(IGaGa.MAE, IGauGa.MAE, RIGGa.MAE, Gamma.MAE, IGaussian.MAE)
Final <- data.frame(Model = c("IGaGa", "IGauGa", "RIGGa", "Gamma", "IGaussian"),
                    Loglikelihood = Loglikelihood %>% round(digits = 4), 
                    AIC = AIC.value %>% round(digits = 4), 
                    BIC = BIC.value %>% round(digits = 4), 
                    MAE = MAE.value %>% round(digits = 4))
Final %>% kable()
