library(knitr)
library(tidyverse)

IGaGa.MLE <- function(alpha, mu, Var.tau, x, ep = 1e-8) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- length(x)
  
  # 确定τ分布的参数λ
  lambda <- 1 / Var.tau + 2
  
  log_f <- function(x, alpha, mu, lambda) {
    f <- function(x) {
      term1 <- gamma(alpha + lambda) * alpha^alpha * (lambda - 1)^lambda * x^(alpha - 1)
      term2 <- mu^alpha * gamma(alpha) * gamma(lambda)
      term3 <- (alpha * x / mu + lambda - 1)^(-alpha - lambda)
      term1 / term2 * term3
    }
    pdf <- sapply(x, function(x) f(x))
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
    
    while (index == 0) {
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
  GIG.m <- function(x, alpha, mu, lambda) {
    compute.b <- function(x) {
      alpha * x / mu + lambda - 1
    }
    a <- alpha + lambda
    b <- sapply(x, function(x) compute.b(x))
    
    result1 <- a / b
    result2 <- log(b) - digamma(a)
    
    
    return(cbind(result1, result2))
  }
  
  alpha_new <- alpha
  mu_new <- mu
  lambda_new <- lambda
  loglikeli_new <- log_f(x, alpha, mu, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    mu <- mu_new
    lambda <- lambda_new
    loglikeli <- loglikeli_new
    
    gig.vals <- GIG.m(x, alpha, mu, lambda)
    b1 <- gig.vals[, 1]
    d1 <- gig.vals[, 2]
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - b1 * x / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    mu_new <- mean(b1 * x)
    cat("Mu:", mu_new, "\n")
    # 更新λ
    c <- 1 - Bb1 - Bd1
    lambda_new <- US.IGa(lambda, c)
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, alpha_new, mu_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Mu = mu_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

IGauGa.MLE <- function(alpha, mu, Var.tau, x, ep = 1e-8) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- length(x)
  
  # 确定τ分布的参数λ
  lambda <- 1 / Var.tau
  
  # 定义对数似然函数
  log_f <- function(x, alpha, mu, lambda) {
    f <- function(x) {
      term1 <- exp(lambda) * sqrt(2 * lambda / pi)
      term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
      term3 <- (lambda / (lambda + 2 * alpha * x / mu))^(alpha / 2 + 1 / 4)
      term4 <- sqrt(lambda * (lambda + 2 * alpha * x / mu))
      term5 <- besselK(term4, alpha + 1 / 2, expon.scaled = FALSE)
      term1 * term2 * term3 * term5
    }
    pdf <- sapply(x, function(x) f(x))
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
  GIG.m <- function(x, alpha, mu, lambda) {
    compute.b <- function(x) {
      2 * alpha * x / mu + lambda
    }
    a <- lambda
    b <- sapply(x, function(x) compute.b(x))
    p <- -alpha - 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  alpha_new <- alpha
  mu_new <- mu
  lambda_new <- lambda
  loglikeli_new <- log_f(x, alpha, mu, lambda)
  
  while (k <= 1000) {
    alpha <- alpha_new
    mu <- mu_new
    lambda <- lambda_new
    loglikeli <- loglikeli_new
    
    gig.vals <- GIG.m(x, alpha, mu, lambda)
    a1 <- gig.vals[, 1]
    b1 <- gig.vals[, 2]
    d1 <- gig.vals[, 3]
    Ba1 <- mean(a1)
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - b1 * x / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    mu_new <- mean(b1 * x)
    cat("Mu:", mu_new, "\n")
    # 更新λ
    lambda_new <- 1 / (Ba1 + Bb1 - 2)
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, alpha_new, mu_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Mu = mu_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

RIGGa.MLE <- function(alpha, mu, Var.tau, x, ep = 1e-8) {
  # 收敛指标与迭代次数
  index <- 0
  k <- 1
  
  # 样本量n
  n <- length(x)
  
  # 确定τ分布的参数λ
  lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
  
  log_f <- function(x, alpha, mu, lambda) {
    f <- function(x) {
      term1 <- exp(1 / (lambda - 1)) * sqrt(2 * lambda / (pi * (lambda - 1)))
      term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
      term3 <- (((lambda * (lambda - 1))^(-1) + 2 * alpha * x / mu) / (lambda / (lambda - 1)))^(-alpha / 2 + 1 / 4)
      term4 <- sqrt(lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * alpha * x / mu))
      term5 <- besselK(term4, -alpha + 1 / 2, expon.scaled = FALSE)
      term1 * term2 * term3 * term5
    }
    pdf <- sapply(x, function(x) f(x))
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
  GIG.m <- function(x, alpha, mu, lambda) {
    compute.b <- function(x) {
      2 * alpha * x / mu + 1 / (lambda * (lambda - 1))
    }
    a <- lambda / (lambda - 1)
    b <- sapply(x, function(x) compute.b(x))
    p <- -alpha + 1 / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(cbind(result1, result2, result3))
  }
  
  alpha_new <- alpha
  mu_new <- mu
  lambda_new <- lambda
  loglikeli_new <- log_f(x, alpha, mu, lambda)
  
  while (k <= 10000) {
    alpha <- alpha_new
    mu <- mu_new
    lambda <- lambda_new
    loglikeli <- loglikeli_new
    
    gig.vals <- GIG.m(x, alpha, mu, lambda)
    a1 <- gig.vals[, 1]
    b1 <- gig.vals[, 2]
    d1 <- gig.vals[, 3]
    Ba1 <- mean(a1)
    Bb1 <- mean(b1)
    Bd1 <- mean(d1)
    # 更新αj
    c <- 1 + mean(log(x) - log(mu) - b1 * x / mu) - Bd1
    alpha_new <- US(alpha, c)
    cat("Alpha:", alpha_new, "\n")
    # 更新βj
    mu_new <- mean(b1 * x)
    cat("Mu:", mu_new, "\n")
    # 更新λ
    lambda_new <- (1 + 2 * Bb1 + sqrt((1 + 2 * Bb1)^2 - 4 * Bb1 * (3 - Ba1))) / (6 - 2 * Ba1)
    cat("Lambda:", lambda_new, "\n")
    # 更新似然函数
    loglikeli_new <- log_f(x, alpha_new, mu_new, lambda_new)
    cat("Loglikelihood:", loglikeli_new, "\n")
    
    if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
      index <- 1
      break
    }
    cat("Iteration number:", k, "\n")
    k <- k + 1
  }
  
  result <- list(Alpha = alpha_new, 
                 Mu = mu_new, 
                 Lambda = lambda_new, 
                 Loglikelihood = loglikeli_new,
                 number = k, 
                 index = index)
  return(result)
}

mle_gamma_optim <- function(x, alpha_init, beta_init) {
  
  neg_loglik <- function(par) {
    alpha <- par[1]
    beta  <- par[2]
    
    if (alpha <= 0 || beta <= 0) return(Inf)
    
    n <- length(x)
    ll <- n * (alpha * log(beta) - lgamma(alpha)) +
      (alpha - 1) * sum(log(x)) -
      beta * sum(x)
    
    return(-ll)  # optim minimizes
  }
  
  fit <- optim(
    par = c(alpha_init, beta_init),
    fn = neg_loglik,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6)
  )
  
  list(
    alpha_hat = fit$par[1],
    beta_hat  = fit$par[2],
    logLik    = -fit$value,
    convergence = fit$convergence
  )
}

mle_ig_optim <- function(x, mu_init, lambda_init) {
  
  neg_loglik <- function(par) {
    mu <- par[1]
    lambda <- par[2]
    
    if (mu <= 0 || lambda <= 0) return(Inf)
    
    n <- length(x)
    ll <- n / 2 * log(lambda) - n / 2 * log(2 * pi) - 3 / 2 * sum(log(x)) - lambda / 2 / mu^2 * sum((x - mu)^2 / x)
    
    return(-ll)
  }
  
  fit <- optim(
    par = c(mu_init, lambda_init),
    fn = neg_loglik,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6)
  )
  
  list(
    mu_hat = fit$par[1],
    lambda_hat = fit$par[2],
    logLik = -fit$value,
    convergence = fit$convergence
  )
}

mle_weibull_optim <- function(x, shape_init, scale_init) {
  
  neg_loglik <- function(par) {
    k <- par[1]       # shape
    lambda <- par[2]  # scale
    
    if (k <= 0 || lambda <= 0) return(Inf)
    
    n <- length(x)
    ll <- n * log(k) -
      n * k * log(lambda) +
      (k - 1) * sum(log(x)) -
      sum((x / lambda)^k)
    
    return(-ll)  # optim minimizes
  }
  
  fit <- optim(
    par = c(shape_init, scale_init),
    fn = neg_loglik,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6)
  )
  
  list(
    shape_hat = fit$par[1],
    scale_hat = fit$par[2],
    logLik = -fit$value,
    convergence = fit$convergence
  )
}

x <- insurance$charges / 1000
n <- length(x)

Res1 <- IGaGa.MLE(1, 1, 1, x)
Res2 <- IGauGa.MLE(1, 1, 1, x)
Res3 <- RIGGa.MLE(1, 1, 1, x)
Res4 <- mle_gamma_optim(x, 1, 1)
Res5 <- mle_ig_optim(x, 1, 1)
Res6 <- mle_weibull_optim(x, 1, 1)

(Mean <- c(Res1$Mu, Res2$Mu, Res3$Mu, Res4$alpha_hat / Res4$beta_hat, Res5$mu_hat, Res6$scale_hat *  gamma(1 + 1 / Res6$shape_hat), mean(x)))

var.IGa <- Res1$Mu^2 / Res1$Alpha + Res1$Mu^2 * (1 + 1 / Res1$Alpha) / (Res1$Lambda - 2)
Var.IGau <- Res2$Mu^2 / Res2$Alpha + Res2$Mu^2 * (1 + 1 / Res2$Alpha) / Res2$Lambda
var.RIG <- Res3$Mu^2 / Res3$Alpha + Res3$Mu^2 * (1 + 1 / Res3$Alpha) * (2 * Res3$Lambda - 1) * (Res3$Lambda - 1) / Res3$Lambda^2
var.Ga <- Res4$alpha_hat^2 / Res4$beta_hat
var.IGaussian <- Res5$mu_hat^3 / Res5$lambda_hat
var.Weibull <- Res6$scale_hat * (gamma(1 + 2 / Res6$shape_hat) - gamma(1 + 1 / Res6$shape_hat)^2)
(Std <- c(sqrt(var.IGa), sqrt(Var.IGau), sqrt(var.RIG), sqrt(var.Ga), sqrt(var.IGaussian), sqrt(var.Weibull), sqrt(var(x))))

(Loglikelihood <- c(Res1$Loglikelihood, Res2$Loglikelihood, Res3$Loglikelihood, Res4$logLik, Res5$logLik, Res6$logLik))
(AIC.value <- c(-2 * Res1$Loglikelihood + 2 * 3, -2 * Res2$Loglikelihood + 2 * 3, -2 * Res3$Loglikelihood + 2 * 3, 
                -2 * Res4$logLik + 2 * 2, -2 * Res5$logLik + 2 * 2, -2 * Res6$logLik + 2 * 2))
(BIC.value <- c(-2 * Res1$Loglikelihood + log(1338) * 3, -2 * Res2$Loglikelihood + log(1338) * 3, -2 * Res3$Loglikelihood + log(1338) * 3, 
                -2 * Res4$logLik + log(1338) * 2, -2 * Res5$logLik + log(1338) * 2, -2 * Res6$logLik + log(1338) * 2))
Final <- data.frame(Mean = Mean[-7], 
                    Std = Std[-7], 
                    Loglikelihood = Loglikelihood, 
                    AIC = AIC.value, 
                    BIC = BIC.value)
kable(Final %>% round(digits = 2))

#=========================================================================================
# 求IGauGa的估计与置信区间

# 渐近正态CI
score_contributions_theta <- function(theta, x) {
  n <- length(x)
  
  # 使用数值梯度计算每个观测的得分贡献
  score_mat <- matrix(0, nrow = n, ncol = 3)
  
  for(i in 1:n) {
    # 为每个观测计算梯度
    grad_i <- numDeriv::grad(
      func = function(th) {
        alpha <- th[1]
        mu <- th[2]
        lambda <- th[3]
        # 计算单个观测的对数似然
        term1 <- exp(lambda) * sqrt(2 * lambda / pi)
        term2 <- alpha^alpha * x[i]^(alpha - 1) / mu^alpha / gamma(alpha)
        term3 <- (lambda / (lambda + 2 * alpha * x[i] / mu))^(alpha / 2 + 1 / 4)
        term4 <- sqrt(lambda * (lambda + 2 * alpha * x[i] / mu))
        term5 <- besselK(term4, alpha + 1 / 2, expon.scaled = FALSE)
        log_val <- log(term1 * term2 * term3 * term5)
        
        return(log_val)
      },
      theta,
      method = "Richardson",
      method.args = list(eps = 1e-10, d = 0.0001, r = 6)
    )
    score_mat[i, ] <- grad_i
  }
  
  return(score_mat)
}

log_f <- function(x, alpha, mu, lambda) {
  f <- function(x) {
    term1 <- exp(lambda) * sqrt(2 * lambda / pi)
    term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
    term3 <- (lambda / (lambda + 2 * alpha * x / mu))^(alpha / 2 + 1 / 4)
    term4 <- sqrt(lambda * (lambda + 2 * alpha * x / mu))
    term5 <- besselK(term4, alpha + 1 / 2, expon.scaled = FALSE)
    term1 * term2 * term3 * term5
  }
  pdf <- sapply(x, function(x) f(x))
  return(sum(log(pdf)))
}

loglikeli_theta <- function(theta, x) {
  alpha <- theta[1]
  mu <- theta[2]
  lambda <- theta[3]
  return(log_f(x, alpha, mu, lambda))
}

# 计算θ的海森矩阵
hessian_theta <- function(theta, x) {
  K <- length(theta)
  
  H <- numDeriv::hessian(
    func = function(th) loglikeli_theta(th, x = x),
    theta,
    method = "Richardson",
    method.args = list(eps = 1e-10, d = 0.0001, r = 6)
  )
  
  return(H)
}

IGauGa.estimation <- c(Res2$Alpha, Res2$Mu, 1 / Res2$Lambda)


theta_hat <- c(Res2$Alpha, Res2$Mu, Res2$Lambda)

S <- score_contributions_theta(theta_hat, x)
B_matrix <- crossprod(S)
H <- -hessian_theta(theta_hat, x)
Sigma_theta <- solve(H) %*% B_matrix %*% solve(H)
se <- diag(Sigma_theta) %>% sqrt()

lower.bound <- theta_hat - qnorm(0.975) * se
upper.bound <- theta_hat + qnorm(0.975) * se

lower.bound[3] <- 1 / Res2$Lambda - qnorm(0.975) * se[3] / Res2$Lambda^2
upper.bound[3] <- 1 / Res2$Lambda + qnorm(0.975) * se[3] / Res2$Lambda^2

IGauGa.estimation %>% round(digits = 4)
lower.bound %>% round(digits = 4)
upper.bound %>% round(digits = 4)


# Bootstrap CI
mean.std.CI <- function(lasample) {
  G <- dim(lasample)[1]
  lasort <- apply(lasample, 2, sort)
  indexx <- floor(c(0.025 * G, 0.975 * G))
  laL <- (lasort[indexx[1], ] + lasort[indexx[1] + 1, ]) / 2
  laU <- (lasort[indexx[2], ] + lasort[indexx[2] + 1, ]) / 2
  results <- c(laL, laU)
  return(results)
}

Boot.IGauGa <- function(n, alpha, mu, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  rIGauGa <- function(n, alpha, mu, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- 1 / Var.tau
    tau <- rinvGauss(n, 1, lambda)
    
    # 计算Gamma分布的速率参数
    rate.para <- alpha / (mu * tau)
    
    # 生成样本x
    x <- rgamma(n, alpha, rate.para)
    return(x)
  }
  
  g <- 1
  G <- 53
  
  boot.estimation <- list(Alpha.boot = numeric(length = G), 
                          Mu.boot = numeric(length = G), 
                          Var.tau.boot = numeric(length = G))
  
  while (g <= G) {
    boot.sample <- rIGauGa(n, alpha, mu, Var.tau)
    Res <- tryCatch({
      IGauGa.MLE(alpha, mu, Var.tau, boot.sample)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      boot.estimation[[1]][g] <- Res$Alpha
      boot.estimation[[2]][g] <- Res$Mu
      boot.estimation[[3]][g] <- 1 / Res$Lambda
      g <- g + 1
    }
  }
  
  return(boot.estimation)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Boot.IGauGa(n, Res2$Alpha, Res2$Mu, 1 / Res2$Lambda)
stopCluster(cl)

Result.alpha <- result_part[[1]][[1]]
Result.mu <- result_part[[1]][[2]]
Result.var.tau <- result_part[[1]][[3]]

for (i in 2:num_cores) {
  Result.alpha <- c(Result.alpha, result_part[[i]][[1]])
  Result.mu <- c(Result.mu, result_part[[i]][[2]])
  Result.var.tau <- c(Result.var.tau, result_part[[i]][[3]])
}

mean.std.CI(Result.alpha %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.mu %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.var.tau %>% as.matrix()) %>% round(digits = 4)
#===============================================================================
# 求RIGauGa的估计与置信区间
# 渐近正态CI
log_f <- function(x, alpha, mu, lambda) {
  f <- function(x) {
    term1 <- exp(1 / (lambda - 1)) * sqrt(2 * lambda / (pi * (lambda - 1)))
    term2 <- alpha^alpha * x^(alpha - 1) / mu^alpha / gamma(alpha)
    term3 <- (((lambda * (lambda - 1))^(-1) + 2 * alpha * x / mu) / (lambda / (lambda - 1)))^(-alpha / 2 + 1 / 4)
    term4 <- sqrt(lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * alpha * x / mu))
    term5 <- besselK(term4, -alpha + 1 / 2, expon.scaled = FALSE)
    term1 * term2 * term3 * term5
  }
  pdf <- sapply(x, function(x) f(x))
  return(sum(log(pdf)))
}

score_contributions_theta <- function(theta, x) {
  n <- length(x)
  
  # 使用数值梯度计算每个观测的得分贡献
  score_mat <- matrix(0, nrow = n, ncol = 3)
  
  for(i in 1:n) {
    # 为每个观测计算梯度
    grad_i <- numDeriv::grad(
      func = function(th) {
        alpha <- th[1]
        mu <- th[2]
        lambda <- th[3]
        # 计算单个观测的对数似然
        term1 <- exp(1 / (lambda - 1)) * sqrt(lambda / (2 * pi * (lambda - 1)))
        term2 <- alpha^alpha * x[i]^(alpha - 1) / mu^alpha / gamma(alpha)
        term3 <- (((lambda * (lambda - 1))^(-1) + 2 * alpha * x[i] / mu) / (lambda / (lambda - 1)))^(-alpha / 2 + 1 / 4)
        term4 <- sqrt(lambda / (lambda - 1) * (1 / (lambda * (lambda - 1)) + 2 * alpha * x[i] / mu))
        term5 <- besselK(term4, -alpha + 1 / 2, expon.scaled = FALSE)
        
        log_val <- log(term1 * term2 * term3 * term5)
        
        return(log_val)
      },
      theta,
      method = "Richardson",
      method.args = list(eps = 1e-10, d = 0.0001, r = 6)
    )
    score_mat[i, ] <- grad_i
  }
  
  return(score_mat)
}

loglikeli_theta <- function(theta, x) {
  alpha <- theta[1]
  mu <- theta[2]
  lambda <- theta[3]
  return(log_f(x, alpha, mu, lambda))
}

# 计算θ的海森矩阵
hessian_theta <- function(theta, x) {
  K <- length(theta)
  
  H <- numDeriv::hessian(
    func = function(th) loglikeli_theta(th, x = x),
    theta,
    method = "Richardson",
    method.args = list(eps = 1e-10, d = 0.0001, r = 6)
  )
  
  return(H)
}

RIGGa.estimation <- c(Res3$Alpha, Res3$Mu, (2 * Res3$Lambda - 1) * (Res3$Lambda - 1) / Res3$Lambda^2)

theta_hat <- c(Res3$Alpha, Res3$Mu, Res3$Lambda)

S <- score_contributions_theta(theta_hat, x)
B_matrix <- crossprod(S)
H <- -hessian_theta(theta_hat, x)
Sigma_theta <- solve(H) %*% B_matrix %*% solve(H)
se <- diag(Sigma_theta) %>% sqrt()

lower.bound <- theta_hat - qnorm(0.975) * se
upper.bound <- theta_hat + qnorm(0.975) * se

V.tau <- (2 * Res3$Lambda - 1) * (Res3$Lambda - 1) / Res3$Lambda^2

lower.bound[3] <- V.tau - qnorm(0.975) * se[3] * abs((3 * Res3$Lambda - 2) / Res3$Lambda^3)
upper.bound[3] <- V.tau + qnorm(0.975) * se[3] * abs((3 * Res3$Lambda - 2) / Res3$Lambda^3)

RIGGa.estimation %>% round(digits = 4)
lower.bound %>% round(digits = 4)
upper.bound %>% round(digits = 4)

# Bootstrap CI
Boot.RIGauGa <- function(n, alpha, mu, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  rRIGGa <- function(n, alpha, mu, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- (3 + sqrt(1 + 4 * Var.tau)) / (4 - 2 * Var.tau)
    tau <- 1 / rinvGauss(n, lambda, lambda / (lambda - 1))
    
    # 计算Gamma分布的速率参数
    rate.para <- alpha / (mu * tau)
    
    # 生成样本x
    x <- rgamma(n, alpha, rate.para)
    return(x)
  }
  
  g <- 1
  G <- 53
  
  boot.estimation <- list(Alpha.boot = numeric(length = G), 
                          Mu.boot = numeric(length = G), 
                          Var.tau.boot = numeric(length = G))
  
  while (g <= G) {
    boot.sample <- rRIGGa(n, alpha, mu, Var.tau)
    Res <- tryCatch({
      RIGGa.MLE(alpha, mu, Var.tau, boot.sample)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      boot.estimation[[1]][g] <- Res$Alpha
      boot.estimation[[2]][g] <- Res$Mu
      boot.estimation[[3]][g] <- (2 * Res$Lambda - 1) * (Res$Lambda - 1) / Res$Lambda^2
      g <- g + 1
    }
  }
  
  return(boot.estimation)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Boot.RIGauGa(n, Res3$Alpha, Res3$Mu, (2 * Res3$Lambda - 1) * (Res3$Lambda - 1) / Res3$Lambda^2)
stopCluster(cl)

Result.alpha <- result_part[[1]][[1]]
Result.mu <- result_part[[1]][[2]]
Result.var.tau <- result_part[[1]][[3]]

for (i in 2:num_cores) {
  Result.alpha <- c(Result.alpha, result_part[[i]][[1]])
  Result.mu <- c(Result.mu, result_part[[i]][[2]])
  Result.var.tau <- c(Result.var.tau, result_part[[i]][[3]])
}

mean.std.CI(Result.alpha %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.mu %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.var.tau %>% as.matrix()) %>% round(digits = 4)

#===============================================================================
# 渐近正态CI
mle_ig_optim <- function(x, mu_init, lambda_init) {
  
  neg_loglik <- function(par) {
    mu <- par[1]
    lambda <- par[2]
    
    if (mu <= 0 || lambda <= 0) return(Inf)
    
    n <- length(x)
    ll <- n / 2 * log(lambda) - n / 2 * log(2 * pi) - 3 / 2 * sum(log(x)) - lambda / 2 / mu^2 * sum((x - mu)^2 / x)
    
    return(-ll)
  }
  
  fit <- optim(
    par = c(mu_init, lambda_init),
    fn = neg_loglik,
    hessian = TRUE,
    method = "BFGS"
  )
  
  theta_hat <- fit$par
  
  # 协方差矩阵
  vcov_hat <- solve(fit$hessian)
  
  # 标准误
  se_hat <- sqrt(diag(vcov_hat))
  
  # 95% Wald CI
  CI_wald <- cbind(
    lower = theta_hat - 1.96 * se_hat,
    upper = theta_hat + 1.96 * se_hat
  )
  
  list(
    mu_hat = fit$par[1],
    lambda_hat = fit$par[2],
    logLik = -fit$value,
    convergence = fit$convergence, 
    CI = CI_wald
  )
}
mle_ig_optim(x, 1, 1)

# Bootstrap CI
Boot.ig <- function(n, mu, lambda) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  mle_ig_optim <- function(x, mu_init, lambda_init) {
    
    neg_loglik <- function(par) {
      mu <- par[1]
      lambda <- par[2]
      
      if (mu <= 0 || lambda <= 0) return(Inf)
      
      n <- length(x)
      ll <- n / 2 * log(lambda) - n / 2 * log(2 * pi) - 3 / 2 * sum(log(x)) - lambda / 2 / mu^2 * sum((x - mu)^2 / x)
      
      return(-ll)
    }
    
    fit <- optim(
      par = c(mu_init, lambda_init),
      fn = neg_loglik,
      method = "L-BFGS-B",
      lower = c(1e-6, 1e-6)
    )
    
    list(
      mu_hat = fit$par[1],
      lambda_hat = fit$par[2],
      logLik = -fit$value,
      convergence = fit$convergence
    )
  }
  
  g <- 1
  G <- 100
  
  boot.estimation <- list(Mu.boot = numeric(length = G), 
                          Lambda.boot = numeric(length = G))
  
  while (g <= G) {
    boot.sample <- rinvGauss(n, mu, lambda)
    Res <- tryCatch({
      mle_ig_optim(boot.sample, mu, lambda)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      boot.estimation[[1]][g] <- Res$mu_hat
      boot.estimation[[2]][g] <- Res$lambda_hat
      g <- g + 1
    }
  }
  
  return(boot.estimation)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Boot.ig(n, Res5$mu_hat, Res5$lambda_hat)
stopCluster(cl)

Result.mu <- result_part[[1]][[1]]
Result.lambda <- result_part[[1]][[2]]

for (i in 2:num_cores) {
  Result.mu <- c(Result.mu, result_part[[i]][[1]])
  Result.lambda <- c(Result.lambda, result_part[[i]][[2]])
}

mean.std.CI(Result.mu %>% as.matrix()) %>% round(digits = 4)
mean.std.CI(Result.lambda %>% as.matrix()) %>% round(digits = 4)
