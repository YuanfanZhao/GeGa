# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

Simulation <- function(number, size, alpha, mu, Var.tau, Ratio) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  
  # 生成GeGa随机数===============================================================
  rIGaGa <- function(n, alpha, mu, Var.tau) {
    # 确定τ分布的参数λ并生成τ分布的随机数
    lambda <- 1 / Var.tau + 2
    tau <- 1 / rgamma(n, lambda, lambda - 1)
    
    # 计算Gamma分布的速率参数
    rate.para <- alpha / (mu * tau)
    
    # 生成样本x
    x <- rgamma(n, alpha, rate.para)
    return(x)
  }
  
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
  
  # NEM算法估计MLE================================================================
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
  
  mom_weibull_optim <- function(mean_x, var_x) {
    
    # ratio of moments
    R <- var_x / mean_x^2
    
    # objective function: squared moment equation
    obj_fun <- function(k) {
      gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1 - R
    }
    
    fit <- uniroot(obj_fun, c(1, 2), extendInt = "yes")
    
    k_hat <- fit$root
    
    # back out lambda
    lambda_hat <- mean_x / gamma(1 + 1 / k_hat)
    
    return(c(k_hat, lambda_hat))
  }
  
  
  
  # 开始Simulation
  k <- 1
  Loglikelihood <- matrix(nrow = number, ncol = 6)
  
  while(k <= number) {
    Mean <- mu
    Variance <- mu^2 * ((1 + 1 / alpha) * Var.tau + 1 / alpha)
    
    alpha.ga <- Mean^2 / Variance
    beta.ga <- Mean / Variance
    mu.igau <- Mean
    lambda.igau <- Mean^3 / Variance
    mom.wei <- mom_weibull_optim(Mean, Variance)
    k.wei <- mom.wei[1]
    lambda.wei <- mom.wei[2]
    
    
    x1 <- rIGaGa(n, alpha, mu, Var.tau)
    x2 <- rIGauGa(n, alpha, mu, Var.tau)
    x3 <- rRIGGa(n, alpha, mu, Var.tau)
    x4 <- rgamma(n, alpha.ga, beta.ga)
    x5 <- rinvGauss(n, mu.igau, lambda.igau)
    x6 <- rweibull(n, k.wei, lambda.wei)
    z <- rmultinom(n, 1, Ratio)
    x <- x1 * z[1, ] + x2 * z[2, ] + x3 * z[3, ] + x4 * z[4, ] + x5 * z[5, ] + x6 * z[6, ]
    
    Res1 <- tryCatch({
      IGaGa.MLE(alpha, mu, Var.tau, x)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Res2 <- tryCatch({
      IGauGa.MLE(alpha, mu, Var.tau, x)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Res3 <- tryCatch({
      RIGGa.MLE(alpha, mu, Var.tau, x)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Res4 <- mle_gamma_optim(x, alpha.ga, beta.ga)
    Res5 <- mle_ig_optim(x, mu.igau, lambda.igau)
    Res6 <- mle_weibull_optim(x, k.wei, lambda.wei)
    
    if (!all(is.na(Res1)) & !all(is.na(Res2)) & !all(is.na(Res3))) {
      if (Res1$index == 1 & Res2$index == 1 & Res3$index == 1) {
        Loglikelihood[k, 1] <- -2 * Res1$Loglikelihood + 2 * 3
        Loglikelihood[k, 2] <- -2 * Res2$Loglikelihood + 2 * 3
        Loglikelihood[k, 3] <- -2 * Res3$Loglikelihood + 2 * 3
        Loglikelihood[k, 4] <- -2 * Res4$logLik + 2 * 2
        Loglikelihood[k, 5] <- -2 * Res5$logLik + 2 * 2
        Loglikelihood[k, 6] <- -2 * Res6$logLik + 2 * 2
        k <- k + 1
      }
    }
  }
  
  return(Loglikelihood)
}

num_cores <- detectCores()
all_times <- num_cores * 11
n <- 800
alpha <- 2
mu <- 3
Var.tau <- 0.5
# Ratio <- c(1, 0, 0, 0, 0, 0)
Ratio <- c(0, 1, 0, 0, 0, 0)
Ratio <- c(0, 0, 1, 0, 0, 0)
Ratio <- c(0.6, 0.2, 0.2, 0, 0, 0)
Ratio <- c(0.2, 0.6, 0.2, 0, 0, 0)
Ratio <- c(0.2, 0.2, 0.6, 0, 0, 0)
# Ratio <- c(1 / 9, 1 / 9, 1 / 9, 2 / 9, 2 / 9, 2 / 9)

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), n, alpha, mu, Var.tau, Ratio)
stopCluster(cl)

Result <- result_part[[1]]
for (i in 2:num_cores) {
  Result <- rbind(Result, result_part[[i]])
}

Rank <- apply(-Result, 1, rank) %>% t()

colMeans(Result, na.rm = TRUE) %>% round(digits = 1)
colMeans(Rank) %>% round(digits = 2)

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)
