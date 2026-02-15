
# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, n, alpha, mu, Var.tau) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(SuppDists)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(numDeriv)))
  suppressWarnings(suppressMessages(library(MASS)))
  
  # 生成IGaGa随机数===============================================================
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
  
  # NEM算法估计MLE================================================================
  IGauGa.MLE <- function(alpha, mu, Var.tau, x, ep = 1e-8) {
    # 收敛指标与迭代次数
    index <- 0
    k <- 1
    
    # 样本量n
    n <- length(x)
    
    # 确定τ分布的参数λ
    lambda <- 1 / Var.tau
    
    # 定义对数似然函数
    
    
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
  
  # 计算θ的得分函数贡献 (每个观测的梯度)
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
  
  # 开始Simulation
  k <- 1
  
  result <- list(Alpha = numeric(length = number), 
                 Mu = numeric(length = number),
                 Var.tau = numeric(length = number), 
                 Loglikelihood = numeric(length = number), 
                 CR = matrix(nrow = number, ncol = 3), 
                 Width = matrix(nrow = number, ncol = 3))
  
  while(k <= number) {
    x <- rIGauGa(n, alpha, mu, Var.tau)
    
    Res <- tryCatch({
      IGauGa.MLE(alpha, mu, Var.tau, x)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    
    if (!all(is.na(Res))) {
      if (Res$index == 1) {
        result[[1]][k] <- Res$Alpha
        result[[2]][k] <- Res$Mu
        result[[3]][k] <- 1 / Res$Lambda
        result[[4]][k] <- Res$Loglikelihood
        
        theta_hat <- c(Res$Alpha, Res$Mu, Res$Lambda)
        
        theta_true <- c(alpha, mu, Var.tau)
        
        S <- score_contributions_theta(theta_hat, x)
        B_matrix <- crossprod(S)
        H <- -hessian_theta(theta_hat, x)
        Sigma_theta <- solve(H) %*% B_matrix %*% solve(H)
        se <- diag(Sigma_theta) %>% sqrt()
        
        lower.bound <- theta_hat - qnorm(0.975) * se
        upper.bound <- theta_hat + qnorm(0.975) * se
        
        lower.bound[3] <- 1 / Res$Lambda - qnorm(0.975) * se[3] / Res$Lambda^2
        upper.bound[3] <- 1 / Res$Lambda + qnorm(0.975) * se[3] / Res$Lambda^2
        
        result[[5]][k, ] <- ifelse(lower.bound <= theta_true & upper.bound >= theta_true, 1, 0)
        result[[6]][k, ] <- upper.bound - lower.bound
        k <- k + 1
      }
    }
  }
  
  return(result)
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 105
n <- 500
alpha <- 1
mu <- 2
Var.tau <- 0.2

# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), n, alpha, mu, Var.tau)
stopCluster(cl)

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)

# ==============================================================================================================================

# 汇总拟合结果并保存
Result_Alpha <- result_part[[1]]$Alpha
Result_Mu <- result_part[[1]]$Mu
Result_Var.tau <- result_part[[1]]$Var.tau
Result_Loglikelihood <- result_part[[1]]$Loglikelihood
Result_CR <- result_part[[1]]$CR
Result_Width <- result_part[[1]]$Width

for (i in 2:num_cores) {
  Result_Alpha <- c(Result_Alpha, result_part[[i]]$Alpha)
  Result_Mu <- c(Result_Mu, result_part[[i]]$Mu)
  Result_Var.tau <- c(Result_Var.tau, result_part[[i]]$Var.tau)
  Result_Loglikelihood <- c(Result_Loglikelihood, result_part[[i]]$Loglikelihood)
  Result_CR <- rbind(Result_CR, result_part[[i]]$CR)
  Result_Width <- rbind(Result_Width, result_part[[i]]$Width)
}

result <- list(Alpha = Result_Alpha, 
               Mu = Result_Mu, 
               Var.tau = Result_Var.tau, 
               Loglikelihood = Result_Loglikelihood, 
               CR = Result_CR, 
               Width = Result_Width)

# ==============================================================================================================================

# 计算A-MLE和A-Moment
MLE.Alpha <- mean(result$Alpha)
MLE.Mu <- mean(result$Mu)
MLE.Var.tau <- mean(result$Var.tau)
MLE.Loglikelihood <- mean(result$Loglikelihood)
CI.CR <- colMeans(Result_CR)
CI.Width <- colMeans(Result_Width)

# ==============================================================================================================================

# 计算MSE
MSE.Alpha <- sum((result$Alpha - alpha)^2 / all_times)
MSE.Mu <- sum((result$Mu - mu)^2 / all_times)
MSE.Var.tau <- sum((result$Var.tau - Var.tau)^2 / all_times)
# ==============================================================================================================================

# 总结最终结果并输出
Final <- data.frame(Alpha.MLE = MLE.Alpha, 
                    Alpha.MSE = MSE.Alpha, 
                    Alpha.CR = CI.CR[1], 
                    Alpha.Width = CI.Width[1], 
                    Mu = MLE.Mu, 
                    Mu.MSE = MSE.Mu,
                    Mu.CR = CI.CR[2], 
                    Mu.Width = CI.Width[2], 
                    Var.tau = MLE.Var.tau, 
                    Var.tau.MSE = MSE.Var.tau, 
                    Var.tau.CR = CI.CR[3], 
                    Var.tau.Width = CI.Width[3], 
                    Loglikelihood = MLE.Loglikelihood)

Final %>% round(digits = 4)
# ==============================================================================================================================

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)
