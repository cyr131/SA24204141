## ----warning = FALSE----------------------------------------------------------
library(MASS)
library(bindata)
library(SA24204141)
  ## settings
  TT <- 10   # 10 or 20
  n <- 100    # 180 or 270
  ID <- 1:n
  N <- n*TT

  ## covariates
  p <- 3   # number of covariates
  CC <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
  zz <- mvrnorm(n*TT, c(0, 0), CC)

  x1 <- matrix(zz[,1], n, TT)
  x2 <- matrix(zz[,2], n, TT)
  X0 <- array(NA, c(n, TT, p))
  X0[,,1] <- 1
  X0[,,2] <- x1
  X0[,,3] <- x2

  ## regression coefficients
  Beta <- c(1,-2,1)   # true regression coefficients


  tMu <- matrix(NA, n, TT)
  for(i in 1:n){
    tMu[i,] <- as.vector(X0[i,,]%*%Beta)
  }
  Prob <- logit(tMu)

  ## true correlation structure (equi-correlation)
  tR <- matrix(0.5, TT, TT)
  diag(tR) <- 1

  ## data generation
  Y0 <- matrix(NA, n, TT)
  for(i in 1:n){
    Y0[i,] <- as.vector(rmvbin(1, margprob=Prob[i,], sigma=tR))
  }

## -----------------------------------------------------------------------------
function(YY, XX, RR, init.beta=NULL, maxit=100){
  # 赋值参数
  th <- 0.0001
  ep <- 10^(-5)
  n <- dim(YY)[1]
  TT <- dim(YY)[2]
  p <- dim(XX)[3]
  if(is.null(init.beta)){ hB <- rep(0, p) }
  else{ hB <- init.beta }
  ## GEE求解
  # Iteration次数
  for(k in 1:maxit){
    hB0 <- hB
    mat1 <- 0
    mat2 <- 0
    # 遍历n个subjects
    for(i in 1:n){
      Mu <- logit(as.vector(XX[i,,]%*%hB))
      S <- YY[i,]-Mu
      vv <- Mu*(1-Mu)
      vv[vv<th] <- th
      A <- diag(vv)
      V <- sqrt(A)%*%RR%*%sqrt(A) + th*diag(TT)
      D <- A%*%XX[i,,]
      mat1 <- mat1+t(D)%*%ginv(V)%*%D
      mat2 <- mat2+t(D)%*%ginv(V)%*%S
    }
    #牛顿法求解
    hB <- hB+as.vector( ginv(mat1)%*%mat2 )
    dd <- sum( abs(hB-hB0) ) / sum( abs(hB0)+0.0001 )
    if( dd<ep ){ break }
  }
  return(hB)
}

## -----------------------------------------------------------------------------
#example
GEE(Y0,X0,tR)

## -----------------------------------------------------------------------------
function(y, x, cor = "EC", init.beta = NULL, maxit = 1000) {
  m <- dim(x)[3]
  # cor == "ID"
  if (cor == "ID") {
    elr_id <-  function(b) {  #b输入为列向量
      m <- dim(x)[3]
      n <- dim(x)[1]    #注意x,y来源
      t <- dim(x)[2]
      g <- matrix(0, n, m)  #n对应i的个数，m对应b的维度3
    #求解经验似然比
      error <- matrix(0,t,n)
      D <- array(0,c(t,3,n))
      A <- array(0,c(t,t,n))
      for (i in 1:n) {
        fitted <- plogis( x[i,,] %*% b )
        error[,i] <- y[i,]-fitted
        D[,,i] <- matrix(rep(fitted*(1-fitted),3),t,3)*x[i,,]
        A[,,i] <- diag(as.vector((fitted*(1-fitted))^(-1)))
      }
      for (i in 1:n) {
        A.half <- A[,,i]^(1/2)
        g[i,] <- crossprod(D[,,i], A.half %*% A.half %*% error[,i])
      }
      g.mu <- rep(0, m)
      el.test(g, g.mu, gradtol = 1e-9)$"-2LLR"
    }
    if(is.null(init.beta)){ hB <- rep(0,m) }
    else{ hB <- init.beta }
    #极大化经验似然比
    result <- optim(par = hB, fn = elr_id, method = "Nelder-Mead",control = list(maxit = maxit))
    return(list(coff=result$par,value=result$value))
  }
  # cor == "EC"
  if (cor == "EC") {
    elr_ec <- function(b) {  # b输入为列向量
      m <- dim(x)[3]
      n <- dim(x)[1]  # 注意x,y格式
      t <- dim(x)[2]
      # n对应i的个数，m对应b的维度3
      g <- matrix(0, n, m)
      theta <- matrix(0, n, t)
      for (i in 1:n) {
        theta[i, ] <- x[i, , ] %*% b
      }
      #求解参数alpha，得到EC结构阵：R
      e <- matrix(0, n, t)
      e <- (y - 1 / (1 + exp(-theta))) / sqrt(1 / (2 + exp(theta) + exp(-theta)))
      alpha <- 0
      for (i in 1:n) {
        for (j in 1:t) {
          for (k in 1:t) {
            if (k != j) alpha <- alpha + e[i, j] * e[i, k]
          }
        }
      }
      alpha <- alpha / (n*t*(t-1)-m)  # 矩估计由b求alpha
      alpha <- min(0.95, alpha)
      alpha <- max(0, alpha)
      R <- exchCpp(t, alpha)
      #求解经验似然比
      error <- matrix(0, t, n)
      D <- array(0, c(t, 3, n))
      A <- array(0, c(t, t, n))
      for (i in 1:n) {
        fitted <- plogis(x[i, , ] %*% b)
        error[, i] <- y[i, ] - fitted
        D[, , i] <- matrix(rep(fitted * (1 - fitted), 3), t, 3) * x[i, , ]
        A[, , i] <- diag(as.vector((fitted * (1 - fitted))^(-1)))
      }
      for (i in 1:n) {
        A.half <- A[, , i]^(1/2)
        g[i, ] <- crossprod(D[, , i], A.half %*% solve(R, A.half %*% error[, i]))
      }

      g.mu <- rep(0, m)
      el.test(g, g.mu, gradtol = 1e-9)$"-2LLR"
    }
    if(is.null(init.beta)){ hB <- numeric(m) }
    else{ hB <- init.beta }
    #极大化经验似然比
    result <- optim(par = hB, fn = elr_ec, method = "Nelder-Mead",control = list(maxit = maxit))
    return(list(coff=result$par,value=result$value))
  }
  # cor == "AR"
  if (cor == "AR") {
    elr_ar <- function(b) {  # b输入为列向量
      m <- dim(x)[3]
      n <- dim(x)[1]    #注意x,y来源
      t <- dim(x)[2]
      #n对应i的个数，m对应b的维度3
      g <- matrix(0, n, m)
      theta<-matrix(0,n,t)
      for (i in 1:n) {
        theta[i, ] <- x[i,,] %*% b
      }
      e<-matrix(0,n,t)
      e <- (y - 1 / (1 + exp(-theta))) / sqrt(1 / (2 + exp(theta) + exp(-theta)))
      #求解参数alpha，得到AR结构阵：R
      alpha<-0
      for (i in 1:n) {
        for (j in 1:(t-1)) {
          alpha<-alpha+e[i,j]*e[i,j+1]
        }}
      alpha<-alpha/((n*(t-1)-m))  #矩估计由b求alpha
      alpha <- min(0.95, alpha)
      alpha <- max(0, alpha)
      R<-ar1Cpp(t,alpha)
      #经验似然比
      error <- matrix(0,t,n)
      D <- array(0,c(t,3,n))
      A <- array(0,c(t,t,n))
      for (i in 1:n) {
        fitted <- plogis( x[i,,] %*% b )
        error[,i] <- y[i,]-fitted
        D[,,i] <- matrix(rep(fitted*(1-fitted),3),t,3)*x[i,,]
        A[,,i] <- diag(as.vector((fitted*(1-fitted))^(-1)))
      }
      for (i in 1:n) {
        A.half <- A[,,i]^(1/2)
        g[i,] <- crossprod(D[,,i], A.half %*% solve(R, A.half %*% error[,i]))
      }
      g.mu <- rep(0, m)
      el.test(g, g.mu, gradtol = 1e-9)$"-2LLR"
    }
    if(is.null(init.beta)){ hB <- numeric(m) }
    else{ hB <- init.beta }
    #极大化经验似然比得到所需参数估计
    result <- optim(par = hB, fn = elr_ar, method = "Nelder-Mead",control = list(maxit = maxit))
    return(list(coff=result$par,value=result$value))
  }
  # cor == "ST"
  if (cor == "ST") {
    elr_st <- function(b) {  # b输入为列向量
      m <- dim(x)[3]
      n <- dim(x)[1]    #注意x,y来源
      t <- dim(x)[2]
      g <- matrix(0, n, m)  #n对应i的个数，m对应b的维度3

      theta<-matrix(0,n,t)
      for (i in 1:n) {
        theta[i, ] <- x[i,,] %*% b
      }
      e<-matrix(0,n,t)
      e <- (y - 1 / (1 + exp(-theta))) / sqrt(1 / (2 + exp(theta) + exp(-theta)))
      # 求解ST矩阵的参数向量alpha
      alpha<-rep(0,t-1)
      for(k in 1:(t-1)){
        for (i in 1:n) {
          for (j in 1:(t-k)) {
            alpha[k]<-alpha[k]+e[i,j]*e[i,j+k]
          }}
        alpha[k]<-alpha[k]/((n*(t-k)-m))  #矩估计由b求alpha
        alpha[k] <- min(0.95, alpha[k])
        alpha[k] <- max(0, alpha[k])
      }
      #得到对应ST结构阵
      R<-statCpp(t,alpha)
      #求解经验似然比
      error <- matrix(0,t,n)
      D <- array(0,c(t,3,n))
      A <- array(0,c(t,t,n))
      for (i in 1:n) {
        fitted <- logit(as.vector(x[i,,] %*% b))
        error[,i] <- y[i,]-fitted
        D[,,i] <- matrix(rep(fitted*(1-fitted),3),t,3)*x[i,,]
        A[,,i] <- diag(as.vector((fitted*(1-fitted))^(-1)))
      }
      for (i in 1:n) {
        A.half <- A[,,i]^(1/2)
        g[i,] <- crossprod(D[,,i], A.half %*% solve(R, A.half %*% error[,i]))
      }
      g.mu <- rep(0, m)
      el.test(g, g.mu, gradtol = 1e-9)$"-2LLR"
    }
    if(is.null(init.beta)){ hB <- numeric(m) }
    else{ hB <- init.beta }
    #极大化经验似然比
    result <- optim(par = hB, fn = elr_st, method = "Nelder-Mead",control = list(maxit = maxit))
    return(list(coff=result$par,value=result$value))
  }
}

## -----------------------------------------------------------------------------
#example
GEE_EL(Y0,X0,"ID")
GEE_EL(Y0,X0,"EC")
GEE_EL(Y0,X0,"AR")
GEE_EL(Y0,X0,"ST")

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix exchCpp(int t, double alpha) {
#    NumericMatrix exch(t, t);
#  
#    // 初始化对角线元素为1
#    for (int i = 0; i < t; ++i) {
#      exch(i, i) = 1.0;
#    }
#  
#    // 设置非对角线元素为alpha
#    for (int i = 0; i < t; ++i) {
#      for (int j = 0; j < t; ++j) {
#        if (i != j) {
#          exch(i, j) = alpha;
#        }
#      }
#    }
#  
#    return exch;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix ar1Cpp(int t, double alpha) {
#     NumericMatrix ar1(t, t);
#     for(int i = 0; i < t; ++i) {
#       for(int j = i; j < t; ++j) {
#         if(i == j) {
#           ar1(i, j) = 1;
#         } else {
#           ar1(i, j) = pow(alpha, j - i);
#           ar1(j, i) = ar1(i, j);
#         }
#       }
#     }
#     return ar1;
#   }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix statCpp(int t, NumericVector a) {
#    int m = a.size();
#    NumericMatrix stat(t, t);
#  
#    // 设置对角线元素为1
#    for(int i = 0; i < t; ++i) {
#      stat(i, i) = 1;
#    }
#  
#    // 设置非对角线元素
#    for(int i = 0; i < t; ++i) {
#      for(int j = 0; j < t; ++j) {
#        if(abs(i - j) > 0 && abs(i - j) <= m) {
#          stat(i, j) = a[abs(i - j) - 1];
#        }
#      }
#    }
#  
#    return stat;
#  }

## -----------------------------------------------------------------------------
function(x){
  #计算logistic函数，允许x为数值，向量或矩阵
  1/(1+exp(-x))}

