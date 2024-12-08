#' @title A package used for illustration.
#' @name Illustration
#' @import MASS
#' @import emplik
#' @import bindata
#' @import ggplot2
#' @import mvtnorm
#' @import coin
#' @import boot
#' @import bootstrap
#' @import DAAG
#' @import lpSolve
#' @import Rcpp
#' @import microbenchmark
NULL

#' @importFrom Rcpp evalCpp
#' @useDynLib SA24204141
NULL

#' @title logistic function
#' @description Convert real numbers to numbers between (0,1) using logistic functions
#' @param x independent variable (Matrix)
#' @return the value of the logistic function corresponding to x
#' @examples
#' \dontrun{
#'  tMu <- matrix(c(0.3,-3,1,-0.5), 2, 2)
#'  logit(tMu)
#' }
#' @export
logit <- function(x){
  #计算logistic函数，允许x为数值，向量或矩阵
  1/(1+exp(-x))}

#' @title GEE function
#' @description Estimates of regression coefficients are obtained by solving generalized estimating equations
#' @param YY Response Variables (Matrix)
#' @param XX Independent variables (3-dimensional arrays)
#' @param RR Designated working correlation matrix (Matrix)
#' @param init.beta initial value of iteration (vector)
#' @param maxit Maximum number of iteration (int)
#' @return final coefficient estimates
#' @examples
#' \dontrun{
#' #data generation
#' CC <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
#' zz <- mvrnorm(1000, c(0, 0), CC)
#' x1 <- matrix(zz[,1], 100, 10)
#' x2 <- matrix(zz[,2], 100, 10)
#' X0 <- array(NA, c(100, 10, 3))
#' X0[,,1] <- 1
#' X0[,,2] <- x1
#' X0[,,3] <- x2
#' Beta <- c(1,-2,1)
#' tMu <- matrix(NA, 100, 10)
#' for(i in 1:100){
#'   tMu[i,] <- as.vector(X0[i,,]%*%Beta)
#' }
#' Prob <- logit(tMu)
#' tR <- matrix(0.5, 10, 10)
#' diag(tR) <- 1
#' Y0 <- matrix(NA, 100, 10)
#' for(i in 1:100){
#'   Y0[i,] <- as.vector(rmvbin(1, margprob=Prob[i,], sigma=tR))
#' }
#' #use GEE function
#' GEE(Y0,X0,tR)
#' }
#' @export
GEE <- function(YY, XX, RR, init.beta=NULL, maxit=100){
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

#' @title GEE_EL function
#' @description Estimates of regression coefficients are obtained by maximizing the empirical likelihood ratio
#' @param y Response Variables (Matrix)
#' @param x Independent variables (3-dimensional arrays)
#' @param cor Designated the type of working correlation;Optional types are:ID、EC、AR、ST
#' @param init.beta initial value of iteration (vector)
#' @param maxit Maximum number of iteration (int)
#' @return final coefficient estimates
#' @examples
#' \dontrun{
#' #data generation
#' CC <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
#' zz <- mvrnorm(1000, c(0, 0), CC)
#' x1 <- matrix(zz[,1], 100, 10)
#' x2 <- matrix(zz[,2], 100, 10)
#' X0 <- array(NA, c(100, 10, 3))
#' X0[,,1] <- 1
#' X0[,,2] <- x1
#' X0[,,3] <- x2
#' Beta <- c(1,-2,1)
#' tMu <- matrix(NA, 100, 10)
#' for(i in 1:100){
#'   tMu[i,] <- as.vector(X0[i,,]%*%Beta)
#' }
#' Prob <- logit(tMu)
#' tR <- matrix(0.5, 10, 10)
#' diag(tR) <- 1
#' Y0 <- matrix(NA, 100, 10)
#' for(i in 1:100){
#'   Y0[i,] <- as.vector(rmvbin(1, margprob=Prob[i,], sigma=tR))
#' }
#' #use GEE function
#' GEE_EL(Y0,X0,cor = "EC")
#' }
#' @importFrom stats plogis optim
#' @import Rcpp
#' @import emplik
#' @export
GEE_EL <- function(y, x, cor = "EC", init.beta = NULL, maxit = 1000) {
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
