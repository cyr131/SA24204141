## -----------------------------------------------------------------------------
R<-function(c,n=1000){
  set.seed(10)
  u<-runif(n)
  x<-sqrt(-2*c^2*log(1-u)) #逆变换生成所需分布
  #画直方图比较
  hist(x, prob = TRUE, ylim = c(0,0.6),main = paste("参数取值为",c))
  v <- seq(0, 10, .01)
  lines(v, v/(c^2)*exp(-v^2/(2*c^2)))
}  #看图可见所生随机数契合分布
R(1)
R(2)
R(1.7)
R(1.2)

## -----------------------------------------------------------------------------
mix<-function(p,n=1000){ #p=p1
  set.seed(10)
  x1<-rnorm(n,0,1)
  x2<-rnorm(n,3,1)
  a<-sample(0:1,n,replace = T,prob = c(1-p,p))
  x<-a*x1+(1-a)*x2
  hist(x,main = paste("p1取值为",p))
}
mix(0.75)#不显双峰

## -----------------------------------------------------------------------------
mix(0.9)#不显双峰
mix(0.7)#不显双峰
mix(0.6)#显双峰
mix(0.5)#显双峰
mix(0.45)#不明显
mix(0.3)#不显双峰

## -----------------------------------------------------------------------------
p_g<-function(lambda,alpha,beta,t=10,n=10000){
  set.seed(100)
  Xt <- numeric(n)
  for (i in 1:n) {
    Nt <- rpois(1, lambda * t) # 生成 N(t)
    Yi <- rgamma(Nt, alpha, beta) # 生成伽马分布的随机变量
    Xt[i] <- sum(Yi) # 计算 X(t)
  } # 模拟复合泊松-伽马过程,生成Xt

# 估计均值和方差
estimated_mean <- mean(Xt)
estimated_var <- var(Xt)
# 理论值
theoretical_mean <- lambda * t * (alpha / beta)
theoretical_var <- lambda * t * ((alpha / beta)^2 + (alpha / beta^2))
# 输出结果
cat("估计的 X(10) 均值:", estimated_mean, "\n")
cat("估计的 X(10) 方差:", estimated_var, "\n")
cat("理论的 X(10) 均值:", theoretical_mean, "\n")
cat("理论的 X(10) 方差:", theoretical_var, "\n")
}  

p_g(2,4,2)#参数选择：lambda:2,gamma:4,beta:2
p_g(1,3,1/2)#参数选择：lambda:1,gamma:3,beta:1/2
p_g(5,2,3)#参数选择：lambda:5,gamma:2,beta:3

## -----------------------------------------------------------------------------
Fbeta<-function(x,n=10000){
  set.seed(123)
  u<-runif(n,0,x)
  return(30*x*mean(u^2*(1-u)^2))
} #估计密度的函数
for (i in 1:9) {
  print(paste0("x=",i/10,"时,","估计值为",Fbeta(i/10),";pbeta值为",pbeta(i/10,3,3)))
}

## -----------------------------------------------------------------------------
f_duiou<-function(sigma,n=1000){
  set.seed(123)
  U<-runif(n/2)
  x1<-sqrt(-2*sigma^2*log(1-U))
  x2<-sqrt(-2*sigma^2*log(U))
  return(c(x1,x2))
}  #对偶生成
f<-function(sigma,n=1000){
  set.seed(123)
  U<-runif(n)
  x<-sqrt(-2*sigma^2*log(1-U))
  return(x)
}  #独立生成

#比较方差
for (i in c(0.5,1,5)) {
  a<-round(var(f_duiou(i)),3)
  b<-round(var(f(i)),3)
  print(paste0("sigma=",i,"时,","对偶生成方差为",a,";独立生成方差值为",b,"相比方差减少:",round(((b-a)/b)*100,3),"%"))
}

## -----------------------------------------------------------------------------
n<-1000
set.seed(123)
a1<-f(1,5000)
a2<-rnorm(15000)
a1<-a1[a1>1][1:n]
a2<-a2[a2>1][1:n]
g1<-a1/sqrt(2*pi*exp(1))
g2<-0.1587*a2^2
paste0("f1为重要函数的估计均值为：",mean(g1))
paste0("f2为重要函数的估计均值为",mean(g2))
paste0("f1为重要函数的估计方差为：",var(g1))
paste0("f2为重要函数的估计方差为",var(g2))


## -----------------------------------------------------------------------------
library(ggplot2)
set.seed(123)
#快排函数
quick <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    pivot <- x[1]
    less <- quick(x[x < pivot])
    greater <- quick(x[x > pivot])
    return(c(less, pivot, greater))
  }
}

# 定义计算排序时间的函数
time <- function(n) {
  numbers <- sample(1:n)
  a<-Sys.time()
  quick(numbers)
  b<-Sys.time()
  return(b-a)
}
# 定义n的值
n <- c(10^4, 2*10^4, 4*10^4, 6*10^4, 8*10^4)

# 初始化存储计算时间的向量
av_time <- c()

# 对每个n值进行100次模拟
for (i in 1:length(n)) {
  n1 <- n[i]
  times <- replicate(100, time(n1))
  av_time[i] <- mean(times)
}

# 计算n log(n)
tn <- n * log(n)

# 绘制散点图和回归线
ggplot(data.frame(tn, av_time), aes(x = tn, y = av_time)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Scatter Plot and Regression Line",
       x = "n log(n)",
       y = "Ave_Time (an)",
       caption = "Regression model: an ~ tn") +
  theme_minimal()

## -----------------------------------------------------------------------------
n<-1000 #样本量
m<-1000 #模拟次数
set.seed(10)
#计算样本偏度
b<-numeric(m)
for (i in 1:m) {
  x<-rnorm(n)
  b[i]<-mean(((x-mean(x))/sd(x))^3)
}
#分位数估计
qvalue<-quantile(b, probs = c(0.025,0.05,0.95,0.975))
qvalue

## -----------------------------------------------------------------------------
# Compute the standard error
q<-c(0.025,0.05,0.95,0.975)
sd<-sqrt(q*(1-q)/(m*(dnorm(quantile(b,probs = q),0,sqrt(6/n)))^2))
print("the standard error of the estimates:")
sd

## -----------------------------------------------------------------------------
# Compare the estimated quantiles with the quantiles of the large sample approximation
print("the estimated quantiles:")
qvalue
print("the quantiles of the large sample approximation:")
qnorm(c(0.025, 0.05, 0.95, 0.975), mean = 0, sd = sqrt(6/n))

## -----------------------------------------------------------------------------
library(MASS)
library(mvtnorm)

## -----------------------------------------------------------------------------
set.seed(10)
# 生成双变量正态样本
n <- 1000 
mu <- c(0, 0)  # 均值
sd <- matrix(c(1, 0.7, 0.7, 1), 2, 2)  # 协方差
data <- mvrnorm(n, mu, sd)
X <- data[, 1]
Y <- data[, 2]
# 相关性检验
pea <- cor.test(X, Y, method = "pearson")
spe <- cor.test(X, Y, method = "spearman")
ken <- cor.test(X, Y, method = "kendall")

print(paste0("Pearson product moment correlation:",round(pea$estimate,3)))
print(paste0("Spearman’s rank correlation coefficient:",round(spe$estimate,3)))
print(paste0("Kendall’s coefficient:",round(ken$estimate,3)))

## -----------------------------------------------------------------------------
df <- 10
mean <- c(0, 0)
Sigma <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
# 生成样本
set.seed(10) 
n <- 10  # 样本大小
data <- rmvt(n, delta = mean, sigma = Sigma, df = df)
X <- data[, 1]
Y <- data[, 2]
# 相关性检验
pea <- cor.test(X, Y, method = "pearson")
spe <- cor.test(X, Y, method = "spearman")
ken <- cor.test(X, Y, method = "kendall")

print(paste0("Pearson product moment correlation:",round(pea$estimate,3)))
print(paste0("Spearman’s rank correlation coefficient:",round(spe$estimate,3)))
print(paste0("Kendall’s coefficient:",round(ken$estimate,3)))



## -----------------------------------------------------------------------------
X <- data[, 1]
Y <- data[, 1]^3+data[,1]
# 相关性检验
pea <- cor.test(X, Y, method = "pearson")
spe <- cor.test(X, Y, method = "spearman")
ken <- cor.test(X, Y, method = "kendall")

print(pea)
print(spe)
print(ken)

## -----------------------------------------------------------------------------
library(coin)
# 代码举例：同时落入两种方法的拒绝域的次数400
da <- matrix(c(400, 251, 276, 73), nrow = 2, byrow = TRUE,
                dimnames = list(Method1 = c("Positive", "Negative"),
                                Method2 = c("Positive", "Negative")))
print(da)
# 进行McNemar检验
mcnemar.test(da)

## -----------------------------------------------------------------------------
set.seed(123)
h0<-numeric(950)
h1<-numeric(50)
m<-10000
FWER1<-FDR1<-TPR1<-numeric(m)
FWER2<-FDR2<-TPR2<-numeric(m)
for (i in 1:m) {
  h0<-runif(950)
  h1<-rbeta(50,0.1,1)
  h<-c(h0,h1)
  h_bon<-p.adjust(h,method='bonferroni')
  h_bh<-p.adjust(h,method='fdr')
  #计算FWER, FDR, and TPR of bonferroni
  h0_bon<-p.adjust(h0,method='bonferroni')
  FWER1[i]<-ifelse(any(h0_bon<0.1),1,0)
  FDR1[i]<-sum(h_bon[1:950]<0.1)/sum(h_bon<0.1)
  TPR1[i]<-sum(h_bon[951:1000]<0.1)/50
  #计算FWER, FDR, and TPR of B-H
  h0_bh<-p.adjust(h0,method='fdr')
  FWER2[i]<-ifelse(any(h0_bh<0.1),1,0)
  FDR2[i]<-sum(h_bh[1:950]<0.1)/sum(h_bh<0.1)  #####
  TPR2[i]<-sum(h_bh[951:1000]<0.1)/50   #####
}

# 将结果整理成3*2表格
results <- matrix(round(c(mean(FWER1),mean(FWER2),mean(FDR1),mean(FDR2),mean(TPR1),mean(TPR2)),4), nrow = 3, ncol = 2,byrow = T)
colnames(results) <- c("Bonferroni校正", "B-H校正")
rownames(results) <- c("FWER", "FDR", "TPR")
print(results)

## -----------------------------------------------------------------------------
library(boot)
set.seed(123)
x<-c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
mle<-function(x,i) 1/mean(x[i])
obj <- boot(data=x,statistic=mle,R=1e4)
print(paste0("lamda的MLE估计为：",round(obj$t0,4)))
print(paste0(" estimate the bias of the estimate：",round(mean(obj$t)-obj$t0,4)))
print(paste0(" estimate the standard error of the estimate：",round(sd(obj$t),4)))

## -----------------------------------------------------------------------------
me<-function(x,i) mean(x[i])
de <- boot(data=x,statistic=me,R=1e4)
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
print("四种CI分别为(区间下界 区间上界)：")
ci$normal[2:3]
ci$basic[4:5]
ci$percent[4:5] 
ci$bca[4:5]

## -----------------------------------------------------------------------------
rm(list = ls())
library(bootstrap)
data1<-data(scor)

## -----------------------------------------------------------------------------
set.seed(123)
#求jackknife的函数
knife<-function(data0){
  n<-nrow(data0)
  ei<-matrix(NA,n,5)
  for (i in 1:n) {
    sigma<-cor(data0[-i,])
    ei[i,]<-eigen(sigma)$values
  } # ei为特征值矩阵
  jacktheta<-ei[,1]/(ei[,1]+ei[,2]+ei[,3]+ei[,4]+ei[,5])
  return(jacktheta)
}

## -----------------------------------------------------------------------------
#theta_hat
sigma<-cor(scor)
ei<-eigen(sigma)$values
theta<-ei[1]/sum(ei)
thetaj<-knife(scor)
#the jackknife estimates of bias 
(nrow(scor)-1)*(mean(thetaj)-theta)

#the jackknife estimates of standard error of θ
sqrt((nrow(scor)-1)*mean((thetaj-theta)^2))


## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits
#Linear
L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)
#Quadratic
L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)
#Exponential
L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)
#cubic
L4 <- lm(magnetic ~ chemical + I(chemical^2)+ I(chemical^3))
plot(log(chemical), log(magnetic), main="Cubic", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] * a +
  L4$coef[3] * a^2+ L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## -----------------------------------------------------------------------------
# cross validation procedure
n <- length(ironslag$magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)

# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- ironslag$magnetic[-k]
  x <- ironslag$chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * ironslag$chemical[k]
  e1[k] <- ironslag$magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * ironslag$chemical[k] + J2$coef[3] * ironslag$chemical[k]^2
  e2[k] <- ironslag$magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * ironslag$chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- ironslag$magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * ironslag$chemical[k] + J4$coef[3] * ironslag$chemical[k]^2 + J4$coef[4] * ironslag$chemical[k]^3
  e4[k] <- ironslag$magnetic[k] - yhat4
}

round(c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)),3)

## -----------------------------------------------------------------------------
#selected according to maximum adjusted R2?
summary(L1)
summary(L2)
summary(L3)
summary(L4)

## -----------------------------------------------------------------------------
# example 8.1
data("chickwts")
x <- sort(as.vector(chickwts[chickwts$feed == "soybean",1]))
m<-length(x)
y <- sort(as.vector(chickwts[chickwts$feed == "linseed",1]))
n<-length(y)

## -----------------------------------------------------------------------------
set.seed(123)
cvm_statistic <- function(x, y) {
  Fn <- ecdf(x)  # ECDF for sample x
  Gm <- ecdf(y)  # ECDF for sample y
  # Cramer-von Mises statistic
  s <- (n * m) / (n + m)^2 * (
    sum((Fn(x) - Gm(x))^2) + 
    sum((Fn(y) - Gm(y))^2)   
  )
  return(s)
}
cvm0<-cvm_statistic(x,y)
M<-1000
cvm<-numeric(M)
for (i in 1:M) {
  z<-sample(c(x,y))
  cvm[i]<-cvm_statistic(z[1:m],z[(m+1):(m+n)])
}
print(paste0("p值为：",round((sum(cvm>=cvm0)+1)/(M+1),3)))

## -----------------------------------------------------------------------------
#直方图
hist(cvm, main = "", xlab = "perm_stats(p=0.423)",freq = FALSE, breaks = "scott")
abline(v=cvm0,col="red") #observed W2

## -----------------------------------------------------------------------------
#生成数据
library(MASS)
# 定义均值向量
mean_vector <- c(0, 0)  # 假设有两个变量，均值为0
# 定义协方差矩阵
cov_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # 假设两个变量之间的相关性为0.5
# 生成20个样本
samples <- mvrnorm(20, mu = mean_vector, Sigma = cov_matrix)


## -----------------------------------------------------------------------------
x<-samples[,1]
y<-samples[,2]
# the observed Spearman correlation
cor0 <- cor(x, y, method = "spearman")
#  the permutation 
perm_test_spearman <- function(x, y, M = 1000) {
  permuted_corrs <- numeric(M)
  combined <- c(x, y)
  n = length(x)
  m = length(y)
  for (i in 1:M) {
    permuted_index <- sample(1:(n+m), n, FALSE)
    permuted_corrs[i] <- cor(combined[permuted_index], combined[-permuted_index], method = "spearman")
  }
  return(permuted_corrs)
}

M <- 1000
permuted_corrs <- perm_test_spearman(x, y, M)

#  p-value
p_value_perm <- mean(abs(permuted_corrs) >= abs(cor0))

# p-value from cor.test
spearman_test <- cor.test(x, y, method = "spearman")

# Results
cat("Observed Spearman correlation:", round(cor0,3), "\n")
cat("Permutation p-value:", round(p_value_perm,4), "\n")
cat("Spearman test p-value:", round(spearman_test$p.value,4), "\n")



## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
set.seed(123)
# M-H生成Cauchy
MH<-function(n,cut=1000){
  x<-numeric(n)
  x[1]<-rnorm(1,0,1)
  U<-runif(n)
  for (i in 2:n) {
    a<-rnorm(1,x[i-1],1)
    ifelse(U[i]<=((1+x[i-1]^2)/(1+a^2)),x[i]<-a,x[i]<-x[i-1])
  }
  x<-x[(cut+1):n]
  return(x)
}
sim_cauchy<-MH(5000)
#比较分位数
quantile(sim_cauchy, probs = seq(0.1, 0.9, by = 0.1)) #deciles of MH
qt(seq(0.1, 0.9, by = 0.1),1) #deciles of the standard Cauchy distribution


## -----------------------------------------------------------------------------
set.seed(123)
gibb<-function(m,cut,n=10,a=1,b=1){
  x<-numeric(m)
  y<-numeric(m)
  x[1]<-rbinom(1,n,1/2)
  y[1]<-rbeta(1,x[1]+a,n-x[1]+b)
  for (i in 2:m) {
    x[i]<-rbinom(1,n,y[i-1])
    y[i]<-rbeta(1,x[i-1]+a,n-x[i-1]+b)
  }
  x<-x[(cut+1):m]
  y<-y[(cut+1):m]
  f<-matrix(c(x,y),nrow = (m-cut),ncol = 2)
}
#. Use the Gibbs sampler to generate a chain with target joint density f(x, y)
# n=10,a=1,b=1
gib<-gibb(5000,1000)

# 绘制马氏链折线图
plot(gib[,1], type = "l", main = "f(x,y)", xlab = "index", ylab = "random")
lines(gib[,2],col="red")
legend("topleft", legend = c("x", "y"), col = c("black", "red"), lty = 1)

## -----------------------------------------------------------------------------
# the Gelman-Rubin
gr <- function(samples) {
  # samples: 一个矩阵，其中每一列代表一个链的样本，每一行代表一个参数
  # φ:mean
  k<-ncol(samples)
  n <- nrow(samples)
  fi<-matrix(NA,n,k)
  for (i in 2:n) {
    fi[i,]<-apply(samples[1:i,],2,mean)
  }  #计算得φ:mean
  fi[1,]<-samples[1,]

  B<-n/(k-1)*sum((apply(fi,2,mean)-mean(apply(fi,2,mean)))^2)
  W<-0
  t<-apply(fi,2,mean)
  for(i in 1:k){
  for (j in 1:n) {
    W<-W+(fi[j,i]-t[i])^2
  }}
  W<-W/(k*(n-1))
  R.hat <- sqrt(((n-1)/n*W+1/n*B) / W)
  return(R.hat)
}
#for 9.3
for (M in seq(4000,6000,100)) {
  set.seed(1)
  A1<-MH(M,0)
  set.seed(2)
  A2<-MH(M,0)
  set.seed(3)
  A3<-MH(M,0)
sam<-cbind(A1,A2,A3)
if(gr(sam)<1.2){
  print(paste0("M=",M,"时收敛；此时R_hat为：",gr(sam)))
  break
}
}



## -----------------------------------------------------------------------------
#对于9.8
for (M in seq(1000,5000,100)) {
  set.seed(1)
  A1<-gibb(M,0)
  set.seed(2)
  A2<-gibb(M,0)
  set.seed(3)
  A3<-gibb(M,0)
sam<-cbind((A1[,1]+A1[,2])/2,(A2[,1]+A2[,2])/2,(A3[,1]+A3[,2])/2)
if(gr(sam)<1.2){
  print(paste0("M=",M,"时收敛；此时R_hat为：",gr(sam)))
  break
}
}

## -----------------------------------------------------------------------------
#见附件

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
#第k项函数
kterm <- function(k, d, a) {
  # 量a的欧几里得范数
  norm_a <- sqrt(sum(a^2))
  # 第k项
  term<-(-1)^k*exp(log(norm_a^(2 * k + 2))+log(gamma(k + 3/2) * gamma((d + 1) / 2))-log(factorial(k) * (2^k))-log((2 * k + 1) * (2 * k + 2))-log(gamma(k + d/2 + 1)))
  
  return(term)
}

#求和函数
n_sum <- function(d, a, n = 100) {
  sum <- 0
  for (k in 0:n) {
    term <- kterm(k, d, a)
    sum <- sum + term
  }
  return(sum)
}

## -----------------------------------------------------------------------------
# Evaluate the sum when a = (1, 2)T
# large n and d
n_sum(d=100, a=c(1,2), n=100)

## -----------------------------------------------------------------------------
#画出更直观的曲线逼近、
d <- 100
a <- c(1, 2)
n_values <- seq(1, 20, by = 1) # 创建一个从1到20的n值向量
sum_values <- sapply(n_values, n_sum, d = d, a = a) # 计算每个n值的级数和
plot(n_values, sum_values, type = "l", xlab = "n", ylab = "Sum", main = "Sum of Series for d = 20 and a = (1, 2)")

## -----------------------------------------------------------------------------
# 计算积分的辅助函数
integrate_function <- function(k, a, c_k) {
  integrand <- function(u) {
    (1 + (u^2) / k)^(-(k+1)/2)
  }
  result <- integrate(integrand, lower = 0, upper = c_k)
  return(result$value)
}

# 计算方程左边的函数
left_side <- function(k, a) {
  c_k <- sqrt(a^2 * (k-1) / (k - a^2))
  integral_value <- integrate_function(k-1, a, c_k)
  (2 * gamma(k/2)) / (sqrt(pi * (k - 1)) * gamma((k - 1)/2)) * integral_value
}

# 计算方程右边的函数
right_side <- function(k, a) {
  c_k <- sqrt(a^2 * k / (k + 1 - a^2))
  integral_value <- integrate_function(k, a, c_k)
  (2 * gamma((k + 1)/2)) / (sqrt(pi * k) * gamma(k/2)) * integral_value
}

#求根
f <- function(a,k) left_side(k,a)-right_side(k,a)
res <- uniroot(f,c(0.5,4),k=25)

unlist(res)[1:3]
unlist(res)[-(1:3)]

## -----------------------------------------------------------------------------
res <- uniroot(f,c(0.5,1.8),k=4)
print(paste0("k=4时："))
unlist(res)[1:3]
res <- uniroot(f,c(0.5,1.8),k=15)
print(paste0("k=15时："))
unlist(res)[1:3]
res <- uniroot(f,c(0.5,4),k=25)
print(paste0("k=25时："))
unlist(res)[1:3]
res <- uniroot(f,c(0.5,5),k=100)
print(paste0("k=100时："))
unlist(res)[1:3]

## -----------------------------------------------------------------------------
# 定义EM算法的函数
myem <- function(Yi, tau=1, max_iter = 100, tol = 1e-6) {
  n <- length(Yi)
  # 初始化参数
  lambda <- 1/mean(Yi)
  
  # EM算法迭代
  for (i in 1:max_iter) {
    # E步：对数似然均值
    Yi[5]<-1+1/lambda
    Yi[6]<-1+1/lambda
    Yi[8]<-1+1/lambda
    
    # M步：更新lambda的估计
    lambda_new <- 1 / mean(Yi)
    
    # 检查收敛性
    if (abs(lambda_new - lambda) < tol) {
      break
    }
    lambda <- lambda_new
  }
  
  return(list(lambda = lambda, iterations = i))
}

# 观察到的数据Yi
Yi <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1

# 运行EM算法
em_result <- myem(Yi, tau)
print(em_result)

## -----------------------------------------------------------------------------
#MLE
#该数据对应的对数似然函数
Y1<-c(0.54, 0.48, 0.33, 0.43, 0.91, 0.21, 0.85)
fu<-function(x) 7*log(x)-sum(Y1)*x-3*x
optimize(fu,lower=0.5,upper=5,maximum=TRUE)


## -----------------------------------------------------------------------------
rm(list=ls())
library(lpSolve)

## -----------------------------------------------------------------------------
f.obj <- c(4, 2, 9)#objective.in
f.con <- matrix (c(2, 1, 1, 1, -1, 3), 
                 nrow=2, byrow=TRUE)#const.mat
f.dir <- c("<=", "<=")#const.dir
f.rhs <- c(2, 3)
#Vector of numeric values for the right-hand sides of the constraints.
print("x,y,z分别取：")
lp ("min", f.obj, f.con, f.dir, f.rhs)$solution
print("此时，目标函数达到最小值：")
lp ("min", f.obj, f.con, f.dir, f.rhs)

## -----------------------------------------------------------------------------
##loops
# mtcars数据集
data(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# 使用for语句循环
for (i in 1:4) {
  model<- lm(formulas[[i]], data = mtcars)
  print(model)
}



## -----------------------------------------------------------------------------
# lapply向量化
fun<-function(formula) {lm(formula, data = mtcars)}
fits <- lapply(formulas, fun)
print(fits)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

# loops
for (i in 1:10) {
  # 拟合模型
  fit <- lm(formulas[[1]], data = bootstraps[[i]])
  # 打印模型摘要
  print(paste0("Fit result for boot ", i, ":"))
  print(fit)
}

## -----------------------------------------------------------------------------
##lapply列表
#without an anonymous function
#...	:optional arguments to FUN
fits2 <- lapply(bootstraps, lm, formula = formulas[[1]])

print(fits2)


## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
## ex3
print("R^2 of ex3 model : ")
lapply(fits,rsq)


## -----------------------------------------------------------------------------
## ex4
print("R^2 of ex4 model : ")
lapply(fits2,rsq)

## -----------------------------------------------------------------------------
set.seed(123)
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)

#Use sapply() and an anonymous function
sapply(trials, function(x) x$p.value)


## -----------------------------------------------------------------------------
#Use sapply() and without anonymous function
#by using [[ directly.
sapply(trials, "[[", "p.value")

## -----------------------------------------------------------------------------
## 编写lapply变形函数--Lappy
#形参：
#FUN--the function to be applied to each element of X
#...--多个X
#FUN.VALUE	--a (generalized) vector; a template for the return value from FUN
Lapply<-function(FUN, ..., FUN.VALUE=numeric(1), sim = TRUE) {
  # 使用Map实现对...中的每个元的并行
  re <- Map(FUN, ...)
  # vapply
  res <- vapply(re, function(x) {
    # 使用FUN.VALUE[[1]]为模板
    if (length(FUN.VALUE) == 1) {
      return(x[[1]])
    } else {
      stop("FUN.VALUE must be a vector of length 1")
    }
  }, FUN.VALUE)
  
  # 如果simplify为TRUE且结果为一维，则简化输出为向量
  if (sim && is.list(res) && length(res) == 1) {
    return(res[[1]])
  }
  
  # 返回结果
  return(res)
}

## -----------------------------------------------------------------------------
##举例使用上述函数
my_fun <- function(x, y) {
  x ^ y
}

# 并行应用myFUN到两个向量上
result <- Lapply(my_fun, 1:3, 4:6, FUN.VALUE = numeric(1),sim = TRUE)

# 打印结果
print(result)
is.vector(result)|is.matrix(result)#检验结果是否向量/矩阵形式

## -----------------------------------------------------------------------------
#编写更快的检验函数
my_chisq.test <- function(x, y) {
  # 默认确保输入是数值型向量且无缺失值
  # 根据卡方检验公式，尽量向量化
  tab <- table(x, y)
  obs <- as.vector(tab)
  ##利用outer外积
  exp <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  # 卡方检验统计量
  result <- sum((obs - exp)^2 / exp)
  return(result)
}

## -----------------------------------------------------------------------------
##举例使用上述函数
set.seed(123)
x <- sample(1:7, 50, replace = TRUE)
y <- sample(1:7, 50, replace = TRUE)

# 计算卡方检验统计量
my_chisq.test(x, y)


## -----------------------------------------------------------------------------
## mytab：
mytab <- function(x, y) {
  
  # 去重
  unique_x <- sort(unique(x))
  unique_y <- sort(unique(y))
  # 计算频数矩阵
  xcol<-function(i){
    sapply(unique_y,function(j) sum(y[x==i]==j))
  } #向量化，得到第i行
  freq<-t(as.matrix(sapply(unique_x,xcol)))
  rownames(freq)<-unique_x
  colnames(freq)<-unique_y
  # 返回频数向量
  return(freq)
}

## -----------------------------------------------------------------------------
##举例使用上述函数
x <- sample(2:7, 50, replace = TRUE)
y <- sample(3:8, 50, replace = TRUE)
mytab(x,y) #行名为x，列为y

## -----------------------------------------------------------------------------
table(x,y)#验证和table结果相同

## -----------------------------------------------------------------------------
#use it to speed up your chi-square test?

my_chisq.test2 <- function(x, y) {
  # 默认确保输入是数值型向量且无缺失值
  # 根据卡方检验公式，尽量向量化
  tab <- mytab(x, y)
  obs <- as.vector(tab)
  ##利用outer外积
  exp <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  # 卡方检验统计量
  result <- sum((obs - exp)^2 / exp)
  return(result)
}

## -----------------------------------------------------------------------------
##举例使用上述函数
set.seed(123)
x <- sample(1:7, 50, replace = TRUE)
y <- sample(1:7, 50, replace = TRUE)

# 计算卡方检验统计量
my_chisq.test2(x, y)#结果是相同的，速度增大了


## ----eval=FALSE---------------------------------------------------------------
#  rm(list=ls())

## ----eval=FALSE---------------------------------------------------------------
#  #the R function
#  gibb<-function(m,cut,n=10,a=1,b=1){
#    x<-numeric(m)
#    y<-numeric(m)
#    x[1]<-rbinom(1,n,1/2)
#    y[1]<-rbeta(1,x[1]+a,n-x[1]+b)
#    for (i in 2:m) {
#      x[i]<-rbinom(1,n,y[i-1])
#      y[i]<-rbeta(1,x[i-1]+a,n-x[i-1]+b)
#    }
#    x<-x[(cut+1):m]
#    y<-y[(cut+1):m]
#    f<-matrix(c(x,y),nrow = (m-cut),ncol = 2)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  #相应的C++版本,具体可见压缩包中的gibbs.cpp
#  library(Rcpp)
#  dir_cpp <- 'D:/Download/RStudio/R-4.2.2/library/Rcpp/myfun/'
#  #路径
#  sourceCpp(paste0(dir_cpp,"gibbs.cpp"))

## ----eval=FALSE---------------------------------------------------------------
#  # 生成1000个数据的MC
#  set.seed(123)
#  g_c<-gibbs(2000,1000) #Rcpp
#  g_r<-gibb(2000,1000)  #R function
#  

## ----eval=FALSE---------------------------------------------------------------
#  #Compare the corresponding generated random numbers
#  #with those by the R function
#  #对x分量标准化以便画qq图
#  g_c[,1]<-(g_c[,1]-min(g_c[,1]))/(max(g_c[,1])-min(g_c[,1]))
#  g_r[,1]<-(g_r[,1]-min(g_r[,1]))/(max(g_r[,1])-min(g_r[,1]))
#  qqplot(g_c[,1],g_r[,1],main="x--qqplot")
#  qqplot(g_c[,2],g_r[,2],main="y--qqplot")

## ----eval=FALSE---------------------------------------------------------------
#  #m=5000
#  library(microbenchmark)
#  ts <- microbenchmark(gC=gibbs(5000,1000),gR=gibb(5000,1000))
#  summary(ts)[,c(1,3,5,6)]

