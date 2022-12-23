#!/usr/bin/env Rscript

#=====
#0
#=====

# define all variables
n = 40
disp = 2.5
sigm = sqrt(disp) 
m = 3
theta = c((-1)^11 * 11, 5, 6, 0.06)

# set error vector with normal distribution
e <- rnorm(n = n, mean = 0, sd = sigm)

# set up x's
#x <- matrix(rep(1, n), nrow=n, ncol=1)
x <- rep(1, n)
for(i in 1:n)
{
    x[i] <- -4 + i * 8/40
}
x_2 <- x*x
x_3 <- x_2*x 
x_4 <- x_3*x

# put all vector and matricies in one expression
data <- c(rep(1, 40), x, x_2, x_3)
x_main <- matrix(data, nrow=40, ncol=4)
y_res <- (x_main%*%theta) + e
#y_res

#=====
#1
#=====

#m = 2   #по итогу theta 5 равно 0, поэтому в гипотезе мы делаем шаг назад и говорим, что число неизвестных(theta) = 4, а число степеней x равно 4-1=3
data_2 <- c(rep(1, 40), x, x_2, x_3, x_4)
x_main_2 <- matrix(data_2, nrow=40, ncol=5)
x_main_2 <- x_main_2[, -5]
#x_main_2

# Tetha' = (X^T * X)^-1 * X^T * Y
# theta with lid 
th <- ((solve(t(x_main_2)%*%x_main_2))%*%t(x_main_2))%*%y_res
th

alpha <- solve(t(x_main_2)%*%x_main_2)
#alpha
al_1_1 <- alpha[4, 4]
#al_1_1

E <- y_res - x_main_2%*%th
#E
norm_E <- sqrt(sum(E**2))
#norm_E

t = th[nrow(th), 1] * sqrt(40 - m - 1) / (sqrt(al_1_1) * norm_E)
#t

left <- qt(.025/2, n - m - 1)
right <- -qt(.025/2, n - m - 1)
#print(paste(left, "<", t, "<", right))

#=====
#2
#=====

for(i in c(0.95, 0.99))
{
    for(j in 1:nrow(th))
    {
        #print(th[j][1] + (qt((1 - i) / 2, n - m - 1) * norm_E) * sqrt(alpha[j, j]))
        l <- th[j][1] + (qt((1 - i) / 2, n - m - 1) * norm_E * sqrt(alpha[j, j])) / sqrt(n - m - 1)
        r <- th[j][1] - (qt((1 - i) / 2, n - m - 1) * norm_E * sqrt(alpha[j, j])) / sqrt(n - m - 1)
        print(paste(l, ", theta=", theta[i], ", ", r))
    }
    print("=============================================")
}

#=====
#3
#=====

sympy("phi_th = th[1][1]")
for(i in nrow(th))
{
    phi_th <- phi_th + i * x**(i + 1)
    print(phi_th)
}