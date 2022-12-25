library(Ryacas)
library(ggplot2)

check_device <- function()
{
    while (!is.null(dev.list())) Sys.sleep(1)
}

# returns m_lid 
fisher_criterion <- function(n, x, y_res)
{
    x_data <- c(rep(1, n), x)
    for(m_lid in 2:10)
    {
        x_data <- c(x_data, x^m_lid)
        X <- matrix(x_data, nrow=n, ncol=m_lid + 1)
        th_lid <- ((solve(t(X)%*%X))%*%t(X))%*%y_res
        alpha_s <- solve(t(X)%*%X)
        al_1_1_s <- alpha_s[nrow(alpha_s), ncol(alpha_s)]
        E_lid <- y_res - X%*%th_lid
        norm_E_lid <- sqrt(sum(E_lid**2))
        t_lid = th_lid[nrow(th_lid), 1] * sqrt(n - m_lid - 1) / (sqrt(al_1_1_s) * norm_E_lid)
        left_lid <- qt(.025/2, n - m_lid - 1)
        right_lid <- -qt(.025/2, n - m_lid - 1)
        if((t_lid > left_lid) & (t_lid < right_lid))
        {
            print("INSIDE INTERVAL!")
            print(paste(left_lid, "<", t_lid, "<", right_lid))
            print(paste("=> solved m with lid = m_lid - 1 = ", m_lid - 1))
            X <- X[, -ncol(X)]
            print(X)
            return (m_lid - 1)
            break
        }
    }
}

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
e <- rnorm(n, 0, sigm)


# set up x's
#x <- matrix(rep(1, n), nrow=n, ncol=1)
x <- rep(1, n)
for(i in 1:n)
{
    x[i] <- -4 + i * 8/40
}
x_2 <- x*x
x_3 <- x_2*x 

# put all vector and matricies in one expression
data <- c(rep(1, n), x, x_2, x_3)
x_main <- matrix(data, nrow=n, ncol=4)
y_res <- (x_main%*%theta) + e
#y_res

#=====
#1
#=====
# ---------------------------
m_lid <- fisher_criterion(n, x, y_res)
print(paste("REAL M = ", m_lid))
x_data <- c(rep(1, n), x)
for(em in 2:m_lid)
{
    x_data <- c(x_data, x^em)
    X <- matrix(x_data, nrow=n, ncol=em + 1)

    # Tetha' = (X^T * X)^-1 * X^T * Y
    # theta with lid 
    th_lid <- ((solve(t(X)%*%X))%*%t(X))%*%y_res
    alpha_s <- solve(t(X)%*%X)
    al_1_1_s <- alpha_s[nrow(alpha_s), ncol(alpha_s)]
    E_lid <- y_res - X%*%th_lid
    norm_E_lid <- sqrt(sum(E_lid**2))
    t_lid = th_lid[nrow(th_lid), 1] * sqrt(n - em - 1) / (sqrt(al_1_1_s) * norm_E_lid)
    left_lid <- qt(.025/2, n - em - 1)
    right_lid <- -qt(.025/2, n - em - 1)
    print(paste(left_lid, "<", t_lid, "<", right_lid))
}
#----------------------------

#=====
#2
#=====
for(i in c(0.95, 0.99))
{
    for(j in 1:nrow(th_lid))
    {
        #print(th[j][1] + (qt((1 - i) / 2, n - m - 1) * norm_E_lid) * sqrt(alpha[j, j]))
        l <- th_lid[j][1] + (qt((1 - i) / 2, n - m_lid - 1) * norm_E_lid * sqrt(alpha_s[j, j])) / sqrt(n - m_lid - 1)
        r <- th_lid[j][1] - (qt((1 - i) / 2, n - m_lid - 1) * norm_E_lid * sqrt(alpha_s[j, j])) / sqrt(n - m_lid - 1)
        print(paste(l, ", theta=", theta[i], ", ", r))
    }
    print("=============================================")
}

#=====
#3
#=====
phi_th <- as.character(th_lid[1])
m_list <- as.list(th_lid[2:nrow(th_lid), 1])
print(m_list)
for(i in seq_along(m_list))
{
    # in every equation: i => i, t = > m_list[[i]]
    phi_th <- paste0(phi_th, " + ", as.character(m_list[[i]]), " * ", "x^", i)
    print(phi_th)
}

real_phi_th <- as.character(theta[1])
m_list <- as.list(theta[2:length(theta)])
for(i in seq_along(m_list))
{
    real_phi_th <- paste0(real_phi_th, " + ", as.character(m_list[[i]]), " * ", "x^", i)
    print(real_phi_th)
}
q_975 <- qt(0.975, n - m_lid - 1)
q_995 <- qt(0.975, n - m_lid - 1)

phi_th <- yac_str(paste0("Simplify(", phi_th, ")"))
phi_th
real_phi_th <-yac_str(paste0("Simplify(", real_phi_th, ")"))
real_phi_th
show <- yac_str(paste0("Plot2D(", phi_th, ", x=-4:4, y=-10:100)"))
show


#TODO

#=====
#4
#=====


#=====
#5
#=====
E_max = max(E_lid) + 0.01
E_min = min(E_lid) - 0.01
l = 7
delta <- (E_max - E_min)/(l-1)
print(paste("E_min = ", E_min, "E_max = ", E_max, "delta = ", delta))

t_interval = seq(E_min, E_max + delta, delta)
#t_interval

fs = c()
for(i in 1:(length(t_interval) - 2))
{
    c <- 0
    for(ep in E_lid)
    {
        if((t_interval[i] <= ep) & (ep < t_interval[i + 1]))
        {
            c <- c + 1
        }
    }
    fs <- c(fs, c / (n * (t_interval[i + 1] - t_interval[i])))
}
print(fs)

data_hist <- data.frame(E_lid)
#data_hist

x <- seq(-5, 5, 0.1)
tmp <- dnorm(x, mean(E_lid), sd(E_lid))
p <- data.frame(tmp)

# histogram + line of normal distribution 
X11() 
ggplot() + geom_histogram(data_hist, mapping = aes(E_lid, after_stat(density)), bins=l) + geom_line(p, mapping = aes(x, tmp))
check_device()

#=====
#6
#=====

# estimate disp
est_sigm <- sqrt(1 / n * norm_E_lid^2)
est_sigm^2

#=====
#7
#=====
p_hist <- c()
for(i in 1:(l-1))
{
    p_hist <- c(p_hist, fs[i] * delta)
}

summ <- 0 
i <- 1
prob <- 0
for(p_lid in p_hist)
{
    p <- pnorm(t_interval[i + 1] / est_sigm) - pnorm(t_interval[i] / est_sigm)
    print(paste("P = ", p))
    prob <- prob + p
    summ <- summ + ((p_lid - p) ^ 2) / p
    i <- i + 1
}
print(summ * 40)

summ <- summ + pnorm(t_interval[1] / est_sigm)
prob <- prob + pnorm(t_interval[1] / est_sigm)
summ <- summ + 1 - pnorm(t_interval[length(t_interval)] / est_sigm)
prob <- prob + 1 - pnorm(t_interval[length(t_interval)] / est_sigm)
summ <- summ * n
print(paste("SUMMM = ", summ))

print(paste("0 < ", summ, " < ", qchisq(.95, l-1)))
