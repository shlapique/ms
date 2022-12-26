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
        print(paste(left_lid, "<", t_lid, "<", right_lid))
        print(th_lid)
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

# TRUE stands for normal
# FALSE stands for uniform
flag <- FALSE
##########

# define all variables
n = 40
disp = 2.5
sigm = sqrt(disp) 
m = 2
if(flag == TRUE)
{
    theta = c((-1)^11 * 11, 5, 6, 0.06)
} else
{

    theta = c((-1)^11 * 11, 5, 6)
}

# set error vector with normal distribution
if(flag == TRUE)
{
    e <- rnorm(n, 0, sigm)
} else
{
    e <- runif(n, -3*sigm, 3*sigm)
}


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
if(flag == TRUE)
{
    data <- c(rep(1, n), x, x_2, x_3)
    x_main <- matrix(data, nrow=n, ncol=4)
} else
{
    data <- c(rep(1, n), x, x_2)
    x_main <- matrix(data, nrow=n, ncol=3)
}
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
        print(paste(l, ", theta=", theta[j], ", ", r))
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
q_995 <- qt(0.995, n - m_lid - 1)

phi_th <- yac_str(paste0("Simplify(", phi_th, ")"))
phi_th
real_phi_th <-yac_str(paste0("Simplify(", real_phi_th, ")"))
real_phi_th

s = m_lid + 1
phi_x <- c()
for(i in 1:s)
{
    phi_x <- c(phi_x, paste0("x^", i-1))
    print(phi_x)
}
phi_x <- yac_str(paste0("Transpose(", ysym(matrix(phi_x)), ")"))
phi_x_T <- yac_str(paste0("Transpose(", phi_x, ")"))
phi_x
phi_x_T
XX <- ysym(alpha_s)
sqrt_alpha_x <- yac_str(paste0("Sqrt(Simplify(", phi_x * XX * phi_x_T, "))"))
sqrt_alpha_x

# q for alpha = 0.95

right_975 <- yac_str(paste0("Simplify(", phi_th, " + ", yac_str(paste0(norm_E_lid*q_975, " * ", sqrt_alpha_x)), " / ", yac_str(paste0("Sqrt(", n - s, ")")), ")"))
left_975 <- yac_str(paste0("Simplify(", phi_th, " - ", yac_str(paste0(norm_E_lid*q_975, " * ", sqrt_alpha_x)), " / ", yac_str(paste0("Sqrt(", n - s, ")")), ")"))
print(paste0(left_975, " < ", "real_phi_th", " < ", right_975))
#
# q for alpha = 0.95

right_995 <- yac_str(paste0("Simplify(", phi_th, " + ", yac_str(paste0(norm_E_lid*q_995, " * ", sqrt_alpha_x)), " / ", yac_str(paste0("Sqrt(", n - s, ")")), ")"))
left_995 <- yac_str(paste0("Simplify(", phi_th, " - ", yac_str(paste0(norm_E_lid*q_995, " * ", sqrt_alpha_x)), " / ", yac_str(paste0("Sqrt(", n - s, ")")), ")"))
print(paste0(left_995, " < ", "real_phi_th", " < ", right_995))


#=====
#4
#=====

data_true <- data.frame("x" = x_main[1:nrow(x_main), 2], "useful" = x_main%*%theta)
data_obs <- data.frame("x" = x_main[1:nrow(x_main), 2], "obs" = y_res)

X11() 
ggplot() + geom_line(data_true, mapping = aes(x=x, y=useful, colour="true useful signal")) + geom_line(data_obs, mapping = aes(x=x, y=obs, colour="observation")) + scale_colour_manual(values=c("observation"="red", "true useful signal"="blue"))
check_device()

x_main
th_lid
if(ncol(x_main) != nrow(th_lid))
{
    print("!!!!!!!!!!!!!!!!!!! M_LID > or < m")
    stop()
} else
{
    es_es <- x_main%*%th_lid
}

data_es <- data.frame("x" = x_main[1:nrow(x_main), 2], "es" = es_es)
X11() 
ggplot() + geom_line(data_true, mapping = aes(x=x, y=useful, colour="true useful signal")) + geom_line(data_es, mapping = aes(x=x, y=es, colour="estimate")) + scale_colour_manual(values=c("estimate"="orange", "true useful signal"="blue"))
check_device()

# alpha = 0.95
print("ALPHAAAAA = 0.95")
u <- theta[1]
for(i in 2:s)
{
    u <- paste0(u, " + ", theta[i], " * ", "x^", i-1)
}
true_sign <- yac_str(paste0("Simplify(", u, ")"))
show_g <- yac_str(paste0("Plot2D({", as_r(left_975), ", ", as_r(right_975), ", ", true_sign, "})"))
#
#
#
# alpha = 0.99
print("ALPHAAAAA = 0.99")
u <- theta[1]
for(i in 2:s)
{
    u <- paste0(u, " + ", theta[i], " * ", "x^", i-1)
}
true_sign <- yac_str(paste0("Simplify(", u, ")"))
show <- yac_str(paste0("Plot2D({", as_r(left_995), ", ", as_r(right_995), ", ", true_sign, "})"))

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
if(flag == TRUE)
{
    ggplot() + geom_histogram(data_hist, mapping = aes(E_lid, after_stat(density)), bins=l) + geom_line(p, mapping = aes(x, tmp))
} else
{
    ggplot() + geom_histogram(data_hist, mapping = aes(E_lid, after_stat(density)), bins=l) 
}
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
print(summ * n)

summ <- summ + pnorm(t_interval[1] / est_sigm)
prob <- prob + pnorm(t_interval[1] / est_sigm)
summ <- summ + 1 - pnorm(t_interval[length(t_interval)] / est_sigm)
prob <- prob + 1 - pnorm(t_interval[length(t_interval)] / est_sigm)
summ <- summ * n
print(paste("SUMMM = ", summ))

print(paste("0 < ", summ, " < ", qchisq(.95, l-1)))
