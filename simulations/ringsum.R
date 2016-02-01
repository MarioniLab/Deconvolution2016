# This justifies why you need to sum in a ring.

N <- 18
out <- list()
for (x in 2:N) {
    y <- integer(N)
    y[x-0:1] <- 1
    out[[x]] <- y
}

y <- integer(N)
y[1] <- y[N] <- 1
out[[1]] <- y
X <- do.call(rbind, out)
X <- rbind(X, diag(N))
w <- rep(c(1, 1e-6), each=N)^0.5

set.seed(10)
err.sum1 <- err.sum2 <- list()
for (it in 1:100) { 
    true.facs <- c(1:10/10, 9:2/10)
    obs <- X %*% true.facs
    obs <- obs * runif(length(obs), 0.8, 1.2)
    err.sum1[[it]] <- as.vector(abs(solve(qr(X*w), obs*w) - true.facs)/true.facs)
    
    true.facs <- c(0.1, 1, 0.2, 0.9, 0.3, 0.8, 0.4, 0.7, 0.5, 0.6, 0.6, 0.5, 0.7, 0.4, 0.8, 0.3, 0.9, 0.2)
    obs <- X %*% true.facs
    obs <- obs * runif(length(obs), 0.8, 1.2)
    err.sum2[[it]] <- as.vector(abs(solve(qr(X*w), obs*w) - true.facs)/true.facs)
}

summary(do.call(rbind, err.sum1)[,1:10])
summary(do.call(rbind, err.sum2)[,c(1,3,5,7,9,2,4,6,8,10)])
