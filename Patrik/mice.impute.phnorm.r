#'@param psi Inverse matrix containing phylogenetic information.
#' This matrix should be passed directly to the mice function.
mice.impute.phnorm <- function(y, ry, x, psi, psiinv, ...) {
	x <- cbind(1, as.matrix(x))
	parm <- .phnorm.draw(y, ry, x, psiinv, ...)

	ymiss <- y[!ry]
	Cinv <- psiinv[!ry, !ry]
	C <- psi[!ry, !ry]
	nmiss <- sum(!ry)
	mu <- vector(mode = "numeric", length = nmiss)
	ch <- vector(mode = "numeric", length = nmiss)

	for (i in 1:nmiss) {
		chh <- C[i,i]
		CC <- C[-i, -i]
		CCinv <- solve(CC)
		yy <- ymiss[-i]
		cih <- C[i, -i]
		mu[i] <- cih %*% CCinv %*% (yy - mean(yy))
		ch[i] <- chh - t(cih) %*% CCinv %*% cih
	}

	return(x[!ry, ] %*% parm$beta + rnorm(nmiss, mean = mu, sd = parm$sigma * ch))
}

.phnorm.draw <- function(y, ry, x, psiinv, ridge = 1e-05, ...) {
	xobs <- x[ry, ]
	yobs <- y[ry]
	psiobs <- psiinv[ry, ry]

	xtx <- t(xobs) %*% psiobs %*% xobs
	pen <- ridge * diag(xtx)
	if (length(pen) == 1) {
		pen <- matrix(pen)
	}

	v <- solve(xtx + diag(pen))
	coef <- v %*% t(xobs) %*% psiobs %*% yobs

	residuals <- yobs - xobs %*% coef
	df <- max(sum(ry) - ncol(x), 1)  # SvB 31/10/2012
	sigma.star <- sqrt(sum((residuals)^2)/rchisq(1, df))  # SvB 01/02/2011
	beta.star <- coef + (t(chol(sym(v))) %*% rnorm(ncol(x))) * sigma.star
	parm <- list(coef, beta.star, sigma.star)  # SvB 10/2/2010
	names(parm) <- c("coef", "beta", "sigma")  # SvB 10/2/2010
	return(parm)
}

# from internal.R
sym <- function(x) {(x + t(x)) / 2}
