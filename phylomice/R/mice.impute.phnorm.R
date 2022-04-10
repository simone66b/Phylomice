#' Imputation by GLM and phylogenetic information
#'
#' Imputes univariate continuous missing data using the generalized least square approach.
#' @aliases mice.impute.phnorm phnorm
#' @param y Incomplete data vector of length \code{n}
#' @param ry Vector of missing data pattern (\code{FALSE}=missing,
#' \code{TRUE}=observed)
#' @param x Matrix (\code{n} x \code{p}) of complete covariates.
#' @param psi Covariance matrix containing phylogenetic information.
#' @param tree Matching Tree
#' @param ... Other named arguments.
#' @return A vector of length \code{nmis} with imputations.
#' @note \code{mice.impute.phnorm} is based on the
#' \code{mice.impute.norm} method from \code{mice} package and the GLM
#' approach by Garland and Ives (2000).
#' @author Patrik Drhlik, Simone P. Blomberg, 2016, 2021
#' @references Van Buuren, S., Groothuis-Oudshoorn, K. (2011). \code{mice}:
#' Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
#' Statistical Software}, \bold{45}(3), 1-67.
#' \url{http://www.jstatsoft.org/v45/i03/}
#'
#' Garland Jr, Theodore, and Anthony R. Ives. Using the past to predict
#' the present: confidence intervals for regression equations in phylogenetic
#' comparative methods. \emph{The American Naturalist}, 155.3 (2000): 346-364.
#' \url{http://www.jstor.org/stable/10.1086/303327}
#' @export
mice.impute.phnorm <- function(y, ry, x, tree, ...) {
    require("ape")
    require("nlme")
    ## x2 <- x
    ## species <- rownames(x2)
    x <- cbind(1, as.matrix(x))
    ## rownames(x) <- species
    
    parm <- .phnorm.draw(y, ry, x, tree)
    species <- tree$tip.label
	ymiss <- which(!ry)
        tree.trimmed <- drop.tip(tree, species[ymiss])
        lambda <- parm$lambda

        ## LAMBDA is a matrix of lambda on the off-diagonals and ones on the diagonal.
        LAMBDA <- 1 + (lambda-1) * (1-diag(length(tree$tip.label)))

        C <- LAMBDA * vcv.phylo(tree, corr=TRUE)
        
	nmiss <- sum(!ry)

                                        # vectorized version
    
	Carr <- lapply(ymiss, function(i) as.matrix(C[-i, -i]))
	cihmat <- lapply(ymiss, function(i) C[i, -i])
	# heights of the trees
	chhvec <- diag(C)
	# matrix of products of cihmat and Carr
	cihCC <- sapply(1:nmiss, function(i) {
            t(cihmat[[i]]) %*% solve(Carr[[i]]) %*% cihmat[[i]]})
        
	# vector of means
	if (!is.matrix(cihCC)) cihCC <- matrix(cihCC)
	mu <- colSums(sapply(1:nmiss, function(i) {
                t(cihmat[[i]]) %*% solve(Carr[[i]]) %*% scale(x[-ymiss[i],], scale=FALSE)}))
    

	# vector of ch
	ch <- sapply(1:nmiss, function(i) {
		chhvec[ymiss[i]] - as.vector(cihmat[[i]] %*% solve(Carr[[i]]) %*% cihmat[[i]])})

	return(x[!ry, ] %*% parm$beta + rnorm(nmiss, mean = mu, sd = parm$sigma * sqrt(ch)))
}

.phnorm.draw <- function (y, ry, x, tree, ...) 
{
    require(ape)
    require(nlme)
    xobs <- x[ry, -1]
    yobs <- y[ry]
    ## names(yobs) <- species[ry]
    tree.obs <- keep.tip(tree, which(ry))
    species.obs <- tree.obs$tip.label
    ## if (is.null(species.obs)) species.obs <- names(xobs)
    datobs <- data.frame(species.obs = species.obs,
                               y = yobs, xobs)
    ## print(datobs)
    ## stop("Stopped")
    form <- formula(paste("y ~", paste(colnames(datobs)[-(1:2)], 
        collapse = " + ")))
    startLambda <- runif(1)
    fit <- try(gls(form,
                   correlation =
                             corPagel(startLambda,
                                      phy = tree.obs, 
                                      form = ~species.obs),
                   data = datobs), silent=TRUE)
    
    if (inherits(fit, "try-error")) {
        fit2 <- lm(form, data = datobs)
        coef <- coef(fit2)
        df <- max(sum(ry) - ncol(x), 1)
        sigma.star <- summary(fit2)$sigma/rchisq(1, df)
        sds <- sqrt(diag(vcov(fit2)))
        beta.star <- rnorm(length(coef), mean = coef, sd = sds)
        names(beta.star) <- names(coef)
        alpha <- 1
        beta <- 1
    } else {
        coef <- coef(fit)
        df <- max(sum(ry) - ncol(x), 1)
        sigma.star <- fit$sigma/rchisq(1, df)
        sds <- sqrt(diag(vcov(fit)))
        beta.star <- rnorm(length(coef), mean = coef, sd = sds)
        names(beta.star) <- names(coef)
        lambda <- fit$modelStruct$corStruct[1]
        lambda.var <- try(fit$apVar[1, 1], silent = TRUE)
        if (inherits(lambda.var, "try-error") || lambda < 0 || lambda > 1) {
            alpha <- 1
            beta <- 1
        } else {
            alpha <- ((1 - lambda)/lambda.var - 1/lambda) * lambda^2
            beta <- alpha * (1/lambda - 1)
        }
    }
    lambda.star <- try(rbeta(1, alpha, beta))
    if (inherits(lambda.star, "try-error")) stop("bad lambda")
    parm <- list(coef, beta.star, sigma.star, lambda.star)
    names(parm) <- c("coef", "beta", "sigma", "lambda")
    ## print(parm)
    return(parm)
}

