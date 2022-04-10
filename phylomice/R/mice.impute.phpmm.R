#' Imputation by GLM and phylogenetic information
#'
#' Imputes univariate continuous missing data using the generalized least square approach.
#' @aliases mice.impute.phpmm phpmm
#' @param y Incomplete data vector of length \code{n}
#' @param ry Vector of missing data pattern (\code{FALSE}=missing,
#' \code{TRUE}=observed)
#' @param x Matrix (\code{n} x \code{p}) of complete covariates.
#' @param ... Other named arguments.
#' @return A vector of length \code{nmis} with imputations.
#' @author Patrik Drhlik, Simone P. Blomberg, 2016, 2021
#' @references Van Buuren, S., Groothuis-Oudshoorn, K. (2011). \code{mice}:
#' Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
#' Statistical Software}, \bold{45}(3), 1-67.
#' \url{http://www.jstatsoft.org/v45/i03/}
#'
#' @export
mice.impute.phpmm <-
    function(y, ry, x, tree, donors = 5L, matchtype=1L, ...)  {
    require(ape)
    require(nlme)

    wy <- !ry
    x <- cbind(1, as.matrix(x))
    rownames(x) <- tree$tip.label
    species <- tree$tip.label
    ynum <- y
    names(y) <- species
 parm <- .phnorm.draw(y, ry, x, tree)
    	ymiss <- which(!ry)
        tree.trimmed <- drop.tip(tree, species[ymiss])
        lambda <- parm$lambda

        ## LAMBDA is a matrix of lambda on the off-diagonals and ones on the diagona.
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

    
    if (is.factor(y)) ynum <- as.integer(y)

    if (matchtype == 0L) {
        yhatobs <- x[ry, , drop = FALSE] %*% parm$coef
        yhatmis <- x[wy, , drop = FALSE] %*% parm$coef
    }
    if (matchtype == 1L) {
        yhatobs <- x[ry, , drop = FALSE] %*% parm$coef
        yhatmis <- x[wy, , drop = FALSE] %*% parm$beta + mu
    }
    if (matchtype == 2L) {
        yhatobs <- x[ry, , drop = FALSE] %*% parm$beta
        yhatmis <- x[wy, , drop = FALSE] %*% parm$beta + mu
    }


    
    idx <- .phpmm.match (yhatobs, yhatmis, tree=tree, donors=donors)
    return(y[ry][idx])
}

.phpmm.match <- function(yhatobs, ry=yhatmis, donors = donors, tree=tree, ...) {
    require(ape)
    mvec <- vector(mode="numeric", length = dim(ry)[1])
    for (i in 1:dim(ry)[1]) {
        z <- ry[i, 1]
    d <- abs(yhatobs-z)
    f <- d > 0
    a1 <- ifelse(any(f), min(d[f]), 1)
    ## print(a1)
    d <- d + runif(length(d), 0, a1/10^10)
    if (donors == 1) 
        return(y[which.min(d)]) ## fail because y doesn't exist!!
    donors <- min(donors, length(d))
    donors <- max(donors, 1)
    ds <- sort(d, index.return=TRUE)
    names(ds$x) <- rownames(d)[ds$ix]
    ds <- ds$x[1:donors]
    
    nmiss <- sum(!ry)
    nall <- length(ry)
    
    if (is.factor(ry)) ry <- as.integer(ry)
    
    probs <- 1-apply(cophenetic(tree), 1, function (x) x/sum(x)) ## columns sum to 1

    ## Values only from observed
    probs[!ry, !ry] <- 0
    diag(probs) <- 0 ## don't choose yourself
    probs <- apply(probs, 1, function(x) x/sum(x)) ## rescale
    probsds <- probs[which(rownames(probs) %in% names(z)),
                     which(rownames(probs) %in% names(ds))]
  ## Resulting vector
  ## Choose one of the values with a certain probability
    donorspecies <- names(ds[1:donors])

    ## print(donorspecies)
    if (any(is.na(donorspecies)))  {
	    m <- sample(rownames(yhatobs), 1) ## m <- rownames(yhatobs)[which.min(d)] 
    } else {
        m <- sample(sort(donorspecies),1, prob=probsds[order(names(probsds))]) 
    }  
        ## print(m)
        mvec[i] <- which(rownames(yhatobs) == m)
        }
  return(mvec)
}


