##' The R version of \code{\link[poetREV]{poet}}, which is written in Cpp.
##'
##' This function is written for the purpose of checking the correctness of
##' \code{\link[poetREV]{poet}}, see the examples in
##' \code{\link[poetREV]{poet}}.
##' @param Y p by n data matrix.
##' @param K the number of latent factors.
##' @param C the shrinkage intensity.
##' @param thres shrinkage method. Default is "soft".
##' @param matrix indicator of shrinkage matrix, "cor" (default) or "vad".
##' @return a list
##' @export
poet2 <- function(Y, K = -Inf, C = -Inf, thres = "soft", matrix = "cor") {
    p <- nrow(Y)
    n <- ncol(Y)
    ## Y is de-meaned by row
    Y <- Y - rowMeans(Y)
    ## Y <- t(scale(t(Y), scale = FALSE))

    if (K == -Inf) {
        K <- floor(mean(poet_Khat(Y))) + 1
        K <- min(K, p)
    }
    if (K > 0) {
        if (K > p) stop("Invalid K argument in poet2(). K must <= p.")
        
        eig <- eigen(t(Y) %*% Y)
        V <- eig$vectors

        ## using the first K pc as a proxy for the space that
        ## spanned by the columns of the factor loading matrix
        F <- V[, 1:K, drop = FALSE] * sqrt(n) ## F is n by K
        LamPCA <- Y %*% F / n
        uhat <- Y - LamPCA %*% t(F) ## uhat is p by n
        Lowrank <- LamPCA %*% t(LamPCA)
        rate <- 1/sqrt(p) + sqrt(log(p)/n)
    } else if (K == 0) {
        ## SigmaY itself is sparse
        uhat <- Y
        rate <- sqrt(log(p)/n)
        Lowrank <- matrix(0, p, p)
    } else {
        stop("Invalid K argument in poet2(). K must >= 0.")
    }

    SuPCA <- uhat %*% t(uhat) / n
    SuDiag <- diag(SuPCA)

    R <- switch(
        matrix,
        "cor" = diag(1/sqrt(SuDiag)) %*% SuPCA %*% diag(1/sqrt(SuDiag)),
        "vad" = SuPCA)
    if (C == -Inf) {
        C <- poet_Cmin(Y, K, thres, matrix) + .1
    }

    ## adaptive threshold
    lambda <- matrix(NA, p, p)
    for (i in 1:p) {
        for (j in 1:i) {
            ## original POET package is wrong for thres = "cor"
            ## See, Fan, J., Liao, Y., & Liu, H. (2016). An overview of the estimation of large covariance and precision matrices. The Econometrics Journal, 19(1), C1-C32. ==> Page 5 Eq.(4)
            ## lambda[i, j] <- sqrt(var(uhat[i, ]*uhat[j, ])*(n-1)/n) * rate * C
            lambda[i, j] <- switch(
                matrix,
                "cor" = rate * C,
                "vad" = sqrt(var(uhat[i, ]*uhat[j, ])*(n-1)/n) * rate * C)
            lambda[j, i] <- lambda[i, j]
        }
    }

    Rthresh <- matrix(NA, p, p)
    switch(thres,
           "soft" = {
               for (i in 1:p) {
                   for (j in 1:i) {
                       if (abs(R[i, j]) <= lambda[i, j] && j < i) {
                           Rthresh[i, j] <- 0
                       } else if (j == i) {
                           Rthresh[i, j] <- R[i, j]
                       } else {
                           Rthresh[i, j] <- sign(R[i, j]) * (abs(R[i, j])-lambda[i, j])
                       }
                       Rthresh[j, i] <- Rthresh[i, j]
                   }
               }
           },
           "hard" = {
               for (i in 1:p) {
                   for (j in 1:i) {
                       if (abs(R[i, j]) <= lambda[i, j] && j < i) {
                           Rthresh[i, j] <- 0
                       } else {
                           Rthresh[i, j] <- R[i, j]
                       }
                       Rthresh[j, i] <- Rthresh[i, j]
                   }
               }
           },
           "scad" = {
               for (i in 1:p) {
                   for (j in 1:i) {
                       if (j == i) {
                           Rthresh[i, j] <- R[i, j]
                       } else if (abs(R[i, j]) <= lambda[i, j]) {
                           Rthresh[i, j] <- 0
                       } else if (abs(R[i, j]) > lambda[i, j] & abs(R[i, j]) < 2*lambda[i, j]) {
                           Rthresh[i, j] <- sign(R[i, j]) * (abs(R[i, j])-lambda[i, j])
                       } else if (abs(R[i, j]) > 2*lambda[i, j] & abs(R[i, j]) < 3.7*lambda[i, j]) {
                           Rthresh[i, j] <- ((3.7-1)*R[i, j]-sign(R[i, j])*3.7*lambda[i, j])/(3.7-2)
                       } else {
                           Rthresh[i, j] <- R[i, j]
                       }
                       Rthresh[j, i] <- Rthresh[i, j]
                   }
               }
           })

    SigmaU <- switch(
        matrix,
        "cor" = diag(sqrt(SuDiag)) %*% Rthresh %*% diag(sqrt(SuDiag)),
        "vad" = Rthresh
    )
    SigmaY <- Lowrank + SigmaU

    ## if (K == 0) {
    ##     return(list(SigmaY = SigmaY, SigmaU = SigmaU))
    ## } else {
    ##     return(list(SigmaU = SigmaU, SigmaY = SigmaY, factors = t(F), loadings = LamPCA))
    ## }
    return(list(SigmaY = SigmaY, SigmaU = SigmaU)) ## coincide with the Cpp version poet
}


##' The R version of \code{\link[poetREV]{poet_Khat}}, which is written in Cpp.
##'
##' This function is written for the purpose of checking the correctness of
##' \code{\link[poetREV]{poet_Khat}}, see the examples in
##' \code{\link[poetREV]{poet_Khat}}.
##' @param Y p by n data matrix.
##' @return a vector with two elements.
##' @export
poet_Khat2 <- function(Y) {
    p <- nrow(Y)
    n <- ncol(Y)
    Y <- Y - rowMeans(Y)

    ##################################################
    ## Hallin and Liska method (JASA, 2007) Page 9

    ## ## for given rmax and c, calculate IC
    ## c <- seq(.05, 5, length.out = 100)
    ## re <- 20 ## re-sample
    ## rmax <- 10
    ## IC <- array(NA, dim = c(2, re, rmax, length(c)))

    ## gT1HL <- rep(NA, re)
    ## gT2HL <- rep(NA, re)
    ## pi <- rep(NA, re)
    ## ni <- rep(NA, re)

    ## for (i in 1:re) {
    ##     ## generate the subsets, "re" of them
    ##     if (i != re) {
    ##         pi[i] <- min(i*floor(p/re)+min(p, 5), p)
    ##         ni[i] <- min(i*floor(n/re)+min(n, 5), n)
    ##     } else {
    ##         pi[i] <- p
    ##         ni[i] <- n
    ##     }
    ##     Yi <- Y[1:pi[i], 1:ni[i]] ## Yi is pi[i] by ni[i]
    ##     frob <- rep(NA, rmax)

    ##     for (k in 1:min(pi[i], ni[i], rmax)) {
    ##         eig <- eigen(t(Yi) %*% Yi)
    ##         V <- eig$vectors
    ##         ## KEY POINT!!! Original code in POET is wrong
    ##         F <- V[, 1:k, drop = FALSE] * sqrt(ni[i]) ## F is ni[i] by k
    ##         LamPCA <- Yi %*% F / ni[i] ## LamPCA is pi[i] by k
    ##         uhat <- Yi - LamPCA %*% t(F) ## uhat is pi[i] by ni[i]

    ##         frob[k] <- sum(diag(uhat %*% t(uhat))) / (pi[i]*ni[i])
    ##         gT1HL[i] <- log((pi[i]*ni[i])/(pi[i]+ni[i])) * (pi[i]+ni[i]) / (pi[i]*ni[i])
    ##         gT2HL[i] <- log(min(pi[i], ni[i])) * (pi[i]+ni[i]) / (pi[i]*ni[i])

    ##         for (l in 1:length(c)) {
    ##             ## only fills in the ICs up to k, which may be < rmax
    ##             IC[1, i, k, l] <- log(frob[k]) + k*c[l]*gT1HL[i]
    ##             IC[2, i, k, l] <- log(frob[k]) + k*c[l]*gT2HL[i]
    ##         }
    ##     }
    ## }
    ## ## WRONG
    ## ## K1HL <- colMeans(apply(IC[1, , , ], 1:2, mean))

    ## ## fixed 1,2 and re, given c, argmin k
    ## rhat <- array(NA, dim = c(2, re, length(c)))
    ## for (i in 1:re) {
    ##     for (l in 1:length(c)) {
    ##         m <- min(pi[i], ni[i], rmax)
    ##         rhat[1, i, l] <- which.min(IC[1, i, 1:m, l])
    ##         rhat[2, i, l] <- which.min(IC[2, i, 1:m, l])
    ##     }
    ## }

    ## Sc1 <- rep(NA, length(c))
    ## Sc2 <- rep(NA, length(c))
    ## for (l in 1:length(c)) {
    ##     Sc1[l] <- sd(rhat[1, , l])
    ##     Sc2[l] <- sd(rhat[2, , l])
    ## }

    ## ## constant we choose in the penalty function
    ## c1 <- c[which(Sc1 == 0)[1]]
    ## ## all Ks in that row are equal
    ## K1HL <- rhat[1, 1, which(Sc1 == 0)[1]]

    ## c2 <- c[which(Sc2 == 0)[1]]
    ## K2HL <- rhat[2, 1, which(Sc2 == 0)[1]]

    ##################################################
    ## Bai and Ng method (Econometrica, 2002) Page 11

    ## penalty corresponds to c = 1
    ## c <- 1
    rmax <- min(15, p)
    IC <- matrix(NA, nrow = 2, ncol = rmax)
    frob <- rep(NA, rmax)

    for (k in 1:rmax) {
        eig <- eigen(t(Y) %*% Y)
        V <- eig$vectors
        ## KEY POINT!!! Original code in POET is wrong
        F <- V[, 1:k, drop = FALSE] * sqrt(n) ## F is n by k
        LamPCA <- Y %*% F / n ## LamPCA is p by k
        uhat <- Y - LamPCA %*% t(F) ## uhat is p by n

        frob[k] <- sum(diag(uhat %*% t(uhat))) / (p*n)
        gT1BN <- log((p*n)/(p+n)) * (p+n) / (p*n)
        gT2BN <- log(min(p, n)) * (p+n) / (p*n)
        IC[1, k] <- log(frob[k]) + k*gT1BN
        IC[2, k] <- log(frob[k]) + k*gT2BN
    }
    K1BN <- which.min(IC[1, ])[1]
    K2BN <- which.min(IC[2, ])[1]

    ## return(c(K1HL = K1HL, K2HL = K2HL, K1BN = K1BN, K2BN = K2BN))
    return(c(K1BN = K1BN, K2BN = K2BN))
}


##' Cmin - Minimum threshold constant
##'
##' (Copy from POET package) This function is for determining the minimum
##' constant in the threshold that guarantees the positive definiteness of the
##' POET estimator.
##'
##' Model: Y_t=Bf_t+u_t, where B, f_t and u_t represent factor loading matrix,
##' common factors and idiosyncratic error respectively. Only Y_t is
##' observable. t=1,...,n. Dimension of Y_t is p. The goal is to estimate the
##' covariance matrices of Y_t and u_t.
##'
##' Note: (1) POET is optimization-free, so no initial value, tolerant, or
##' maximum iterations need to be specified as inputs.
##'
##' (2) We can apply the adaptive thresholding (Cai and Liu 2011, JASA) on
##' either the correlation matrix or the covariance matrix, specified by the
##' option 'matrix'.
##'
##' (3) If no factor structure is assumed, i.e., no common factors exist and
##' var(Y_t) itself is sparse, set K=0.
##' 
##' @param Y p by n matrix of raw data, where p is the dimensionality, n is the
##'     sample size. It is recommended that Y is de-meaned, i.e., each row has
##'     zero mean.
##' @param K number of factors. K is pre-determined by the users. Suggestions on
##'     choosing K:
##'
##' (1) A simple way of determining K is to count the number of very spiked
##' (much larger than others) eigenvalues of the p by p sample covariance matrix
##' of Y.
##'
##' (2) A formal data-driven way of determining K is described in Bai and Ng
##' (2002):"Determining the number of factors in approximate factor models",
##' Econometrica, 70, 191-221. This procedure requires a one-dimensional
##' optimization.
##'
##' (3) POET is very robust to over-estimating K. But under-estimating K can
##' result to VERY BAD performance. Therefore we strongly recommend choosing a
##' relatively large K (normally less than 8) to avoid missing any important
##' common factor.
##'
##' (4) K=0 corresponds to threshoding the sample covariance directly.
##' @param thres choice of thresholding. Users can choose from three
##'     thresholding methods:
##'
##' 'soft': soft thresholding
##'
##' 'hard': hard thresholding
##'
##' 'scad': scad thresholding
##'
##' 'alasso': adaptive lasso thresholding
##'
##' Details are found in Rothman et al. (2009):
##' "Generalized thresholding of large covariance matrices." JASA, 104, 177-186
##' @param matrix the option of thresholding either correlation or covairance
##'     matrix. Users can choose from:
##'
##' 'cor': threshold the error correlation matrix then transform back to
##' covariance matrix
##'
##' 'vad': threshold the error covariance matrix directly.
##' @return a scalar, the minimum threshold constant
##' @examples
##' n <- 50; p <- 100; k <- 3
##' set.seed(1)
##' B <- matrix(rnorm(100 * 3), 100, 3)
##' F <- MASS::mvrnorm(n, rep(0, k), diag(1, k))
##' U <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
##' Y <- t(F %*% t(B) + U)
##' poetREV::poet_Cmin(Y, 3, "soft", "cor")
##' ## [1] 0.2605228
##' @export
poet_Cmin <- function(Y, K, thres, matrix) {
    min_eig <- function(Y, K, C, thres, matrix) {
        SigmaU <- poet(Y, K, C, thres, matrix)$SigmaU
        min(eigen(SigmaU)$values)
    }
    f <- function(x) min_eig(Y, K, x, thres, matrix)
    if (f(50)*f(-50) < 0) {
        r <- uniroot(f, c(-50, 50), tol = 1e-3)
        res <- max(0, r$root)
    } else {
        res <- 0
    }
    return(res)
}


##' Cross-Validation of POET (parallel version)
##'
##' Cross-Validation criterion is the Kullback–Leibler loss (KLL) function. Note
##' that the n observations are partitioned in order (not randomly), i.e. if n =
##' 100, 5-folds cv will use the partitions as follows: [1, ..., 20], [21, ...,
##' 40], [41, ..., 60], [61, ..., 80], [81, ..., 100].
##'
##' If you want to split the data in a random way, you can simply modify this
##' function to meet your requirement.
##' 
##' @param Y p by n data matrix.
##' @param K_seq sequence of K for cv. K is the number of latent factor.
##' @param C_seq sequence of C for cv. C is the shrinkage intensity.
##' @param thres shrinkage method.
##' @param matrix indicator of shrinkage matrix, "cor" (default) or "vad".
##' @param kfolds kfolds-cv.
##' @return A list with two elements
##' \item{params}{params expanded by K_seq and C_seq.}
##' \item{cri}{the calculated criteria across the params using kfolds-cv.}
##' \item{K_opt}{the optimal K over the K_seq.}
##' \item{C_opt}{the optimal C over the C_seq.}
##' 
##' @references Fan, J., Liao, Y., & Mincheva, M. (2013). Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society Series B: Statistical Methodology, 75(4), 603-680.
##'
##' @examples
##' n <- 50; p <- 100; k <- 3
##' set.seed(1)
##' B <- matrix(rnorm(100 * 3), 100, 3)
##' F <- MASS::mvrnorm(n, rep(0, k), diag(1, k))
##' U <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
##' Y <- t(F %*% t(B) + U)
##' res_cv1 <- poetREV::cv_poet(Y, 1:5, seq(.5, 1.5, .1), "soft", "cor", kfolds = 10)
##' 
##' library(doParallel, quietly = TRUE)
##' ## cl <- makeCluster(detectCores()-1)
##' cl <- makeCluster(2) ## So funny... https://stackoverflow.com/questions/15648772
##' registerDoParallel(cl)
##' res_cv2 <- poetREV::cv_poet_par(Y, 1:5, seq(.5, 1.5, .1), "soft", "cor", kfolds = 10)
##' stopCluster(cl)
##' 
##' res_cv1$params - res_cv2$params
##' res_cv1$cri - res_cv2$cri
##' 
##' ## We can obtain K by Information Criteria first, then cv the C_seq
##' ## To reduce the computational burden
##' K_opt <- floor(mean(poetREV::poet_Khat(Y)))
##' res_cv3 <- poetREV::cv_poet(Y, K_opt, seq(.5, 1.5, .1), "soft", "cor", kfolds = 10)
##' 
##' @seealso \code{\link[poetREV]{cv_poet}}
##' @export
cv_poet_par <- function(Y, K_seq, C_seq, thres = "soft", matrix = "cor", kfolds = 5) {
    Y <- as.matrix(Y)
    n <- ncol(Y)
    index <- c(rep(1:kfolds, each = n%/%kfolds), rep(kfolds, n%%kfolds))
    params <- expand.grid(K = K_seq, C = C_seq)
    ## require(doParallel, quietly = TRUE)
    i <- NULL
    cri <- foreach(i = seq_len(nrow(params)), .combine = "c", .packages = "poetREV") %dopar% {
        subcri <- rep(NA, kfolds)
        for (j in seq_len(kfolds)) {
            Y_train <- Y[, index != j]
            Y_test <- Y[, index == j]
            S_test <- cov(t(Y_test))
            K <- params[i, 1]
            C <- params[i, 2]
            res_poet <- poetREV::poet(Y_train, K, C, thres, matrix)
            SigmaY <- res_poet$SigmaY
            evals <- eigen(SigmaY, only.values = TRUE)$values
            ## IMPORTANT: NA treatment is necessary
            if (all(evals > 0)) {
                subcri[j] <- sum(log(evals)) + sum(diag(S_test %*% solve(SigmaY)))
            } else {
                subcri[j] <- NA
            }
        }
        if (all(is.na(subcri))) Inf else mean(subcri, na.rm = TRUE)
    }
    loc <- which.min(cri)
    param_opt <- params[loc, ]
    return(list(params = params, cri = cri, K_opt = param_opt[1], C_opt = param_opt[2]))
}
