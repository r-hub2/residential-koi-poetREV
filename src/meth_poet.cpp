#include <limits>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace arma;
using namespace Rcpp;

//' Determine the number of latent factors in POET (approximate factor model)
//'
//' This function is for calculating the optimal number of factors in an
//' approximate factor model.
//'
//' (Copy from POET package) This method was proposed by Bai & Ng (2002) and
//' Hallin & Liska (2007). They propose two penalty functions and in turn
//' minimize the corresponding information criteria. Notice that this method may
//' underestimate K. POET is very robust to over-estimating K. But
//' under-estimating K can result to VERY BAD performance. Therefore we strongly
//' recommend choosing a relatively large K (normally less than 8) to avoid
//' missing any important common factor.
//'
//' IMPORTANT: the calculation of Hallin & Liska method (2007) is wrong in the
//' original POET package, here we only calculate the criteria from Bai & Ng
//' (2002).
//' 
//' @param Y p by n matrix of raw data, where p is the dimensionality, n is the
//'   sample size. It is recommended that Y is de-meaned, i.e., each row has
//'   zero mean.
//' @param kmax the given upper bound of number of factors. The optimal k is
//'   searched from 1 to kmax.
//' @return A vector with two elements
//' \item{K1BN}{(The 1st element) estimated number of factors based on the first infomation criterion using Bai & Ng method}
//' \item{K2BN}{(The 2nd element) estimated number of factors based on the second infomation criterion using Bai & Ng method}
//'
//' @references Hallin, M., & Liška, R. (2007). Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603-617.
//' @references Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.
//'
//' @examples
//' n <- 50; p <- 100; k <- 3
//' set.seed(1)
//' B <- matrix(rnorm(100 * 3), 100, 3)
//' F <- MASS::mvrnorm(n, rep(0, k), diag(1, k))
//' U <- MASS::mvrnorm(n, rep(0, p), diag(1, p))
//' Y <- t(F %*% t(B) + U)
//' poetREV::poet_Khat(Y)
//' poetREV::poet_Khat2(Y)
// [[Rcpp::export]]
arma::rowvec poet_Khat(arma::mat Y, int kmax = 15) {
    int p = Y.n_rows, n = Y.n_cols;
    kmax = min(kmax, p);
    for (int j = 0; j < p; j++) {
        Y.row(j) -= mean(Y.row(j));
    }
    vec d; mat V;
    eig_sym(d, V, Y.t() * Y, "std");
    d = reverse(d); V = fliplr(V);
    mat F, LamPCA, uhat, IC(2, kmax);
    double frob;
    double gT1BN = log((p*n)/(double)(p+n)) * (p+n) / (p*n);
    double gT2BN = log(min(p, n)) * (p+n) / (p*n);
    for (int k = 0; k < kmax; k++) {
        F = V.cols(0, k) * sqrt(n);
        LamPCA = Y * F / n;
        uhat = Y - LamPCA * F.t();
        frob = sum(diagvec(uhat * uhat.t())) / (p*n);
        IC(0, k) = log(frob) + (k+1)*gT1BN;
        IC(1, k) = log(frob) + (k+1)*gT2BN;
    }
    rowvec out(2);
    out(0) = IC.row(0).index_min() + 1;
    out(1) = IC.row(1).index_min() + 1;
    return out;
}

double sign(double x) {
    if (x > 0.0) return 1;
    else if (x < 0.0) return -1;
    else return 0.0;
}

double soft_thresh(double x, double a) {
    // NOTE: This version is not symmetric if a < 0 (further used in poet_Cmin)
    // if (x > a) return x - a;
    // else if (x < -a) return x + a;
    // else return 0.0;
    return sign(x) * max(abs(x) - a, 0.0);
}

//' POET Revised by Wenliang Ding
//'
//' Estimates large covariance matrices in approximate factor models by
//' thresholding principal orthogonal complements.
//'
//' (Copy from POET package) This function is for POET, proposed by Fan, Liao
//' and Mincheva (2012) 'Large Covariance Estimation by Thresholding Principal
//' Orthogonal Complements', manuscript of Princeton University
//'
//' Model: Y_t=Bf_t+u_t, where B, f_t and u_t represent factor loading matrix,
//' common factors and idiosyncratic error respectively. Only Y_t is
//' observable. t=1,...,n. Dimension of Y_t is p. The goal is to estimate the
//' covariance matrices of Y_t and u_t.
//'
//' Note: (1) POET is optimization-free, so no initial value, tolerant, or
//' maximum iterations need to be specified as inputs.
//'
//' (2) We can apply the adaptive thresholding (Cai and Liu 2011, JASA) on
//' either the correlation matrix or the covariance matrix, specified by the
//' option 'matrix'.
//'
//' (3) If no factor structure is assumed, i.e., no common factors exist and
//' var(Y_t) itself is sparse, set K=0.
//'
//' IMPORTANT: the calculation of Hallin & Liska method (2007) is wrong in the
//' original POET package, here we only calculate the criteria from Bai & Ng
//' (2002).
//' 
//' @param Y p by n matrix of raw data, where p is the dimensionality, n is the
//'   sample size. It is recommended that Y is de-meaned, i.e., each row has
//'   zero mean.
//' @param K number of factors. K is pre-determined by the users. Default value
//'   is set at the average value obtained from the Hallin&Liska and Bai&Ng
//'   methods. Suggestions on choosing K:
//'
//' A simple way of determining K is to count the number of very spiked (much
//' larger than others) eigenvalues of the p by p sample covariance matrix of Y.
//'
//' A formal data-driven way of determining K is described in Bai and Ng
//' (2002):"Determining the number of factors in approximate factor models",
//' Econometrica, 70, 191-221. This procedure requires a one-dimensional
//' optimization.
//'
//' POET is very robust to over-estimating K. But under-estimating K can result
//' to VERY BAD performance. Therefore we strongly recommend choosing a
//' relatively large K (normally less than 8) to avoid missing any important
//' common factor.
//' 
//' K=0 corresponds to threshoding the sample covariance directly.
//' @param C the positive constant for thresholding, user-specified. Default
//'   value is set at C=0.5 Our experience shows that C=0.5 performs quite well
//'   for soft thresholding.
//' @param thres choice of thresholding. Users can choose from three
//'   thresholding methods:
//' 
//' 'soft': soft thresholding;
//' 
//' 'hard' hard thresholding;
//' 
//' 'scad': scad thresholding;
//' 
//' 'alasso': adaptive lasso thresholding;
//' 
//' Default value is set at thres='soft'.
//' 
//' Details are found in Rothman et al. (2009): "Generalized thresholding of large covariance matrices." JASA, 104, 177-186
//' @param matrix the option of thresholding either correlation or covairance matrix. Users can choose from:
//' 
//' 'cor': threshold the error correlation matrix then transform back to covariance matrix
//' 
//' 'vad': threshold the error covariance matrix directly.
//' 
//' Default value is set at matrix='cor'.
//' 
//' @return A list with two elements
//' \item{SigmaY}{estimated p by p covariance matrix of y_t}
//' \item{SigmaU}{estimated p by p covariance matrix of u_t}
//'
//' @references Fan, J., Liao, Y., & Mincheva, M. (2013). Large covariance estimation by thresholding principal orthogonal complements. Journal of the Royal Statistical Society Series B: Statistical Methodology, 75(4), 603-680.
//'
//' @examples
//' n <- 50; p <- 100
//' set.seed(1)
//' Y <- t(MASS::mvrnorm(n, rep(0, p), diag(1, p)))
//' res_poet <- poetREV::poet(Y, NULL, .5, "soft", "cor")
//' ## res_poet$SigmaY
//' ## res_poet$SigmaU
//' res_poet2 <- poetREV::poet2(Y, -Inf, .5, "soft", "cor")
//' ## res_poet2$SigmaY
//' ## res_poet2$SigmaU
//' mean(abs(res_poet$SigmaY - res_poet2$SigmaY)) ## [1] 5.15726e-16
//' mean(abs(res_poet$SigmaU - res_poet2$SigmaU)) ## [1] 1.015261e-16
//' 
// [[Rcpp::export]]
List poet(arma::mat Y, Nullable<NumericVector> K = R_NilValue, double C = .5, std::string thres = "soft", std::string matrix = "cor") {
    int p = Y.n_rows, n = Y.n_cols;
    // Y is de-meaned by row
    for (int i = 0; i < p; i++) {
        Y.row(i) -= mean(Y.row(i));
    }    
    // Set the number of factors
    int K_;
    if (K.isNotNull()) {
        NumericVector getk(K.get());
        K_ = getk[0];
        if (K_ > p) stop("Invalid K argument in poet(). K must <= p.");
    } else {
        K_ = floor(mean(poet_Khat(Y))) + 1;
        K_ = min(K_, p);
    }
    mat uhat(p, n), Lowrank(p, p);
    double rate;
    if (K_ > 0) {
        vec d; mat V;
        eig_sym(d, V, Y.t() * Y, "std");
        d = reverse(d); V = fliplr(V);
        // using the first K pc as a proxy for the space that
        // spanned by the columns of the factor loading matrix
        mat F = V.cols(0, K_-1) * sqrt(n);
        mat LamPCA = Y * F / n;
        uhat = Y - LamPCA * F.t(); // uhat is p by n
        Lowrank = LamPCA * LamPCA.t();
        rate = 1/sqrt(p) + sqrt(log(p)/n);
    } else if (K_ == 0) {
        // SigmaY itself is sparse
        uhat = Y;
        Lowrank.fill(0.0);
        rate = sqrt(log(p)/n);
    } else {
        stop("Invalid K argument in poet(). K must >= 0.");
    }
    mat SuPCA = uhat * uhat.t() / n;
    vec SuDiag = SuPCA.diag();
    mat R(p, p);
    if (matrix == "cor") {
        R = diagmat(1/sqrt(SuDiag)) * SuPCA * diagmat(1/sqrt(SuDiag));
    } else if (matrix == "vad") {
        R = SuPCA;
    } else {
        stop("Invalid matrix argument in poet().");
    }
    // adaptive threshold
    mat lambda(p, p);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j <= i; j++) {
            // original POET package is wrong for thres = "cor"
            // See, An Overview on the Estimation of Large Covariance and Precision Matrices, Page 5
            if (matrix == "cor") {lambda(i, j) = rate * C;}
            if (matrix == "vad") {lambda(i, j) = sqrt(var(uhat.row(i) % uhat.row(j))*(n-1)/(double)n) * rate * C;}
            lambda(j, i) = lambda(i, j);
        }
    }
    mat Rthresh(p, p);
    if (thres == "soft") {
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                if (j < i) Rthresh(i, j) = soft_thresh(R(i, j), lambda(i, j));
                else Rthresh(i, j) = R(i, j);
                Rthresh(j, i) = Rthresh(i, j);
            }
        }
    } else if (thres == "hard") {
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                if (abs(R(i, j)) <= lambda(i, j) && j < i) {
                    Rthresh(i, j) = 0.0;
                } else {
                    Rthresh(i, j) = R(i, j);
                }
                Rthresh(j, i) = Rthresh(i, j);
            }
        }
    } else if (thres == "scad") {
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                if (j == i) {
                    Rthresh(i, j) = R(i, j);
                } else if (abs(R(i, j)) <= lambda(i, j)) {
                    Rthresh(i, j) = 0.0;
                } else if (abs(R(i, j)) > lambda(i, j) && abs(R(i, j)) <= 2*lambda(i, j)) {
                    Rthresh(i, j) = sign(R(i, j)) * (abs(R(i, j)) - lambda(i, j));
                } else if (abs(R(i, j)) > 2*lambda(i, j) && abs(R(i, j)) <= 3.7*lambda(i, j)) {
                    Rthresh(i, j) = (2.7*R(i, j) - sign(R(i, j))*3.7*lambda(i, j)) / 1.7;
                } else {
                    Rthresh(i, j) = R(i, j);
                }
                Rthresh(j, i) = Rthresh(i, j);
            }
        }
    }
    // calculate SigmaU
    mat SigmaU(p, p);
    if (matrix == "cor") SigmaU = diagmat(sqrt(SuDiag)) * Rthresh * diagmat(sqrt(SuDiag));
    if (matrix == "vad") SigmaU = Rthresh;
    // calculate SigmaY
    mat SigmaY = Lowrank + SigmaU;
    return List::create(
        Named("SigmaY") = SigmaY,
        Named("SigmaU") = SigmaU
        );
}


//' Cross-Validation of POET (non-parallel, Cpp version)
//'
//' Cross-Validation criterion is the Kullback–Leibler loss (KLL) function. Note
//' that the n observations are partitioned in order (not randomly), i.e. if n =
//' 100, 5-folds cv will use the partitions as follows: [1, ..., 20], [21, ...,
//' 40], [41, ..., 60], [61, ..., 80], [81, ..., 100].
//'
//' Cpp version of \code{\link[poetREV]{cv_poet_par}}, see
//' \code{\link[poetREV]{cv_poet_par}} for examples.
//' 
//' @param Y p by n data matrix.
//' @param K_seq sequence of K for cv. K is the number of latent factor.
//' @param C_seq sequence of C for cv. C is the shrinkage intensity.
//' @param thres shrinkage method.
//' @param matrix indicator of shrinkage matrix, "cor" (default) or "vad".
//' @param kfolds kfolds-cv.
//' @return A list with two elements
//' \item{params}{params expanded by K_seq and C_seq.}
//' \item{cri}{the calculated criteria across the params using kfolds-cv.}
//' \item{K_opt}{the optimal K over the K_seq.}
//' \item{C_opt}{the optimal C over the C_seq.}
//' @seealso \code{\link[poetREV]{cv_poet_par}}
// [[Rcpp::export]]
List cv_poet(arma::mat Y, arma::vec K_seq, arma::vec C_seq, std::string thres = "soft", std::string matrix = "cor", int kfolds = 5) {
    // obtain index of cv data
    int n = Y.n_cols, subn = floor(n/(double)kfolds);
    vec index(n); index.fill(kfolds);
    for (int j = 0; j < kfolds-1; j++) {
        for (int i = 0; i < subn; i++) {
            index(i+j*subn) = j+1;
        }
    }
    // obtain params for cv
    int K_len = K_seq.size(), C_len = C_seq.size();
    int n_param = K_len * C_len;
    mat params(n_param, 2);
    for (int i = 0; i < C_len; i++) {
        for (int j = 0; j < K_len; j++) {
            params(j + K_len*i, 0) = K_seq(j);
            params(j + K_len*i, 1) = C_seq(i);
        }
    }
    // perform cv
    vec cri(n_param);
    Rcout << "Perform cv for POET..." << endl;
    Progress prog(n_param, true);
    for (int i = 0; i < n_param; i++) {
        vec subcri(kfolds); subcri.fill(0.0);
        ivec flag(kfolds); flag.fill(1);
        for (int j = 1; j <= kfolds; j++) {
            // extract train and test samples
            uvec idx = find(index != j);
            mat Y_train = Y.cols(idx);
            idx = find(index == j);
            mat Y_test = Y.cols(idx);
            // perform POET
            NumericVector K = NumericVector::create(params(i, 0));
            double C = params(i, 1);
            List res_poet = poet(Y_train, K, C, thres, matrix);
            mat SigmaY = res_poet["SigmaY"];
            mat S_test = cov(Y_test.t());
            vec d; mat V;
            eig_sym(d, V, SigmaY, "std");
            d = reverse(d); V = fliplr(V);
            // IMPORTANT: NA treatment is necessary
            if (all(d > 0.0)) {
                subcri(j-1) = sum(log(d)) + sum(diagvec(S_test * SigmaY.i()));
            } else {
                flag(j-1) = 0;
            }
        }
        uvec idx = find(flag == 1);
        if (idx.size() == 0) {
            cri(i) = numeric_limits<double>::infinity();
        } else {
            cri(i) = mean(subcri(idx));
        }
        // Rcout << i+1 << "th finished, total:" << n_param << endl;
        prog.increment(); // update progress
    }
    uvec idx = find(cri == min(cri));
    mat param_opt = params.rows(idx);
    return List::create(
        Named("params") = params,
        Named("cri") = cri,
        Named("K_opt") = param_opt(0, 0),
        Named("C_opt") = param_opt(0, 1)
        );
}
