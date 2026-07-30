#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;               	// 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

// [[Rcpp::export]]
VectorXd bivariateNormalEigen(Map<VectorXd> traitmeans, Map<VectorXd> traitsds, Map<VectorXd> corr, Map<VectorXd> randnorm) {
  VectorXd Mu(2); // create a 2 x 1 matrix of zeros to hold the means 
  MatrixXd Sigma(2,2); // create a 2 x 2 matrix of zeros to hold the covariance matrix
  // To produce a distribution with mean equal to traitmean
  // and variance equal to traitsd^2, the parameters of
  // the corresponding lognormal distribution are:
  Mu(0) = log(pow(traitmeans(0),2)/sqrt(pow(traitsds(0),2)+pow(traitmeans(0),2)));
  Mu(1) = log(pow(traitmeans(1),2)/sqrt(pow(traitsds(1),2)+pow(traitmeans(1),2)));
  Sigma(0,0) = sqrt(log(1+pow(traitsds(0),2)/pow(traitmeans(0),2)));
  Sigma(0,1) = 0.0;
  Sigma(1,0) = 0.0;
  Sigma(1,1) = sqrt(log(1+pow(traitsds(1),2)/pow(traitmeans(1),2)));
  MatrixXd Corr(2,2);
  Corr(0,1) = corr(0);
  Corr(1,0) = corr(0);
  Corr(0,0) = 1;
  Corr(1,1) = 1;
  MatrixXd Cov = Sigma * Corr * Sigma; // generate the covariance matrix
  SelfAdjointEigenSolver<MatrixXd> es(Cov);
  VectorXd Eigenvalues = es.eigenvalues();
  MatrixXd Eigenvectors = es.eigenvectors();
  MatrixXd DiagSqrtEval(2,2);
  DiagSqrtEval(0,0) = sqrt(Eigenvalues(0));
  DiagSqrtEval(1,1) = sqrt(Eigenvalues(1));
  DiagSqrtEval(0,1) = 0.0;
  DiagSqrtEval(1,0) = 0.0;
  MatrixXd RootCov = Eigenvectors * DiagSqrtEval * Eigenvectors.transpose();
  VectorXd X(1,2);
  X(0) = randnorm(0);
  X(1) = randnorm(1);
  VectorXd X2 = RootCov * X;
  VectorXd X3(2); 
  X3(0) = exp(Mu(0) + X2(0));
  X3(1) = exp(Mu(1) + X2(1));
  return X3;
}
