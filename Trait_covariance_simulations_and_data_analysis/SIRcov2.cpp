#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericVector bivariateNormal(NumericVector means, NumericVector sds, double correl) {
  NumericVector Mu(2); // create a 2 x 1 matrix of zeros to hold the means 
  NumericMatrix Sigma(2,2); // create a 2 x 2 matrix of zeros to hold the covariance matrix
  // To produce a distribution with mean equal to traitmean
  // and variance equal to traitsd^2, the parameters of
  // the corresponding lognormal distribution are:
  Mu(0) = log(pow(means(0),2)/sqrt(pow(sds(0),2)+pow(means(0),2)));
  Mu(1) = log(pow(means(1),2)/sqrt(pow(sds(1),2)+pow(means(1),2)));
  Sigma(0,0) = sqrt(log(1+pow(sds(0),2)/pow(means(0),2)));
  Sigma(0,1) = 0.0;
  Sigma(1,0) = 0.0;
  Sigma(1,1) = sqrt(log(1+pow(sds(1),2)/pow(means(1),2)));
  NumericMatrix Corr(2,2);
  Corr(0,1) = correl;
  Corr(1,0) = correl;
  Corr(0,0) = 1;
  Corr(1,1) = 1;
  Eigen::MatrixXd Cov(2,2);
  // generate the covariance matrix
  Cov(0,0) = Sigma(0,0)*(Corr(0,0)*Sigma(0,0)+Corr(1,0)*Sigma(0,1)) + (Corr(0,1)*Sigma(0,0) + Corr(1,1)*Sigma(0,1))*Sigma(1,0);
  Cov(0,1) = Sigma(0,1)*(Corr(0,0)*Sigma(0,0)+Corr(1,0)*Sigma(0,1)) + (Corr(0,1)*Sigma(0,0) + Corr(1,1)*Sigma(0,1))*Sigma(1,1);
  Cov(1,0) = Sigma(0,0)*(Corr(0,0)*Sigma(1,0)+Corr(1,0)*Sigma(1,1)) + (Corr(0,1)*Sigma(1,0) + Corr(1,1)*Sigma(1,1))*Sigma(1,0);
  Cov(1,1) = Sigma(0,1)*(Corr(0,0)*Sigma(1,0)+Corr(1,0)*Sigma(1,1)) + (Corr(0,1)*Sigma(1,0) + Corr(1,1)*Sigma(1,1))*Sigma(1,1);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Cov);
  Eigen::VectorXd Eigenvalues = es.eigenvalues();
  Eigen::MatrixXd Eigenvectors = es.eigenvectors();
  Eigen::MatrixXd DiagSqrtEval(2,2);
  DiagSqrtEval(0,0) = sqrt(Eigenvalues(0));
  DiagSqrtEval(1,1) = sqrt(Eigenvalues(1));
  DiagSqrtEval(0,1) = 0.0;
  DiagSqrtEval(1,0) = 0.0;
  Eigen::MatrixXd RootCov = Eigenvectors * DiagSqrtEval * Eigenvectors.transpose();
  NumericVector X = rnorm(2);
  NumericVector X2(2);
  X2(0) = RootCov(0,0)*X(0) + RootCov(0,1)*X(1);
  X2(1) = RootCov(1,0)*X(0) + RootCov(1,1)*X(1);
  NumericVector X3(2); 
  X3(0) = exp(Mu(0) + X2(0));
  X3(1) = exp(Mu(1) + X2(1));
  return X3;
}

// [[Rcpp::export]]
List SIRcovCS(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double cSD = params(4); // SD in contact rate
  double sSD = params(5); // SD in shedding rate
  double corr = params(6); // correlation between c and s
  NumericVector traitmeans = {c,s};
  NumericVector traitsds = {cSD,sSD};
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 5); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: c value; column 2: s value; column 3: no. of 2 infections; column 4: alive or not? (alive = 0; not alive=Inf)
  for(int i=0; i<100000; i++) {
    storage(i,4) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    NumericVector newTraits = bivariateNormal(traitmeans, traitsds, corr);
    storage(i,1) = newTraits(0);
    storage(i,2) = newTraits(1);
    storage(i,3) = 0;
    storage(i,4) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, cVec, sVec, aliveVec; 
  NumericVector idAlive, cAlive, sAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) { // simulate until epidemic ends, however long that takes
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    cVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all c values
    sVec = subStorage(_, 2); //[ Range(0, ID-1), 2]; // all s values
    aliveVec = subStorage(_, 4); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    cAlive = cVec[!is_na(aliveVec)];
    sAlive = sVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate = cAlive * sAlive/(1 + sAlive) * S; // secondary infection rate for each infected individual
    NumericVector gRate(I, g); // rate of recovery for each infected individual
    NumericVector dRateI(I, a); // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
     
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 3) += 1;
      // uncomment the next few lines if you want to allow heritability => evolution
      // get the traits of the individual causing the infection
      // iTraits = {cAlive(event), sAlive(event)};
      // nTraits = bivariateNormal(iTraits, traitsds, corr);
      // get the traits of the newly infected individual
      nTraits = bivariateNormal(traitmeans, traitsds, corr);
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = nTraits(0);
      storage(ID, 2) = nTraits(1);
      storage(ID, 3) = 0;
      storage(ID, 4) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
     
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}

// [[Rcpp::export]]
List SIRcovCA(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double cSD = params(4); // SD in contact rate
  double aSD = params(5); // SD in virulence
  double corr = params(6); // correlation between c and a
  NumericVector traitmeans = {c,a};
  NumericVector traitsds = {cSD,aSD};
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 5); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: c value; column 2: s value; column 3: no. of 2 infections; column 4: alive or not? (alive = 0; not alive=Inf)
  for(int i=0; i<100000; i++) {
    storage(i,4) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    NumericVector newTraits = bivariateNormal(traitmeans, traitsds, corr);
    storage(i,1) = newTraits(0);
    storage(i,2) = newTraits(1);
    storage(i,3) = 0;
    storage(i,4) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, cVec, aVec, aliveVec; 
  NumericVector idAlive, cAlive, aAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) { // simulate until epidemic ends, however long that takes
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    cVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all c values
    aVec = subStorage(_, 2); //[ Range(0, ID-1), 2]; // all a values
    aliveVec = subStorage(_, 4); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    cAlive = cVec[!is_na(aliveVec)];
    aAlive = aVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate = cAlive * s/(1 + s) * S; // secondary infection rate for each infected individual
    NumericVector gRate(I, g); // rate of recovery for each infected individual
    NumericVector dRateI = aAlive; // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 3) += 1;
      // uncomment the next few lines if you want to allow heritability => evolution
      // get the traits of the individual causing the infection
      // iTraits = {cAlive(event), sAlive(event)};
      // nTraits = bivariateNormal(iTraits, traitsds, corr);
      // get the traits of the newly infected individual
      nTraits = bivariateNormal(traitmeans, traitsds, corr);
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = nTraits(0);
      storage(ID, 2) = nTraits(1);
      storage(ID, 3) = 0;
      storage(ID, 4) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}


// [[Rcpp::export]]
List SIRcovCG(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double cSD = params(4); // SD in contact rate
  double gSD = params(5); // SD in recovery rate
  double corr = params(6); // correlation between c and g
  NumericVector traitmeans = {c,g};
  NumericVector traitsds = {cSD,gSD};
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 5); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: c value; column 2: s value; column 3: no. of 2 infections; column 4: alive or not? (alive = 0; not alive=Inf)
  for(int i=0; i<100000; i++) {
    storage(i,4) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    NumericVector newTraits = bivariateNormal(traitmeans, traitsds, corr);
    storage(i,1) = newTraits(0);
    storage(i,2) = newTraits(1);
    storage(i,3) = 0;
    storage(i,4) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, cVec, gVec, aliveVec; 
  NumericVector idAlive, cAlive, gAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) { // simulate until epidemic ends, however long that takes
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    cVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all c values
    gVec = subStorage(_, 2); //[ Range(0, ID-1), 2]; // all s values
    aliveVec = subStorage(_, 4); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    cAlive = cVec[!is_na(aliveVec)];
    gAlive = gVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate = cAlive * s/(1 + s) * S; // secondary infection rate for each infected individual
    NumericVector gRate = gAlive; // rate of recovery for each infected individual
    NumericVector dRateI(I, a); // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 3) += 1;
      // uncomment the next few lines if you want to allow heritability => evolution
      // get the traits of the individual causing the infection
      // iTraits = {cAlive(event), sAlive(event)};
      // nTraits = bivariateNormal(iTraits, traitsds, corr);
      // get the traits of the newly infected individual
      nTraits = bivariateNormal(traitmeans, traitsds, corr);
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = nTraits(0);
      storage(ID, 2) = nTraits(1);
      storage(ID, 3) = 0;
      storage(ID, 4) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}

// [[Rcpp::export]]
List SIRcovSA(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double sSD = params(4); // SD in shedding rate
  double aSD = params(5); // SD in virulence
  double corr = params(6); // correlation between s and a
  NumericVector traitmeans = {s,a};
  NumericVector traitsds = {sSD,aSD};
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 5); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: c value; column 2: s value; column 3: no. of 2 infections; column 4: alive or not? (alive = 0; not alive=Inf)
  for(int i=0; i<100000; i++) {
    storage(i,4) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    NumericVector newTraits = bivariateNormal(traitmeans, traitsds, corr);
    storage(i,1) = newTraits(0);
    storage(i,2) = newTraits(1);
    storage(i,3) = 0;
    storage(i,4) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, sVec, aVec, aliveVec; 
  NumericVector idAlive, sAlive, aAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) { // simulate until epidemic ends, however long that takes
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    sVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all s values
    aVec = subStorage(_, 2); //[ Range(0, ID-1), 2]; // all a values
    aliveVec = subStorage(_, 4); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    sAlive = sVec[!is_na(aliveVec)];
    aAlive = aVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate = c * sAlive/(1 + sAlive) * S; // secondary infection rate for each infected individual
    NumericVector gRate(I, g); // rate of recovery for each infected individual
    NumericVector dRateI = aAlive; // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 3) += 1;
      // uncomment the next few lines if you want to allow heritability => evolution
      // get the traits of the individual causing the infection
      // iTraits = {cAlive(event), sAlive(event)};
      // nTraits = bivariateNormal(iTraits, traitsds, corr);
      // get the traits of the newly infected individual
      nTraits = bivariateNormal(traitmeans, traitsds, corr);
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = nTraits(0);
      storage(ID, 2) = nTraits(1);
      storage(ID, 3) = 0;
      storage(ID, 4) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}

// [[Rcpp::export]]
List SIRcovSG(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double sSD = params(4); // SD in shedding rate
  double gSD = params(5); // SD in recovery rate
  double corr = params(6); // correlation between s and g
  NumericVector traitmeans = {s,g};
  NumericVector traitsds = {sSD,gSD};
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 5); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: c value; column 2: s value; column 3: no. of 2 infections; column 4: alive or not? (alive = 0; not alive=Inf)
  for(int i=0; i<100000; i++) {
    storage(i,4) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    NumericVector newTraits = bivariateNormal(traitmeans, traitsds, corr);
    storage(i,1) = newTraits(0);
    storage(i,2) = newTraits(1);
    storage(i,3) = 0;
    storage(i,4) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, sVec, gVec, aliveVec; 
  NumericVector idAlive, sAlive, gAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) { // simulate until epidemic ends, however long that takes
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    sVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all s values
    gVec = subStorage(_, 2); //[ Range(0, ID-1), 2]; // all a values
    aliveVec = subStorage(_, 4); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    sAlive = sVec[!is_na(aliveVec)];
    gAlive = gVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate = c * sAlive/(1 + sAlive) * S; // secondary infection rate for each infected individual
    NumericVector gRate = gAlive; // rate of recovery for each infected individual
    NumericVector dRateI(I, a); // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 3) += 1;
      // uncomment the next few lines if you want to allow heritability => evolution
      // get the traits of the individual causing the infection
      // iTraits = {cAlive(event), sAlive(event)};
      // nTraits = bivariateNormal(iTraits, traitsds, corr);
      // get the traits of the newly infected individual
      nTraits = bivariateNormal(traitmeans, traitsds, corr);
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = nTraits(0);
      storage(ID, 2) = nTraits(1);
      storage(ID, 3) = 0;
      storage(ID, 4) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}


// [[Rcpp::export]]
List SIRcovAG(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double aSD = params(4); // SD in virulence
  double gSD = params(5); // SD in recovery rate
  double corr = params(6); // correlation between a and g
  NumericVector traitmeans = {a,g};
  NumericVector traitsds = {aSD,gSD};
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 5); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: c value; column 2: s value; column 3: no. of 2 infections; column 4: alive or not? (alive = 0; not alive=Inf)
  for(int i=0; i<100000; i++) {
    storage(i,4) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    NumericVector newTraits = bivariateNormal(traitmeans, traitsds, corr);
    storage(i,1) = newTraits(0);
    storage(i,2) = newTraits(1);
    storage(i,3) = 0;
    storage(i,4) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, aVec, gVec, aliveVec; 
  NumericVector idAlive, aAlive, gAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) { // simulate until epidemic ends, however long that takes
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    aVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all s values
    gVec = subStorage(_, 2); //[ Range(0, ID-1), 2]; // all a values
    aliveVec = subStorage(_, 4); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    aAlive = aVec[!is_na(aliveVec)];
    gAlive = gVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate(I, c * s/(1 + s) * S); // secondary infection rate for each infected individual
    NumericVector gRate = gAlive; // rate of recovery for each infected individual
    NumericVector dRateI = aAlive; // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 3) += 1;
      // uncomment the next few lines if you want to allow heritability => evolution
      // get the traits of the individual causing the infection
      // iTraits = {cAlive(event), sAlive(event)};
      // nTraits = bivariateNormal(iTraits, traitsds, corr);
      // get the traits of the newly infected individual
      nTraits = bivariateNormal(traitmeans, traitsds, corr);
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = nTraits(0);
      storage(ID, 2) = nTraits(1);
      storage(ID, 3) = 0;
      storage(ID, 4) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 4) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}



// [[Rcpp::export]]
List SIRvarG(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double gSD = params(4); // SD in recovery rate
  
  // for computing new individuals 
  double mu = log(pow(g,2.0)/sqrt(pow(gSD,2.0)+pow(g,2.0)));
  double sigma = sqrt(log(pow(gSD,2.0)/pow(g,2.0) + 1));
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 4); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: a value; column 2: no. of 2 infections; column 3: alive or not? (alive = 0; not alive=NA)
  for(int i=0; i<100000; i++) {
    storage(i,3) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    storage(i,1) = rlnorm(1, mu, sigma)[0];
    storage(i,2) = 0;
    storage(i,3) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, gVec, aliveVec; 
  NumericVector idAlive, gAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) {
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    gVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all a values
    aliveVec = subStorage(_, 3); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    gAlive = gVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate(I, c * s/(1 + s) * S); // secondary infection rate for each infected individual
    NumericVector gRate = gAlive; // rate of recovery for each infected individual
    NumericVector dRateI(I, a); // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 2) += 1;
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = rlnorm(1, mu, sigma)[0];
      storage(ID, 2) = 0;
      storage(ID, 3) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}

// [[Rcpp::export]]
List SIRvarA(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double aSD = params(4); // SD in virulence
  
  // for computing new individuals 
  double mu = log(pow(a,2.0)/sqrt(pow(aSD,2.0)+pow(a,2.0)));
  double sigma = sqrt(log(pow(aSD,2.0)/pow(a,2.0) + 1));
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 4); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: a value; column 2: no. of 2 infections; column 3: alive or not? (alive = 0; not alive=NA)
  for(int i=0; i<100000; i++) {
    storage(i,3) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    storage(i,1) = rlnorm(1, mu, sigma)[0];
    storage(i,2) = 0;
    storage(i,3) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, aVec, aliveVec; 
  NumericVector idAlive, aAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) {
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    aVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all a values
    aliveVec = subStorage(_, 3); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    aAlive = aVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate(I, c * s/(1 + s) * S); // secondary infection rate for each infected individual
    NumericVector gRate(I, g);  // rate of recovery for each infected individual
    NumericVector dRateI = aAlive; // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 2) += 1;
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = rlnorm(1, mu, sigma)[0];
      storage(ID, 2) = 0;
      storage(ID, 3) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}


// [[Rcpp::export]]
List SIRvarS(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double sSD = params(4); // SD in shedding rate
  
  // for computing new individuals 
  double mu = log(pow(s,2.0)/sqrt(pow(sSD,2.0)+pow(s,2.0)));
  double sigma = sqrt(log(pow(sSD,2.0)/pow(s,2.0) + 1));
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 4); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: a value; column 2: no. of 2 infections; column 3: alive or not? (alive = 0; not alive=NA)
  for(int i=0; i<100000; i++) {
    storage(i,3) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    storage(i,1) = rlnorm(1, mu, sigma)[0];
    storage(i,2) = 0;
    storage(i,3) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, sVec, aliveVec; 
  NumericVector idAlive, sAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) {
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    sVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all a values
    aliveVec = subStorage(_, 3); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    sAlive = sVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate = c * sAlive/(1 + sAlive) * S; // secondary infection rate for each infected individual
    NumericVector gRate(I, g);  // rate of recovery for each infected individual
    NumericVector dRateI(I, a); // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 2) += 1;
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = rlnorm(1, mu, sigma)[0];
      storage(ID, 2) = 0;
      storage(ID, 3) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}


// [[Rcpp::export]]
List SIRvarC(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  double cSD = params(4); // SD in contact rate
  
  // for computing new individuals 
  double mu = log(pow(c,2.0)/sqrt(pow(cSD,2.0)+pow(c,2.0)));
  double sigma = sqrt(log(pow(cSD,2.0)/pow(c,2.0) + 1));
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 4); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: a value; column 2: no. of 2 infections; column 3: alive or not? (alive = 0; not alive=NA)
  for(int i=0; i<100000; i++) {
    storage(i,3) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    storage(i,1) = rlnorm(1, mu, sigma)[0];
    storage(i,2) = 0;
    storage(i,3) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, cVec, aliveVec; 
  NumericVector idAlive, cAlive; 
  NumericVector iTraits, nTraits; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) {
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    //Rcout << "The value of idVec : " << idVec << "\n";
    cVec = subStorage(_, 1); //[ Range(0, ID-1), 1]; // all a values
    aliveVec = subStorage(_, 3); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    //Rcout << "The value of idAlive : " << idAlive << "\n";
    cAlive = cVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate = cAlive * s/(1 + s) * S; // secondary infection rate for each infected individual
    NumericVector gRate(I, g);  // rate of recovery for each infected individual
    NumericVector dRateI(I, a); // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 2) += 1;
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = rlnorm(1, mu, sigma)[0];
      storage(ID, 2) = 0;
      storage(ID, 3) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 3) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}



// [[Rcpp::export]]
List SIRnovar(NumericVector params, IntegerVector x) {
  double c = params(0); // contact rate
  double s = params(1); // shedding rate/infectiousness
  double a = params(2); // virulence
  double g = params(3); // recovery rate
  
  int S=x(0); // initial size of the susceptible population
  int I=x(1); // initial size of the infected population
  int R=x(2); // initial size of the recovered population
  
  // keep track of time
  double t = 0.0;
  NumericVector tVec(1000000); 
  tVec(0) = t;
  
  // set up storage
  IntegerMatrix Popn(1000000,3); // store population sizes and times
  NumericMatrix storage(100000, 3); // store individual traits and fitnesses
  // store initial population sizes
  Popn(0,0) = S; 
  Popn(0,1) = I; 
  Popn(0,2) = R; 
  
  // store initial individuals
  // column 0: ID; column 1: no. of 2 infections; column 2: alive or not? (alive = 0; not alive=NA)
  for(int i=0; i<100000; i++) {
    storage(i,2) = NA_REAL;
  }
  for (int i=0; i < I; i++) {
    storage(i,0) = i;
    storage(i,1) = 0;
    storage(i,2) = 0.0;
  }
  
  // ID to assign to the next infection
  int ID = I;
  // new row to store population sizes
  int nextRow = 1;
  
  NumericMatrix subStorage; 
  NumericVector idVec, aliveVec, idAlive; 
  double totalRate, rand;
  IntegerVector whichEvent;
  int event, thisID;
  
  // GO!
  while(I > 0) {
    //Rcout << "The value of I : " << I << "\n"; 
    subStorage = storage( Range(0, ID-1), _); 
    
    idVec = subStorage(_, 0); // all IDs
    aliveVec = subStorage(_, 2); //[ Range(0, ID-1), 3]; // all alive statuses
    // get the IDs, c, and s values for only the living individuals
    idAlive = idVec[!is_na(aliveVec)];
    
    // compute all rates
    NumericVector iRate(I, c * s/(1 + s) * S); // secondary infection rate for each infected individual
    NumericVector gRate(I, g); // rate of recovery for each infected individual
    NumericVector dRateI(I, a); // mortality rate for each infected individual
    // combine all rates into a single vector
    NumericVector rateVec(iRate.length()+gRate.length()+dRateI.length());
    for(int i=0; i < I; i++) {
      rateVec(i) = iRate(i);
      rateVec(I+i) = gRate(i);
      rateVec(2*I+i) = dRateI(i);
    }
    
    // compute the total rate of events across all patches
    totalRate = std::accumulate(rateVec.begin(), rateVec.end(), 0.0);
    // compute fractional rates
    NumericVector fracRates = rateVec / totalRate;
    // compute cumulative fraction rates
    NumericVector cumsumRates(fracRates.length());
    std::partial_sum(fracRates.begin(), fracRates.end(), cumsumRates.begin());
    
    // generate a random uniform
    rand = runif(1)[0]; 
    // identify which event is happening
    whichEvent = ifelse(rand > cumsumRates, 1, 0);
    event = std::accumulate(whichEvent.begin(), whichEvent.end(), 0);
    
    // increment time
    t += rexp(1, totalRate)[0];
    
    // a new infection occurs
    if (event < I) {
      // get the ID (=row) of the individual causing the infection
      thisID = idAlive(event);
      // increment its number of secondary infections
      storage(thisID, 1) += 1;
      // add the new individual to storage
      storage(ID, 0) = ID;
      storage(ID, 1) = 0;
      storage(ID, 2) = 0.0;
      I += 1; // add one more infected individual to the population
      S -= 1; // remove one susceptible individual
      ID += 1; // increment the next infected ID
    }
    else if (event < 2*I) { // an infected individual recovers
      // get the ID (=row) of the individual recovering
      thisID = idAlive(event-I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 2) = NA_REAL;
      I -= 1; // remove one infected individual from the population
      R += 1; // add one more recovered individual to the population
    }
    else  { // an infected individual dies
      // get the ID (=row) of the individual dying
      thisID = idAlive(event-2*I);
      // set its infection status back to NA so it will no longer be among the living infecteds
      storage(thisID, 2) = NA_REAL;
      I -= 1; // remove one infected individual from the population
    }
    
    // store the population state and increment nextRow
    tVec(nextRow) = t;
    Popn(nextRow,0) = S; 
    Popn(nextRow,1) = I; 
    Popn(nextRow,2) = R; 
    nextRow += 1;
  }
  NumericVector finalT = tVec[ Range(0, nextRow-1) ];
  IntegerMatrix finalPopn = Popn( Range(0,nextRow-1), _);
  NumericMatrix finalInds = storage( Range(0, ID-1), _);
  List finalState(3);
  finalState(0) = finalT;
  finalState(1) = finalPopn;
  finalState(2) = finalInds;
  return finalState;
}




