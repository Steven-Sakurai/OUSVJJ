// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <random>
#include <cmath>

using namespace std;
typedef vector< double > stdVec;

// [[Rcpp::export]]
List SimOUSVJJ(NumericVector par, int len, double y0, double v0, double Delta, 
    unsigned int seed1, unsigned int seed2, unsigned int seed3, unsigned int seed4, unsigned int seed5) {
	NumericVector y(len);
	stdVec v(len);

	default_random_engine gen1{seed1};
	default_random_engine gen2{seed2};
	default_random_engine gen3{seed3};
	default_random_engine gen4{seed4};
    default_random_engine gen5{seed5};

	double mu = par[0], kappa = par[1], theta = par[2], alpha = par[3], rho = par[4], 
		l0 = par[5], l1 = par[6], rho_j = par[7],
        mu_s = par[8], sd_s = sqrt(par[9]), mu_v = par[10], sd_v = sqrt(par[11]);

	normal_distribution<double> norm{0,1};
        
    y[0] = y0; v[0] = v0;
	for(int i = 1; i < len; ++i) {
		double rnorm1 = norm(gen1);
		double rnorm2 = rho*rnorm1 + sqrt(1-rho*rho)*norm(gen2);

        v[i] = v[i-1] + kappa*(theta - v[i-1])*Delta + alpha*sqrt(Delta)*rnorm1;

        y[i] = y[i-1] + (mu - (0.5)*v[i-1]*v[i-1])*Delta + v[i-1]*sqrt(Delta)*rnorm2;    

		poisson_distribution<int> pois((l0 + l1*v[i-1])*Delta);
        int njumps = pois(gen3);
        
        while(njumps != 0) {
            double z1 = norm(gen4);
            double z2 = norm(gen5); 
            y[i] += z1*sd_s + mu_s;
            v[i] += ( rho_j*z1 + sqrt(1 - (rho_j*rho_j))*z2 ) * sd_v + mu_v;
            --njumps;
        }    
	}
    List ret;
    ret["y"] = y;
    ret["v"] = v;
    return ret;
}