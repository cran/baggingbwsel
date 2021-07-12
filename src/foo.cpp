#include <Rcpp.h>
#include <Rmath.h>
#include <Rinternals.h>
#define DELMAX 1000
#define _USE_MATH_DEFINES
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double bw_boot(NumericVector x, double h, double g)
{
	int n = x.size();
	double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, u = 0.0;
	double RK = 0.5 / M_SQRT_PI;
	double c0 = 1.0 / (sqrt(2.0) * M_SQRT_PI);
	double c1 = 1.0 / sqrt(2.0 * h * h + 2.0 * g * g);
	double c2 = 1.0 / sqrt(h * h + 2.0 * g * g);
	double c3 = 1.0 / (sqrt(2.0) * g);
	for(int i = 0; i < n; i++)
	{
		for(int j = 1; j < i; j++)
		{
			double dif = x[i] - x[j];
			double delta1 = dif * c1;
			double delta2 = dif * c2;
			double delta3 = dif * c3;
			delta1 *= delta1;
			delta2 *= delta2;
			delta3 *= delta3;
			sum1 += exp(- delta1 / 2.0);
			sum2 += exp(- delta2 / 2.0);
			sum3 += exp(- delta3 / 2.0);
		}
	}
	u = - 2.0 * c1 * sum1 / (n * n) + 2.0 * (c1 * sum1 - 2.0 * c2 * sum2 + c3 * sum3) / n + c1 * (1.0 - 1.0 / n) - 2.0 * c2 + c3;
	u = (RK / h + c0 * u) / n;
	return u;
}


// [[Rcpp::export]]
double Cbw_boot(int n, double d, NumericVector cnt, double h, double g)
{
	double sum1 = 0.0, sum2 = 0.0, term1 = 0.0, term2 = 0.0, u = 0.0;
	int nbin = cnt.size();

	double c1 = 1.0 / sqrt(2.0 * h * h + 2.0 * g * g);
	double c2 = 1.0 / sqrt(h * h + 2.0 * g * g);
	double invn = 1.0 / n;
	double invnn = invn * invn;
	for(int i = 0; i < nbin; i++)
	{
		double delta = i * d;
		double delta1 = delta * c1;
		double delta2 = delta * c2;
		delta1 *= delta1;
		delta2 *= delta2;
		if(delta1 >= DELMAX && delta2 >= DELMAX) break;
		term1 = exp(- delta1 / 2.);
		term2 = exp(- delta2 / 2.);
		sum1 += term1 * cnt[i];
		sum2 += term2 * cnt[i];
	}
	u = M_SQRT1_2 / h - 2.0 * c1 * sum1 * invnn + 2.0 * (c1 * sum1 - 2.0 * c2 * sum2) * invn + (1.0 - 1.0 / n) * c1 - 2.0 * c2;
	return u;
}


// [[Rcpp::export]]
double Cbw_ucv_nb(NumericVector x, double h)
{
	double sum = 0.0, u = 0.0;
	int n = x.size();
	double invn = 1.0 / n;
	double cte = invn / (h * M_SQRT_PI);
	for(int i = 0; i < n; i++)
	{
		for(int j = 1; j < i; j++)
		{
			double delta = (x[i] - x[j]) / h;
			delta *= delta;
			if(delta >= DELMAX) break;
			sum += exp(-delta / 4.) - sqrt(8.0) * exp(-delta / 2.);
		}
	}
	u = (0.5 + sum * invn) * cte;
	return u;
}


// [[Rcpp::export]]
double Cbw_ucv(int n, double d, NumericVector cnt, double h)
{
	double sum = 0.0, term, u;
	double invn = 1.0 / n, cte = n * h, sqrt8 = sqrt(8.0);
	int nbin = cnt.size();
	for(int i = 0; i < nbin; i++)
	{
		double delta = i * d / h;
		delta *= delta;
		if(delta >= DELMAX) break;
		term = exp(-delta / 4.) - sqrt8 * exp(-delta / 2.);
		sum += term * cnt[i];
	}
	u = (0.5 + sum * invn) / cte;
	return u;
}


//// [[Rcpp::export]]
//double Cbw_ucv(int n, double d, NumericVector cnt, double h)
//{
//	double sum = 0.0, term, u;
//	int nbin = cnt.size();
//	for(int i = 0; i < nbin; i++)
//	{
//		double delta = i * d / h;
//		delta *= delta;
//		if(delta >= DELMAX) break;
//		term = exp(-delta / 4.) - sqrt(8.0) * exp(-delta / 2.);
//		sum += term * cnt[i];
//	}
//	u = (0.5 + sum/n) / (n * h * M_SQRT_PI);
//	return u;
//}


// [[Rcpp::export]]
List bw_den(int nbin, NumericVector x)
{
	int nb = nbin, n = x.size();
	double xmin, xmax, rang, dd;
	xmin = R_PosInf; xmax = R_NegInf;
	for(int i = 0; i < n; i++)
	{
		if(x[i] < xmin) xmin = x[i];
		if(x[i] > xmax) xmax = x[i];
	}
	rang = (xmax - xmin) * 1.01;
	dd = rang / nb;
	NumericVector cnt(nb);
	for(int i = 1; i < n; i++)
	{
		int ii = (int)(x[i] / dd);
		for(int j = 0; j < i; j++)
		{
			int jj = (int)(x[j] / dd);
			cnt[abs(ii - jj)] += 1.0;
		}
	}
	return List::create(dd, cnt);
}

// [[Rcpp::export]]
NumericVector bw_den_binned(IntegerVector sx)
{
	int nb = sx.size();
	NumericVector cnt(nb);

	for(int ii = 0; ii < nb; ii++)
	{
		int w = sx[ii];
		cnt[0] += w*(w-1.);
		for(int jj = 0; jj < ii; jj++) cnt[ii - jj] += w * sx[jj];
	}
	cnt[0] *= 0.5;
	return cnt;
}




// [[Rcpp::export]]
NumericVector prnw(const NumericVector& x, const NumericVector& y, double h, double x0)
{
	int n = x.size();
	double sum1 = 0.0, sum2 = 0.0;
	for(int i = 0; i < n; i++)
	{
		double delta = (x0 - x[i]) / h;
		delta *= delta;
		double term = exp(- delta / 2.);
		sum1 += term;
		sum2 += term * y[i];
	}
	double nwx0 = sum2 / sum1;
	double prx0 = sum1;
	NumericVector out = NumericVector::create(prx0, nwx0);
	return out;
}


// [[Rcpp::export]]
double nw_cv(const NumericVector& x, const NumericVector& y, double h)
{
	int n = x.size();
	double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
	for(int j = 0; j < n; j++)
	{
		double x0 = x[j];
		for(int i = 0; i < n; i++)
		{
			double delta = (x0 - x[i]) / h;
			delta *= delta;
			double term = exp(- delta / 2.);
			sum1 += term;
			sum2 += term * y[i];
		}
		double aux = (y[j] - sum2 / sum1) / (1.0 - 1.0 / sum1);
		sum3 += aux * aux;
		sum1 = 0.0;
		sum2 = 0.0;
	}
	return sum3;
}

// [[Rcpp::export]]
double nw_cv_binning(const NumericVector& x, const NumericVector& y, int nb, double d, double h)
{
    double dd = d / h, cvscore = 0.0;
    for(int i = 0; i < nb; i++)
    {
        if(x[i] == 0) continue;
        double sum1 = 0.0, sum2 = 0.0;
        for(int j = 0; j < nb; j++)
        {
            if(x[j] == 0) continue;
            double delta = (i - j) * dd;
            delta *= delta;
            double kereval = exp(- 0.5 * delta);
            sum1 += kereval * x[j];
            sum2 +=  kereval * y[j];
        }
        double ymean = y[i] / (double)x[i];
        double aux = (ymean - sum2 / sum1) / (1.0 - 1.0 / sum1);

        cvscore += (double)x[i] * aux * aux;
    }
    return cvscore;
}


// [[Rcpp::export]]
double new_nw_cv_binning(const NumericVector& x, const NumericVector& y, const IntegerVector& ind, int nb, double d, double h)
{
    double dd = d / h, cvscore = 0.0;
    for(int i = 0; i < nb; i++)
    {
        double sum1 = 0.0, sum2 = 0.0;
        int ii = ind[i];
        for(int j = 0; j < nb; j++)
        {
            int jj = ind[j];
            double delta = (ii - jj) * dd;
            delta *= delta;
            double kereval = exp(- 0.5 * delta);
            sum1 += kereval * x[j];
            sum2 +=  kereval * y[j];
        }
        double ymean = y[i] / x[i];
        double aux = (ymean - sum2 / sum1) / (1.0 - 1.0 / sum1);

        cvscore += x[i] * aux * aux;
    }
    return cvscore;
}

// [[Rcpp::export]]
double nw_binning(int k, const NumericVector& x, const NumericVector& y, int nb, double d, double h)
{
    double dd = d / h;
        double sum1 = 0.0, sum2 = 0.0;
        for(int j = 0; j < nb; j++)
        {
            double delta = (k - 1 - j) * dd;
            delta *= delta;
            double kereval = exp(- delta / 2.0);
            sum1 += kereval * x[j];
            sum2 +=  kereval * y[j];
        }
    return sum2/sum1;
}
