#include <stdio.h>
#include <iostream>  
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <fstream>
#include <fstream>
#include "matrix.h"
#include "mex.h"


# define M_PIL       3.141592653589793238462643383279502884L /* pi */  
//using Eigen::MatrixXd;  
using namespace Eigen;  
using namespace Eigen::internal;  
using namespace Eigen::Architecture;  
using namespace std;  


double logdet(MatrixXd A)
{
	MatrixXd L=(A).llt().matrixL();
	double y=0;
	for(int i=0;i<L.rows();++i)
		y+=log(L(i,i));
	y=2*y;
	return y;
}

MatrixXd logmvnpdf(MatrixXd data, MatrixXd mu, MatrixXd cov)
{
	int D=(data).cols();
	double const0=-0.5*D*log(2*M_PIL);
	MatrixXd xc = (data).rowwise()-(mu).colwise().mean();
	MatrixXd term1=-0.5 * ((xc *(cov).inverse()).array() * xc.array()).rowwise().sum();
	double term2=const0-0.5*logdet(cov);
	MatrixXd logp=term1.transpose().array()+term2;
	return logp;
}

MatrixXd mylog(MatrixXd A)
{
	MatrixXd logA=A;
	for(int i=0;i<A.rows();++i)
		for(int j=0;j<A.cols();++j)
			logA(i,j)=log(A(i,j));
	return logA;

}

MatrixXd myrow(MatrixXd A, int n)
{
	MatrixXd row(1,A.cols());
	for (int i=0;i<A.cols();++i)
		row(0,i)=A(n,i);
	return row;
}

double newProbability1(MatrixXd XX, MatrixXd Z, MatrixXd Mu, MatrixXd Cov, int K)
{
	int dim=XX.cols();
	MatrixXd nat1=MatrixXd::Zero(dim,dim);
	MatrixXd nat2=MatrixXd::Zero(1,dim);	
	MatrixXd covr;
	RowVectorXd mu; 
	MatrixXd cov;
	MatrixXd covrow;
				
	double newP=0.0;
	

	  for(int n=0;n<XX.rows();++n)
	  {
		nat1=MatrixXd::Zero(dim,dim);
		nat2=MatrixXd::Zero(1,dim);
	    for(int k=0;k<K;++k)
		{
			covrow=myrow(Cov,k);
			covrow.resize(dim,dim);
			nat1+=Z(n,k)*(covrow.inverse());
		    nat2+=Z(n,k)*(myrow(Mu,k))*(covrow.inverse());
		}
		covr=nat1.inverse();
	    mu=nat2*covr;
        newP+=logmvnpdf(myrow(XX,n), mu,covr).sum();
  
	  }
	  return newP;
}
  


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    
	double *data,*membership,*mean,*cov;
	int N, Dim,K,n,dim,k;
	double prob;
	data=mxGetPr(prhs[0]);
	membership=mxGetPr(prhs[1]);
	mean=mxGetPr(prhs[2]);
	cov=mxGetPr(prhs[3]);

	N=mxGetM(prhs[0]);
	Dim=mxGetN(prhs[0]);
	K=mxGetN(prhs[1]);

    MatrixXd XX=MatrixXd::Zero(N,Dim);
	MatrixXd ZZ=MatrixXd::Zero(N,K);
	MatrixXd Mu=MatrixXd::Zero(K,Dim);
	MatrixXd Cov=MatrixXd::Zero(K,Dim*Dim);

	for (n=0;n<N;++n)
		{
			for(dim=0;dim<Dim;++dim)
			   XX(n,dim)=data[n+dim*N];
			for(k=0;k<K;++k)
			   ZZ(n,k)=membership[n+k*N];
	    }
	for(k=0;k<K;++k)
	{
		for(dim=0;dim<Dim;++dim)
			Mu(k,dim)=mean[k+dim*K];
		for(dim=0;dim<Dim*Dim;++dim)
			Cov(k,dim)=cov[k+dim*K];
	}
	prob=newProbability1(XX, ZZ, Mu, Cov, K);
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=prob;
}