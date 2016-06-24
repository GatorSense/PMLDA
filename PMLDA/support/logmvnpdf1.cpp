
#include <stdio.h>
#include <iostream>
#include "mex.h"
#include <Eigen/Dense>
# define M_PIL       3.141592653589793238462643383279502884L /* pi */
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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    
	double *data,*mean,*cov;
	int N, Dim,n,dim,dim1;
	double prob;
	data=mxGetPr(prhs[0]);
	mean=mxGetPr(prhs[1]);
	cov=mxGetPr(prhs[2]);

	N=mxGetM(prhs[0]);
	Dim=mxGetN(prhs[0]);


    MatrixXd XX=MatrixXd::Zero(N,Dim);
	MatrixXd Mu=MatrixXd::Zero(1,Dim);
	MatrixXd Cov=MatrixXd::Zero(Dim,Dim);

	for (n=0;n<N;++n)
		for(dim=0;dim<Dim;++dim)
			XX(n,dim)=data[n+dim*N];

	for(dim=0;dim<Dim;++dim)
			Mu(0,dim)=mean[dim];

	for(dim1=0;dim1<Dim;++dim1)
		for(dim=0;dim<Dim;++dim)
			Cov(dim1,dim)=cov[dim1+dim*Dim];

	MatrixXd result=logmvnpdf(XX,Mu,Cov);
	double results=result.sum();
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=results;
}
