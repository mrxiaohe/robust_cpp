#include <R.h>
#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

IntegerVector order(NumericVector x);

template <class RandomAccessIterator, class StrictWeakOrdering>
void sort(RandomAccessIterator first, RandomAccessIterator last, StrictWeakOrdering comp);

struct val_order{
		int order;
		double value;
};

bool compare(const val_order & a, const val_order & b){return (a.value<b.value);}


// arma in; arma out
arma::mat sumbatC1(arma::mat X, arma::colvec T) {
    mat y = X.rows(find(T == 1));
    return y;
}

// rcpp in; arma out
arma::mat submatrix(NumericMatrix X, NumericVector T) {
    mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
    colvec tIdx(T.begin(), T.size(), false); 
    mat y = Xmat.rows(find(tIdx == 1));
    return y;
}

// rcpp in; rcpp out
NumericMatrix sumbatC3(NumericMatrix X, LogicalVector condition) { 
    int n=X.nrow(), k=X.ncol();
    NumericMatrix out(sum(condition),k);
    for (int i = 0, j = 0; i < n; i++) {
        if(condition[i]) {
            out(j,_) = X(i,_);
            j = j+1;
        }
    }
    return(out);
}


NumericVector eigenval(arma::mat M, Function eigen) {
    List output=eigen(M);
    
    return output(0);
}

double gvar_C(arma::mat M){
	arma::mat temp=arma::cov(M);
	arma::vec values=arma::eig_sym(temp);
	double output=arma::as_scalar(prod(values));
	return output;
}

double sum(NumericVector v){
	int n=v.size();
	double Sum=0.0;
	for(int i=0;i<n;i++)
		Sum+=v(i);
		
	return Sum;
}

// [[Rcpp::export]]
//NumericVector mgvar_while(NumericVector flag, NumericMatrix m){
RcppExport SEXP mgvar_while(SEXP X, SEXP M){
	NumericVector flag(X);
	NumericMatrix m(M);
 	int nr=m.nrow();
 	NumericVector f=flag;
 	arma::mat m_sub;
 	NumericVector varvec(nr, 0.0);
	while(sum(f)>0.0){
		int ic=0;
	 	NumericVector chk;
	 	IntegerVector remi;
		for(int i=0; i<nr; i++){
 			if(f(i)==1.0) {
 				ic+=1;
	 			NumericVector f2=ifelse(f==1.0, 0.0, 1.0);
 				m_sub=submatrix(m, f2);
 				arma::rowvec m_i=m(i,_);
 				m_sub.insert_rows(m_sub.n_rows, m_i);
	 			chk.push_back(gvar_C(m_sub));
 				remi.push_back(i);
 			}
		}
		IntegerVector sor=order(chk);
 		int k=remi(sor(0));
		varvec(k)=chk(sor(0));
		f(k)=0.0;
		R_CheckUserInterrupt();
	}
	return varvec;
}

IntegerVector order(NumericVector x){
	int n=x.size();
	std::vector<int> output(n);
	std::vector<val_order> index(n);
	for(int i=0;i<x.size();i++){
		index[i].value=x(i);
		index[i].order=i;
	}
	std::sort(index.begin(), index.end(), compare); 
	for(int i=0;i<x.size();i++){
		output[i]=index[i].order;
	}
	return wrap(output);
}

double median(NumericVector x) {
    double output;
  	int n=x.size();
  	std::vector<double> xcopy(n);
    std::copy(x.begin(), x.end(), xcopy.begin());
    std::sort(xcopy.begin(), xcopy.end());
    if(n%2==0) {
	   output=((xcopy[n/2] + xcopy[n/2-1])/2);
    } else {
      output=xcopy[n/2];
   }
    return output;
}

double mad(NumericVector x, double constant=1.4826){
	NumericVector xcopy=Rcpp::clone(x);
	int n=xcopy.size();
	double center=median(xcopy);
	NumericVector diff=xcopy-center;
	NumericVector diff_abs=ifelse(diff<0, -diff, diff);
	return median(diff_abs)*constant;
}

NumericVector outer_pos(NumericVector x, NumericVector y){
	std::vector<double> output;
	int n=x.size();
	double temp;
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			temp=x[j]-x[i];
			if(temp>0){
				output.push_back((y[j]-y[i])/temp);
			}
		}
	} 
	return Rcpp::wrap(output);
}

double hd_C(NumericVector x, NumericVector q){
	int n=x.size();
	double m1=(n+1.0)*q(0), m2=(n+1.0)*(1.0-q(0));
	double output=0.0;
	IntegerVector i=seq_len(n);
	NumericVector i2(n);
	i2=i*1.0;
	NumericVector w=pbeta(i2*1.0/n, m1, m2) - pbeta((i2-1.0)/n, m1, m2);
	NumericVector xcopy=Rcpp::clone<Rcpp::NumericVector>(x);
	std::sort(xcopy.begin(), xcopy.end());
	for(int j=0; j<n; j++){
		output+=(xcopy(j)*w(j));
	}
	//Rprintf("%lf", output);
	return(output);
}

// [[Rcpp::export]]
RcppExport SEXP tshd_C(SEXP X, SEXP Y, SEXP hd){
	NumericVector x(X), y(Y);
	IntegerVector HD(hd);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	if(sum(HD)==0){
		coef(1)=median(v1v2);
	}else{
		coef(1)=hd_C(v1v2,wrap(0.5));
	}
	NumericVector res=y-coef(1)*x;
	coef(0)=hd_C(res,wrap(0.5));
	return coef;
}

NumericVector tshd_C(NumericVector x, NumericVector y, IntegerVector HD){
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	if(sum(HD)==0){
		coef(1)=median(v1v2);
	}else{
		coef(1)=hd_C(v1v2,wrap(0.5));
	}
	NumericVector res=y-coef(1)*x;
	coef(0)=hd_C(res,wrap(0.5));
	return coef;
}


// [[Rcpp::export]]
RcppExport SEXP pbvar_C(SEXP X, SEXP BETA){
	NumericVector x(X);
	NumericVector beta(BETA);
	int n=x.size();
	double Median=median(x);
	NumericVector w=ifelse((x-Median)<=0, -(x-Median), (x-Median));
	double pbvar;
	std::sort(w.begin(), w.end());
	double omega=w[floor((1.0-beta(0))*(n+0.5))-1];
	if(omega>0.0){
		double len=0;
		NumericVector z(n);
		for(int i=0;i<n;i++){
			if((x(i)-Median)/omega<= -1.0)
				z(i)=-1.0;
			if ((x(i)-Median)/omega>=1.0)
				z(i)=1.0;
			if(fabs((x(i)-Median)/omega) < 1.0){
				z(i)=(x(i)-Median)/omega;
				len+=1.0;
			}
		}
		double z_sq_sum=0;
		for(int i=0;i<n;i++){
			z_sq_sum+=pow(z(i), 2.0);
		}
		pbvar=(double)n*pow(omega, 2)*z_sq_sum/pow(len, 2);
		return wrap(pbvar);
	} else
		return wrap(0.0);
}

// [[Rcpp::export]]
RcppExport SEXP tsp1reg_C(SEXP X, SEXP Y, SEXP hd){
	NumericVector x(X), y(Y);
	IntegerVector HD(hd);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	coef(1)=median(v1v2);
	if(sum(HD)==0){
		coef(0)=median(y)-coef(1)*median(x);
	}else{
		NumericVector temp1(1), temp2(1);
		temp1(0)=hd_C(y,wrap(0.5));
		temp2(0)=hd_C(x,wrap(0.5));
		coef(0)=temp1(0)-coef(1)*temp2(0);
	}
	NumericVector res=y-coef(1)*x-coef(0);
	return List::create(_["coef"]=coef, _["res"]=res);
}


// [[Rcpp::export]]

RcppExport SEXP stsregp1_C(SEXP X, SEXP Y){
	NumericVector x(X), y(Y);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector allvar(n_s);
	for(int i=0; i<n_s;i++){
		NumericVector temp(pbvar_C(wrap(y-v1v2(i)*x), wrap(0.2)));
		allvar(i)=temp(0);
		R_CheckUserInterrupt();
	}
	IntegerVector b1_id=order(allvar);
	NumericVector coef(2);
	coef(1)=v1v2[b1_id(0)];
	coef(0)=median(y)-coef[1]*median(x);
	NumericVector res=y-coef(1)*x-coef(0);
	return List::create(_["coef"]=coef, _["res"]=res);
}


NumericVector stsregp1_coef(SEXP X, SEXP Y){
	NumericVector x(X), y(Y);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector allvar(n_s);
	for(int i=0; i<n_s;i++){
		NumericVector temp(pbvar_C(wrap(y-v1v2(i)*x), wrap(0.2)));
		allvar(i)=temp(0);
		R_CheckUserInterrupt();
	}
	IntegerVector b1_id=order(allvar);
	NumericVector coef(2);
	coef(1)=v1v2[b1_id(0)];
	coef(0)=median(y)-coef[1]*median(x);
	return coef;
}


NumericVector tsp1reg_C(NumericVector x, NumericVector y, IntegerVector HD){
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	coef(1)=median(v1v2);
	if(sum(HD)==0){
		coef(0)=median(y)-coef(1)*median(x);
	}else{
		NumericVector temp1(1), temp2(1);
		temp1(0)=hd_C(y,wrap(0.5));
		temp2(0)=hd_C(x,wrap(0.5));
		coef(0)=temp1(0)-coef(1)*temp2(0);
	}
	NumericVector res=y-coef(1)*x-coef(0);
	return coef;
}

RcppExport SEXP stsreg_for(SEXP X, SEXP Y, SEXP IT){
	NumericMatrix x(X);
	NumericVector y(Y);
	IntegerVector it(IT);
	int ncols=x.ncol();
	int nrows=x.nrow();
	NumericVector temp(ncols);
	for(int i=0; i<ncols;i++){
		NumericVector tempcoef=tsp1reg_C(x(_,i), y, wrap(1));
		temp(i)=tempcoef(1);
	}
	arma::colvec res=(Rcpp::as<arma::colvec>(y))
	              -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp));
	double alpha=median((Rcpp::as<NumericVector>(wrap(res))));
//	NumericVector tempold(temp);
	arma::colvec r(nrows);
	for(int i=0; i<it(0); i++){
		for(int j=0; j<ncols; j++){
			NumericVector tempcol=x(_,j);
			r=(Rcpp::as<arma::colvec>(y))
			  -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
		      -alpha + temp(j)*(Rcpp::as<arma::colvec>(tempcol));
		    temp(j)=stsregp1_coef(wrap(tempcol), wrap(r))(1);
		}
		alpha=median(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))));
		R_CheckUserInterrupt();
//		NumericVector tempold(temp);
	}
	res=(Rcpp::as<arma::colvec>(y))
	    -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
	    -alpha;
	return List::create(_["alpha"]=alpha, _["beta"]=temp, _["res"]=res);
}


RcppExport SEXP tshdreg_for(SEXP X, SEXP Y, SEXP IT, SEXP TOL){
	NumericMatrix x(X);
	NumericVector y(Y);
	NumericVector tol(TOL);
	IntegerVector it(IT);
	int ncols=x.ncol();
	int nrows=x.nrow();
	IntegerVector HD(1);
	HD(0)=1;
	NumericVector temp(ncols), tempold(ncols);
	for(int i=0; i<ncols;i++){
		temp(i)=tshd_C(x(_,i), y, HD)(1);
		tempold(i)=temp(i);
	}
	arma::colvec res=(Rcpp::as<arma::colvec>(y))
	                 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp));
	double alpha=hd_C((Rcpp::as<NumericVector>(wrap(res))), wrap(0.5));
	std::vector<double> diff(ncols);
	arma::colvec r(nrows);
	for(int i=0; i<it(0); i++){
		for(int j=0; j<ncols; j++){
			NumericVector tempcol=x(_,j);
			r=(Rcpp::as<arma::colvec>(y))
			  -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
		      -alpha + temp(j)*(Rcpp::as<arma::colvec>(tempcol));
		    temp(j)=tshd_C(tempcol, Rcpp::as<Rcpp::NumericVector>(wrap(r)), HD)(1);
			diff[j]=temp(j)-tempold(j);
			if(diff[j]<0.0) diff[j]=-1.0*diff[j];
			tempold(j)=temp(j);
		}
		R_CheckUserInterrupt();
		std::sort(diff.begin(), diff.end());
		if(diff[diff.size()-1]<tol(0)) break;
		alpha=hd_C(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))), wrap(0.5));
		NumericVector tempold(clone(temp));
	}
	res=(Rcpp::as<arma::colvec>(y))
	    -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
	    -alpha;
	return List::create(_["alpha"]=alpha, _["beta"]=temp, _["res"]=res);
}




RcppExport SEXP tsreg_for(SEXP X, SEXP Y, SEXP IT, SEXP hd){
	NumericMatrix x(X);
	NumericVector y(Y);
	IntegerVector it(IT);
	int ncols=x.ncol();
	int nrows=x.nrow();
	IntegerVector HD(hd);
	NumericVector temp(ncols), tempold(ncols);
	for(int i=0; i<ncols;i++){
		temp(i)=tsp1reg_C(x(_,i), y, HD)(1);
		tempold(i)=temp(i);
	}
	arma::colvec res=(Rcpp::as<arma::colvec>(y))
	                 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp));
	double alpha;
	if(sum(HD)==0){
		alpha=median(res);
	} else {
		alpha=hd_C((Rcpp::as<NumericVector>(wrap(res))), wrap(0.5));
	}
	std::vector<double> diff(ncols);
	arma::colvec r(nrows);
	for(int i=0; i<it(0); i++){
		for(int j=0; j<ncols; j++){
			NumericVector tempcol=x(_, j);
			r=(Rcpp::as<arma::colvec>(y))
			  -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
		      -alpha + temp(j)*(Rcpp::as<arma::colvec>(tempcol));
		    temp(j)=tsp1reg_C(x(_,j), Rcpp::as<Rcpp::NumericVector>(wrap(r)), HD)(1);
			tempold(j)=temp(j);
		}
		R_CheckUserInterrupt();
		if(sum(HD)==0){
			alpha=median(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))));
		} else {
			alpha=hd_C(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))), wrap(0.5));
		}
	}
	res=(Rcpp::as<arma::colvec>(y))
	    -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
	    -alpha;
	return List::create(_["alpha"]=alpha, _["beta"]=temp, _["res"]=res);
}




double sign(double x)
{
	 double output;
	 if (std::isnan(x)) {
	    	output=NA_REAL;
	} else {
    		output=(((x > 0) ? 1 : ((x == 0)? 0 : -1)));
    	}
    return output;
}


NumericVector unidepth1(NumericVector x, NumericVector y){
	int nx=x.size();
	int ny=y.size();
	NumericVector output(ny);
	double pup, pdown;
	for(int i=0; i<ny; i++){
		double temp1=0, temp2=0;
		for(int j=0; j<nx; j++){
			if(y(i)<=x(j)) temp1+=1;
			if(y(i)<x(j)) temp2+=1;
		}
		pup=temp1/((double)nx);
		pdown=1-temp2/((double)nx);
		if(pup<pdown)
		{	output(i)=pup;  }
		else output(i)=pdown;
	} 
	return output;
}

NumericVector subsetVec(NumericVector A, int start, int end) {
  NumericVector B(end-start+1) ;
  std::copy(A.begin() + start, A.begin() + end+1, B.begin() ) ;
  return B;
}


RcppExport SEXP fdepthv2_for(SEXP M, SEXP PTS){
	NumericMatrix pts(PTS);
	NumericMatrix m(M);
	int nrows=m.nrow();
	//IntegerVector mdep_dim(MDEP_DIM);
	//int mdep_nr=mdep_dim(1), mdep_nc=mdep_dim(0); 
	int mdep_nc=(nrows*nrows-nrows)/2, mdep_nr=nrows;
	if(!R_IsNA(pts(0,0))) mdep_nr=pts.nrow();
	NumericMatrix mdep(mdep_nr, mdep_nc);
	NumericVector dis_temp;
	int ic=0;
	for(int iall=0; iall<nrows; iall++){
		for(int i=0; i<nrows; i++){
			R_CheckUserInterrupt();
			if(iall<i){
				ic+=1;
				NumericVector B=m(i,_)-m(iall,_);
				NumericVector dis;
				NumericVector BB=B*B;
				double bot=sum(BB);
				if(bot!=0.0){
					if(R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							NumericVector A=m(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
							//dis_temp.push_back(sum(temp));
						}
					}
					if(!R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							NumericVector A=m(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
						for(int j=0; j<pts.nrow(); j++){
							NumericVector A=pts(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
					}
					//
					// For ic_th projection, store depths of
					// points in mdep[ic,]
					//
					if(R_IsNA(pts(0,0))){
						NumericVector unidepth_out=unidepth1(dis, dis);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					} else if(!R_IsNA(pts(0,0))){
						//Rprintf("%lf ", sum(dis));
						NumericVector dis_1=subsetVec(dis, 0, nrows-1);
						NumericVector dis_2=subsetVec(dis, nrows, nrows+pts.nrow()-1);
						//Rprintf("%d \n", 2*nrows-1);
						NumericVector unidepth_out=unidepth1(dis_1, dis_2);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					}
				}
				if(bot==0.0) {
					for(int k=0; k<mdep_nr;k++)
						mdep(k,ic-1)=0.0;					
				}
			}
		}
	}
	return(mdep);
}
	
RcppExport SEXP fdepthv2_for2(SEXP M, SEXP PTS){
	NumericMatrix pts(PTS);
	NumericMatrix m(M);
	int nrows=m.nrow();
	//IntegerVector mdep_dim(MDEP_DIM);
	//int mdep_nr=mdep_dim(1), mdep_nc=mdep_dim(0); 
	int mdep_nc=(nrows*nrows-nrows)/2, mdep_nr=nrows;
	if(!R_IsNA(pts(0,0))) mdep_nr=pts.nrow();
	NumericMatrix mdep(mdep_nr, mdep_nc);
	NumericVector dis_temp;
	int ic=0;
	for(int iall=0; iall<nrows; iall++){
		for(int i=0; i<nrows; i++){
			R_CheckUserInterrupt();
			if(iall<i){
				ic+=1;
				NumericVector B=m(i,_)-m(iall,_);
				NumericVector dis;
				NumericVector BB=B*B;
				double bot=sum(BB);
				if(bot!=0.0){
					if(R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							NumericVector A=m(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
							//dis_temp.push_back(sum(temp));
						}
					}
					if(!R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							NumericVector A=m(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
						for(int j=0; j<pts.nrow(); j++){
							NumericVector A=pts(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
					}
					//
					// For ic_th projection, store depths of
					// points in mdep[ic,]
					//
					if(R_IsNA(pts(0,0))){
						NumericVector unidepth_out=unidepth1(dis, dis);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					} else if(!R_IsNA(pts(0,0))){
						//Rprintf("%lf ", sum(dis));
						NumericVector dis_1=subsetVec(dis, 0, nrows-1);
						NumericVector dis_2=subsetVec(dis, nrows, nrows+pts.nrow()-1);
						//Rprintf("%d \n", 2*nrows-1);
						NumericVector unidepth_out=unidepth1(dis_1, dis_2);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					}
				}
				if(bot==0.0) {
					for(int k=0; k<mdep_nr;k++)
						mdep(k,ic-1)=0.0;					
				}
			}
		}
	}
	return(mdep);
}

NumericVector idealf(NumericVector x){
	NumericVector xcopy=Rcpp::clone(x), output(2);
	int n=xcopy.size();
	int j=floor((double)n/4.0 + 5.0/12.0)-1;
	std::sort(xcopy.begin(), xcopy.end());
	double g=n/4.0-j-1+5.0/12.0;
	int k=n-j-1;
	//Rprintf("%lf  %d  %d", g, j, k);
	output(0)=(1.0-g)*xcopy(j)+g*xcopy(j+1);
	output(1)=(1.0-g)*xcopy(k)+g*xcopy(k-1);
	return output;
}

/*
library(RcppArmadillo)
source("Rallfun-v22.R")
source("/Users/praguewatermelon/Dropbox/robustmethod_CPP.R")
dyn.load("robustmethods_CPP.so")
set.seed(1)
x=rmul(60, 4)
rmdzeroG(m, est=skip_C)
*/

RcppExport SEXP outpro_for(SEXP M, SEXP GVAL, SEXP CENTER, SEXP MM){
	NumericMatrix m(M);
	int nrows=m.nrow();
	NumericVector gval(GVAL), center(CENTER);
	LogicalVector mm(MM);
	NumericVector outid; 
	LogicalVector flag(nrows, false);

	for(int i=0; i<nrows; i++){
		NumericVector B=m(i,_)-center, BB=B*B, dis(nrows);
		double bot=sum(BB); 
		if(bot!=0){
			for(int j=0;j<nrows;j++){
				NumericVector A=m(j,_)-center, temp=sum(A*B)*B/bot;
				dis(j)=pow(sum(temp*temp), 0.5);
			}
			NumericVector temp=idealf(dis);
			double cu;
			if(!mm(0))
				cu=median(dis)+gval(0)*(temp(1)-temp(0));
			else if(mm(0))
				cu=median(dis)+gval(0)*mad(dis);
			//Rprintf("%lf\n", sum(dis));
			for(int k=0; k<nrows; k++){
				if(dis(k)>cu)
					flag(k)=true;
			}
		}
	}
	return flag;
}


/*

IntegerVector outpro_for(NumericMatrix m, NumericVector gval, NumericVector center, 
				LogicalVector mm){
	int nrows=m.nrow();
	IntegerVector outid; 
	IntegerVector flag(nrows, 0);
	for(int i=0; i<nrows; i++){
		NumericVector B=m(i,_)-center, BB=B*B, dis(nrows);
		double bot=sum(BB); 
		if(bot!=0){
			for(int j=0;j<nrows;j++){
				NumericVector A=m(j,_)-center, temp=sum(A*B)*B/bot;
				dis(j)=pow(sum(temp*temp), 0.5);
			}
			NumericVector temp=idealf(dis);
			double cu;
			if(!mm(0))
				cu=median(dis)+gval(0)*(temp(1)-temp(0));
			else if(mm(0))
				cu=median(dis)+gval(0)*mad(dis);
			//Rprintf("%lf\n", sum(dis));
			for(int k=0; k<nrows; k++){
				if(dis(k)>cu)
					flag(k)=1;
			}
		}
	}
	if(sum(flag)>0){
		for(int i=0;i<nrows; i++){
			if(flag(i)==1)
				outid.push_back(i+1);
		}
	}
	return outid;
}

NumericVector outpro_keep(NumericMatrix m, NumericVector gval, NumericVector center, 
				LogicalVector mm){
	int nrows=m.nrow();
	IntegerVector keepid; 
	NumericVector flag(nrows, 1.0);
	for(int i=0; i<nrows; i++){
		NumericVector B=m(i,_)-center, BB=B*B, dis(nrows);
		double bot=sum(BB); 
		if(bot!=0){
			for(int j=0;j<nrows;j++){
				NumericVector A=m(j,_)-center, temp=sum(A*B)*B/bot;
				dis(j)=pow(sum(temp*temp), 0.5);
			}
			NumericVector temp=idealf(dis);
			double cu;
			if(!mm(0))
				cu=median(dis)+gval(0)*(temp(1)-temp(0));
			else if(mm(0))
				cu=median(dis)+gval(0)*mad(dis);
			//Rprintf("%lf\n", sum(dis));
			for(int k=0; k<nrows; k++){
				if(dis(k)>cu)
					flag(k)=0.0;
			}
		}
	}
	return flag;
}



RcppExport SEXP outpro_loop(SEXP M, SEXP GVAL, SEXP CENTER, SEXP MM, SEXP NBOOT, 
							SEXP N, SEXP NCOL){
	List mlist(M), centerlist(CENTER);
	NumericVector gval(GVAL);
	LogicalVector mm(MM);
	IntegerVector nboot(NBOOT), n(N), tempncols(NCOL);
	int temp_ncols=tempncols(0);
	NumericMatrix colmeans(nboot(0), temp_ncols);

	for(int i=0; i<nboot(0); i++){
		NumericMatrix tempmat=mlist(i);
		IntegerVector temp=outpro_keep(tempmat, gval, centerlist(i), mm);
		int temp_n=temp.size();
		for(int j=0; j<temp_ncols; j++){
			NumericVector tempcol(temp_n);
			for(int k=0; k<temp_n; k++){
					tempcol(k)=tempmat(temp(k), j);
			}
			colmeans(i, j)=sum(tempcol)/temp_n;
		}
	}
	return colmeans;
}



RcppExport SEXP outpro_loop2(SEXP M, SEXP GVAL, SEXP CENTER, SEXP MM, SEXP NBOOT, 
							SEXP N, SEXP NCOL){
	List mlist(M), centerlist(CENTER);
	NumericVector gval(GVAL);
	LogicalVector mm(MM);
	IntegerVector nboot(NBOOT), n(N), tempncols(NCOL);
	int temp_ncols=tempncols(0);
	NumericMatrix colmeans(nboot(0), temp_ncols);
	arma::mat m_sub;
	for(int i=0; i<nboot(0); i++){
		NumericMatrix tempmat=mlist(i);
		NumericVector temp=outpro_keep(tempmat, gval, centerlist(i), mm);
		m_sub=submatrix(tempmat, temp);
		//Rprintf("%d %d, ", m_sub.n_cols, m_sub.n_rows);
		for(int j=0; j<temp_ncols; j++){
			colmeans(i, j)=arma::mean(m_sub.col(j));
		}
	}
	return colmeans;
}


rmdzeroG_vers2<-function(x,est=skip,grp=NA,nboot=500,SEED=TRUE,...){
#
#   Do ANOVA on dependent groups
#   using #   depth of zero among  bootstrap values
#   based on difference scores.
#
#   The data are assumed to be stored in x in list mode
#   or in a matrix. In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, columns correspond to groups.
#
#   grp is used to specify some subset of the groups, if desired.
#   By default, all J groups are used.
#
#   The default number of bootstrap samples is nboot=500
#
	argslist<-formals()
	#require("parallel")
	if(!is.list(x) && !is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
	if(is.list(x)){
		# put the data in an n by J matrix
		mat<-matrix(0,length(x[[1]]),length(x))
		for (j in 1:length(x))mat[,j]<-x[[j]]
	}
	if(is.matrix(x))mat<-x
	if(!is.na(grp[1])){
		mat<-mat[,grp]
	}
	mat<-elimna(mat) # Remove rows with missing values.
	J<-ncol(mat)
	jp<-0
	Jall<-(J^2-J)/2
	dif<-matrix(NA,nrow=nrow(mat),ncol=Jall)
	ic<-0
	#browser()
	for(j in 1:J){
		for(k in 1:J){
			if(j<k){
				ic<-ic+1
				dif[,ic]<-mat[,j]-mat[,k]; 
			}
		}
	}
	dif<-as.matrix(dif)
	if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
	data <- as.data.frame(t(matrix(sample(nrow(mat), size = nrow(mat) * nboot, replace = T),
                nrow = nboot)))
	if(identical(est, skip_C)){
		if(!("MM" %in% names(argslist))) MM=FALSE
		else if(!is.logical(argslist$MM)) stop("MM needs to be a boolean")
		output<-mapply(function(index) outpro_preloop(dif[index,],...), data, SIMPLIFY=F)
		m<-lapply(output, "[[", "m")
		center<-lapply(output, "[[", "center")
		gval<-outpro_gval(dif,...)
		bvec<-.Call("outpro_loop2", M=m, GVAL=gval, CENTER=center, MM=MM, NBOOT=nboot, 
					N=nrow(mat), NCOL=Jall)
		#bvec<-as.data.frame(t(matrix(bvec, nrow=nboot, byrow=TRUE)))
	} else {
		bvec<-mapply(function(index) est(dif[index,],...)$center, data)
		bvec<-t(bvec)
		#browser()
	}
	center<-est(dif,...)$center
	bcen<-colMeans(bvec)
	cmat<-var(bvec-bcen+center)
	zvec<-rep(0,Jall)
	m1<-rbind(bvec,zvec)
	bplus<-nboot+1
	discen<-mahalanobis(m1,center,cmat)
	sig.level<-sum(discen[bplus]<=discen)/bplus
	list(p.value=sig.level,center=center)
}

#data<-matrix(sample(60, 10*60, T), nrow=10)

#output<-mapply(function(index) outpro_preloop(m[index,]), as.data.frame(t(data)), SIMPLIFY=FALSE)
#m<-lapply(output, "[[", "m")
#center<-lapply(output, "[[", "center")
#.Call("outpro_loop", M=m, GVAL=3.057516, CENTER=center, MM=FALSE, NBOOT=10)

outpro_preloop<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
xlab="VAR 1",ylab="VAR 2",STAND=FALSE,tr=.2,q=.5,pr=TRUE,...){
	#
# Detect outliers using a modification of the
# Stahel-Donoho  projection method.
#
# Determine center of data cloud, for each point,
# connect it with center, project points onto this line
# and use distances between projected points to detect
# outliers. A boxplot method is used on the
# projected distances.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data.
#
# op=T
# means the .5 depth contour is plotted
# based on data with outliers removed.
#
# op=F
# means .5 depth contour is plotted without removing outliers.
#
#  MM=F  Use interquatile range when checking for outliers
#  MM=T  uses MAD.
#
#  If value for center is not specified,
#  there are four options for computing the center of the
#  cloud of points when computing projections:
#
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)#  cop=7 uses the spatial (L1) median
#
#  args q and tr having are not used by this function. They are included to deal
#  with situations where smoothers have optional arguments for q and tr
#
#  When using cop=2, 3 or 4, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  STAND=T means that marginal distributions are standardized before
#  checking for outliers.
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
	m<-as.matrix(m)
	if(pr){
		if(!STAND){
			#if(ncol(m)>1)cat("STAND=FALSE. If measures are on different scales,", 
			#					"might want to use STAND=TRUE\n")
		}
	}
	library(MASS)
	m=elimna(m)
	m<-as.matrix(m)
	nv=nrow(m)
	if(ncol(m)==1){
		dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
		dis<-sqrt(dis)
		dis[is.na(dis)]=0
		crit<-sqrt(qchisq(.975,1))
		chk<-ifelse(dis>crit,1,0)
		vec<-c(1:nrow(m))
		outid<-vec[chk==1]
		keep<-vec[chk==0]
	}
	if(ncol(m)>1){
		if(STAND)m=standm(m,est=median,scat=mad)
		if(cop==1 && is.na(center[1])){
			if(ncol(m)>2)center<-dmean(m,tr=.5,cop=1)
			if(ncol(m)==2){
				tempd<-NA
				for(i in 1:nrow(m))
				tempd[i]<-depth(m[i,1],m[i,2],m)
				mdep<-max(tempd)
				flag<-(tempd==mdep)
				if(sum(flag)==1)center<-m[flag,]
				if(sum(flag)>1)center<-apply(m[flag,],2,mean)
			}
		}
		if(cop==2 && is.na(center[1])){
			center<-cov.mcd(m)$center
		}
		if(cop==4 && is.na(center[1])){
			center<-cov.mve(m)$center
		}
		if(cop==3 && is.na(center[1])){
			center<-apply(m,2,median)
		}
		if(cop==5 && is.na(center[1])){
			center<-tbs(m)$center
		}
		if(cop==6 && is.na(center[1])){
			center<-rmba(m)$center
		}
		if(cop==7 && is.na(center[1])){
			center<-spat(m)
		}
	}
	list(m=m, center=center)
}

outpro_gval<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
xlab="VAR 1",ylab="VAR 2",STAND=FALSE,tr=.2,q=.5,pr=TRUE,...){
	#
# Detect outliers using a modification of the
# Stahel-Donoho  projection method.
#
# Determine center of data cloud, for each point,
# connect it with center, project points onto this line
# and use distances between projected points to detect
# outliers. A boxplot method is used on the
# projected distances.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data.
#
# op=T
# means the .5 depth contour is plotted
# based on data with outliers removed.
#
# op=F
# means .5 depth contour is plotted without removing outliers.
#
#  MM=F  Use interquatile range when checking for outliers
#  MM=T  uses MAD.
#
#  If value for center is not specified,
#  there are four options for computing the center of the
#  cloud of points when computing projections:
#
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)#  cop=7 uses the spatial (L1) median
#
#  args q and tr having are not used by this function. They are included to deal
#  with situations where smoothers have optional arguments for q and tr
#
#  When using cop=2, 3 or 4, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  STAND=T means that marginal distributions are standardized before
#  checking for outliers.
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
	m<-as.matrix(m)
	if(pr){
		if(!STAND){
			#if(ncol(m)>1)cat("STAND=FALSE. If measures are on different scales,", 
			#					"might want to use STAND=TRUE\n")
		}
	}
	library(MASS)
	m=elimna(m)
	m<-as.matrix(m)
	nv=nrow(m)
	if(ncol(m)==1){
		dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
		dis<-sqrt(dis)
		dis[is.na(dis)]=0
		crit<-sqrt(qchisq(.975,1))
		chk<-ifelse(dis>crit,1,0)
		vec<-c(1:nrow(m))
		outid<-vec[chk==1]
		keep<-vec[chk==0]
	}
	if(ncol(m)>1){
		if(STAND)m=standm(m,est=median,scat=mad)
		if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
		if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
	}
	gval
}
	



rmdzeroG<-function(x,est=skip,grp=NA,nboot=500,SEED=TRUE,...){
#
#   Do ANOVA on dependent groups
#   using #   depth of zero among  bootstrap values
#   based on difference scores.
#
#   The data are assumed to be stored in x in list mode
#   or in a matrix. In the first case
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J.
#   If stored in a matrix, columns correspond to groups.
#
#   grp is used to specify some subset of the groups, if desired.
#   By default, all J groups are used.
#
#   The default number of bootstrap samples is nboot=500
#
if(!is.list(x) && !is.matrix(x))stop('Data must be stored in a matrix or in list mode.')
if(is.list(x)){
# put the data in an n by J matrix
mat<-matrix(0,length(x[[1]]),length(x))
for (j in 1:length(x))mat[,j]<-x[[j]]
}
if(is.matrix(x))mat<-x
if(!is.na(grp[1])){
mat<-mat[,grp]
}
mat<-elimna(mat) # Remove rows with missing values.
J<-ncol(mat)
jp<-0
Jall<-(J^2-J)/2
dif<-matrix(NA,nrow=nrow(mat),ncol=Jall)
ic<-0
for(j in 1:J){
for(k in 1:J){
if(j<k){
ic<-ic+1
dif[,ic]<-mat[,j]-mat[,k]
}}}
dif<-as.matrix(dif)
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
data <- matrix(sample(nrow(mat), size = nrow(mat) * nboot, replace = T),
                nrow = nboot)
bvec <- matrix(NA, ncol = ncol(dif), nrow = nboot)
        for(i in 1:nboot) {
                bvec[i, ] <- est(dif[data[i,],],...)$center
        }  #bvec is an nboot by Jm matrix
center<-est(dif,...)$center
bcen<-apply(bvec,2,mean)
cmat<-var(bvec-bcen+center)
zvec<-rep(0,Jall)
m1<-rbind(bvec,zvec)
bplus<-nboot+1
discen<-mahalanobis(m1,center,cmat)
sig.level<-sum(discen[bplus]<=discen)/bplus
list(p.value=sig.level,center=center)
}
*/