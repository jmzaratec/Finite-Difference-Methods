/*
* Program structure created by Paul Johnson. May be edited and 
* handed in as part of the coursework. Leave this here if this 
* is the case.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


#include <algorithm>
#include <cstdlib>
#include <ctime>
using namespace std;


class BlackScholes {

public:
	double normal_distribution(double x){
	// declare the number of steps
	int n;
	// declare Pi, h, minus_infty,normal_dist
	double Pi, h, minus_infty,normal_dist;
	// assign constants
	n=10000;
	Pi = 3.141592653589793;
	minus_infty = -10.;
	h = (x - minus_infty)/2./n;
	// summation and the current value of x
	double sum,x_i;
	// add the two ends y_0 and y_2n
	sum = exp(-0.5*minus_infty*minus_infty) + 
		exp(-0.5*x*x);
	// find x_i at position 1
	x_i = minus_infty + h;
	// sum = sum + 4y_1
	sum = sum + 4.*exp(-0.5*x_i*x_i);
	// now loop over i=1 to n-1 to add on all other terms
	for(int i=1;i<=n-1;i=i+1)
	{
	  // add on terms y_2, y_4, y_6 ,..., y_{n-2}
	  x_i = minus_infty + 2.*i*h;
	  sum = sum + 2.*exp(-0.5*x_i*x_i);
	  // add on terms y_3, y_5, ..., y_{n-1}
	  x_i = x_i + h;	
	  sum = sum + 4.*exp(-0.5*x_i*x_i);
	}
return 1/sqrt(2.*Pi)*h/3.*sum; }

// A function to compute the price of a European Call Option using BlackSholes formula
double call(double S0,double T, double t, double r, double q, double sigma, double X){
	
	
	double d1,d2, call;
	call=0;
	if(S0==0.){	call= max(S0-X,0.);
	return call;}
	else if (t>=T){call= max(S0-X,0.);
	return call;}
	else {
	//Black Scholes Formula
		d1 = (log(S0/X) +(r-q)*(T-t) + (sigma*sigma*0.5)*(T-t))/( sigma*sqrt(T-t));
		d2 = d1-(sigma*sqrt(T-t));
	//cout<<d1<<" "<<d2<<"\n";
		call  = S0*exp((-q)*(T-t))*normal_distribution(d1)-X*exp((-r)*(T-t))*normal_distribution(d2);
	//To avoid error in case t=>T
	return call;}
}





};


int main(){
	//Declaration 
	BlackScholes o1;
	double r=0.069;
	double T=1.4;
	double X=1600;
	double d1=0.045;
	double S0=X;
	float secsElapsed;
	//Sce1
	//*
	double Beta=1.;
	double k=0;
	double sigma=0.32;//*/
	//Sce2
	/*
	double Beta=0.28;
	double k=0;
	double sigma=63;//*/
	//Sce3
	/*
	double Beta=0.28;
	double k=0.39;
	double sigma=63;//*/

	// output
	std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/12_finiteDifferenceProject/output.csv");
	 output. precision(15);
	if(!output.is_open()) { std::cout << " File not opened \n";    
	 throw;  }

	output<<"FINETE dt  dS, Stock, V(S-t=0), Analytical\n";

	for (int jMax= 350; jMax<=351;jMax=jMax+1){
	cout<<jMax<<"\n";
	clock_t startTime = clock();
	
	int iMax=20000;
	double SMax=S0*2 ;
	double dS= SMax/jMax;
	vector<double> stock (jMax+1);

	// set Stock
	for(int j=0; j<=jMax;j++){	stock[j]=j*dS;	}
	double dt= T/iMax;

	//Make a vector to access option values at any time
	vector<vector<double> > V(iMax+1,vector<double>(jMax+1));

	// Boundaries
	for(int j=0; j<=jMax;j++){	V[iMax][j]=max(stock[j]-X,0.);}

	// Boundaries
	for(int i=iMax-1;i>=0; i--)	{			
	V[0][0]=0.; 	
	V[i][jMax]=stock[jMax]*exp(-d1*(T-i*dt)) - X*exp(-r*(T-i*dt));
	}
	// Inside the grid
for(int i=iMax-1;i>=0; i--){
	for (int j=1; j<jMax; j++){
		// Calculate values of A, B, C 
		double A, B, C;
		A=(pow(dS,2*Beta-2)*exp(-k*dt*i)*0.5*sigma*sigma*pow(j,2*Beta) + 0.5*(r-d1)*j)*dt;
		B=1.- pow(dS,2*Beta-2)*exp(-k*dt*i)*sigma*sigma *pow(j,2*Beta)*dt;
		C=(pow(dS,2*Beta-2)*exp(-k*dt*i)*0.5*sigma*sigma*pow(j,2*Beta) - 0.5*(r-d1)*j)*dt;
	  	V[i][j]= (A*V[i+1][j+1]+ B*V[i+1][j] + C*V[i+1][j-1])/(1. + r*dt);
	}
	}
/*// output GRID
//Writing headings
	output<<"stock[j]"<<"\n";
	for(int i=0; i<=iMax;i++){	output<<", "<<dt*i;	}
	output<<"\n";
	for(int j=0; j<=jMax;j++){
	output<<stock[j]<<",";
		for(int i=0; i<=iMax;i++){
			output<<V[i][j]<<",";
	//cout<<"V("<<stock[j]<<","<<i*dt<<")=";
	//cout<<vNew[j]<<"\n";
	}
	output<<"\n";
	}//*/

	//output
for(int j=1; j<=jMax;j++){
	output<<iMax<<"  "<<jMax<<",";	
	if(j==jMax){
	secsElapsed = (float)(clock() - startTime)/CLOCKS_PER_SEC;	
	output<<stock[j]<<","<<V[0][j]<<","<<o1.call(stock[j],T,0.,r,d1,sigma,X)<<","<<secsElapsed<<"\n";
	}
	else{	output<<stock[j]<<","<<V[0][j]<<","<<o1.call(stock[j],T,0.,r,d1,sigma,X)<<"\n";	}
	
	}

	}//end of loop
	  std::cout << " File write successful \n";
	// CLOSE file
  output.close();
	system("pause");
}