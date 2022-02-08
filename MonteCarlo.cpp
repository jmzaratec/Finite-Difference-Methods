
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
using namespace std;


class MonteCarlo {

//public: MonteCarlo(void);

public:
double uniformRandom()
{
  
	//We cannot get 1
  return ((double)(rand())+1)/((double)(RAND_MAX)+1);
}



// a function to return a sample from the normal distribution
double normalRandom()
{

	double a=0;
	double x1=uniformRandom();
	double x2=uniformRandom();

	double pi= atan(1.)* 4;

	a=(cos(2.*pi*x2))*(sqrt(-2.*log(x1)));

	return a;

  
}
double monteCarloCallOption(double S, double T, double t, double r, double q, double sigma, double X, int N){

// terminal stock value
	double S_T;
	//store the running sum
	double sum=0;
	

	for (int i=0; i<N; i++){// loop simulation each payoff
		
		// get noral distributed number
		
		double phi = normalRandom();
		// calculate_S_T
		S_T=S*exp((r-q-0.5*sigma*sigma)*(T-t) + sigma* phi*sqrt((T-t)));
		// add in the payoff sum = sum + 
			sum = sum + max(S_T-X,0.);
	}

	return sum*exp(-r*(T-t))/N;

}

double monteCarloPutOption(double S, double T, double t, double r, double q, double sigma, double X, int N){

// terminal stock value
	double S_T;
	//store the running sum
	double sum=0;
	

	for (int i=0; i<N; i++){// loop simulation each payoff
		
		// get noral distributed number
		
		double phi = normalRandom();
		// calculate_S_T
		S_T=S*exp((r-q-0.5*sigma*sigma)*(T-t) + sigma* phi*sqrt((T-t)));
		// add in the payoff sum = sum + 
			sum = sum + max(-S_T+X,0.);
	}

	return sum*exp(-r*(T-t))/N;

}

double binaryCallOption(double S, double T, double t, double r, double q, double sigma, double X, int N){

	// terminal stock value
	double S_T;
	//store the running sum
	double sum=0;
	double I=0;
	
	for (int i=0; i<N; i++){// loop simulation each payoff
		
		// get noral distributed number
		
		double phi = normalRandom();
		// calculate_S_T
		S_T=S*exp((r-q-0.5*sigma*sigma)*(T-t) + sigma* phi*sqrt((T-t)));
		// add in the payoff sum = sum + 
		if(S_T>X){	I=1;}
		if(S_T<X){	I=0;}
		//if(S_T=X){	I=0.5;}// To be fair
		sum = sum + I;//*X;
		//sum= sum+I;

	}

	return (sum*exp(-r*(T-t)))/N;

}

double AsianCallOKOK(double S, double T, double t, double r, double q, double sigma, double X, int N, int M){

if(M==0.){cout<<"Please insert a value between [1,T]";
	return 0;
}
/*else  if(M==1.){
	return 0;
}*/

else {
double sum=0;
double S_T=0;
	for (int i=0; i<N; i++){
	double sumti=0;
	double S_Ti=S;
		for (int j=1; j<=M; j++){
		double ti=  ((T-t)/M); /// remove t
		//cout<<"ti"<<ti<<" "<<T<<" ";
		S_Ti =S_Ti*exp((r-q-0.5*sigma*sigma)*(ti) + sigma* normalRandom()*sqrt(ti));
		
		//cout<<i<<" "<<S_Ti<<"\n ";
		sumti= sumti+ S_Ti;
		}
		sumti= sumti/M;
		//cout<<"ave: "<<sumti<<" ";
		S_T=S_Ti;
		//cout<<"S_T: "<<S_T<<"\n ";
		sum = sum + max(sumti-S_T,0.);
		//sum = sum + max(sumti-X,0.);
	}
		//cout<<"\n payoff T: "<<sum/N<<" ";
		return (sum*exp(-r*(T-t)))/N;
	}
}

double BarrierPut(double S, double T, double t, double r, double q, double sigma, double X, int N, double B){

if(B<S){cout<<"Please insert a value Higher or equal to S";
	return 0;
}
/*else  if(M==1.){
	return 0;
}*/

else {
double sum=0;
double S_T=0;
double payoff=0;

for (int i=0; i<N; i++){// loop simulation each payoff
		
		// get noral distributed number
		
		double phi = normalRandom();
		// calculate_S_T
		S_T=S*exp((r-q-0.5*sigma*sigma)*(T-t) + sigma* phi*sqrt((T-t)));
		if(S_T>B){	payoff=0.;	}
		else{payoff=max(-S_T+X,0.);}
		
		// add in the payoff sum = sum + 
			sum = sum + payoff;
	}

	return sum*exp(-r*(T-t))/N;
	}
}
};


class DBonds{
public:
	double NB, F, NS, S;
	double maturity,sigma,interestRate, S0, X;
// auxiliry function
double h(double N, double x){
	double h;
	
	h= 0.5+sqrt((0.25-(0.25*exp(-1*pow(x/(N+(0.3333333333333)),2)*(N+(0.16666666666666666))))));
	

	return h;
}
double f(double S,double X, double F)
{
	
	if(S<X)   { return F;}
	else {return 0.;}
}
double qstar(double maturity, double sigma, double r, double S0, double X, int N, double NB, double F){
double d1,d2;
double qstar;
 
d1 = (log(S0/X) +(r)*(maturity) + (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));
d2 = (log(S0/X) +(r)*(maturity) - (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));

qstar = h(N,d1);


return qstar;


}
double qnormal(double maturity, double sigma, double r, double S0, double X, int N, double NB, double F){
double d1,d2;
double qnormal;
 
d1 = (log(S0/X) +(r)*(maturity) + (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));
d2 = (log(S0/X) +(r)*(maturity) - (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));

qnormal=h(N,d2);

return qnormal;
}
double f_Put_heston(  double S, double X, double dt, double sigma){
	//double a;
	//a= max(S-X,0.);

	  // if nearest to node use different formula
	//which of my node is closest to strike price
	//cout<<"Value: "<<log(S)-log(X)<<" Com "<<sigma*sqrt(dt);
	if(log(S)-log(X)>(sigma*sqrt(dt))){   return 0;}
	else if(log(S)-log(X)<(-sigma*sqrt(dt))) {return  (X-(S*((exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt)))/(2*sigma*sqrt(dt)))));}
	else{return (((X*(sigma*sqrt(dt)-log((S/X))))/(2*sigma*sqrt(dt)))+(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))-(X/S)));}


}
//no sirve
double f_Put_hestonConvertible(  double S, double NB, double dt, double sigma, double F, double Z, double Cp){
	//double a;
	//a= max(S-X,0.);

	  // if nearest to node use different formula
	//which of my node is closest to strike price
	//cout<<"Value: "<<log(S)-log(X)<<" Com "<<sigma*sqrt(dt);
	  //min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	
	  if(log(S)-log(F/Z)>(sigma*sqrt(dt))){   return  (-F/Z+(S*((exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt)))/(2*sigma*sqrt(dt))))) ;}
	  else if(log(S)-log(F/Z)<(sigma*sqrt(dt))&&log(S)-log(F/Z)>(-sigma*sqrt(dt))) {return  (((-F/Z*(sigma*sqrt(dt)+log((S/NB*F))))/(2*sigma*sqrt(dt)))-(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))+(NB*F/S)));}
	  else if(log(S)-log(F/Z)<(-sigma*sqrt(dt)) && log(S)-log(NB*F)>(sigma*sqrt(dt))) {return F ;}

	
	else if(log(S)-log(NB*F)>(sigma*sqrt(dt))&& log(S)-log(F/Z)<(-sigma*sqrt(dt))) {   return F;}
	else if(log(S)-log(NB*F)<(sigma*sqrt(dt))&& log(S)-log(F/Z)>(-sigma*sqrt(dt))) {return (((NB*F*(sigma*sqrt(dt)-log((S/NB*F))))/(2*sigma*sqrt(dt)))+(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))-(NB*F/S)));}
	else {return  (NB*F-(S*((exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt)))/(2*sigma*sqrt(dt)))));}
	



}

// Extrafunction
double binomialPut(double maturity, double sigma, double interestRate, double S0, double X, int N){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]);
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];
}
double amePut(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/(N); //dt=maturity/(N-1); in order to match results on internet
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp((interestRate-div2)*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		  //tree[i][j] = max(tree[i][j],Z*tree[i][j]);// For dBonds
		 
		//  cout << "S_{"<<i<<"," << j << "} =";
		 // cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
 
  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  double tem1=max(X-tree[i][j],0.);
			//cout<<i<<" "<<j<<"    "<<tem1<<"       "<<tree[i][j]<<"\n";
		  valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),tem1);
		  
		  
		  //check values
	// cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }
  }

   
	
return valueTree[0][0];
}
double amePutBarrier(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2, double Barrier){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/(N); //dt=maturity/(N-1); in order to match results on internet
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp((interestRate-div2)*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		  //tree[i][j] = max(tree[i][j],Z*tree[i][j]);// For dBonds
		 
		//  cout << "S_{"<<i<<"," << j << "} =";
		 // cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
 
  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  double tem1=max(X-tree[i][j],0.);
			//cout<<i<<" "<<j<<"    "<<tem1<<"       "<<tree[i][j]<<"\n";
		  valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),tem1);
		  if(tree[i][j]>Barrier){valueTree[i][j]=0.;}
		 // 	cout<<i<<" "<<j<<"    "<<valueTree[i][j]<<"       "<<tree[i][j]<<"\n";
		  //check values
	// cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }
  }

   
	
return valueTree[0][0];
}

};

int main(){

	double r=0.1, sigma=0.2;
	double T=0.5;
	double t=0.;
	double X1=40, X2=10.7; // model  parameters
	double q=0.;
	double S=42;//1*exp(-q);// variable initial stock price*/
	double B=44;
	int N=5000000;

	MonteCarlo Object1;
	DBonds Object2;


std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/12_finiteDifferenceProject/bin.csv");
	 output. precision(10);
	if(!output.is_open()) { std::cout << " File not opened \n";    
	 throw;  }


	for (int N1=1000;N1<=2000;N1++){
	//cout<<Object1.BarrierPut(S,  T,t,  r, q,  sigma,  X1,  N,B)<<" \n ";
	//cout<<Object2.amePut(T,sigma,r,S,X1,2001,q)<<" \n ";
	output<<Object2.amePutBarrier(T,sigma,r,S,X1,N1,q,B)<<" \n ";
	cout<<N1<<"\n";
	
	
	}
	
	  std::cout << " File write successful \n";//
	// CLOSE file
  output.close();//*/


	system("pause");

}