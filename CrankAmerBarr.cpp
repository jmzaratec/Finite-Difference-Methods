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

using namespace std;
 double EuropeanCall(){
 /*
  *   Initialise values and setup storage
  */
  

  // declare parameters for the option
  //double S0=2.,strike=2.,sigma=0.4,interestRate=0.05,dividend=0.,SMax=4.,maturity=1.;

  // Model parameters
	double r=0.1, sigma=0.2;
	double T=0.5;
	double X=40;
	double d1=0.;
	double S0=42;
	double Beta=1.;
	double k=0;

  // declare number of points
  int iMax=2000,jMax=100;
  double SMax=S0*2 ;

  // local parameters
  double dS=SMax/jMax,dt=T/iMax;
  // storage for stock and option value
  vector<double> S(jMax+1),vold(jMax+1),vnew(jMax+1);
  vector<vector<double> > V(iMax+1,vector<double>(jMax+1));
  
  // set values of vnew and vold at i==iMax
  for(int j=0;j<=jMax;j++)
  {
    S[j]=j*dS;
    // update option values as the payoff
    //vold[j] = max(-X + S[j],0.);
    //vnew[j] = max(-X + S[j],0.);
	V[iMax][j]=max(-X + S[j],0.);
	//cout << "V("<<S[j]<<","<<T<<")="<<vnew[j]<<endl;
  }
  /*
  *   Timestep through domain
  */
  for(int i=iMax-1;i>=0;i--)
  {
    // declare storage for matrix diagonals and rhs
    vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
    // MATRIX SETUP
    /*
    *   Set boundary conditions at j=0
		American so this changes V>max(X-S,0)
    */ 
    a[0] = 0.;
	b[0] = 1.;
	c[0] = 0.;
	d[0] = 0;
    /*
    *   loop through equations for 0<j<jMax
    */
    for(int j=1;j<jMax;j++)
    {
      a[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) - 0.25*(r-d1)*j;//0.25*sigma*sigma*j*j - 0.25*(r-d1)*j;
	  b[j] = -0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)  - 0.5*r - 1./dt;
	  c[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) + 0.25*(r-d1)*j;
	// d[j] = - a[j]*vold[j-1] - (1./dt-0.5*sigma*sigma*j*j - 0.5*r)*vold[j] - c[j]*vold[j+1];
	  d[j] = - a[j]* V[i+1][j-1] + ((-1./dt)+(0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)) + 0.5*r)* V[i+1][j] - c[j]* V[i+1][j+1];
	 
    }
    
    /*
    *   Set boundary conditions at j=jmax
    */
	a[jMax] = 0.;
	b[jMax] = 1.;
	c[jMax] = 0.;
	d[jMax] = S[jMax]*exp(-d1*(T-i*dt)) - X*exp(-r*(T-i*dt));
	// check values of a, b, c and d
	//for(int j=0;j<=jMax;j++)
	//{
	//	cout << j << " " << a[j] << " " << b[j] << " " << c[j] << " " << d[j] << endl;
	//}
    
    // FINISHED MATRIX SETUP
    
    // MATRIX SOLVE
    // sor variables
    int iterMax=100;
    double error,tol=1.e-8,omega=1.;
    // sor loop
    for(int sor=0;sor<iterMax;sor++)
    {
      // write equation for boundary value at j=0
	 // vnew[0] = (d[0] - c[0]*vnew[1]) / b[0];
	  V[i][0] = (d[0] - c[0]*V[i][1]) / b[0];
	  //Make it american	  vnew[0] = max(vnew[0],-X+S[0]);

	  // all points in the middle j=1, j<jMax
	  for(int j=1;j<jMax;j++)
	  {
		 // vnew[j] = (d[j] - a[j]*vnew[j-1] - c[j]*vnew[j+1])/b[j];
		  V[i][j] = (d[j] - a[j]*V[i][j-1] - c[j]*V[i][j+1])/b[j];
		  //make it american  vnew[j] = max(vnew[j],X-S[j]);
	  }
	  // boundary value at j=jMax
	 // vnew[jMax] = (d[jMax] - a[jMax]*vnew[jMax-1]) / b[jMax];
	  V[i][jMax] = (d[jMax] - a[jMax]*V[i][jMax-1]) / b[jMax];

	  //make it american  vnew[jMax] = max(vnew[jMax],X-S[jMax]);
    }// finish sor

    // FINISH MATRIX SOLVE
    /*
    *   don't forget to set old values equal to new
    */ 

	for(int j=0;j<=jMax;j++)
	{
		//cout << i << " " << S[j] << " " << vold[j] << " " << vnew[j] << endl;
	//cout << "V("<<S[j]<<","<<i*dt<<")="<<vnew[j]<<endl;
	}

	vold = vnew;
    
  } // Finish Timestep

  //cout << "V("<<S[jMax/2]<<","<<0.*dt<<")="<<vnew[jMax/2]<<endl;
 // cout << "V("<<S[jMax/2]<<","<<0.*dt<<")="<<V[0][jMax/2]<<endl;
    	

/*/ output		
 std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/12_finiteDifferenceProject/output2.csv");
	 output. precision(10);
	if(!output.is_open()) { std::cout << " File not opened \n";    
	 throw;  }
//Writing headings
	output<<"stock[j]"<<"\n";
	for(int i=0; i<=iMax;i++){	output<<", "<<dt*i;	}

	output<<"\n";


// check results
	for(int j=0; j<=jMax;j++){
	output<<S[j]<<",";
		for(int i=0; i<=iMax;i++){
			output<<V[i][j]<<",";
	//cout<<"V("<<stock[j]<<","<<i*dt<<")=";
	//cout<<vNew[j]<<"\n";
	}
	output<<"\n";
	}

  std::cout << " File write successful \n";//
	// CLOSE file
  output.close();

  system("pause");//*/
  return 0;
 }
 double AmericanPut(){
  /*
  *   Initialise values and setup storage
  */
  

  // declare parameters for the option
  //double S0=2.,strike=2.,sigma=0.4,interestRate=0.05,dividend=0.,SMax=4.,maturity=1.;

  // Model parameters
	double r=0.1, sigma=0.2;
	double T=0.5;
	double X=40;
	double d1=0.;
	double S0=42;
	double Beta=1.;
	double k=0;

  // declare number of points
  int iMax=2000,jMax=100;
  double SMax=S0*2 ;

  // local parameters
  double dS=SMax/jMax,dt=T/iMax;
  // storage for stock and option value
  vector<double> S(jMax+1),vold(jMax+1),vnew(jMax+1);
  vector<vector<double> > V(iMax+1,vector<double>(jMax+1));
  
  // set values of vnew and vold at i==iMax
  for(int j=0;j<=jMax;j++)
  {
    S[j]=j*dS;
    // update option values as the payoff
    //vold[j] = max(-X + S[j],0.);
    //vnew[j] = max(-X + S[j],0.);
	V[iMax][j]=max(X - S[j],0.);
	//cout << "V("<<S[j]<<","<<T<<")="<<vnew[j]<<endl;
  }
  /*
  *   Timestep through domain
  */
  for(int i=iMax-1;i>=0;i--)
  {
    // declare storage for matrix diagonals and rhs
    vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
    // MATRIX SETUP
    /*
    *   Set boundary conditions at j=0
		American so this changes V>max(X-S,0)
    */ 
    a[0] = 0.;
	b[0] = 1.;
	c[0] = 0.;
	d[0] = X*exp(-r*(T-i*dt));;
    /*
    *   loop through equations for 0<j<jMax
    */
    for(int j=1;j<jMax;j++)
    {
      a[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) - 0.25*(r-d1)*j;//0.25*sigma*sigma*j*j - 0.25*(r-d1)*j;
	  b[j] = -0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)  - 0.5*r - 1./dt;
	  c[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) + 0.25*(r-d1)*j;
	// d[j] = - a[j]*vold[j-1] - (1./dt-0.5*sigma*sigma*j*j - 0.5*r)*vold[j] - c[j]*vold[j+1];
	  d[j] = - a[j]* V[i+1][j-1] + ((-1./dt)+(0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)) + 0.5*r)* V[i+1][j] - c[j]* V[i+1][j+1];
	 
    }
    
    /*
    *   Set boundary conditions at j=jmax
    */
	a[jMax] = 0.;
	b[jMax] = 1.;
	c[jMax] = 0.;
	d[jMax] = 0.;
	// check values of a, b, c and d
	//for(int j=0;j<=jMax;j++)
	//{
	//	cout << j << " " << a[j] << " " << b[j] << " " << c[j] << " " << d[j] << endl;
	//}
    
    // FINISHED MATRIX SETUP
    
    // MATRIX SOLVE
    // sor variables
    int iterMax=100;
    double error,tol=1.e-8,omega=1.;
    // sor loop
    for(int sor=0;sor<iterMax;sor++)
    {
      // write equation for boundary value at j=0
	 // vnew[0] = (d[0] - c[0]*vnew[1]) / b[0];
	  V[i][0] = (d[0] - c[0]*V[i][1]) / b[0];
	  //Make it american	  vnew[0] = max(vnew[0],-X+S[0]);
	  V[i][0]=max(V[i][0],X-S[0]);
	  // all points in the middle j=1, j<jMax
	  for(int j=1;j<jMax;j++)
	  {
		 // vnew[j] = (d[j] - a[j]*vnew[j-1] - c[j]*vnew[j+1])/b[j];
		  V[i][j] = (d[j] - a[j]*V[i][j-1] - c[j]*V[i][j+1])/b[j];
		  //make it american  vnew[j] = max(vnew[j],X-S[j]);
		    V[i][j]=max(V[i][j],X-S[j]);
	  }
	  // boundary value at j=jMax
	 // vnew[jMax] = (d[jMax] - a[jMax]*vnew[jMax-1]) / b[jMax];
	  V[i][jMax] = (d[jMax] - a[jMax]*V[i][jMax-1]) / b[jMax];
	  //make it american  vnew[jMax] = max(vnew[jMax],X-S[jMax]);
	      V[i][jMax]=max(V[i][jMax],X-S[jMax]);
    }// finish sor

    // FINISH MATRIX SOLVE
    /*
    *   don't forget to set old values equal to new
    */ 

	for(int j=0;j<=jMax;j++)
	{
		//cout << i << " " << S[j] << " " << vold[j] << " " << vnew[j] << endl;
	//cout << "V("<<S[j]<<","<<i*dt<<")="<<vnew[j]<<endl;
	}

	vold = vnew;
    
  } // Finish Timestep

  //cout << "V("<<S[jMax/2]<<","<<0.*dt<<")="<<vnew[jMax/2]<<endl;
  cout << "V("<<S[jMax/2]<<","<<0.*dt<<")="<<V[0][jMax/2]<<endl;
    	

/*/ output		
 std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/12_finiteDifferenceProject/output2.csv");
	 output. precision(10);
	if(!output.is_open()) { std::cout << " File not opened \n";    
	 throw;  }
//Writing headings
	output<<"stock[j]"<<"\n";
	for(int i=0; i<=iMax;i++){	output<<", "<<dt*i;	}

	output<<"\n";


// check results
	for(int j=0; j<=jMax;j++){
	output<<S[j]<<",";
		for(int i=0; i<=iMax;i++){
			output<<V[i][j]<<",";
	//cout<<"V("<<stock[j]<<","<<i*dt<<")=";
	//cout<<vNew[j]<<"\n";
	}
	output<<"\n";
	}

  std::cout << " File write successful \n";//
	// CLOSE file
  output.close();

  system("pause");//*/
  return 0;

 }
double ABPO(double r, double sigma,	double T,double X,double d1,double S0,double Beta,double k, double Barrier)
{
  // declare number of points
  int iMax=1000,jMax=100;//2720;
  double SMax=S0*2 ;
    // local parameters
   double dS=SMax/jMax;
	double dt=T/iMax;
  // storage for stock and option value
  vector<double> S(jMax+1),vold(jMax+1),vnew(jMax+1);
  vector<vector<double> > V(iMax+1,vector<double>(jMax+1));
  
  // initial set up
  for(int j=0;j<=jMax;j++)  {
	  S[j]=j*dS;
	V[iMax][j]=max(X - S[j],0.);
	if(S[j]>Barrier){V[iMax][j]=0;} //barrier
	  }
//LOOP
  for(int i=iMax-1;i>=0;i--)
  {
       // declare storage for matrix diagonals and rhs
    vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
    // MATRIX SETUP
    /*
    *   Set boundary conditions at j=0
		American so this changes V>max(X-S,0)
    */ 
    a[0] = 0.;
	b[0] = 1.;
	c[0] = 0.;
	d[0] = X*exp(-r*(T-i*dt))-S[0];
    for(int j=1;j<jMax;j++)    {
      a[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) - 0.25*(r-d1)*j;
	  b[j] = -0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)  - 0.5*r - 1./dt;
	  c[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) + 0.25*(r-d1)*j;
	  d[j] = - a[j]* V[i+1][j-1] + ((-1./dt)+(0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)) + 0.5*r)* V[i+1][j] - c[j]* V[i+1][j+1];
	    }
        //boundaries
	a[jMax] = 0.;
	b[jMax] = 1.;
	c[jMax] = 0.;
	d[jMax] = 0.;
    // FINISHED MATRIX SETUP
    
    // MATRIX SOLVE
    // sor 
    int iterMax=100;
    double error,tol=1.e-8,omega=1.;
    // sor loop
    for(int sor=0;sor<iterMax;sor++)    {
      
		V[i][0] = (d[0] - c[0]*V[i][1]) / b[0];
	  //Make it american	  vnew[0] = max(vnew[0],-X+S[0]);
	  V[i][0]=max(V[i][0],X-S[0]);
	   if(S[0]>Barrier){V[i][0]=0;}
	  for(int j=1;j<jMax;j++)	  {
		  V[i][j] = (d[j] - a[j]*V[i][j-1] - c[j]*V[i][j+1])/b[j];
		  //make it american  vnew[j] = max(vnew[j],X-S[j]);
		    V[i][j]=max(V[i][j],X-S[j]);
			if(S[j]>Barrier){V[i][j]=0;}//barrier
		  }
	  // boundary value at j=jMax
	 /  V[i][jMax] = (d[jMax] - a[jMax]*V[i][jMax-1]) / b[jMax];
	  //make it american  vnew[jMax] = max(vnew[jMax],X-S[jMax]);
	      V[i][jMax]=max(V[i][jMax],X-S[jMax]);
		  	if(S[jMax]>Barrier){V[i][jMax]=0;}
    }// finish sor

    // FINISH MATRIX SOLVE
  
  } // Finish Timestep

 
return V[0][jMax/2];
}
 double interpolateV(vector<double>& x,vector<vector<double> >& y,double a,int degree){
  int istar;
  double dx = x[1]-x[0];
istar = a/dx;
//cout<<a<<"  "<<istar<<"  "<<istar*dx<<"  "<<(istar+1)*dx;
	// interpolate in here
  double sum=0.;
  double l_j;
  int i;
  //for eah j=istar up to istar+degree
 for (int j = istar; j<istar+degree;j++){
  l_j=1.;
  //for each i= istar up to istar +degree in degree, here istar and istar +1 
	for (int i=istar;i<istar+degree;i++){ 
 		if (i==j) continue;  //if j==i do nothing, so i++ //
  		l_j=l_j *(a-x[i])/(x[j]-x[i]);
		}
		sum = sum + l_j*y[0][j]; 
	}

	  return sum;
	}

int main(){
   cout. precision(15);
	
	double r=0.069;
	double T=1.4;
	double X=1600;
	double d1=0.045;
	double S0=1570.7;//X;
	double Barrier=1700;
	//*Sce3
	double Beta=0.28;
	double k=0.39;
	double sigma=63;//*/

  // declare number of points
  int iMax=20000,jMax=2720;
  double SMax=2720 ;

  // local parameters
   double dS=SMax/jMax; 
	double dt=T/iMax;
  // storage for stock and option value
  vector<double> S(jMax+1),vold(jMax+1),vnew(jMax+1);
  vector<vector<double> > V(iMax+1,vector<double>(jMax+1));
  
  // set values of vnew and vold at i==iMax
  for(int j=0;j<=jMax;j++)  {
   S[j]=j*dS;
  	V[iMax][j]=max(X - S[j],0.);
	if(S[j]>Barrier){V[iMax][j]=0;}
	  }
    for(int i=iMax-1;i>=0;i--)  {
     // cout<<i<<"\n";
	  // declare storage for matrix diagonals and rhs
    vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
    // MATRIX SETUP
    /*
    *   Set boundary conditions at j=0
		American so this changes V>max(X-S,0)
    */ 
    a[0] = 0.;
	b[0] = 1.;
	c[0] = 0.;
	d[0] = X*exp(-r*(T-i*dt)); //ABPO(r,sigma,dt*(iMax-i),X,d1,S[0],Beta,k,Barrier);//Change this to adjuste the boundaries to a Smin
      for(int j=1;j<jMax;j++)
    {
      a[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) - 0.25*(r-d1)*j;
	  b[j] = -0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)  - 0.5*r - 1./dt;
	  c[j] = 0.25*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i) + 0.25*(r-d1)*j;
	  d[j] = - a[j]* V[i+1][j-1] + ((-1./dt)+(0.5*sigma*sigma*pow(j,2*Beta)*pow(dS,2*Beta-2)*exp(-k*dt*i)) + 0.5*r)* V[i+1][j] - c[j]* V[i+1][j+1];
	 
    }
    
  //Boundaries
	a[jMax] = 0.;
	b[jMax] = 1.;
	c[jMax] = 0.;
	d[jMax] = 0.;

    // FINISHED MATRIX SETUP
    
    // MATRIX SOLVE
    // sor variables
    int iterMax=100;
    double error,tol=1.e-8,omega=1.;
    // sor loop
    for(int sor=0;sor<iterMax;sor++)    {
    		V[i][0] = (d[0] - c[0]*V[i][1]) / b[0];
	  //Make it american	  vnew[0] = max(vnew[0],-X+S[0]);
	  V[i][0]=max(V[i][0],X-S[0]);
	  if(S[0]>Barrier){V[i][0]=0;}//Barrier
	   for(int j=1;j<jMax;j++)  {
		   V[i][j] = (d[j] - a[j]*V[i][j-1] - c[j]*V[i][j+1])/b[j];
		  //make it american  vnew[j] = max(vnew[j],X-S[j]);
		    V[i][j]=max(V[i][j],X-S[j]);
			if(S[j]>Barrier){V[i][j]=0;}
			  }
	  // boundary value at j=jMax
	   V[i][jMax] = (d[jMax] - a[jMax]*V[i][jMax-1]) / b[jMax];
	  //make it american  vnew[jMax] = max(vnew[jMax],X-S[jMax]);
	      V[i][jMax]=max(V[i][jMax],X-S[jMax]);
		  	if(S[jMax]>Barrier){V[i][jMax]=0;}
    }// finish sor

    // FINISH MATRIX SOLVE

  } // Finish Timestep


// output GRID!!		
 std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/12_finiteDifferenceProject/output3.csv");
	 output. precision(15);
	if(!output.is_open()) { std::cout << " File not opened \n";    
	 throw;  }
//Writing headings
	output<<"stock[j]"<<"\n";
	for(int i=0; i<=iMax;i++){	output<<", "<<dt*i;	}
	output<<"\n";
// check results
	for(int j=0; j<=jMax;j++){
	output<<S[j]<<",";
		for(int i=0; i<=iMax;i++){
			output<<V[i][j]<<",";
	}
	output<<"\n";
	}

  std::cout << " File write successful \n";//
	// CLOSE file
  //output.close();//*/


  /*/output DELTA		
 std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/12_finiteDifferenceProject/outputAme.csv");
	 output. precision(15);
	if(!output.is_open()) { std::cout << " File not opened \n";    
	 throw;  }
//Writing headings
	output<<"stock[j], V(S-t=0),Delta"<<"\n";
	
	output<<S[0]<<",";
		output<<V[0][0]<<",";
		//output<<(1/(4*dS))*(V[0][j+1]-V[0][j-1]+V[1][j+1]-V[1][j-1]);
		output<<"\n";
	for(int j=1; j<jMax;j++){
		output<<S[j]<<",";
		output<<V[0][j]<<",";
		output<<(1/(4*dS))*(V[0][j+1]-V[0][j-1]+V[1][j+1]-V[1][j-1]);
		output<<"\n";	
	//cout<<"V("<<stock[j]<<","<<i*dt<<")=";
	//cout<<vNew[j]<<"\n";
	}
		output<<S[jMax]<<",";
		output<<V[0][jMax]<<",";
		//output<<(1/(4*dS))*(V[0][j+1]-V[0][j-1]+V[1][j+1]-V[1][j-1]);
		output<<"\n";

	std::cout << " File write successful \n";//
	// CLOSE file
  output.close();//*/
  double a=1570.7;
 output << " Value at x=,"<<a<<", is ,"<<interpolateV(S,V,a,2) << "\n";
 output << " Value at x=,"<<a<<", is ,"<<interpolateV(S,V,a,3) << "\n";
 output << " Value at x=,"<<a<<", is ,"<<interpolateV(S,V,a,5) << "\n";
 output << " Value at x=,"<<a<<", is ,"<<interpolateV(S,V,a,10) << "\n";
 output << " Value at x=,"<<a<<", is ,"<<interpolateV(S,V,a,30) << "\n";
 output << " Value at x=,"<<a<<", is ,"<<interpolateV(S,V,a,50) << "\n";
 output << " Value at x=,"<<a<<", is ,"<<interpolateV(S,V,a,70) << "\n";


  output.close();//*/

  system("pause");

}