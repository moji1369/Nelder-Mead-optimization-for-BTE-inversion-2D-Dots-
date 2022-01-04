#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <time.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <iomanip>
#include <limits>

using namespace std;

//this is the function that solves the forward  problem (BTE) using adjoint (backward) Monte Carlo (MC) method:
double BTE_bkwrd(double x[]);

//******************************************************************************
//******************************************************************************
//******************************************************************************

int main ()
{ 
    ofstream result;   //the file that the final optimized parameters will be written to
    int n={8};         //size of the optimized parameters (8=2*3+2: 1 function (same for all branches) with a piecewise linear functions with 3 lines,
	               //where each line is parameterized with 2 parameters; and another 2 for the line parameterizing the transmissivity
    int maxfun=10000;  //max number of function calculations allowed in the optimization
    int maxiter=10000; //max number of iterations allowed in the optimization
    double tolf=1e-30; //function tolerance
    double tolx=1e-30; //parameter tolerance
    double x[8]={4.5,11,1.4e13,9.7,5.5e13,11,0.5,0.5};      //initial condition

    double IC[n]; 
    for (int i=0;i<n;++i) {x[i]=x[i]*1.0; IC[i]=x[i];}      //initial simplex of NM
    
    //NM parameters:
    double rho={1}; 
    double chi={2}; 
    double psi={0.5}; 
    double sigma={0.95};

    //initialization of the simplex:
    double onesn[n];
    for (int i=0;i<n;++i) {onesn[i] = 1.0;}
    double two2np1[n];
    for (int i=0;i<n;++i) {two2np1[i] = i+2.0;}
    double one2n[n];
    for (int i=0;i<n;++i) {one2n[i] = i+1.0;}
    double xin[n];
    for (int i=0;i<n;++i) {xin[i]=x[i];}
    double v[n][n+1];
    for (int i=0;i<n;++i) {for (int j=0;j<n+1;++j) {v[i][j]=0.0;}}
    double fv[n+1];
    for (int i=0;i<n+1;++i) {fv[i] = 0.0;}
    for (int i=0;i<n;++i) {v[i][0]=xin[i];}
    fv[0]=BTE_bkwrd(x);
    int func_evals={1};
    int itercount={0};
    cout<<"Nelder-Mead algorithm with Lagarias correction. Created by Mojtaba Forghani \n\n"; 
    cout<<setw(5)<<"iteration"<<setw(15)<<"function eval";
    cout<<setw(15)<<"function val"<<"\n";
    cout<<setw(5)<<itercount<<setw(15)<<func_evals<<setw(15)<<fv[0]<<"\n";
    double usual_delta={0.1};
    double zero_term_delta ={0.00025};

    //construction of the initial simplex:
    for (int j=0;j<n;++j)
    {
        double y[n];
        for (int i=0;i<n;++i) {y[i]= xin[i];} 
        if (y[j] != 0) { y[j]=(1.0+usual_delta)*y[j];}
        else { y[j]=zero_term_delta;}
        for (int i=0;i<n;++i) { v[i][j+1]=y[i]; }
        for (int i=0;i<n;++i) {x[i]=y[i];}
        double f;
        f=BTE_bkwrd(x);
        fv[j+1]=f;
    }

    //sorting the simplex:
    int J[n+1]; for (int i=0;i<(n+1);++i) { J[i]=i; }
    for (int i=(n);i>=0;i--) { for (int j=1;j<=i;j++) { if (fv[j-1]>fv[j])
    {
        double a=fv[j-1]; int b=J[j-1];
        fv[j-1]=fv[j]; J[j-1]=J[j];
        fv[j]=a; J[j]=b; } } 
    }
    double v_temp[n][n+1];
    for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) {v_temp[j][i]=v[j][J[i]];}}
    for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) {v[j][i]=v_temp[j][i];}}
    string how = "initial simplex";
    itercount=itercount+1;
    func_evals=n+1;
    cout<<setw(5)<<itercount<<setw(15)<<func_evals;
    cout<<setw(15)<<fv[0]<<setw(25)<<how<<"\n";

    //iterations of the NM algorithm:
    while (func_evals<maxfun && itercount<maxiter) 
    {
        //sorting of simplex values at each iteration:
        double max_fv=abs(fv[0]-fv[1]);
        for (int i=1;i<(n+1);++i) 
        { 
            if (abs(fv[0]-fv[i])<max_fv) {max_fv=abs(fv[0]-fv[i]);}
        }
        double max_v=abs(v[0][0]-v[0][1]);
        for (int i=1;i<(n+1);++i) {for (int j=0;j<n;++j) 
        { 
            if (abs(v[j][0]-v[j][i])<max_v) {max_v=abs(v[j][0]-v[j][i]);}}
        }

        //stopping criteria:
        if (max_fv <= tolf && max_v<=tolx) {break;}
        double xbar[n];
        for (int i=0;i<n;++i) {xbar[i]=0.0;}
        for (int i=0;i<n;++i) {for (int j=0;j<n;++j) 
        {
            xbar[i]=xbar[i]+v[i][j]/n;}
        }
        
        //check if "expand", "reflect", "shrink" or "contract" step should take place:
        double xr[n];
        for (int i=0;i<n;++i) {xr[i]=(1.0+rho)*xbar[i]-rho*v[i][n];}
        for (int i=0;i<n;++i) {x[i]=xr[i];}
        double fxr;
        fxr=BTE_bkwrd(x);
        func_evals=func_evals+1;
        if (fxr < fv[0])
        {
            double xe[n];
            for (int i=0;i<n;++i) {xe[i]=(1.0+rho*chi)*xbar[i]-rho*chi*v[i][n];}
            for (int i=0;i<n;++i) {x[i]=xe[i];}
            double fxe;
            fxe=BTE_bkwrd(x);
            func_evals=func_evals+1;
            if (fxe<fxr) 
            { 
                for (int i=0;i<n;++i) {v[i][n]=xe[i];}
                fv[n]=fxe;
                how="expand";
            }
            else
            {
                for (int i=0;i<n;++i) {v[i][n]=xr[i];}
                fv[n]=fxr;
                how="reflect";
            } 
        }
        else
        { 
            if (fxr<fv[n-1]) 
            { 
                for (int i=0;i<n;++i) {v[i][n]=xr[i];}
                fv[n]=fxr;
                how="reflect";
            }
            else
            {
                if (fxr<fv[n])
                {
                    double xc[n];
                    for (int i=0;i<n;++i) 
                    {
                        xc[i]=(1.0+psi*rho)*xbar[i]-psi*rho*v[i][n];
                    }
                    for (int i=0;i<n;++i) {x[i]=xc[i];}
                    double fxc;
                    fxc=BTE_bkwrd(x);
                    func_evals=func_evals+1;
                    if (fxc<=fxr)
                    {
                        for (int i=0;i<n;++i) {v[i][n]=xc[i];}
                        fv[n]=fxc;
                        how="contract outside";
                    }
                    else {how="shrink";}
                }
                else
                {
                    double xcc[n];
                    for (int i=0;i<n;++i) 
                    {
                        xcc[i]=(1.0-psi)*xbar[i]+psi*v[i][n];
                    }
                    for (int i=0;i<n;++i) {x[i]=xcc[i];}
                    double fxcc;
                    fxcc=BTE_bkwrd(x);
                    func_evals=func_evals+1;
                    if (fxcc<fv[n])
                    {
                        for (int i=0;i<n;++i) {v[i][n]=xcc[i];}
                        fv[n]=fxcc;
                        how="contract inside";
                    }
                    else {how="shrink";}
                }
                if (how=="shrink")
                { 
                    for (int j=1;j<(n+1);++j)
                    {
                        for (int i=0;i<n;++i) 
                        {
                            v[i][j]=v[i][0]+sigma*(v[i][j]-v[i][0]);
                        }
                        for (int i=0;i<n;++i) {x[i]=v[i][j];}
                        fv[j]=BTE_bkwrd(x);
                    }
                    func_evals=func_evals+n;
                }
            }
        }

        //sorting of simplex values:
        int J[n+1]; for (int i=0;i<(n+1);++i) { J[i]=i; }
        for (int i=(n);i>=0;i--) { for (int j=1;j<=i;j++) { if (fv[j-1]>fv[j])
        {
            double a=fv[j-1]; int b=J[j-1];
            fv[j-1]=fv[j]; J[j-1]=J[j];
            fv[j]=a; J[j]=b; } } 
        }
        double v_temp[n][n+1];
        for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) 
        {
            v_temp[j][i]=v[j][J[i]];}
        }
        for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) 
        {
            v[j][i]=v_temp[j][i];}
        }

        //updating and outputing each iteration result:
        itercount=itercount+1;
        cout<<setw(5)<<itercount<<setw(15)<<func_evals;
        cout<<setw(15)<<fv[0]<<setw(25)<<how<<"\n";
        for (int i=0;i<n;++i) { cout<<v[i][0]<<" "; }
		cout<<"\n";
        result.open("result_5n_ord_3_1_2_2e4_loglog_init_4_52_11_18_1_38e13_9_71_5_51e13_10_92_0_67_0_57_alpha_0_1_del_0_1.txt");
        for (int i=0;i<n;++i) {result<<IC[i]<<" ";} result<<"\n";
        result<<fv[0]<<"\n"; for (int i=0;i<n;++i) { result<<v[i][0]<<" "; }
        result<<"\n";
        result.close();
    }
    for (int i=0;i<n;++i) {x[i]=v[i][0];}
    cout<<"final x is: \n";
    for (int i=0;i<n;++i) {cout<<"   "<<setw(15)<<x[i]<<"\n";}
    result.close();
    cin>>n;
}

//******************************************************************************
//******************************************************************************

double BTE_bkwrd ( double x[8] ) 
{
    double y={0.0};      //objective function value
    
    //random generation setup:
    srand (time(NULL));
    typedef boost::mt19937 RNGType;
    RNGType rng( time(0) );
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator< RNGType, boost::uniform_real<> > uni(rng, uni_dist);

    //physical (model) parameters:
    double PI = 3.14159265;
    double hbar=1.054517e-34;                            
    double boltz=1.38065e-23;                              
    double Teq = 300;    //base temperature
    int Nmode_Si_first=1;
    int Nmode_Si_LA_TA=1000;
    int Nmode_Si_last=1399;
    int Nmode_Al_first=1;
    int Nmode_Al_last=1498;
    int Ntt=101;         //number of measurements

    //load and setup Al and Si material data:
    ifstream file_Si("BVK_Si.txt");
    double BVK_Si[4][1399];
    if(file_Si.is_open())
    {
        for(int i = 0; i < 1399; ++i)
        {
        for(int j = 0; j < 4; ++j)
        {
            file_Si >> BVK_Si[j][i];
        }
        }
        file_Si.close();
    }
    ifstream file_Al("BVK_Al.txt");
    double BVK_Al[3][1498];
    if(file_Al.is_open())
    {
        for(int i = 0; i < 1498; ++i)
        {
        for(int j = 0; j < 3; ++j)
        {
                file_Al >> BVK_Al[j][i];
        }
        }
        file_Al.close();
    }
    double SD_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        SD_Si[i] = BVK_Si[1][i+Nmode_Si_first-1];
    }
    double SD_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        SD_Al[i] = BVK_Al[1][i+Nmode_Al_first-1];
    }
    double Dom_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        Dom_Si[i] = BVK_Si[3][i+Nmode_Si_first-1];
    }
    double Dom_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        Dom_Al[i]=Dom_Si[0]; 
    }
    double V_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        V_Si[i] = BVK_Si[2][i+Nmode_Si_first-1];
    }
    double V_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        V_Al[i] = BVK_Al[2][i+Nmode_Al_first-1];
    }
    double F_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        F_Si[i] = BVK_Si[0][i+Nmode_Si_first-1];
    }
    double F_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        F_Al[i] = BVK_Al[0][i+Nmode_Al_first-1];
    }
    double tau_inv_Al[Nmode_Al_last-Nmode_Al_first+1];
    double tau_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i< Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        tau_inv_Al[i]=1e11;
        tau_Al[i]=1.0/tau_inv_Al[i];
    }
    double BL=1/5.32e18/0.9/2/PI/2/PI;
    double BT=1/5.07e18/1.6/4/pow(PI,2);
    double wb=1.2e6;
    double A=3e-45;
    double tau_inv_Si[Nmode_Si_last-Nmode_Si_first+1];
    double tau_Si[Nmode_Si_last-Nmode_Si_first+1];

    //**************************************************************************
    //**************************************************************************
    //relaxation times function parameterization:
    
    //ORDER 0 TAU:
    double Del=0.25e13;
    double tau_inv_Si_temp[Nmode_Si_last-Nmode_Si_first+1];
    /*for (int i=0; i < Nmode_Si_LA_TA-Nmode_Si_first+1; ++i)
    {
        tau_inv_Si_temp[i]=2.464268144426204e+10;
        tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
        tau_Si[i]=1.0/tau_inv_Si[i];
    }
    for (int i=(Nmode_Si_LA_TA-Nmode_Si_first+1);
	    i < (Nmode_Si_last-Nmode_Si_first+1); ++i)
    {
    	tau_inv_Si_temp[i]=2.464268144426204e+10;
        tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
        tau_Si[i]=1.0/tau_inv_Si[i];
    }*/
        
    //************************************************************************** 
    //**************************************************************************
    
    //ORDER 1 TAU:
    /*for (int i=0; i < Nmode_Si_LA_TA-Nmode_Si_first+1; ++i)
    {
        tau_inv_Si_temp[i]=pow(10.0,(x[1]-x[0])/(log10(F_Si[Nmode_Si_LA_TA-Nmode_Si_first])-
            log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
        tau_Si[i]=1.0/tau_inv_Si[i];
    }
    for (int i=(Nmode_Si_LA_TA-Nmode_Si_first+1);
	    i < (Nmode_Si_last-Nmode_Si_first+1); ++i)
    {
    	tau_inv_Si_temp[i]=pow(10.0,(x[1]-x[0])/(log10(F_Si[Nmode_Si_LA_TA-Nmode_Si_first])-
		    log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
        tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
        tau_Si[i]=1.0/tau_inv_Si[i];
    }*/
    
    //************************************************************************** 
    //************************************************************************** 
    
    //ORDER 2 TAU:    
    /*double a_1=(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]));
    double b_1=-log10(F_Si[0])*a_1+x[0];
    double a_2=(x[1]-x[3])/(log10(F_Si[999])-log10(x[2]));
    double b_2=-log10(x[2])*a_2+x[3];
    double X_1=log10(x[2]-Del);
    double Y_1=x[0]+x[3]-x[0]*log10((x[2]-Del)/F_Si[0])/log10(x[2]/F_Si[0]);
    double X_2=log10(x[2]+Del);
    double Y_2=x[3]+(x[1]-x[3])*log10(1.0+Del/x[2])/log10(F_Si[999]/x[2]);
    double a=-((x[1]-x[3])/log10(F_Si[999]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
	         pow(log10((x[2]+Del)/(x[2]-Del)),3.0)*log10(1.0-(Del/x[2])*(Del/x[2]));
    double b=0.5*(((x[1]-x[3])/log10(F_Si[999]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
             log10((x[2]+Del)/(x[2]-Del))-3*a*log10(x[2]*x[2]-Del*Del));
    double c=(x[3]-x[0])/log10(x[2]/F_Si[0])-3*a*log10(x[2]-Del)*log10(x[2]-Del)-2*b*log10(x[2]-Del);
    double d=x[0]-log10(F_Si[0])*(x[3]-x[0])/log10(x[2]/F_Si[0])+0.5*((x[1]-x[3])/
	         log10(F_Si[999]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/log10((x[2]+Del)/
			 (x[2]-Del))*log10(x[2]-Del)*log10(x[2]-Del)-0.5*a*log10(x[2]-Del)*log10(x[2]-
			 Del)*log10((x[2]+Del)*(x[2]+Del)*(x[2]+Del)/(x[2]-Del));
    for (int i=0;i<1000;++i)
    {
	    if (log10(F_Si[i])<X_1)
	    {
            tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else if (log10(F_Si[i])<X_2)
        {
            tau_inv_Si_temp[i]=pow(10,a*pow(log10(F_Si[i]),3.0)+b*pow(log10(F_Si[i]),2.0)+c*log10(F_Si[i])+d);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else
        {
            tau_inv_Si_temp[i]=pow(10,(x[1]-x[3])/(log10(F_Si[999])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
    }
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || log10(x[2])>log10(F_Si[999]))
    {
    	y=y+1e5;
    }
    //**************************************************************************
    a_1=(x[7]-x[4])/(log10(x[6])-log10(F_Si[0]));
    b_1=-log10(F_Si[0])*a_1+x[4];
    a_2=(x[5]-x[7])/(log10(F_Si[999])-log10(x[6]));
    b_2=-log10(x[6])*a_2+x[7];
    X_1=log10(x[6]-Del);
    Y_1=x[4]+x[7]-x[4]*log10((x[6]-Del)/F_Si[0])/log10(x[6]/F_Si[0]);
    X_2=log10(x[6]+Del);
    Y_2=x[7]+(x[5]-x[7])*log10(1.0+Del/x[6])/log10(F_Si[999]/x[6]);
    a=-((x[5]-x[7])/log10(F_Si[999]/x[6])-(x[7]-x[4])/log10(x[6]/F_Si[0]))/
	         pow(log10((x[6]+Del)/(x[6]-Del)),3.0)*log10(1.0-(Del/x[6])*(Del/x[6]));
    b=0.5*(((x[5]-x[7])/log10(F_Si[999]/x[6])-(x[7]-x[4])/log10(x[6]/F_Si[0]))/
             log10((x[6]+Del)/(x[6]-Del))-3*a*log10(x[6]*x[6]-Del*Del));
    c=(x[7]-x[4])/log10(x[6]/F_Si[0])-3*a*log10(x[6]-Del)*log10(x[6]-Del)-2*b*log10(x[6]-Del);
    d=x[4]-log10(F_Si[0])*(x[7]-x[4])/log10(x[6]/F_Si[0])+0.5*((x[5]-x[7])/
	         log10(F_Si[999]/x[6])-(x[7]-x[4])/log10(x[6]/F_Si[0]))/log10((x[6]+Del)/
			 (x[6]-Del))*log10(x[6]-Del)*log10(x[6]-Del)-0.5*a*log10(x[6]-Del)*log10(x[6]-
			 Del)*log10((x[6]+Del)*(x[6]+Del)*(x[6]+Del)/(x[6]-Del));
    for (int i=1000;i<(1399);++i)
    {
	    if (log10(F_Si[i])<X_1)
	    {
            tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else if (log10(F_Si[i])<X_2)
        {
            tau_inv_Si_temp[i]=pow(10,a*pow(log10(F_Si[i]),3.0)+b*pow(log10(F_Si[i]),2.0)+c*log10(F_Si[i])+d);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else
        {
            tau_inv_Si_temp[i]=pow(10,(x[1]-x[3])/(log10(F_Si[999])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
    }
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || log10(x[2])>log10(F_Si[999]))
    {
    	y=y+1e5;
	}
	for (int i=1000;i<(1399);++i)
    {
	    if (log10(F_Si[i])<X_1)
	    {
            tau_inv_Si_temp[i]=pow(10,(x[7]-x[4])/(log10(x[6])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[4]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else if (log10(F_Si[i])<X_2)
        {
            tau_inv_Si_temp[i]=pow(10,a*pow(log10(F_Si[i]),3.0)+b*pow(log10(F_Si[i]),2.0)+c*log10(F_Si[i])+d);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else
        {
            tau_inv_Si_temp[i]=pow(10,(x[5]-x[7])/(log10(F_Si[999])-log10(x[6]))*(log10(F_Si[i])-log10(x[6]))+x[7]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
    }
    if (log10(F_Si[0])>x[6] || X_1>log10(x[6]) || log10(x[6])>X_2 || log10(x[6])>log10(F_Si[1398]))
    {
    	y=y+1e5;
    }*/
	
    //**************************************************************************
    //**************************************************************************
	
    //ORDRER 3 TAU:
    double a_1=(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]));
    double b_1=-log10(F_Si[0])*a_1+x[0];
    double a_2=(x[5]-x[3])/(log10(x[4])-log10(x[2]));
    double b_2=-log10(x[2])*a_2+x[3];
    double a_3=(x[1]-x[5])/(log10(F_Si[999])-log10(x[4]));
    double b_3=-log10(x[4])*a_3+x[5];
    double X_1=log10(x[2]-Del);
    double Y_1=x[0]+x[3]-x[0]*log10((x[2]-Del)/F_Si[0])/log10(x[2]/F_Si[0]);
    double X_2=log10(x[2]+Del);
    double Y_2=x[3]+(x[5]-x[3])*log10(1.0+Del/x[2])/log10(x[4]/x[2]);
    double X_3=log10(x[4]-Del);
    double Y_3=x[2]+x[5]-x[3]*log10((x[4]-Del)/x[2])/log10(x[4]/x[2]);
    double X_4=log10(x[4]+Del);
    double Y_4=x[5]+(x[1]-x[5])*log10(1.0+Del/x[4])/log10(F_Si[999]/x[4]);
    double a1=-((x[5]-x[3])/log10(x[4]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
             pow(log10((x[2]+Del)/(x[2]-Del)),3.0)*log10(1.0-(Del/x[2])*(Del/x[2]));
    double b1=0.5*(((x[5]-x[3])/log10(x[4]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
             log10((x[2]+Del)/(x[2]-Del))-3*a1*log10(x[2]*x[2]-Del*Del));
    double c1=(x[3]-x[0])/log10(x[2]/F_Si[0])-3*a1*log10(x[2]-Del)*log10(x[2]-Del)-2*b1*log10(x[2]-Del);
    double d1=x[0]-log10(F_Si[0])*(x[3]-x[0])/log10(x[2]/F_Si[0])+0.5*((x[5]-x[3])/log10(x[4]/x[2])-
	         (x[3]-x[0])/log10(x[2]/F_Si[0]))/log10((x[2]+Del)/(x[2]-Del))*log10(x[2]-Del)*log10(x[2]-Del)
			  -0.5*a1*log10(x[2]-Del)*log10(x[2]-Del)*log10((x[2]+Del)*(x[2]+Del)*(x[2]+Del)/(x[2]-Del));
    double a2=-((x[1]-x[5])/log10(F_Si[999]/x[4])-(x[5]-x[3])/log10(x[4]/x[2]))/
              pow(log10((x[4]+Del)/(x[4]-Del)),3.0)*log10(1.0-(Del/x[4])*(Del/x[4]));
    double b2=0.5*(((x[1]-x[5])/log10(F_Si[999]/x[4])-(x[5]-x[3])/log10(x[4]/x[2]))/
              log10((x[4]+Del)/(x[4]-Del))-3*a2*log10(x[4]*x[4]-Del*Del));
    double c2=(x[5]-x[3])/log10(x[4]/x[2])-3*a2*log10(x[4]-Del)*log10(x[4]-Del)-2*b2*log10(x[4]-Del);
    double d2=x[3]-log10(x[2])*(x[5]-x[3])/log10(x[4]/x[2])+0.5*((x[1]-x[5])/log10(F_Si[999]/x[4])-
	          (x[5]-x[3])/log10(x[4]/x[2]))/log10((x[4]+Del)/(x[4]-Del))*log10(x[4]-Del)*log10(x[4]-Del)
	          -0.5*a2*log10(x[4]-Del)*log10(x[4]-Del)*log10((x[4]+Del)*(x[4]+Del)*(x[4]+Del)/(x[4]-Del));
	for (int i=0;i<1000;++i)
    {
    	if (log10(F_Si[i])<X_1)
    	{
    		tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
    		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_2)
        {
        	tau_inv_Si_temp[i]=pow(10,a1*pow(log10(F_Si[i]),3.0)+b1*pow(log10(F_Si[i]),2.0)+c1*log10(F_Si[i])+d1);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_3)
        {
        	tau_inv_Si_temp[i]=pow(10,(x[5]-x[3])/(log10(x[4])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_4)
        {
        	tau_inv_Si_temp[i]=pow(10,a2*pow(log10(F_Si[i]),3.0)+b2*pow(log10(F_Si[i]),2.0)+c2*log10(F_Si[i])+d2);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else
        {
        	tau_inv_Si_temp[i]=pow(10,(x[1]-x[5])/(log10(F_Si[999])-log10(x[4]))*(log10(F_Si[i])-log10(x[4]))+x[5]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
	}
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || X_2>X_3 || X_3>log10(x[4]) || log10(x[4])>X_4 || log10(x[4])>log10(F_Si[999]))
    {
    	y=y+1e5;
	}
	//**************************************************************************
    /*a_1=(x[9]-x[6])/(log10(x[8])-log10(F_Si[0]));
    b_1=-log10(F_Si[0])*a_1+x[6];
    a_2=(x[11]-x[9])/(log10(x[10])-log10(x[8]));
    b_2=-log10(x[8])*a_2+x[9];
    a_3=(x[7]-x[11])/(log10(F_Si[999])-log10(x[10]));
    b_3=-log10(x[10])*a_3+x[11];
    X_1=log10(x[8]-Del);
    Y_1=x[6]+x[9]-x[6]*log10((x[8]-Del)/F_Si[0])/log10(x[8]/F_Si[0]);
    X_2=log10(x[8]+Del);
    Y_2=x[9]+(x[11]-x[9])*log10(1.0+Del/x[8])/log10(x[10]/x[8]);
    X_3=log10(x[10]-Del);
    Y_3=x[8]+x[11]-x[9]*log10((x[10]-Del)/x[8])/log10(x[10]/x[8]);
    X_4=log10(x[10]+Del);
    Y_4=x[11]+(x[7]-x[11])*log10(1.0+Del/x[10])/log10(F_Si[999]/x[10]);
    a1=-((x[11]-x[9])/log10(x[10]/x[8])-(x[9]-x[6])/log10(x[8]/F_Si[0]))/
        pow(log10((x[8]+Del)/(x[8]-Del)),3.0)*log10(1.0-(Del/x[8])*(Del/x[8]));
    b1=0.5*(((x[11]-x[9])/log10(x[10]/x[8])-(x[9]-x[6])/log10(x[8]/F_Si[0]))/
        log10((x[8]+Del)/(x[8]-Del))-3*a1*log10(x[8]*x[8]-Del*Del));
    c1=(x[9]-x[6])/log10(x[8]/F_Si[0])-3*a1*log10(x[8]-Del)*log10(x[8]-Del)-2*b1*log10(x[8]-Del);
    d1=x[6]-log10(F_Si[0])*(x[9]-x[6])/log10(x[8]/F_Si[0])+0.5*((x[11]-x[9])/log10(x[10]/x[8])-
	    (x[9]-x[6])/log10(x[8]/F_Si[0]))/log10((x[8]+Del)/(x[8]-Del))*log10(x[8]-Del)*log10(x[8]-Del)
		 -0.5*a1*log10(x[8]-Del)*log10(x[8]-Del)*log10((x[8]+Del)*(x[8]+Del)*(x[8]+Del)/(x[8]-Del));
    a2=-((x[7]-x[11])/log10(F_Si[999]/x[10])-(x[11]-x[9])/log10(x[10]/x[8]))/
        pow(log10((x[10]+Del)/(x[10]-Del)),3.0)*log10(1.0-(Del/x[10])*(Del/x[10]));
    b2=0.5*(((x[7]-x[11])/log10(F_Si[999]/x[10])-(x[11]-x[9])/log10(x[10]/x[8]))/
        log10((x[10]+Del)/(x[10]-Del))-3*a2*log10(x[10]*x[10]-Del*Del));
    c2=(x[11]-x[9])/log10(x[10]/x[8])-3*a2*log10(x[10]-Del)*log10(x[10]-Del)-2*b2*log10(x[10]-Del);
    d2=x[9]-log10(x[8])*(x[11]-x[9])/log10(x[10]/x[8])+0.5*((x[7]-x[11])/log10(F_Si[999]/x[10])-
        (x[11]-x[9])/log10(x[10]/x[8]))/log10((x[10]+Del)/(x[10]-Del))*log10(x[10]-Del)*log10(x[10]-Del)
        -0.5*a2*log10(x[10]-Del)*log10(x[10]-Del)*log10((x[10]+Del)*(x[10]+Del)*(x[10]+Del)/(x[10]-Del));*/
	for (int i=1000;i<(1399);++i)
    {
    	if (log10(F_Si[i])<X_1)
    	{
    		tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
    		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_2)
        {
        	tau_inv_Si_temp[i]=pow(10,a1*pow(log10(F_Si[i]),3.0)+b1*pow(log10(F_Si[i]),2.0)+c1*log10(F_Si[i])+d1);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_3)
        {
        	tau_inv_Si_temp[i]=pow(10,(x[5]-x[3])/(log10(x[4])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_4)
        {
        	tau_inv_Si_temp[i]=pow(10,a2*pow(log10(F_Si[i]),3.0)+b2*pow(log10(F_Si[i]),2.0)+c2*log10(F_Si[i])+d2);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else
        {
        	tau_inv_Si_temp[i]=pow(10,(x[1]-x[5])/(log10(F_Si[999])-log10(x[4]))*(log10(F_Si[i])-log10(x[4]))+x[5]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
	}
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || X_2>X_3 || X_3>log10(x[4]) || log10(x[4])>X_4 || log10(x[4])>log10(F_Si[999]))
    {
    	y=y+1e5;
	}
	/*for (int i=1000;i<(1399);++i)
    {
    	if (log10(F_Si[i])<X_1)
    	{
    		tau_inv_Si_temp[i]=pow(10,(x[9]-x[0])/(log10(x[8])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[6]);
    		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_2)
        {
        	tau_inv_Si_temp[i]=pow(10,a1*pow(log10(F_Si[i]),3.0)+b1*pow(log10(F_Si[i]),2.0)+c1*log10(F_Si[i])+d1);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_3)
        {
        	tau_inv_Si_temp[i]=pow(10,(x[11]-x[9])/(log10(x[10])-log10(x[8]))*(log10(F_Si[i])-log10(x[8]))+x[9]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_4)
        {
        	tau_inv_Si_temp[i]=pow(10,a2*pow(log10(F_Si[i]),3.0)+b2*pow(log10(F_Si[i]),2.0)+c2*log10(F_Si[i])+d2);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else
        {
        	tau_inv_Si_temp[i]=pow(10,(x[7]-x[11])/(log10(F_Si[999])-log10(x[10]))*(log10(F_Si[i])-log10(x[10]))+x[11]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
	}
    if (log10(F_Si[0])>x[8] || X_1>log10(x[8]) || log10(x[8])>X_2 || X_2>X_3 || X_3>log10(x[10]) || log10(x[10])>X_4 || log10(x[10])>log10(F_Si[1398]))
    {
    	y=y+1e5;
    }*/

    //**************************************************************************
    //**************************************************************************
	
    double de_dT_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0;i <Nmode_Si_last-Nmode_Si_first+1;++i)
    {
        de_dT_Si[i] = pow(hbar*F_Si[i]/Teq,2)/boltz*exp(hbar*F_Si[i]/boltz/Teq)/pow(exp(hbar*F_Si[i]/boltz/Teq)-1,2);
    }
    double de_dT_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0;i <Nmode_Al_last-Nmode_Al_first+1;++i)
    {
        de_dT_Al[i] = pow(hbar*F_Al[i]/Teq,2)/boltz*exp(hbar*F_Al[i]/boltz/Teq)/pow(exp(hbar*F_Al[i]/boltz/Teq)-1,2);
    }
	double Lambda_Si[Nmode_Si_last-Nmode_Si_first+1];
    double Lambda_Si_max=0.0;
    for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    { 
        Lambda_Si[i]=V_Si[i]/tau_inv_Si[i];
        if (Lambda_Si[i]>Lambda_Si_max)
        {
        Lambda_Si_max=Lambda_Si[i];
        }
    }
    double Lambda_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    { 
        Lambda_Al[i]=V_Al[i]/tau_inv_Al[i];
    }
    double C_Si[Nmode_Si_last-Nmode_Si_first+1];
    double C_Si_ave=0.0;
    for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    { 
        C_Si[i]=SD_Si[i]*de_dT_Si[i]*Dom_Si[i];
        C_Si_ave=C_Si_ave+C_Si[i];
    }
    double C_Al[Nmode_Al_last-Nmode_Al_first+1];
    double C_Al_ave=0.0;
    for (int i=0;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    { 
        C_Al[i]=SD_Al[i]*de_dT_Al[i]*Dom_Al[i];
        C_Al_ave=C_Al_ave+C_Al[i];
    }
    double C_Al_Si_ratio;
    C_Al_Si_ratio=C_Al_ave/C_Si_ave;
    //------------------------------------------------------------------------------
    //calculating the objective function: 

    int L_num={8};    //number of length-scales in the TTG experiments:
    double **tt;
    tt=new double*[L_num];
    for(int i=0;i<L_num;i++) {tt[i] = new double[Ntt];}
    for(int j=0;j<Ntt;j++)
    { 
        tt[0][j]=5.0*j*1e-11;
        tt[1][j]=5.0*j*1e-11;
        tt[2][j]=5.0*j*1e-11;
        tt[3][j]=5.0*j*1e-11;
        tt[4][j]=5.0*j*1e-11;
        tt[5][j]=5.0*j*1e-11;
        tt[6][j]=5.0*j*1e-11;
        tt[7][j]=5.0*j*1e-11;
    }

    //length-scales:
    double D_val[L_num];
    D_val[0]=10e-9; D_val[1]=50e-9; D_val[2]=100e-9; 
    D_val[3]=500e-9; D_val[4]=1e-6; D_val[5]=5e-6; 
    D_val[6]=10e-6; D_val[7]=50e-6; //D_val[8]=100e-6;
    double L_val[L_num];
    for (int i=0;i<L_num;++i) {L_val[i]=D_val[i]*2.0; }
    double L_Al = 100.0e-9;
    double L_Si = 200e9;

    //number of MC particles for each length-scale:
    int N[L_num];
    N[0]=20000; N[1]=20000; N[2]=20000; N[3]=20000; N[4]=20000;
    N[5]=20000; N[6]=20000; N[7]=20000; //N[8]=10000;

    //forward BTE simulations using adjoint MC method:
    double T_eff[L_num];
    for (int i=0;i<L_num;i++) { T_eff[i]= (double) 1.0/N[i];}
    double k_Si=0.0;
    for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    {
        k_Si = k_Si+SD_Si[i]*Dom_Si[i]*V_Si[i]*V_Si[i]/tau_inv_Si[i]*de_dT_Si[i]/3;
    }
    double k_Al=0.0;
    for (int i=0;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    {
        k_Al = k_Al+SD_Al[i]*Dom_Al[i]*V_Al[i]*V_Al[i]/tau_inv_Al[i]*de_dT_Al[i]/3;
    }
    double cumul_base_Si[Nmode_Si_last-Nmode_Si_first+1];
    cumul_base_Si[0] = SD_Si[0]*de_dT_Si[0]*Dom_Si[0];
    double cumul_base_Al[Nmode_Al_last-Nmode_Al_first+1];
    cumul_base_Al[0] = SD_Al[0]*de_dT_Al[0]*Dom_Al[0];
    double cumul_coll_Si[Nmode_Si_last-Nmode_Si_first+1];
    cumul_coll_Si[0] = SD_Si[0]*de_dT_Si[0]*Dom_Si[0]*tau_inv_Si[0];
    double cumul_coll_Al[Nmode_Al_last-Nmode_Al_first+1];
    cumul_coll_Al[0] = SD_Al[0]*de_dT_Al[0]*Dom_Al[0]*tau_inv_Al[0];
    for (int i=1;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    {
        cumul_base_Al[i] = cumul_base_Al[i-1]+SD_Al[i]*de_dT_Al[i]*Dom_Al[i];
        cumul_coll_Al[i] = cumul_coll_Al[i-1]+SD_Al[i]*de_dT_Al[i]*tau_inv_Al[i]*Dom_Al[i];
    }
    for (int i=1;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    {
        cumul_base_Si[i] = cumul_base_Si[i-1]+SD_Si[i]*de_dT_Si[i]*Dom_Si[i];
        cumul_coll_Si[i] = cumul_coll_Si[i-1]+SD_Si[i]*de_dT_Si[i]*tau_inv_Si[i]*Dom_Si[i];
    }
    double RAND1=uni();
    double T21[Nmode_Si_last-Nmode_Si_first+1];
    int ind_shift=6;
    //**************************************************************************
    //**************************************************************************
    //interface transmissivity functions:
    
    //ORDER 0 T21:
    double T21_temp[Nmode_Si_last-Nmode_Si_first+1];
    /*for (int i=0; i < Nmode_Si_LA_TA-Nmode_Si_first+1; ++i)
    {
        T21_temp[i]=x[0+ind_shift];
	    T21[i]=T21_temp[i]; 
		if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
    }
    for (int i=(Nmode_Si_LA_TA-Nmode_Si_first+1);
	    i < (Nmode_Si_last-Nmode_Si_first+1); ++i)
    {
    	T21_temp[i]=x[0+ind_shift];
	    T21[i]=T21_temp[i]; 
		if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
    }*/
    
    //**************************************************************************
    //**************************************************************************
    
    //ORDER 1 T21:
    for (int i=0; i < Nmode_Si_LA_TA-Nmode_Si_first+1; ++i)
    {
        T21_temp[i]=(x[1+ind_shift]-x[0+ind_shift])/(F_Si[Nmode_Si_LA_TA-Nmode_Si_first]-
            F_Si[0])*(F_Si[i]-F_Si[0])+x[0+ind_shift];
	    T21[i]=T21_temp[i]; 
		if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
    }
    for (int i=(Nmode_Si_LA_TA-Nmode_Si_first+1);
	    i < (Nmode_Si_last-Nmode_Si_first+1); ++i)
    {
    	T21_temp[i]=(x[1+ind_shift]-x[0+ind_shift])/(F_Si[Nmode_Si_LA_TA-Nmode_Si_first]-
            F_Si[0])*(F_Si[i]-F_Si[0])+x[0+ind_shift];
	    T21[i]=T21_temp[i]; 
		if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
    }
    
    //**************************************************************************
    //**************************************************************************
    
    //ORDER 2 T21:    
    /*double a_1T=(x[3+ind_shift]-x[0+ind_shift])/(x[2+ind_shift]-F_Si[0]);
    double b_1T=-F_Si[0]*a_1T+x[0+ind_shift];
    double a_2T=(x[1+ind_shift]-x[3+ind_shift])/(F_Si[999]-x[2+ind_shift]);
    double b_2T=-x[2+ind_shift]*a_2T+x[3+ind_shift];
    double X_1T=x[2+ind_shift]-Del;
    double Y_1T=X_1T*a_1T+b_1T;
    double X_2T=x[2+ind_shift]+Del;
    double Y_2T=X_2T*a_2T+b_2T;
    double aT=((a_1T-a_2T)*x[2+ind_shift]+b_1T-b_2T)/4.0/Del/Del/Del;
    double bT=(a_2T-a_1T)/4.0/Del-3.0*aT*x[2+ind_shift];
    double cT=a_1T-3.0*aT*(x[2+ind_shift]-Del)*(x[2+ind_shift]-Del)-2.0*bT*(x[2+ind_shift]-Del);
    double dT=a_2T*(x[2+ind_shift]+Del)+b_2T-aT*(x[2+ind_shift]+Del)*(x[2+ind_shift]+Del)*
	(x[2+ind_shift]+Del)-bT*(x[2+ind_shift]+Del)*(x[2+ind_shift]+Del)-cT*(x[2+ind_shift]+Del);
    for (int i=0; i < Nmode_Si_LA_TA-Nmode_Si_first+1; ++i)
    {
    	if (F_Si[i]<X_1T)
	    {
            T21_temp[i]=(x[3+ind_shift]-x[0+ind_shift])/(x[2+ind_shift]-
			F_Si[0])*(F_Si[i]-F_Si[0])+x[0+ind_shift];
	        T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
        else if (F_Si[i]<X_2T)
        {
            T21_temp[i]=aT*F_Si[i]*F_Si[i]*F_Si[i]+bT*F_Si[i]*F_Si[i]+cT*F_Si[i]+dT;
            T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
        else
        {
            T21_temp[i]=(x[1+ind_shift]-x[3+ind_shift])/(F_Si[999]-
			x[2+ind_shift])*(F_Si[i]-x[2+ind_shift])+x[3+ind_shift];
            T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
    }
    if (F_Si[0]>x[2+ind_shift] || X_1T>x[2+ind_shift] || 
	x[2+ind_shift]>X_2T || x[2+ind_shift]>F_Si[999])
    {
    	y=y+1e5;
    }
    //**************************************************************************
    /*double a_1T=(x[7+ind_shift]-x[4+ind_shift])/(x[6+ind_shift]-F_Si[0]);
    double b_1T=-F_Si[0]*a_1T+x[4+ind_shift];
    double a_2T=(x[5+ind_shift]-x[7+ind_shift])/(F_Si[999]-x[6+ind_shift]);
    double b_2T=-x[6+ind_shift]*a_2T+x[7+ind_shift];
    double X_1T=x[6+ind_shift]-Del;
    double Y_1T=X_1T*a_1T+b_1T;
    double X_2T=x[6+ind_shift]+Del;
    double Y_2T=X_2T*a_2T+b_2T;
    double aT=((a_1T-a_2T)*x[6+ind_shift]+b_1T-b_2T)/4.0/Del/Del/Del;
    double bT=(a_2T-a_1T)/4.0/Del-3.0*aT*x[6+ind_shift];
    double cT=a_1T-3.0*aT*(x[6+ind_shift]-Del)*(x[6+ind_shift]-Del)-
	2.0*bT*(x[6+ind_shift]-Del);
    double dT=a_2T*(x[6+ind_shift]+Del)+b_2T-aT*(x[6+ind_shift]+Del)*(x[6+ind_shift]+Del)*
	(x[6+ind_shift]+Del)-bT*(x[6+ind_shift]+Del)*(x[6+ind_shift]+Del)-cT*(x[6+ind_shift]+Del);*/
    /*for (int i=(Nmode_Si_LA_TA-Nmode_Si_first+1);
	    i < (Nmode_Si_last-Nmode_Si_first+1); ++i)
    {
    	if (F_Si[i]<X_1T)
	    {
            T21_temp[i]=(x[3+ind_shift]-x[0+ind_shift])/(x[2+ind_shift]-
			F_Si[0])*(F_Si[i]-F_Si[0])+x[0+ind_shift];
	        T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
        else if (F_Si[i]<X_2T)
        {
            T21_temp[i]=aT*F_Si[i]*F_Si[i]*F_Si[i]+bT*F_Si[i]*F_Si[i]+cT*F_Si[i]+dT;
            T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
        else
        {
            T21_temp[i]=(x[1+ind_shift]-x[3+ind_shift])/(F_Si[999]-
			x[2+ind_shift])*(F_Si[i]-x[2+ind_shift])+x[3+ind_shift];
            T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
    }
    if (F_Si[0]>x[2+ind_shift] || X_1T>x[2+ind_shift] || 
	x[2+ind_shift]>X_2T || x[2+ind_shift]>F_Si[999])
    {
    	y=y+1e5;
	}*/
	/*for (int i=(Nmode_Si_LA_TA-Nmode_Si_first+1);
	    i < (Nmode_Si_last-Nmode_Si_first+1); ++i)
    {
    	if (F_Si[i]<X_1T)
	    {
            T21_temp[i]=(x[7+ind_shift]-x[4+ind_shift])/(x[6+ind_shift]-
			F_Si[0])*(F_Si[i]-F_Si[0])+x[4+ind_shift];
	        T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
        else if (F_Si[i]<X_2T)
        {
            T21_temp[i]=aT*F_Si[i]*F_Si[i]*F_Si[i]+bT*F_Si[i]*F_Si[i]+cT*F_Si[i]+dT;
            T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
        else
        {
            T21_temp[i]=(x[5+ind_shift]-x[7+ind_shift])/(F_Si[999]-
			x[6+ind_shift])*(F_Si[i]-x[6+ind_shift])+x[7+ind_shift];
            T21[i]=T21_temp[i]; 
		    if (T21_temp[i]<0){T21[i]=0;} if (T21_temp[i]>1.0){T21[i]=1.0;}
        }
    }
    if (F_Si[0]>x[6+ind_shift] || X_1T>x[6+ind_shift] || x[6+ind_shift]>X_2T ||
	 x[6+ind_shift]>F_Si[1398])
    {
    	y=y+1e5;
    }*/
    
    //************************************************************************** 
    //************************************************************************** 
    int pos_num=1;
    double pos_x[1]={0.0e-9};
    double pos_y[1]={0.0e-9};
    double pos_z[1]={0.0e-9};
    double **T_grid;
    T_grid=new double*[L_num];
    for(int i=0;i<L_num;i++) {T_grid[i]=new double[Ntt];}
    for(int i=0;i<L_num;i++) {for (int j=0;j<Ntt;j++) {T_grid[i][j]=0.0;}}
    for (int L_ind=0;L_ind<L_num;++L_ind)
    {
    	for (int pos_ind=0;pos_ind<pos_num;++pos_ind)
        {
            for (int counter=1;counter<=N[L_ind];++counter)
            {
                //double x0=pos_x[pos_ind];
                //double y0=pos_y[pos_ind];
                double x0=-D_val[L_ind]/2.0+D_val[L_ind]*uni();
                double y0=-D_val[L_ind]/2.0+D_val[L_ind]*uni();
                double z0=pos_z[pos_ind];
                //traject<<x0<<" "<<y0<<" "<<z0<<"\n";
                double t0 = 0.0;  
                //--------------------------------------------------------------
	        double Rx=uni();
                int i1 = 0;
                int i3;
                int i2;
                if (z0<=L_Al)                                  
                {         
                    i3 = Nmode_Al_last-Nmode_Al_first+1;
                    i2 = (int) floor((i1+i3)/2.0);
                    while (i3-i1>1)
                    {
                    if (Rx<cumul_base_Al[i2-1]/cumul_base_Al[Nmode_Al_last-Nmode_Al_first])
                    {
                    i3  = i2;
                    i2 = (int) floor((i1+i3)/2.0);
                    }
                    else
                    {
                    i1 = i2;
                    i2 = (int) floor((i1+i3)/2.0);
                    }
                    } 
                }
                else
                {         
                    i3 = Nmode_Si_last-Nmode_Si_first+1;
                    i2 = (int) floor((i1+i3)/2.0);
                    while (i3-i1>1)
                    {
                    if (Rx<cumul_base_Si[i2-1]/cumul_base_Si[Nmode_Si_last-Nmode_Si_first])
                    {
                    i3  = i2;
                    i2 = (int) floor((i1+i3)/2.0);
                    }
                    else
                    {
                    i1 = i2;
                    i2 = (int) floor((i1+i3)/2.0);
                    }
                    } 
                }
                int ind_mod=i3-1;                                 
                //double R_angle=PI*uni();
                double R_angle=acos(1.0-2.0*uni());
			    double phi_angle=2.0*PI*uni();
                double v_x;
			    double v_y;
			    double v_z;
                if (z0<=L_Al)                         
                {
			    v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
			    v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
			    v_z=V_Al[ind_mod]*cos(R_angle); 
			    }
                else
                {
			    v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
			    v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
			    v_z=V_Si[ind_mod]*cos(R_angle); 
		     	}
                int finished = 0;
                int im =0;
                double Delta_t;   
                if (z0<=L_Al)         
                { Delta_t = -tau_Al[ind_mod]*log(1.0-uni()); }
                else
                { Delta_t = -tau_Si[ind_mod]*log(1.0-uni()); }
                int part_state;
                if (z0<=L_Al)
                { part_state=11; }
                else
                { part_state=22; }
                while (finished==0)
                {
            	   double x1=x0+Delta_t*v_x;
            	   double y1=y0+Delta_t*v_y;
                    double z1=z0+Delta_t*v_z;
                    double t1=t0+Delta_t;
                    double p10x=-z0*(x1-x0)/(z1-z0)+x0;
                    double p10y=-z0*(y1-y0)/(z1-z0)+y0;
                    double p12x=(L_Al-z0)*(x1-x0)/(z1-z0)+x0;
                    double p12y=(L_Al-z0)*(y1-y0)/(z1-z0)+y0;
                    double p13y=(D_val[L_ind]/2.0-x0)*(y1-y0)/(x1-x0)+y0;
                    double p13z=(D_val[L_ind]/2.0-x0)*(z1-z0)/(x1-x0)+z0;
                    double p14y=(-D_val[L_ind]/2.0-x0)*(y1-y0)/(x1-x0)+y0;
                    double p14z=(-D_val[L_ind]/2.0-x0)*(z1-z0)/(x1-x0)+z0;
                    double p15x=(D_val[L_ind]/2.0-y0)*(x1-x0)/(y1-y0)+x0;
                    double p15z=(D_val[L_ind]/2.0-y0)*(z1-z0)/(y1-y0)+z0;
                    double p16x=(-D_val[L_ind]/2.0-y0)*(x1-x0)/(y1-y0)+x0;
                    double p16z=(-D_val[L_ind]/2.0-y0)*(z1-z0)/(y1-y0)+z0;
                    double p23y=(L_val[L_ind]/2.0-x0)*(y1-y0)/(x1-x0)+y0;
                    double p23z=(L_val[L_ind]/2.0-x0)*(z1-z0)/(x1-x0)+z0;
                    double p24y=(-L_val[L_ind]/2.0-x0)*(y1-y0)/(x1-x0)+y0;
                    double p24z=(-L_val[L_ind]/2.0-x0)*(z1-z0)/(x1-x0)+z0;
                    double p25x=(L_val[L_ind]/2.0-y0)*(x1-x0)/(y1-y0)+x0;
                    double p25z=(L_val[L_ind]/2.0-y0)*(z1-z0)/(y1-y0)+z0;
                    double p26x=(-L_val[L_ind]/2.0-y0)*(x1-x0)/(y1-y0)+x0;
                    double p26z=(-L_val[L_ind]/2.0-y0)*(z1-z0)/(y1-y0)+z0;
                    double p27x=(L_Al+L_Si-z0)*(x1-x0)/(z1-z0)+x0;
                    double p27y=(L_Al+L_Si-z0)*(y1-y0)/(z1-z0)+y0;
                    if (part_state==11 && z0>=0.0 && z1<0.0 && p10x>=-D_val[L_ind]/2.0 && 
                            p10x<=D_val[L_ind]/2.0 && p10y>=-D_val[L_ind]/2.0 && p10y<=D_val[L_ind]/2.0)
                    {
                   	    z1=0.0;
                	    x1=p10x;
                	    y1=p10y;
                        t1=(z1-z0)/v_z+t0;
                        part_state=10;
                    }
                    else if (part_state==11 && z0<=L_Al && z1>L_Al && p12x>=-D_val[L_ind]/2.0 &&
					        p12x<=D_val[L_ind]/2.0 && p12y>=-D_val[L_ind]/2.0 && p12y<=D_val[L_ind]/2.0)
                    {
                	    z1=L_Al;
                	    x1=p12x;
                	    y1=p12y;
                        t1=(z1-z0)/v_z+t0;
                        part_state=12; 
                    }
                    else if (part_state==11 && x0<=D_val[L_ind]/2.0 && x1>D_val[L_ind]/2.0 &&
					        p13y>=-D_val[L_ind]/2.0 && p13y<=D_val[L_ind]/2.0 && p13z>=0.0 && p13z<=L_Al)
                    {
                	    x1=D_val[L_ind]/2.0;
                	    y1=p13y;
                	    z1=p13z;
                        t1=(x1-x0)/v_x+t0;
                        part_state=13; 
                    }
                    else if (part_state==11 && x0>=-D_val[L_ind]/2.0 && x1<-D_val[L_ind]/2.0 &&
					        p14y>=-D_val[L_ind]/2.0 && p14y<=D_val[L_ind]/2.0 && p14z>=0.0 && p14z<=L_Al)
                    {
                    	x1=-D_val[L_ind]/2.0;
                    	y1=p14y;
                	    z1=p14z;
                        t1=(x1-x0)/v_x+t0;
                        part_state=14; 
                    }
                    else if (part_state==11 && y0<=D_val[L_ind]/2.0 && y1>D_val[L_ind]/2.0 &&
					        p15x>=-D_val[L_ind]/2.0 && p15x<=D_val[L_ind]/2.0 && p15z>=0.0 && p15z<=L_Al)
                    {
                	    y1=D_val[L_ind]/2.0;
                	    x1=p15x;
                	    z1=p15z;
                        t1=(y1-y0)/v_y+t0;
                        part_state=15; 
                    }
                    else if (part_state==11 && y0>=-D_val[L_ind]/2.0 && y1<-D_val[L_ind]/2.0 &&
					        p16x>=-D_val[L_ind]/2.0 && p16x<=D_val[L_ind]/2.0 && p16z>=0.0 && p16z<=L_Al)
                    {
                    	y1=-D_val[L_ind]/2.0;
                    	x1=p16x;
                    	z1=p16z;
                        t1=(y1-y0)/v_y+t0;
                        part_state=16; 
                    }
                    else if (part_state==22 && z0>=L_Al && z1<L_Al && 
					        p12x>=-D_val[L_ind]/2.0 && p12x<=D_val[L_ind]/2.0 && p12y>=-D_val[L_ind]/2.0 && p12y<=D_val[L_ind]/2.0)
                    {
                    	z1=L_Al;
                    	x1=p12x;
                    	y1=p12y;
                        t1=(z1-z0)/v_z+t0;
                        part_state=21; 
                    }
                    else if (part_state==22 && z0>=L_Al && z1<L_Al && 
					        p12x>=-L_val[L_ind]/2.0 && p12x<=L_val[L_ind]/2.0 && p12y>=-L_val[L_ind]/2.0 && p12y<=L_val[L_ind]/2.0
				            && (p12x<-D_val[L_ind]/2.0 || p12x>D_val[L_ind]/2.0 || p12y<-D_val[L_ind]/2.0 || p12y>D_val[L_ind]/2.0))
                    {
                	    z1=L_Al;
                	    x1=p12x;
                	    y1=p12y;
                        t1=(z1-z0)/v_z+t0;
                        part_state=20; 
                    }
                    else if (part_state==22 && x0<=L_val[L_ind]/2.0 && x1>L_val[L_ind]/2.0 && 
					        p23y>=-L_val[L_ind]/2.0 && p23y<=L_val[L_ind]/2.0 && p23z>=L_Al && p23z<=(L_Al+L_Si))
                    {
                    	x1=L_val[L_ind]/2.0;
                    	y1=p23y;
                	    z1=p23z;
                        t1=(x1-x0)/v_x+t0;
                        part_state=23; 
                    }
                    else if (part_state==22 && x0>=-L_val[L_ind]/2.0 && x1<-L_val[L_ind]/2.0 && 
					        p24y>=-L_val[L_ind]/2.0 && p24y<=L_val[L_ind]/2.0 && p24z>=L_Al && p24z<=(L_Al+L_Si))
                    {
                    	x1=-L_val[L_ind]/2.0;
                	    y1=p24y;
                	    z1=p24z;
                        t1=(x1-x0)/v_x+t0;
                        part_state=24; 
                    }
                    else if (part_state==22 && y0<=L_val[L_ind]/2.0 && y1>L_val[L_ind]/2.0 && 
					        p25x>=-L_val[L_ind]/2.0 && p25x<=L_val[L_ind]/2.0 && p25z>=L_Al && p25z<=(L_Al+L_Si))
                    {
                	    y1=L_val[L_ind]/2.0;
                	    x1=p25x;
                	    z1=p25z;
                        t1=(y1-y0)/v_y+t0;
                        part_state=25; 
                    }
                    else if (part_state==22 && y0>=-L_val[L_ind]/2.0 && y1<-L_val[L_ind]/2.0 &&
					        p26x>=-L_val[L_ind]/2.0 && p26x<=L_val[L_ind]/2.0 && p26z>=L_Al && p26z<=(L_Al+L_Si))
                    {
                    	y1=-L_val[L_ind]/2.0;
                	    x1=p26x;
                	    z1=p26z;
                        t1=(y1-y0)/v_y+t0;
                        part_state=26; 
                    }
                    /*else if (part_state==22 && z0<=(L_Al+L_Si) && z1>(L_Al+L_Si) && p27x>=-L/2.0 && p27x<=L/2.0 && p27y>=-L/2.0 && p27y<=L/2.0)
                    {
                    	z1=L_Al+L_Si;
                	    x1=p27x;
                	    y1=p27y;
                        t1=(z1-z0)/v_z+t0;
                        part_state=27; 
                    }*/
                    //**********************************************************
                    int traj_exists=0;
                    int start_im=im;
                    while ((im+1)<=Ntt && tt[L_ind][im] < t1)
                    {
                        im = im+1;
                        traj_exists=1;
                    }
                    if (traj_exists==1)
                    {
                        double Zpos[im-start_im];
                        for (int j=0;j<im-start_im;++j)
                        {
                            Zpos[j] = z0 + v_z*(tt[L_ind][j+start_im]-t0);
                            if (Zpos[j]<=L_Al && part_state!=21 && part_state!=20 && part_state!=22)
                            {  
                            T_grid[L_ind][j+start_im]= T_grid[L_ind][j+start_im]+ (double)T_eff[L_ind];
                            }
                        }
                    }
                    //**********************************************************
                    if (part_state==11)
                    {
                        double Rx=uni();
                        int i1 = 0;
                        int i3 = Nmode_Al_last-Nmode_Al_first+1;
                        int i2 = (int) floor((i1+i3)/2.0);
                        while (i3-i1>1)
                        {
                        if (Rx<cumul_coll_Al[i2-1]/cumul_coll_Al[Nmode_Al_last-Nmode_Al_first])
                        {
                            i3  = i2;
                            i2 = (int) floor((i1+i3)/2.0);
                        }
                        else
                        {
                            i1 = i2;
                            i2 = (int) floor((i1+i3)/2.0);
                        }
                        }
                        ind_mod=i3-1;
                        double R_angle=acos(1.0-2.0*uni());
                        double phi_angle=2.0*PI*uni();
                        v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
			            v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
			            v_z=V_Al[ind_mod]*cos(R_angle);
					    x0 = x1;
					    y0 = y1; 
		                z0 = z1;                             
                        t0 = t1;                                
                        Delta_t = -tau_Al[ind_mod]*log(1.0-uni());
                    }
                    else if (part_state==22)
                    {
                        double Rx=uni();
                        int i1 = 0;
                        int i3 = Nmode_Si_last-Nmode_Si_first+1;
                        int i2 = (int) floor((i1+i3)/2.0);
                        while (i3-i1>1)
                        {
                        if (Rx<cumul_coll_Si[i2-1]/cumul_coll_Si[Nmode_Si_last-Nmode_Si_first])
                        {
                            i3  = i2;
                            i2 = (int) floor((i1+i3)/2.0);
                        }
                        else
                        {
                            i1 = i2;
                            i2 = (int) floor((i1+i3)/2.0);
                        }
                        }
                        ind_mod=i3-1;
                        double R_angle=acos(1.0-2.0*uni());
                        double phi_angle=2.0*PI*uni();
                        v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
			            v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
			            v_z=V_Si[ind_mod]*cos(R_angle);
					    x0 = x1;
					    y0 = y1; 
		                z0 = z1;                             
                        t0 = t1;                                
                        Delta_t = -tau_Si[ind_mod]*log(1-uni());
                    }
                    if (part_state==12)
                       {
                        if(ind_mod<832-Nmode_Si_first+1 || (ind_mod>=1000-Nmode_Si_first+1 && ind_mod<1399-Nmode_Si_first+1))
                        {
                            if (uni() < (SD_Si[ind_mod]*V_Si[ind_mod])/(SD_Al[ind_mod]*V_Al[ind_mod])*T21[ind_mod])
                            {
                                double R_angle=asin(sqrt(uni()));
                                double phi_angle=2.0*PI*uni();
                                v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
			                    v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
			                    v_z=V_Si[ind_mod]*cos(R_angle);
			                    x0 = x1;
					            y0 = y1; 
		                        z0 = z1;
                                t0 = t1; 
                                Delta_t = -tau_Si[ind_mod]*log(1.0-uni());
                                part_state=22;
                            }
                            else
                            {
                        	    double R_angle=asin(sqrt(uni()));
                                double phi_angle=2.0*PI*uni();
                                v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
			                    v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
			                    v_z=-V_Al[ind_mod]*cos(R_angle);
			                    Delta_t=Delta_t-(t1-t0);
			                    x0 = x1;
					            y0 = y1; 
		                        z0 = z1;
                                t0 = t1;  
                                part_state=11;
                            }
                        }
                        else
                        {
                    	    double R_angle=asin(sqrt(uni()));
                            double phi_angle=2.0*PI*uni();
                            v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                            v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                            v_z=-V_Al[ind_mod]*cos(R_angle);
						    Delta_t=Delta_t-(t1-t0);
                            x0 = x1;
                            y0 = y1; 
                            z0 = z1;
                            t0 = t1;  
                            part_state=11;
                        }
                    }
                    else if (part_state==21)
                    {
                        if (ind_mod<832-Nmode_Si_first+1 || (ind_mod>=1000-Nmode_Si_first+1 && ind_mod<1399-Nmode_Si_first+1))
                        {
                            if (uni() < T21[ind_mod])
                            { 
                                double R_angle=asin(sqrt(uni()));
                                double phi_angle=2.0*PI*uni();
                                v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
			                    v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
			                    v_z=-V_Al[ind_mod]*cos(R_angle);
			                    x0 = x1;
					            y0 = y1; 
		                        z0 = z1;
                                t0 = t1; 
                                Delta_t = -tau_Al[ind_mod]*log(1.0-uni());
                                part_state=11;
                            }
                            else
                            {
                        	    double R_angle=asin(sqrt(uni()));
							    double phi_angle=2.0*PI*uni();
                                v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
			                    v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
			                    v_z=V_Si[ind_mod]*cos(R_angle);
			                    Delta_t=Delta_t-(t1-t0);
			                    x0 = x1;
					            y0 = y1; 
		                        z0 = z1;
                                t0 = t1; 
                                part_state=22;
                            }   
                        }
                        else
                        {
                    	    double R_angle=asin(sqrt(uni()));
						    double phi_angle=2.0*PI*uni();
                            v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
                            v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
                            v_z=V_Si[ind_mod]*cos(R_angle);
                            Delta_t=Delta_t-(t1-t0);
                            x0 = x1;
                            y0 = y1; 
                            z0 = z1;
                            t0 = t1;  
                            part_state=22;
                        }
                    }
                    else if (part_state==10)
                    {
                    	double R_angle=asin(sqrt(uni()));
                        double phi_angle=2.0*PI*uni();
                        v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        v_z=V_Al[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=11;
                    }
                    else if (part_state==13)
                    {
                    	double R_angle=asin(sqrt(uni()));
				    	double phi_angle=2.0*PI*uni();
                        v_z=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        v_x=-V_Al[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=11;
                    }
                    else if (part_state==14)
                    {
                    	double R_angle=asin(sqrt(uni()));
                        double phi_angle=2.0*PI*uni();
                        v_z=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        v_x=V_Al[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=11;
                    }
                    else if (part_state==15)
                    {
                    	double R_angle=asin(sqrt(uni()));
                        double phi_angle=2.0*PI*uni();
                        v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        v_z=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        v_y=-V_Al[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=11;
                    }
                    else if (part_state==16)
                    {
                    	double R_angle=asin(sqrt(uni()));
                        double phi_angle=2.0*PI*uni();
				        v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        v_z=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        v_y=V_Al[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=11;
                    }
                    else if (part_state==20)
                    {
                    	double R_angle=asin(sqrt(uni()));
                        double phi_angle=2.0*PI*uni();
                        v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
                        v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
                        v_z=V_Si[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=22;
                    }
                    else if (part_state==23)
                    {
                    	//double R_angle=PI*uni();
                        //double R_angle=asin(sqrt(uni()));
				    	//double phi_angle=PI/2.0+PI*uni();
					    //double phi_angle=2.0*PI*uni();
                        //v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_z=V_Al[ind_mod]*cos(R_angle);
                        //v_z=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_x=-V_Si[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1-L_val[L_ind];
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=22;
                    }
                    else if (part_state==24)
                    {
                    	//double R_angle=PI*uni();
                     	//double R_angle=asin(sqrt(uni()));
                        //double phi_angle=-PI/2.0+PI*uni();
                        //double phi_angle=2.0*PI*uni();
                        //v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_z=V_Al[ind_mod]*cos(R_angle);
                        //v_z=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_x=V_Si[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1+L_val[L_ind];
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=22;
                    }
                    else if (part_state==25)
                    {
                     	//double R_angle=PI*uni();
                	    //double R_angle=asin(sqrt(uni()));
                        //double phi_angle=PI+PI*uni();
                        //double phi_angle=2.0*PI*uni();
                        //v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_z=V_Al[ind_mod]*cos(R_angle);
                        //v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_z=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_y=-V_Si[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1-L_val[L_ind]; 
                        z0 = z1;
                        t0 = t1;
                        part_state=22;
                    }
                    else if (part_state==26)
                    {
                    	//double R_angle=PI*uni();
                    	//double R_angle=asin(sqrt(uni()));
                        //double phi_angle=PI*uni();
                        //double phi_angle=2.0*PI*uni();
					    //v_x=V_Al[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_y=V_Al[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_z=V_Al[ind_mod]*cos(R_angle);
                        //v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
                        //v_z=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_y=V_Si[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1+L_val[L_ind]; 
                        z0 = z1;
                        t0 = t1;
                        part_state=22;
                    }
                    /*else if (part_state==27)
                    {
                    	//double R_angle=PI-PI*uni()/2.0;
                        double R_angle=asin(sqrt(uni()));
					    double phi_angle=2.0*PI*uni();
                        v_x=V_Si[ind_mod]*sin(R_angle)*cos(phi_angle);
                        v_y=V_Si[ind_mod]*sin(R_angle)*sin(phi_angle);
                        //v_z=V_Al[ind_mod]*cos(R_angle);
                        v_z=-V_Si[ind_mod]*cos(R_angle);
                        Delta_t=Delta_t-(t1-t0);
                        x0 = x1;
                        y0 = y1; 
                        z0 = z1;
                        t0 = t1;
                        part_state=22;
                    }*/
                    if (t0>tt[L_ind][Ntt-1])
                    {
                    finished = 1;
                    }      
                    //traject<<x0<<" "<<y0<<" "<<z0<<"\n";
                } 
            }    
        }
    }
    //reading temperatures for each length-scale:
    ifstream file_res_10n("Al_Si_bkwrd_10n20n100n_SURF0n_5n_periodic_noisy.txt");
    ifstream file_res_50n("Al_Si_bkwrd_50n100n100n_SURF0n_5n_periodic_noisy.txt");
    ifstream file_res_100n("Al_Si_bkwrd_100n200n100n_SURF0n_5n_periodic_noisy.txt");
    ifstream file_res_500n("Al_Si_bkwrd_500n1m100n_SURF0n_5n_periodic_noisy.txt");
    ifstream file_res_1m("Al_Si_bkwrd_1m2m100n_SURF0n_5n_periodic_noisy.txt");
    ifstream file_res_5m("Al_Si_bkwrd_5m10m100n_SURF0n_5n_periodic_noisy.txt");
    ifstream file_res_10m("Al_Si_bkwrd_10m20m100n_SURF0n_5n_periodic_noisy.txt");
    ifstream file_res_50m("Al_Si_bkwrd_50m100m100n_SURF0n_5n_periodic_noisy.txt");

    double res[L_num][Ntt+3][2];
    if(file_res_10n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_10n >> res[0][j][i];}} file_res_10n.close();
    } 
    if(file_res_50n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_50n >> res[1][j][i];}} file_res_50n.close();
    }
    if(file_res_100n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_100n >> res[2][j][i];}} file_res_100n.close();
    }
    if(file_res_500n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_500n >> res[3][j][i];}} file_res_500n.close();
    }
    if(file_res_1m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_1m >> res[4][j][i];}} file_res_1m.close();
    }
    if(file_res_5m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_5m >> res[5][j][i];}} file_res_5m.close();
    }
    if(file_res_10m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_10m >> res[6][j][i];}} file_res_10m.close();
    }
    if(file_res_50m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<(Ntt+3);++j)
    {
        file_res_50m >> res[7][j][i];}} file_res_50m.close();
    }
    
    //calculating the objective function:
    for (int i=0;i<L_num;i++) {for (int j=0;j<Ntt;++j)
    {
        y=y+(double) abs(T_grid[i][j]-res[i][j+3][1])/800.0;
    }}
    y=y+0.1*abs(1.438416377372402E+02-k_Si)/1.438416377372402E+02;
    for (int i=0;i<L_num;++i) {delete [] tt[i];} delete [] tt;
    for (int i=0;i<L_num;++i) {delete [] T_grid[i];} delete [] T_grid;
    for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;i++) 
	{ if (tau_inv_Si_temp[i]<=0.0) {y=y+1e5;}}
	//for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;i++) 
	//{ if (T21_temp[i]<0.0 || T21_temp[i]>1.0) {y=y+1e5;}}
    return y;
}
