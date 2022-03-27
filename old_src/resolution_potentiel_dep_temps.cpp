#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <bibli_fonctions.h>
#include <Python.h>
#include <unistd.h>
#include <fstream>
#define K 350
using namespace std;


//constantes physiques//

complex<double> I(0.,1.);
double h_barre = 1.054e-34;
double V_0=-1.6e-9;
double tau=1e-15;
double a_0=1e-15;
double m=2000*9.11e-31;
double k_0=sqrt(2*m*1.9e-9)/h_barre;
double x_0=30e-15;
double x_min = -10e-15;
double x_max = 10e-15;

//données du problème//


int n = 100; //nombre d'itération en temps

double delta_t = 1e-18;
double delta_x = (x_max-x_min)/K; //de l'ordre de 1e-18; 

complex<double> alpha[1]={0. + I*( h_barre*0.5/m)};
complex<double>  beta[1] ={0. +I*( h_barre*0.25/(m*2*delta_x*delta_x))};

//données en fonction//

double normaliser(complex<double> phi[K]) {
  double norme= 0;
  for (int i=0;i<K;i++)
    norme+=(real(phi[i])*real(phi[i]) + imag(phi[i])*imag(phi[i]))*delta_x;
  return norme;
}

double*** v_fun(int k, int n) {  //pour l'instant considéré réel !!!
  double*** matrice = (double***)malloc(k*sizeof(double**));
  for (int i =0;i<k;i++) {
    matrice[i] = (double**)malloc(n*sizeof(double*));
    for (int j = 0;j<n;j++) {
      matrice[i][j]= (double*)malloc(2*sizeof(double));
      //on entre ici partie réelle et imaginaire de V/(ih_bare)
      matrice[i][j][0] = 0;  
      matrice[i][j][1] = 0;
    }
  }
  return matrice;
}

double*** lambda_fun (int k,int n) {
  double*** pot = v_fun(k,n);
  double*** matrice = (double***)malloc(k*sizeof(double**));
  for (int i =0;i<k;i++) {
     matrice[i] = (double**)malloc(n*sizeof(double*));
      for (int j = 0;j<n;j++) {
       matrice[i][j]= (double*)malloc(2*sizeof(double));
       matrice[i][j][0] = 1/delta_t-0.5*pot[i][j][0];
       matrice[i][j][1] = 2*imag(beta[0])-0.5*pot[i][j][1];
     }
  }
  return matrice;
}

double*** mu_fun (int k,int n) {
  double*** pot = v_fun(k,n);
  double*** matrice = (double***)malloc(k*sizeof(double**));
  for (int i =0;i<k;i++) {
     matrice[i] = (double**)malloc(n*sizeof(double*));
      for (int j = 0;j<n;j++) {
       matrice[i][j]= (double*)malloc(2*sizeof(double));
       matrice[i][j][0] = 1/delta_t+0.5*pot[i][j][0];
       matrice[i][j][1] = -2*imag(beta[0])+0.5*pot[i][j][1];
     }
  }
  return matrice;
}



int main () {
  complex<double> I(0.,1.);
  double*** mu =mu_fun(K,n);
  double*** v = v_fun(K,n);
  double* proba =(double*)malloc(K*sizeof(double));
  double norm;
  fstream fich;

  
  complex<double> phi_0[K];
  for (int i=0;i<K;i++)
    phi_0[i] =  pow(M_PI/2.,0.25)*(1/sqrt(a_0))*cos(k_0*(x_min+i*delta_x))*exp(-pow((x_0-(x_min+i*delta_x))/a_0,2)) + I*(pow( (M_PI/2.),0.25)*(1/sqrt(a_0))*exp(-pow((x_0-(x_min+i*delta_x))/a_0,2))*sin(k_0*(x_min+i*delta_x)));

  norm=normaliser(phi_0);
   for(int i =0;i<K;i++)
     phi_0[i]=phi_0[i]/sqrt(norm);

   for (int p=0;p<K;p++)
      proba[p]= (real(phi_0[p])*real(phi_0[p]) + imag(phi_0[p])*imag(phi_0[p]));


    fich.open("densite_tps.dat", ios::out);
    for (int q=0;q<K;q++)
      fich << proba[q] << " " << real(phi_0[q]) << " " << imag(phi_0[q])<< " " << q*delta_x << endl;
    fich.close();


  complex<double> a[K],c[K];
  for (int i =0;i<K;i++) {
    a[i]= -beta[1];
    c[i]= -beta[1];
  }
      Py_Initialize();
      PyRun_SimpleString("from numpy import *");
      PyRun_SimpleString("from matplotlib.pyplot import *");
      PyRun_SimpleString("fig=figure()");

    for (int i=0;i<n;i++) {

        //AFFICHAGE DE RE_PHI(X) EN FONCTION DU TEMPS
    
    PyRun_SimpleString("A=loadtxt('densite_tps.dat')");
    PyRun_SimpleString("real=A[:,1]");
    PyRun_SimpleString("imag=A[:,2]");
    PyRun_SimpleString("proba=A[:,0]");
    PyRun_SimpleString("x=A[:,3]");
    PyRun_SimpleString("plot(x,real)");
    PyRun_SimpleString("plot(x,imag)");
    PyRun_SimpleString("plot(x,proba)");
    PyRun_SimpleString("draw()");
    PyRun_SimpleString("pause(0.0001)");
    PyRun_SimpleString("fig.clear()");
      
      complex<double> b[K];
      for (int j=0;j<K;j++)
	b[j]= 1/delta_t -0.5*v[j][i][0] + I*(2*imag(beta[0] -0.5*v[j][i][1]));
      
   
      complex<double> y[K];
      for (int j=1;j<K-1;j++)
	y[j]= -imag(beta[0])*imag(phi_0[j-1]) +mu[j][i][0]*real(phi_0[j])-mu[j][i][1]*imag(phi_0[j]) - imag(beta[0])*imag(phi_0[j+1]) +I*(imag(beta[0])*real(phi_0[j-1]) +mu[j][i][1]*real(phi_0[j])+mu[j][i][0]*imag(phi_0[j]) + imag(beta[0])*real(phi_0[j+1]));
      
      y[0] = mu[0][i][0]*real(phi_0[0]) - mu[0][i][1]*imag(phi_0[0]) - imag(beta[0])*imag(phi_0[1])+ I*(mu[0][i][1]*real(phi_0[0]) +mu[0][i][0]*imag(phi_0[0])+imag(beta[0])*real(phi_0[1]));
      y[K-1] =mu[K-1][i][0]*real(phi_0[K-1])- mu[K-1][i][1]*imag(phi_0[K-1])- imag(beta[0])*imag(phi_0[K-2])+ I*(mu[K-1][i][1]*real(phi_0[K-1]) + mu[K-1][i][0]*imag(phi_0[K-1])+imag(beta[0])*real(phi_0[K-2]));
    complex<double> x[K];
    
    choleski_complex(a,b,c,x,y,K);
    

    norm=normaliser(x);

      for (int q=0;q<K;q++)
      phi_0[q]=x[q]/sqrt(norm);
    
    for (int p=0;p<K;p++)
      proba[p]= (real(phi_0[p])*real(phi_0[p]) + imag(phi_0[p])*imag(phi_0[p]));

     fich.open("densite_tps.dat", ios::out);
    for (int q=0;q<K;q++)
      fich << proba[q] << " " << real(phi_0[q]) << " " << imag(phi_0[q])<< " " << q*delta_x << endl;
    fich.close();
									  

    }
											 Py_Finalize();  
 
    free(proba);
    for (int q=0;q<K ;q++) {
      for(int i=0;i<n;i++) {
	free(v[q][i]);
	free(mu[q][i]);
      }
      free(v[q]);
      free(mu[q]);
    }
    free(v);
    free(mu);
  
  return 0;
}
