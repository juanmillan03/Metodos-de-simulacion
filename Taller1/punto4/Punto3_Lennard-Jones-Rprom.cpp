#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"
#include "Random64.h"
using namespace std;

//Constantes del problema físico
double Lx=60, Ly=60; 
const int Nx=5, Ny=5;
const int N=Nx*Ny;
const double G=1.0;
const double KHertz=1.0e4;
const double Epsilon = 1.0;
const double Sigma = 10.0;
const  double Rpared = 50;

//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);

//--------------- Declarar las clases-----------
class Cuerpo;
class Colisionador;

//--------- Declarar las interfases de las clases---------
class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};// Inline
  void SumeFuerza(vector3D dF){F+=dF;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  friend class Colisionador;
};
class Colisionador{
private:
  double Rgiro;
public:
  void CalculeTodasLasFuerzas(Cuerpo * Molecula);
  void CalculeFuerzaPared(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void CalculeRgiro(Cuerpo * Molecula);
  double GetRgiro(void){return Rgiro;}; // Inline
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}

//------- Funciones de la clase Colisionador --------


void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i,j;
  //Borro las fuerzas de todos las moléculas
  for(i=0;i<N+1;i++)
    Molecula[i]. BorreFuerza();
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=i+1;j<N+1;j++)
    if (j>=N){
      CalculeFuerzaPared(Molecula[i],Molecula[j]);
    }
    else{
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
    }
}

void Colisionador::CalculeFuerzaPared(Cuerpo & Molecula1,Cuerpo & Molecula2){
  vector3D r21=Molecula2.r-Molecula1.r; double d=r21.norm();
  double s=Molecula1.R-Molecula2.R+d;
  if(s>0){ //Si hay colisión
    //Calcular el vector normal
    vector3D n=r21*(1.0/d);
    //Calculo la fuerza
    vector3D F2=n*(KHertz*pow(s,1.5));
    //Las sumo a los granos
    Molecula2.SumeFuerza((-1)*F2);  Molecula1.SumeFuerza(F2);
  }

}


void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2){
  //Determinar si hay colision
  vector3D r21=Molecula2.r-Molecula1.r; double r=r21.norm();
  double s=(Molecula1.R+Molecula2.R)-r;
  //if(s>0){ //Si hay colisión
    //Calcular el vector normal
    vector3D n=r21*(1.0/r);
    //Calculo la fuerza
    vector3D F1=n*12*Epsilon*((pow((Sigma/r),12)) - (pow((Sigma/r),6)))/r;
    //Las sumo a las moléculas
    Molecula2.SumeFuerza(F1);  Molecula1.SumeFuerza(F1*(-1));
  //}
}

void Colisionador::CalculeRgiro(Cuerpo * Molecula){
  vector3D sum;
  sum.load(0,0,0);
  double giro = 0.0;

  for(int i=0;i<N;i++){
    sum = sum + Molecula[i].r;
  }

  vector3D Rcm = sum/N;

  for(int i=0;i<N;i++){
    vector3D temp = Molecula[i].r - Rcm;
    giro = giro + temp.norm2();
  }

  Rgiro = sqrt(giro/N);

}

//----------- Funciones Globales -----------

int main(){
  Cuerpo Molecula[N+1];
  Colisionador Newton;
  Crandom ran64(1);
  int i,ix,iy;
  //Parametros de la simulación
  double m0=1.0; double R0=3.0;
  double kT=0.05; 
  //Variables auxiliares para la condición inicial
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);
  double theta; double V0=sqrt(kT/m0);
  double x0,y0,Vx0,Vy0;

  double Mpared=100*m0;
  //Variables auxiliares para correr la simulacion
  double t,dt=5e-4,tmax=100; 
  


  //Inicializar las paredes
  //------------------(  x0,       y0,z0,Vx0,Vy0,Vz0,    m0,    R0)
  Molecula[N].Inicie(0,0, 0,  0,  0,  0,Mpared,Rpared);   //Pared circular  
  
  //INICIO
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      theta=2*M_PI*ran64.r();
      x0=(ix-2)*dx; y0=(iy-2)*dy; Vx0=V0*cos(theta); Vy0=V0*sin(theta);
      //----------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
      Molecula[iy*Nx+ix].Inicie(x0,y0, 0,Vx0,Vy0,  0,m0,R0);
      
    }
      

  //CORRO
  
  double sum = 0.0;
  double tot = 0.0;
  for(t=0;t<tmax;t+=dt){

    Newton.CalculeRgiro(Molecula);
    
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,xi);    
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Um2chiplusxi);
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Molecula); for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,xi);
    if (t>= 20.0){
      sum = sum + Newton.GetRgiro();
      tot = tot + 1.0;
    }
  }
  cout << sum/tot<< endl;
  
  return 0;
}
