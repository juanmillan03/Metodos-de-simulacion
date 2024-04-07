#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"
using namespace std;

//Constantes del problema físico
const double epsilon=1.0;
const double r_0 = 10.0;

//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);
//Deaclaración de la clase
class Cuerpo;

//Deaclaración de la interfase

class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
};

//Implementación de las funciones
void Cuerpo::Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(void){
  double norma=r.norm();
  double F1=12*epsilon*((pow((r_0/norma),12)) - (pow((r_0/norma),6)))/(norma*norma);
  F=r*F1;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}

//----------- Funciones Globales -----------


int main(){
  double x0=10.0,m0=1.0, R0=3.0;
  double kT=1.0; 
  double V0=sqrt(kT/m0);
  double t,dt=0.1,ttotal=100;
  Cuerpo Part;

  
  //----------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)
  Part.Inicie(x0, 0, 0,  V0, 0,  0, m0, R0);

  std::ofstream outfile;
  outfile.open("One_lennard.txt");

  for(t=0;t<ttotal;t+=dt){

    outfile<< t <<" "<< Part.Getx() <<endl;
    
    Part.Mueva_r(dt,xi);
    Part.CalculeFuerza(); Part.Mueva_V(dt,Um2lambdau2);
    Part.Mueva_r(dt,chi);
    Part.CalculeFuerza(); Part.Mueva_V(dt,lambda);
    Part.Mueva_r(dt,Um2chiplusxi);
    Part.CalculeFuerza(); Part.Mueva_V(dt,lambda);
    Part.Mueva_r(dt,chi);
    Part.CalculeFuerza(); Part.Mueva_V(dt,Um2lambdau2);
    Part.Mueva_r(dt,xi);
    
  }
  outfile.close();
  return 0;
}
