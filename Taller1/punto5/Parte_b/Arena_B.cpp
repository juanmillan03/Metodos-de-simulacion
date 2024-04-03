/*
Fuerza de hertz deformacion entre particulas, que depende de la geometria y 
interpenetracíon 
Velocidad relativa de contacto, V traslacion y rotacional

*/
#include<iostream>
#include<fstream>
#include<cmath>
#include "vector.h"
#include "Random64.h"
#include <string>


const double Ly=60;
const int N=200,Ns=80,Ntot=N+Ns+3;
const double g=9.8;
const double Khertz=1.0e4;
const double Gamma=150, Kcundall=500,mu=0.4;


const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1.0-2*lambda)/2;
const double Um2chiplusxi=1.0-2*(chi*xi);

class molecule; 
class colisionador;
// ---------------------------------------
void InicieAnimacion(double Lx);
void InicieCuadro(double Lx);
void TermineCuadro(void);
// --------------------------------------

// -------------------------------------
class molecule{ 
    private: 
    vector3D r,V,F;
    double m,R;
    double theta,omega,tau,I;
    public:
    void inicie(double x0,double y0,double Vx0,double Vy0,
    double theta0,double omega0, double m0,double R0);
    void BorreFuerza(void){F.load(0,0,0);tau=0;}
    void SumeFuerza(vector3D dF,double dtau){F+=dF;tau+=dtau;};
    void Mueva_r(double dt,double coeficiente);
    void Mueva_V(double dt,double coeficiente);
    void Dibujese(void);
    double getx(void){return r.x();};
    double gety(void){return r.y();};
    friend class colisionador;
};
// -------------------------------------
class colisionador{
private:
    double xCundall[Ntot][Ntot],sold[Ntot][Ntot];

public:
    void Inicie(void);
    void CalculeTodasLasFuerzas(molecule*moleculas,double dt,int Nlive);   
    void CalculeFuerzaEntre( molecule & molecula1, molecule & molecula2,
    double & xCundall,double & sold,double dt);

};
int main(int argc, char **argv){
    if (argc!=2){
        std::cout<<"Usage: Meter el valor Lx "<<std::endl;
        return 1;
    }
    const double Lx = std::atof(argv[1]);
    std::ofstream outfile("mu_experimental.dat", std::ios::app);
    if (!outfile) {
        std::cerr << "No se pudo abrir el archivo: " << "mu_experimental.dat" << std::endl;
        return 1;
    }
    molecule moleculas[Ntot];
    colisionador hertz;
    Crandom ran64(26);
    
    //Parametros de la simulación
    double m0=1.0, R0=2.0;
    double kT=10;
    double theta,V0=std::sqrt(kT/m0);
    int j,Nlive;
    double Ncuadros=5,t,tdibujo,dt=1e-3,
        tmax=0.25*Ncuadros*std::sqrt(Ly/g),Tcuadro=190*tmax/(10*Ncuadros);
    double Omega0,OmegaMax=8.0;     
    double Rpared=100*Lx,Mpared=100*m0,Rs=Lx/(2*Ns);
 

    InicieAnimacion(Lx);
    for(int i=0;i<Ns;i++){
//   Grano[N+1].Inicie(Lx/2,  -Rpared,  0,  0,     0,     0,Mpared,Rpared); //Pared abajo
        moleculas[N+i].inicie(Rs*(2*i+1),0,0,0,0,0,Mpared,Rs);
    }
    //                  (x0,y0,vx0,vy0,theta0,omega0,m0,r0)
    moleculas[N+Ns].inicie(Lx/2,Ly+Rpared,0,0,0,0,Mpared,Rpared);
    moleculas[N+Ns+1].inicie(Lx+Rpared,Ly/2,0,0,0,0,Mpared,Rpared);
    moleculas[N+Ns+2].inicie(-Rpared,Ly/2,0,0,0,0,Mpared,Rpared);


    // inicializar los granos
    Nlive=1;
    double tlive;
    double x0=Lx/2,y0=Ly-2*R0;
    Omega0= OmegaMax*(2*ran64.r()-1);
    moleculas[0].inicie(x0,y0,0,0,0,Omega0,m0,1.6+0.8*ran64.r());
    for (t=tdibujo=0;t<tmax*210;t+=dt,tdibujo+=dt){
        if( t>tmax*Nlive&& Nlive<=N-1){
            Omega0= OmegaMax*(2*ran64.r()-1);
            double R=1.6+0.8*ran64.r();
            moleculas[Nlive].inicie(x0,y0,0,0,0,Omega0,m0,R);
            ++Nlive;
        }
        if(tdibujo>tmax*205){
            InicieCuadro(Lx);
            for (int i = 0; i<N+Ns; i++)moleculas[i].Dibujese();    
            TermineCuadro();
            tdibujo=0; 
        }
        // std::cout<<moleculas[1].getx()<<"  "<< moleculas[1].gety()<<"   "<<moleculas[2].getx()<<"    "<< moleculas[2].gety()<<std::endl;
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_r(dt,xi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt,Nlive); 
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_V(dt,Um2lambdau2);
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_r(dt,chi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt,Nlive);
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_V(dt,lambda);
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_r(dt,Um2chiplusxi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt,Nlive);
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_V(dt,lambda);
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_r(dt,chi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt,Nlive);
        for (int i = 0; i < Nlive; i++) moleculas[i].Mueva_V(dt,Um2lambdau2);
        for (int i = 0; i < Nlive; i++)moleculas[i].Mueva_r(dt,xi); 
        
    }
    double xmin,ymax,xmax,xymax,tang;
    xmin=moleculas[0].getx();
    xmax=moleculas[0].getx();
    ymax=moleculas[0].gety();
    for (int i = 1; i < Nlive; i++)
    {
        if (xmin>moleculas[i].getx()){xmin=moleculas[i].getx();}
        if (xmax<moleculas[i].getx()){xmax=moleculas[i].getx();}
        if (ymax<moleculas[i].gety()){ymax=moleculas[i].gety();xymax=moleculas[i].getx();}
    }
    tang=(ymax/std::fabs(Lx/2.0-xmin)+ymax/std::fabs(Lx/2.0-xmax))/2;
    outfile<<Lx<<" "<<tang<<" "<<Nlive<<std::endl;
    return 0;
}
// -----------------------------FUncion molecule
void molecule::inicie(double x0,double y0,double Vx0,double Vy0,
        double theta0,double omega0, double m0,double R0){
    r.load(x0,y0,0);V.load(Vx0,Vy0,0);m=m0;R=R0;
    theta= theta0;omega=omega0;I=2.0/5.0*m*R*R;
}
void molecule::Mueva_r(double dt,double coeficiente){
    r+=V*(coeficiente*dt);theta+=omega*(dt*coeficiente);
}
void molecule::Mueva_V(double dt,double coeficiente){
    V+=F*(coeficiente*dt/m);omega+=tau*(coeficiente*dt/I);
}
void molecule::Dibujese(void){
    std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

// ---------------------funciones colisionador

void colisionador::Inicie(void){
    int i,j; //j>i
    for(i=0;i<Ntot;i++)
        for(j=0;j<Ntot;j++)
            xCundall[i][j]=sold[i][j]=0;
}


void colisionador::CalculeTodasLasFuerzas(molecule*moleculas,double dt,int Nlive){
    for(int i=N;i<N+Ns+3;i++)moleculas[i].BorreFuerza();
    for (int i = 0; i<Nlive; i++){
        moleculas[i].BorreFuerza();
    }
    vector3D Fg;
    for(int i=0;i<Nlive;i++){
        Fg.load(0,-moleculas[i].m*g,0);
        moleculas[i].SumeFuerza(Fg,0);
    }
    for(int i=0;i<N;i++){
        for (int j = i+1; j<N+Ns+3; j++)
            if (i<=Nlive){
                CalculeFuerzaEntre(moleculas[i],moleculas[j],xCundall[i][j],sold[i][j],dt);  
            }
            else continue; 
    }           
            
            
    
}
void colisionador::CalculeFuerzaEntre(molecule & molecula1, molecule & molecula2,
        double & xCundall,double & sold,double dt){

    // calcular el vector normal 
    vector3D dr=molecula2.r-molecula1.r;double d=dr.norm();
    double R1=molecula1.R,R2=molecula2.R;
    double s=R1+R2-d;

    // existe contacto ???
    if (s>0){
        //Vectores unitarios
        vector3D n=dr*(1.0)/d,t,k; t.load(n.y(),-n.x(),0); k.load(0,0,1);
        
        // Velocidad de contacto
        vector3D Rw; Rw.load(0,0,R2*molecula2.omega+R1*molecula1.omega);
        vector3D Vc=(molecula2.V-molecula1.V)-(Rw^n);
        double Vcn=(Vc*n), Vct=Vc*t;
        // feu  rza herzt

        // fuerza plasta Kuramoto-kano 
        double m1=molecula1.m, m2=molecula2.m; double M=m1*m2/(m1+m2);
        double Fn=-Gamma*std::sqrt(s)*M*Vcn+Khertz*std::pow(s,1.5);
        // Fuerza tangencial de Cundall

        xCundall+=Vct*dt; double Ft=-Kcundall*xCundall; double Ftmax=mu*std::fabs(Fn);
        if(std::fabs(Ft)>Ftmax) Ft=Ft/std::fabs(Ft)*Ftmax;
        
        vector3D F1,F2,tau1,tau2;
        F2=n*Fn+Ft*t;tau2=((n*(-R2))^F2);F1=F2*(-1);tau1=((n*R1)^F1);
        molecula1.SumeFuerza(F2*(-1),tau1*k); molecula2.SumeFuerza(F2,tau2*k);

        }
        if(sold>=0 && s<0)xCundall=0;
        sold=s;
}
void InicieAnimacion(double Lx){
    // std::cout<<"set terminal gif animate"<<std::endl; 
    std::cout<<"set terminal pdf"<<std::endl; 
    std::cout<<"set output 'gif/Arena"<<Lx<<".pdf'"<<std::endl;   
    std::cout<<"unset key"<<std::endl;
    std::cout<<"set grid"<<std::endl;
    std::cout<<"set xrange[-10:"<<Lx+10<<"]"<<std::endl;
    std::cout<<"set yrange[-10:"<<Ly+10<<"]"<<std::endl;
    std::cout<<"set size ratio -1"<<std::endl;
    std::cout<<"set parametric"<<std::endl;
    std::cout<<"set trange [0:7]"<<std::endl;
    std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(double Lx){
    std::cout<<"plot 0,0 ";
    std::cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    std::cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    std::cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    std::cout<<std::endl;
}

