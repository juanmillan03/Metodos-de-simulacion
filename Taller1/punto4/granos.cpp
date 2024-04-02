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

const double Lx=60,Ly=60;
const int Nx=5,Ny=5;
const int N=Nx*Ny,Ntot=N+4;
const double g=9.8;
const double Khertz=1.0e4;
const double Gamma=15, Kcundall=500,mu=0.4;


const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1.0-2*lambda)/2;
const double Um2chiplusxi=1.0-2*(chi*xi);
// ---------------------------------------
void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);
// --------------------------------------
class molecule; 
class colisionador;
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
    void CalculeTodasLasFuerzas(molecule*moleculas,double dt);   
    void CalculeFuerzaEntre( molecule & molecula1, molecule & molecula2,
    double & xCundall,double & sold,double dt);
    
};


int main(){
    //std::ofstream outfile;
    //outfile.open("granos_hertz.dat");

    molecule moleculas[Ntot];
    colisionador hertz;
    Crandom ran64(26);
    
    //Parametros de la simulación
    double m0=1.0, R0=2.0;
    double kT=10;
    //Variables auxiliares para la condición inicial
    double dx=Lx/(Nx+1),dy=Ly/(Ny+1);   
    double theta,V0=std::sqrt(kT/m0);
    double x0,y0,Vx0,Vy0;
    double omega0=3*M_PI;
    //Variables auxiliares para correr la simulacion
    int Ncuadros=1000;
    double t,tdibujo,dt=1e-3,tmax=Lx/V0,Tcuadro=tmax/Ncuadros;


    InicieAnimacion();
    // inicializar las paredes
    double Rpared=100*Lx,Mpared=100*m0;
    //                  (x0,y0,vx0,vy0,theta0,omega0,m0,r0)
    moleculas[N].inicie(Lx/2,Ly+Rpared,0,0,0,0,Mpared,Rpared);
    moleculas[N+1].inicie(Lx/2,-Rpared,0,0,0,0,Mpared,Rpared);
    moleculas[N+2].inicie(Lx+Rpared,Ly/2,0,0,0,0,Mpared,Rpared);
    moleculas[N+3].inicie(-Rpared,Ly/2,0,0,0,0,Mpared,Rpared);


    // inicializar los granos
    for (int i=0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            theta=2*M_PI*ran64.r();
            x0=(i+1)*dx;y0=(j+1)*dy;Vx0=V0*cos(theta);Vy0=V0*sin(theta);
            moleculas[j*Nx+i].inicie(x0,y0,Vx0,Vy0,0,omega0,m0,R0);
        }
        
    }

    for ( t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
        if(tdibujo>Tcuadro){
            InicieCuadro();
            for (int i = 0; i < N; i++)
            {
                moleculas[i].Dibujese();
            }
            TermineCuadro();
            tdibujo=0; 
        }
        // std::cout<<moleculas[1].getx()<<"  "<< moleculas[1].gety()<<"   "<<moleculas[2].getx()<<"    "<< moleculas[2].gety()<<std::endl;
        for (int i = 0; i < N; i++)moleculas[i].Mueva_r(dt,xi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt); 
        for (int i = 0; i < N; i++)moleculas[i].Mueva_V(dt,Um2lambdau2);
        for (int i = 0; i < N; i++)moleculas[i].Mueva_r(dt,chi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt);
        for (int i = 0; i < N; i++)moleculas[i].Mueva_V(dt,lambda);
        for (int i = 0; i < N; i++)moleculas[i].Mueva_r(dt,Um2chiplusxi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt);
        for (int i = 0; i < N; i++)moleculas[i].Mueva_V(dt,lambda);
        for (int i = 0; i < N; i++)moleculas[i].Mueva_r(dt,chi);
        hertz.CalculeTodasLasFuerzas(moleculas,dt);
        for (int i = 0; i < N; i++) moleculas[i].Mueva_V(dt,Um2lambdau2);
        for (int i = 0; i < N; i++)moleculas[i].Mueva_r(dt,xi); 
        
    }
    //outfile.close();
    return 0;
}




// -----------------------------FUncion molecule
void molecule::inicie(double x0,double y0,double Vx0,double Vy0,
        double theta0,double omega0, double m0,double R0){
    r.load(x0,y0,0);V.load(Vx0,Vy0,0);m=m0;R=R0;
    theta= theta0;omega=omega0;I=2.0/5.0*m*R*R;
}
void molecule::Mueva_r(double dt,double coeficiente){
    r+=V*(coeficiente*dt);theta+=omega*dt*coeficiente;
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
    for(i=0;i<Ntot;i++){
        for(j=0;j<Ntot;j++)
        xCundall[i][j]=sold[i][j]=0;
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
        // defino variables 
        vector3D F1,F2,tao1,tao2;
        //Vectores unitarios
        vector3D normal=dr/d,t,k; t.load(normal.y(),-normal.x(),0); k.load(0,0,1);
        
        // Velocidad de contacto
        vector3D Rw; Rw.load(0,0,R2*molecula2.omega+R1*molecula1.omega);
        vector3D Vc=(molecula2.V-molecula1.V)-(Rw^normal);
        double Vcn=(Vc*normal), Vct=Vc*t;
        // feurza herzt

        double Fn;
        // fuerza plasta Kuramoto-kano 
        double m1=molecula1.m,m2=molecula2.m; double M=m1*m2/(m1+m2);
        Fn=(-Gamma*std::sqrt(s)*M)*Vcn+Khertz*std::pow(s,1.5);
        // Fuerza tangencial de Cundall

        xCundall+=Vct*dt; double Ft=-Kcundall*xCundall; double Ftmax=mu*fabs(Fn);
        if(fabs(Ft)>Ftmax)Ft=Ft/fabs(Ft)*Ftmax;

        F2=normal*Fn+Ft*normal;tao2=((normal*(-R2))^F2);F1=F2*(-1);tao1=((normal*R1)^F1);
        molecula1.SumeFuerza(F2*(-1),tao1*k); molecula2.SumeFuerza(F2,tao2*k);

        }
        if(sold>=0 && s<0)xCundall=0;
        sold=s;
}

void colisionador::CalculeTodasLasFuerzas(molecule*moleculas,double dt){
    for(int i=0;i<N+4;i++){
        moleculas[i].BorreFuerza();
    }
    vector3D Fg;
    for(int i=0;i<N;i++){
        Fg.load(0,-moleculas[i].m*g,0);
        moleculas[i].SumeFuerza(Fg,0);
    }
    for(int i=0;i<N;i++){
        for (int j = i+1; j <N+4; j++)
        {
            CalculeFuerzaEntre(moleculas[i],moleculas[j],xCundall[i][j],sold[i][j],dt);

        }  
    }
}
void InicieAnimacion(void){
    // std::cout<<"set terminal gif animate"<<std::endl; 
    // std::cout<<"set output 'particulas.gif'"<<std::endl;
    std::cout<<"unset key"<<std::endl;
    std::cout<<"set grid"<<std::endl;
    std::cout<<"set xrange[-10:"<<Lx+10<<"]"<<std::endl;
    std::cout<<"set yrange[-10:"<<Ly+10<<"]"<<std::endl;
    std::cout<<"set size ratio -1"<<std::endl;
    std::cout<<"set parametric"<<std::endl;
    std::cout<<"set trange [0:7]"<<std::endl;
    std::cout<<"set isosamples 12"<<std::endl;  
}
void InicieCuadro(void){
    std::cout<<"plot 0,0 ";
    std::cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    std::cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    std::cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    std::cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    std::cout<<std::endl;
}

