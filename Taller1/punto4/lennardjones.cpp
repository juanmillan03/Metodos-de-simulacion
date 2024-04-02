#include<iostream>
#include<fstream>
#include<cmath>
#include "vector.h"
#include "Random64.h"

const double Lx=60,Ly=20;
const int Nx=1,Ny=1;
const int N=Nx*Ny,Ntot=N;
const double epsilon=1.0;
const double poso=10.0;


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
    public:
    void inicie(double x0,double y0,double Vx0,double Vy0, double m0,double R0);
    void BorreFuerza(void){F.load(0,0,0);};
    void SumeFuerza(vector3D dF){F+=dF;};
    void Mueva_r(double dt,double coeficiente);
    void Mueva_V(double dt,double coeficiente);
    void Dibujese(void);
    double getx(void){return r.x();};
    double gety(void){return r.y();};
    friend class colisionador;
};
// -------------------------------------
class colisionador{
public:
    void CalculeTodasLasFuerzas(molecule*moleculas,double dt);   
    void CalculeFuerzaEntre( molecule & molecula1,double dt);
    
};


int main(){
    //std::ofstream outfile;
    //outfile.open("granos_hertz.dat");

    molecule moleculas[Ntot];
    colisionador hertz;
    //Parametros de la simulación
    double m0=1.0, R0=3.0;
    double kT=1.0;
    //Variables auxiliares para la condición inicial
    double V0=std::sqrt(kT/m0);
    double x0,y0,Vx0,Vy0;
    //Variables auxiliares para correr la simulacion
    int Ncuadros=1000;
    double t,tdibujo,dt=1e-3,tmax=Lx/V0,Tcuadro=tmax/Ncuadros;


    InicieAnimacion();
    // inicializar los granos
    moleculas[0].inicie(20.0,0,0,0,m0,R0);


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
void molecule::inicie(double x0,double y0,double Vx0,double Vy0, double m0,double R0){
    r.load(x0,y0,0);V.load(Vx0,Vy0,0);m=m0;R=R0;
}
void molecule::Mueva_r(double dt,double coeficiente){
    r+=V*(coeficiente*dt);
}
void molecule::Mueva_V(double dt,double coeficiente){
    V+=F*(coeficiente*dt/m);
}
void molecule::Dibujese(void){
    std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

// ---------------------funciones colisionador


void colisionador::CalculeFuerzaEntre(molecule & molecula1,double dt){
    // calcular el vector normal 
    vector3D dr= molecula1.r;double d=dr.norm()+10;
    vector3D n=dr/d;
    vector3D F1=(12*epsilon/d)*(std::pow(12,poso/d)-std::pow(6,poso/d))*n;
    molecula1.SumeFuerza(F1*(-1));
}

void colisionador::CalculeTodasLasFuerzas(molecule*moleculas,double dt){
    for(int i=0;i<N;i++){
        moleculas[i].BorreFuerza();
        CalculeFuerzaEntre(moleculas[i],dt);
    }
}
void InicieAnimacion(void){
    // std::cout<<"set terminal gif animate"<<std::endl; 
    // std::cout<<"set output 'particulas.gif'"<<std::endl;
    std::cout<<"unset key"<<std::endl;
    std::cout<<"set grid"<<std::endl;
    std::cout<<"set xrange["<<0<<":"<<Lx<<"]"<<std::endl;
    std::cout<<"set yrange["<<-Ly<<":"<<Ly<<"]"<<std::endl;
    std::cout<<"set size ratio -1"<<std::endl;
    std::cout<<"set parametric"<<std::endl;
    std::cout<<"set trange [0:7]"<<std::endl;
    std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
    std::cout<<"plot 0,0";
    std::cout<<",[h=0.1:60] h,[h=0.1:60] 12/h*((10/h)**(12)-(10/h)**(6))";
    

//     std::cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
//     std::cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
//     std::cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
//     std::cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    std::cout<<std::endl;
}

