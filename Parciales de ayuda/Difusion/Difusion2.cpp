#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

const int Lx=256;
const int Ly=64;

const int Q=9;

const double tau=0.55;
const double Utau=1.0/tau;
const double UmUtau= 1-Utau;

const double mu=Lx/2;
const double sigma=Ly/9;




//-------------------------------------------------Clase LatticeBoltzman-------------------------------------------------
class LatticeBoltzman{
private:
    double w[Q];      //pesos
    int Vx[Q], Vy[Q]; //Vectores velocidad
    double *f, *fnew; //Funciones de distribucion
public:
    LatticeBoltzman(void);
    ~LatticeBoltzman(void);
    int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew);
    double Jy(int ix, int iy, bool UseNew);
    double feq(double rho0, double Ux0, double Uy0, int i);
    void Inicie(double mu, double sigma, double Ux0, double Uy0);
    void Colision(void);
    void Adveccion(void);
    void Print(const char * NameFile);
    double Sigma2(void);
};

LatticeBoltzman::LatticeBoltzman(void){
    //Cargar los pesos
    w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
    //Cargar los vectores velocidad
    Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;
    Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;
             Vx[5]=1; Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
             Vy[5]=1; Vy[6]=1; Vy[7]=-1; Vy[8]=-1;
             
    //Crear los arreglos dinamicos
    int ArraySize=Lx*Ly*Q;
    f=new double [ArraySize]; fnew=new double [ArraySize];
}

LatticeBoltzman::~LatticeBoltzman(void){
    delete[] f; delete[] fnew;
}

double LatticeBoltzman::rho(int ix, int iy, bool UseNew){
    double sum; int i, n0;
    for(sum=0, i=0;i<Q;i++){
        n0=n(ix,iy,i);
        if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
    }
    return sum;
}

double LatticeBoltzman::Jx(int ix, int iy, bool UseNew){
    
    return 0;
}

double LatticeBoltzman::Jy(int ix, int iy, bool UseNew){
    
    return 0;
}

double LatticeBoltzman::feq(double rho0, double Ux0, double Uy0, int i){
    double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
    return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

void LatticeBoltzman::Inicie(double mu, double sigma, double Ux0, double Uy0){
    int ix, iy, i, n0; double rho0;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++)
            for(i=0;i<Q;i++){ //En cada direccion
                rho0=(1/(sigma*sqrt(2*M_PI)))*exp(-0.5*((ix-mu)/sigma)*((ix-mu)/sigma));;
                n0=n(ix,iy,i);
                f[n0]=feq(rho0,Ux0,Uy0,i);
            }  
}

void LatticeBoltzman::Colision(void){
    int ix, iy, i, n0; double rho0, Ux0, Uy0;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
            for(i=0;i<Q;i++){ //En cada direccion
                n0=n(ix,iy,i);
                fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
            }
        }
}


void LatticeBoltzman::Adveccion(void){
    int ix, iy, i, ixnext, iynext, n0, n0next;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++)
            for(i=0;i<Q;i++){ //En cada direccion
                ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
                n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
                f[n0next]=fnew[n0];
            } 
}

double LatticeBoltzman::Sigma2(void){
    int ix, iy; double rho0, xAux, sAux; double xprom=0; double N=0; double sigma2=0;
    for(ix=0;ix<Lx;ix++)      //Calcule N
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,true);
            N+=rho0;
        }
    
    for(ix=0;ix<Lx;ix++)      //Calcule x promedio
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,true);
            xAux=(rho0*ix)/N;
            xprom+=xAux;
        }

    for(ix=0;ix<Lx;ix++)      //Calcule la varianza
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,true);
            sAux=(rho0*(ix-xprom)*(ix-xprom))/N;
            sigma2+=sAux;
        }
    
    return sigma2;
    

    
}

void LatticeBoltzman::Print(const char * NameFile){
    ofstream MyFile(NameFile); double rho0; int ix, iy;
    
    for(ix=0;ix<Lx;ix++){
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,true); 
            MyFile<<ix<<" "<<Sigma2()<<endl;
        }
        MyFile<<endl;
    }
    MyFile.close();
}

int main(void){
    LatticeBoltzman Aire;
    int t, tmax=1000;
    double rho0=1.0;

    ofstream MyFile("Varianza.dat");

    //INICIE
    Aire.Inicie(mu,sigma,0,0);
    
    //CORRA
    for(t=0;t<tmax;t++){
        Aire.Colision();
        Aire.Adveccion();

        MyFile<<t<<" "<<Aire.Sigma2()<<endl;
    }
    return 0;
}