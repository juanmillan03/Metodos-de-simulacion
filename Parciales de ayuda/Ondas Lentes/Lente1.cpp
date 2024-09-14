#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

const int Lx=512;
const int Ly=128;

const int Q=5;
const double W0=1.0/3;

const double C=0.5;
const double C2=C*C;
const double AUX0=1-3*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau= 1-Utau;




//-------------------------------------------------Clase LatticeBoltzman-------------------------------------------------
class LatticeBoltzman{
private:
    double w[Q];      //Pesos
    int Vx[Q], Vy[Q]; //Vectores velocidad
    double *f, *fnew; //Funciones de distribucion
public:
    LatticeBoltzman(void);
    ~LatticeBoltzman(void);
    int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew);
    double Jy(int ix, int iy, bool UseNew);
    double feq(double rho0, double Jx0, double Jy0, int i);
    void Inicie(double rho0, double Jx0, double Jy0);
    void Colision(void);
    void ImponerCampos(int t);
    void Adveccion(void);
    void Print(const char * NameFile);
};

LatticeBoltzman::LatticeBoltzman(void){
    //Cargar los pesos
    w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
    //Cargar los vectores velocidad
    Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;
    Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;
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
    double sum; int i, n0;
    for(sum=0, i=0;i<Q;i++){
        n0=n(ix,iy,i);
        if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
    }
    return sum;
}

double LatticeBoltzman::Jy(int ix, int iy, bool UseNew){
    double sum; int i, n0;
    for(sum=0, i=0;i<Q;i++){
        n0=n(ix,iy,i);
        if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
    }
    return sum;
}

double LatticeBoltzman::feq(double rho0, double Jx0, double Jy0, int i){
    if(i>0)
        return 3*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
    else
        return rho0*AUX0;

}

void LatticeBoltzman::Inicie(double rho0, double Jx0, double Jy0){
    int ix, iy, i, n0;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++)
            for(i=0;i<Q;i++){ //En cada direccion
                n0=n(ix,iy,i);
                f[n0]=feq(rho0,Jx0,Jy0,i);
            }  
}

void LatticeBoltzman::Colision(void){
    int ix, iy, i, n0; double rho0, Jx0, Jy0;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
            for(i=0;i<Q;i++){ //En cada direccion
                n0=n(ix,iy,i);
                fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i);
            }
        }
}

void LatticeBoltzman::ImponerCampos(int t){
    int i, ix, iy, n0;
    double lambda, A, omega, rho0, Jx0, Jy0; A=10; lambda=10; omega=2*M_PI/lambda*C;
    //Una fuente oscilante en el medio
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++){
            if(ix==0){
                rho0=A*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
            for(i=0;i<Q;i++){
                n0=n(ix,iy,i);
                fnew[n0]=feq(rho0,Jx0,Jy0,i);
            }
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

void LatticeBoltzman::Print(const char * NameFile){
    ofstream MyFile(NameFile); double rho0; int ix, iy;
    for(ix=0;ix<Lx;ix++){
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,true);
            MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
        }
        MyFile<<endl;
    }
    MyFile.close();
}

int main(void){
    LatticeBoltzman Ondas;
    int t, tmax=512;
    double rho0=0, Jx0=0, Jy0=0;

    //INICIE
    Ondas.Inicie(rho0,Jx0,Jy0);
    
    //CORRA
    for(t=0;t<tmax;t++){
        Ondas.Colision();
        Ondas.ImponerCampos(t);
        Ondas.Adveccion();
    }
    //Print
    Ondas.Print("Lente1.dat");
    return 0;
}