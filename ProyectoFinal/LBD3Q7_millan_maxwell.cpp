#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <omp.h>
using namespace std;

const int Lx=320;
const int Ly=320;
const int Lz=320;



const int Q=7;
const double W0=1.0/4;

const double C=0.5;
const double C2=C*C;
const double AUX0=1-4*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau= 1-Utau;

const double D = 0.7;




//-------------------------------------------------Clase LatticeBoltzman-------------------------------------------------
class LatticeBoltzman{
private:
    double w[Q];      //Pesos
    int Vx[Q], Vy[Q], Vz[Q]; //Vectores velocidad
    double *f, *fnew; //Funciones de distribucion
public:
    LatticeBoltzman(void);
    ~LatticeBoltzman(void);
    int n(int ix, int iy, int iz, int i){return (ix*Ly*Lz+iy*Lz+iz)*Q+i;};
    double rho(int ix, int iy, int iz, bool UseNew);
    double Jx(int ix, int iy, int iz, bool UseNew);
    double Jy(int ix, int iy, int iz, bool UseNew);
    double Jz(int ix, int iy, int iz, bool UseNew);
    double feq(double rho0, double Jx0, double Jy0, double Jz0, int i);
    void Inicie(double rho0, double Jx0, double Jy0, double Jz0);
    void Colision(void);
    void ImponerCampos(int t);
    void Adveccion(void);
    void Print(const char * NameFile,int z);
};
LatticeBoltzman::LatticeBoltzman(void){
    //Cargar los pesos
    w[0]=W0; w[1]=w[2]=w[3]=w[4]=w[5]=w[6]=1.0/8;
    //Cargar los vectores velocidad
    Vx[0]=0; Vx[1]=1; Vx[2]=-1; Vx[3]=0; Vx[4]=0;  Vx[5]=0; Vx[6]=0;
    Vy[0]=0; Vy[1]=0; Vy[2]=0;  Vy[3]=1; Vy[4]=-1; Vy[5]=0; Vy[6]=0;
    Vz[0]=0; Vz[1]=0; Vz[2]=0;  Vz[3]=0; Vz[4]=0;  Vz[5]=1; Vz[6]=-1;
    //Crear los arreglos dinamicos
    int ArraySize=Lx*Ly*Lz*Q;
    f=new double [ArraySize]; fnew=new double [ArraySize];
}
LatticeBoltzman::~LatticeBoltzman(void){
    delete[] f; delete[] fnew;
}
double LatticeBoltzman::rho(int ix, int iy, int iz, bool UseNew){
    double sum; int i, n0;
    for(sum=0, i=0;i<Q;i++){
        n0=n(ix,iy,iz,i);
        if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
    }
    return sum;
}
double LatticeBoltzman::Jx(int ix, int iy, int iz, bool UseNew){
    double sum; int i, n0;
    for(sum=0, i=0;i<Q;i++){
        n0=n(ix,iy,iz,i);
        if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
    }
    return sum;
}
double LatticeBoltzman::Jy(int ix, int iy, int iz, bool UseNew){
    double sum; int i, n0;
    for(sum=0, i=0;i<Q;i++){
        n0=n(ix,iy,iz,i);
        if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
    }
    return sum;
}
double LatticeBoltzman::Jz(int ix, int iy, int iz, bool UseNew){
    double sum; int i, n0;
    for(sum=0, i=0;i<Q;i++){
        n0=n(ix,iy,iz,i);
        if(UseNew) sum+=Vz[i]*fnew[n0]; else sum+=Vz[i]*f[n0];
    }
    return sum;
}
double LatticeBoltzman::feq(double rho0, double Jx0, double Jy0, double Jz0, int i){
    if(i>0)
        return 4*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0+Vz[i]*Jz0);
    else
        return rho0*AUX0;

}
void LatticeBoltzman::Inicie(double rho0, double Jx0, double Jy0, double Jz0){
    int ix, iy, iz, i, n0;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++)
            for(iz=0;iz<Lz;iz++)
                for(i=0;i<Q;i++){ //En cada direccion
                    n0=n(ix,iy,iz,i);
                    f[n0]=feq(rho0,Jx0,Jy0,Jz0,i);
                }  
}
void LatticeBoltzman::Colision(void){
    int ix, iy, iz, i, n0, n1, n2, n3, n4, n5, n6; double rho0, Jx0, Jy0, Jz0;
    #pragma omp parallel for collapse(3) private(ix, iy, iz, i, n0, n1, n2, n3, n4, n5, n6, rho0, Jx0, Jy0, Jz0) 
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++)
            for(iz=0;iz<Lz;iz++){
                rho0=rho(ix,iy,iz,false); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false); 
                for(i=0;i<Q;i++){ //En cada direccion
                    //Bounce-Back
                if (ix==Lx-2 || ix==1 || iy==Ly-2 || iy==1 || iz==Lz-2 || iz==1){
                    n0 = n(ix, iy, iz, 0);
                    n1 = n(ix, iy, iz, 1);
                    n3 = n(ix, iy, iz, 3);
                    n2 = n(ix, iy, iz, 2);
                    n4 = n(ix, iy, iz, 4);
                    n5 = n(ix, iy, iz, 5);
                    n6 = n(ix, iy, iz, 6);

                    fnew[n0] = D*f[n0];
                    fnew[n1] = D*f[n2];
                    fnew[n2] = D*f[n1];
                    fnew[n3] = D*f[n4];
                    fnew[n4] = D*f[n3];
                    fnew[n5] = D*f[n6];
                    fnew[n6] = D*f[n5];
                }
                if (ix==Lx-1 || ix==0 || iy==Ly-1 || iy==0 || iz==Lz-1 || iz==0){
                    n0 = n(ix, iy, iz, 0);
                    n1 = n(ix, iy, iz, 1);
                    n3 = n(ix, iy, iz, 3);
                    n2 = n(ix, iy, iz, 2);
                    n4 = n(ix, iy, iz, 4);
                    n5 = n(ix, iy, iz, 5);
                    n6 = n(ix, iy, iz, 6);
                    fnew[n0] = 0;
                    fnew[n1] = 0;
                    fnew[n2] = 0;
                    fnew[n3] = 0;
                    fnew[n4] = 0;
                    fnew[n5] = 0;
                    fnew[n6] = 0;
                }
                else{
                    for(i=0;i<Q;i++){ //En cada direccion
                    n0=n(ix,iy,iz,i);
                    fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);
                    }
                }
                }
            }
}
void LatticeBoltzman::ImponerCampos(int t){
    int i, ix, iy, iz, n0,n1,n2;
    double lambda, omega, rho0, Jx0, Jy0, Jz0,rho1, Jx1, Jy1, Jz1,rho2, Jx2, Jy2, Jz2; lambda=13; omega=2*M_PI/lambda*C;
    //Una fuente oscilante en el medio
    ix=Lx/2; iy=Ly/2; iz=Lz/2;
    rho0=20*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
    rho1=20*sin(omega*t+M_PI/2); Jx1=Jx(ix,iy-iy/2,iz+iz/2,false); Jy1=Jy(ix,iy-iy/2,iz+iz/2,false); Jz1=Jz(ix,iy-iy/2,iz+iz/2,false);
    rho2=20*sin(omega*t-M_PI/2); Jx2=Jx(ix,iy+iy/2,iz-iz/2,false); Jy2=Jy(ix,iy+iy/2,iz-iz/2,false); Jz2=Jz(ix,iy+iy/2,iz-iz/2,false);
    for(i=0;i<Q;i++){
        n0=n(ix,iy,iz,i);
        fnew[n0]=feq(rho0,Jx0,Jy0,Jz0,i);
        n1=n(ix,iy-iy/2,iz+iz/2,i);
        fnew[n1]=feq(rho1,Jx1,Jy1,Jz1,i);
        n2=n(ix,iy+iy/2,iz-iz/2,i);
        fnew[n2]=feq(rho2,Jx2,Jy2,Jz2,i);
    }
}
void LatticeBoltzman::Adveccion(void){
    int ix, iy, iz, i, ixnext, iynext, iznext, n0, n0next;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++)
            for(iz=0;iz<Lz;iz++)
                for(i=0;i<Q;i++){ //En cada direccion
                    ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly; iznext=(iz+Vz[i]+Lz)%Lz;
                    n0=n(ix,iy,iz,i); n0next=n(ixnext,iynext,iznext,i);
                    f[n0next]=fnew[n0];
                } 
}
void LatticeBoltzman::Print(const char * NameFile,int z){
    ofstream MyFile(NameFile); double rho0; int ix, iy;
    int iz = z;
    for(ix=0;ix<Lx;ix++){
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,iz,true);
            MyFile<<(float)ix/10<<" "<<(float)iy/10<<" "<<rho0<<endl;
        }
        MyFile<<endl;
    }
    MyFile.close();
}
int main(void){
    LatticeBoltzman Ondas;
    int t, tmax=1000;
    double rho0=0, Jx0=0, Jy0=0, Jz0=0;

    //INICIE
    Ondas.Inicie(rho0,Jx0,Jy0,Jz0);
    
    //Correr
    // ...

    auto start = std::chrono::high_resolution_clock::now(); // Start timer

    for(t=0;t<tmax;t++){
        Ondas.Colision();
        Ondas.ImponerCampos(t);
        Ondas.Adveccion();
        if(t%20==0){
            #pragma omp for
            for(int z=Lz/4;z<Lz;z=z+Lz/4){
                char filename[30];
                sprintf(filename, "D3/%d/Ondas_t%d.txt", z,t);
                Ondas.Print(filename,z);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now(); // End timer
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // Calculate duration in milliseconds
    std::cout << "Total time for all iterations: " << duration.count()/1000 << " seconds" << std::endl;
    //Print
    
    return 0;
}