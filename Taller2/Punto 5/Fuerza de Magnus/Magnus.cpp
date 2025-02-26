#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

const int Lx=512;
const int Ly=64;
const int ixc=128;
const int iyc=32;
const double N=100;

const int Q=9;






//-------------------------------------------------Clase LatticeBoltzman-------------------------------------------------
class LatticeBoltzman{
private:
    double w[Q];      //pesos
    int Vx[Q], Vy[Q]; //Vectores velocidad
    double *f, *fnew; //Funciones de distribucion
    double tau;
    double Utau=1.0/tau;
    double UmUtau= 1-Utau;
    double nu=(1.0/3.0)*(tau-0.5);
public:
    LatticeBoltzman(double TAU);
    ~LatticeBoltzman(void);
    int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew);
    double Jy(int ix, int iy, bool UseNew);
    double feq(double rho0, double Ux0, double Uy0, int i);
    void Inicie(double rho0, double Ux0, double Uy0);
    void Colision(void);
    void ImponerCampos(double Ufan, double Omega);
    void Adveccion(void);
    double dUx(int ix, int iy);
    double dUy(int ix, int iy);
    double dUxy(int ix, int iy);
    double dUyx(int ix, int iy);
    double sigmaxx(int ix,int iy);
    double sigmayy(int ix,int iy);
    double sigmaxy(int ix,int iy);
    std::vector<double> Calcule_dF(double Px, double Py, double dAx, double dAy);
    std::vector<double> CalculeFuerza(double Ufan, double Omega);
    void Print(const char * NameFile, double Ufan);
};

LatticeBoltzman::LatticeBoltzman(double TAU){
    //Cargar los pesos
    tau=TAU;
    Utau=1.0/tau;
    UmUtau= 1-Utau;
    nu=(1.0/3.0)*(tau-0.5);
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

double LatticeBoltzman::feq(double rho0, double Ux0, double Uy0, int i){
    double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
    return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

void LatticeBoltzman::Inicie(double rho0, double Ux0, double Uy0){
    int ix, iy, i, n0;
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++)
            for(i=0;i<Q;i++){ //En cada direccion
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

void LatticeBoltzman::ImponerCampos(double Ufan, double Omega){
    int i, ix, iy, n0;
    double rho0; int R=8; double R2=R*R;
    double Ux0, Uy0;
    //Ir por todas las celdas y mirar si es ventilador o obstaculo
    for(ix=0;ix<Lx;ix++)      //Para cada celda
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,false);
            //Ventilador
            if(ix==0)
                for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ufan,0,i);}
            //Obstaculo
            else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2){
                Ux0=-Omega*(iy-iyc); Uy0=Omega*(ix-ixc);
                for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ux0,Uy0,i);}
            }
            // Un punto extra
            else if(ix==ixc && iy==iyc+R+1){
                Ux0=-Omega*(iy-iyc); Uy0=Omega*(ix-ixc);
                for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ux0,Uy0,i);}
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

double LatticeBoltzman::dUx(int ix, int iy){
    int i, ixnext, iynext; double rho0, Ux0next;
    double sum=0;
    
    for(i=0;i<Q;i++){
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
        rho0=rho(ixnext,iynext,false);
        Ux0next=Jx(ixnext,iynext,false)/rho0; 
        sum+=w[i]*Vx[i]*Ux0next;
    }
    return 3*sum;
}

double LatticeBoltzman::dUy(int ix, int iy){
    int i, ixnext, iynext; double rho0, Uy0next;
    double sum=0;
    
    for(i=0;i<Q;i++){
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
        rho0=rho(ixnext,iynext,false);
        Uy0next=Jy(ixnext,iynext,false)/rho0; 
        sum+=w[i]*Vy[i]*Uy0next;
    }
    return 3*sum;
}

double LatticeBoltzman::dUxy(int ix, int iy){
    int i, ixnext, iynext; double rho0, Ux0next;
    double sum=0;
    
    for(i=0;i<Q;i++){
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
        rho0=rho(ixnext,iynext,false);
        Ux0next=Jx(ixnext,iynext,false)/rho0; 
        sum+=w[i]*Vy[i]*Ux0next;
    }
    return 3*sum;
}

double LatticeBoltzman::dUyx(int ix, int iy){
    int i, ixnext, iynext; double rho0, Uy0next;
    double sum=0;
    
    for(i=0;i<Q;i++){
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
        rho0=rho(ixnext,iynext,false);
        Uy0next=Jy(ixnext,iynext,false)/rho0; 
        sum+=w[i]*Vx[i]*Uy0next;
    }
    return 3*sum;
}

double LatticeBoltzman::sigmaxx(int ix,int iy){
    double rho0, p; 
    
    rho0=rho(ix,iy,false); 
    p=rho0/3.0;

    return 2*(nu*rho0)*dUx(ix,iy)-p;
}

double LatticeBoltzman::sigmayy(int ix,int iy){
    double rho0, p; 
    
    rho0=rho(ix,iy,false); 
    p=rho0/3.0;

    return 2*(nu*rho0)*dUy(ix,iy)-p;
}

double LatticeBoltzman::sigmaxy(int ix,int iy){
    double rho0, p; 
    
    rho0=rho(ix,iy,false); 
    p=rho0/3.0;

    return (nu*rho0)*(dUxy(ix,iy)+dUyx(ix,iy));
}

std::vector<double> LatticeBoltzman::Calcule_dF(double Px, double Py, double dAx, double dAy){

    int i, ix, iy, ixn, iyn;
    ix= int(Px);
    iy= int(Py);
    double dFx, dFy, u=Px-ix, v=Py-iy;
    double sxx_interp, syy_interp, sxy_interp;
    std::vector<double> dF(2, 0.0);
    double sxx[4], syy[4], sxy[4];

    for (i=0; i<4; i++){
        ixn = ix+(i % 2);
        iyn = iy+(i / 2);
        sxx[i] = sigmaxx(ixn, iyn);
        syy[i] = sigmayy(ixn, iyn);
        sxy[i] = sigmaxy(ixn, iyn);
    }

    sxx_interp = sxx[0]*(1-u)*(1-v) + sxx[1]*u*(1-v) + sxx[2]*(1-u)*v + sxx[3]*u*v;
    syy_interp = syy[0]*(1-u)*(1-v) + syy[1]*u*(1-v) + syy[2]*(1-u)*v + syy[3]*u*v;
    sxy_interp = sxy[0]*(1-u)*(1-v) + sxy[1]*u*(1-v) + sxy[2]*(1-u)*v + sxy[3]*u*v;

    dFx = sxx_interp*dAx + sxy_interp*dAy; dFy = sxy_interp*dAx + syy_interp*dAy; 

    return dF={dFx, dFy};
}

std::vector<double> LatticeBoltzman::CalculeFuerza(double Ufan, double Omega) {
    std::vector<double> F(2, 0.0); // Fx, Fy
    int ix, iy;
    double dA, dAx, dAy, x, y;
    double R = 8, R2=R*R; // radio del cilindro
    double theta = 2 * M_PI / N; // ángulo entre cada elemento
    double Re = Ufan*2*R/nu;
    std::vector<double> dF;
       
    for (int i = 0; i < N; i++) {
    x = R * cos(i * theta+theta/2)+ixc;
    y = R * sin(i * theta+theta/2)+iyc;
    dA = R * theta; 
    dAx = dA*cos(i * theta+theta/2); 
    dAy = dA*sin(i * theta+theta/2); 

    dF = Calcule_dF(x, y, dAx, dAy);
    F[0] += dF[0]; // sumar la componente x de la fuerza
    F[1] += dF[1]; // sumar la componente y de la fuerza
    //cout<<x<<' '<<y<<' '<<dF[0]<<' '<<dF[1]<<endl; //--------------> Archivo FuerzaVec.dat
    //cout<<x<<' '<<y<<' '<<dAx<<' '<<dAy<<endl; //--------------> Archivo VecArea.dat
    }
    cout<<Omega<<' '<<abs(F[1])<<endl;
    return F;
}



void LatticeBoltzman::Print(const char * NameFile, double Ufan){
    ofstream MyFile(NameFile); double rho0, Ux0, Uy0; int ix, iy;
    for(ix=0;ix<Lx;ix+=4){
        for(iy=0;iy<Ly;iy+=4){
            rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix,iy,true)/rho0;
            MyFile<<ix<<" "<<iy<<" "<<Ux0/Ufan*4<<" "<<Uy0/Ufan*4<<endl;
        }
        MyFile<<endl;
    }
    MyFile.close();
}

int main(int argc, char **argv){
    double tau = 0.7;
    double Omega0 = std::atof(argv[1]);
    LatticeBoltzman Aire(tau);
    int t, tmax=1000;
    double rho0=1.0, Ufan0=0.1; 
    double R=8;
    std::vector<double> Fuerza(2, 0.0);

    //INICIE
    Aire.Inicie(rho0,Ufan0,0);
    
    //CORRA
    for(t=0;t<tmax;t++){
        Aire.Colision();
        Aire.ImponerCampos(Ufan0,Omega0);
        Aire.Adveccion();
        

    }
    Aire.CalculeFuerza(Ufan0,Omega0);
    
    //Print
    Aire.Print("Magnus.dat", Ufan0);
    return 0;
}