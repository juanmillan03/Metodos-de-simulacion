#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <omp.h>
#include "Random64.h"
#include <sys/stat.h>   // Para la función mkdir
#include <sys/types.h>  // Para el tipo de datos mode_t


const double deltax=0.4;//metro por celda
//--------------------Dimensiones reales del recinto-----------
const double Lx_real=19.7;
const double Ly_real=26.5;
const double LzPequeño_real=5;
const double b_real=3;
const double a_real=Lx_real/2;
const double Lz_real=LzPequeño_real+b_real;
//--------------------Dimensiones simulacion del recinto------
const int LzPequeño=LzPequeño_real/deltax+2;
const int Lx=Lx_real/deltax+2;
const int Ly=Ly_real/deltax+2;
const int Lz=Lz_real/deltax+3;
const int a=a_real/deltax;
const int b=b_real/deltax;
const double deltaT=0.5*deltax/300.0;//segundo por click 
//--------------------Dimensiones reales del resto de elementos-----------
const double h_silla_real=0.4;
const double L_silla_real=0.4;
const double x_sillaI_real=3.0;
const double y_sillaI_real=4.0;
const double deltax_silla_real=0.7;
const double deltay_silla_real=2*deltax_silla_real;
const double h_mesa_real=0.65;
const double L_mesa_real=0.6;
const double x_mesaI_real=3.0+0.3;//Se le suma el radio de la mesa a la posicion inicial de la silla
const double y_mesaI_real=4.0+0.3;//Se le suma el radio de la mesa a la posicion inicial de la silla
const double deltax_mesa_real=2*deltax_silla_real;
const double deltay_mesa_real=2*deltax_silla_real;
const double h_sentado_real=1.3;
const double h_parado_real=1.7;
//--------------------Dimensiones simulacion del resto de elementos-----------
const int h_silla=h_silla_real/deltax+2;
const int L_silla=L_silla_real/deltax;
const int x_sillaI=x_sillaI_real/deltax+2;
const int y_sillaI=y_sillaI_real/deltax+2;
const int deltax_silla=deltax_silla_real/deltax;
const int deltay_silla=deltay_silla_real/deltax;
const int h_mesa=h_mesa_real/deltax+2;
const int L_mesa=L_mesa_real/deltax;
const int x_mesaI=x_mesaI_real/deltax+2;
const int y_mesaI=y_mesaI_real/deltax+2;
const int deltax_mesa=deltax_mesa_real/deltax;
const int deltay_mesa=deltay_mesa_real/deltax;
const int h_sentado=h_sentado_real/deltax+2;
const int h_parado=h_parado_real/deltax+2;



const int Q=7;
const double W0=1.0/4;

const double C=0.5;
const double C2=C*C;
const double AUX0=1-4*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau= 1-Utau;

const double D_paredes = 0.987;
const double D_acero = 0.96;




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
    bool Pared(int ix, int iy, int iz);
    bool Techo(int ix, int iy, int iz);
    bool Silla(int ix, int iy, int iz);
    bool Mesa(int ix, int iy, int iz);
    void Rebote(int ix, int iy, int iz, double D);
    void Frontera(int ix, int iy, int iz);
    double feq(double rho0, double Jx0, double Jy0, double Jz0, int i);
    void Inicie(double rho0, double Jx0, double Jy0, double Jz0);
    void Colision(void);
    void ImponerCampos(int t);
    void Adveccion(void);
    void Print(const char * NameFile,int z, bool UseNew);
    friend class Fuentes; // Declarar a la clase Fuentes como amiga
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

bool LatticeBoltzman::Pared(int ix, int iy, int iz){
    if (ix==Lx-2 || ix==1 || iy==Ly-2 || iy==1 || iz==Lz-2 || iz==1){
        return true;
    }
    else{
        return false;
    }
}

bool LatticeBoltzman::Techo(int ix, int iy, int iz){
    int ixc=Lx/2, izc=LzPequeño;

    if((ix-ixc)*(ix-ixc)*b*b+(iz-izc)*(iz-izc)*a*a>=a*a*b*b && iz>=LzPequeño){
        return true;
    }
    else{
        return false;
    }
}

bool LatticeBoltzman::Silla(int ix, int iy, int iz){
//Para que se vean a baja resolucion usar r2=Lx/4, h_silla=Lz/2, deltax_silla=10, deltay_silla=20

    double r2 = Lx/4/*(L_silla / 2) * (L_silla / 2)*/;

    // Verificamos si está en el rango correcto para iz
    if (iz > 1 && iz <= Ly/2/*h_silla*/) {
        // Recorremos los desplazamientos en x e y para evaluar las posiciones
        for (int dx = 0; dx <= 4; ++dx) {
            for (int dy = 0; dy <= 2; ++dy) {
                double x_shift = x_sillaI + dx * 10/*deltax_silla*/;
                double y_shift = y_sillaI + dy * 20/*deltay_silla*/;

                // Comprobamos si la posición (ix, iy) está dentro del radio r2
                if ((ix - x_shift) * (ix - x_shift) + (iy - y_shift) * (iy - y_shift) <= r2) {
                    return true; 
                }
            }
        }
    }

    return false; 
}

bool LatticeBoltzman::Mesa(int ix, int iy, int iz){
//Para que se vean a baja resolucion usar r2=Lx/3, h_silla=Lz/2, deltax_silla=deltay_silla=20

    double r2 = Lx/3/*(L_mesa / 2) * (L_mesa / 2)*/;

    // Verificamos si está en el rango correcto para iz
    if (iz > 1 && iz <= Ly/2/*h_mesa*/) {
        // Recorremos los desplazamientos en x e y para evaluar las posiciones
        for (int dx = 0; dx <= 2; ++dx) {
            for (int dy = 0; dy <= 2; ++dy) {
                double x_shift = 15/*x_mesaI*/ + dx * 20/*deltax_mesa*/;
                double y_shift = y_mesaI + dy * 20/*deltay_mesa*/;

                // Comprobamos si la posición (ix, iy) está dentro del radio r2
                if ((ix - x_shift) * (ix - x_shift) + (iy - y_shift) * (iy - y_shift) <= r2) {
                    return true; 
                }
            }
        }
    }

    return false; 
}



void LatticeBoltzman::Rebote(int ix, int iy, int iz, double D){
    int n0, n1, n2, n3, n4, n5, n6;

    n0 = n(ix, iy, iz, 0);n1 = n(ix, iy, iz, 1);n3 = n(ix, iy, iz, 3);
    n2 = n(ix, iy, iz, 2);n4 = n(ix, iy, iz, 4);
    n5 = n(ix, iy, iz, 5);n6 = n(ix, iy, iz, 6);
    fnew[n0] = D*f[n0];fnew[n1] = D*f[n2];fnew[n2] = D*f[n1];
    fnew[n3] = D*f[n4];fnew[n4] = D*f[n3];
    fnew[n5] = D*f[n6];fnew[n6] = D*f[n5];
}

void LatticeBoltzman::Frontera(int ix, int iy, int iz){
    int n0, n1, n2, n3, n4, n5, n6;

    n0 = n(ix, iy, iz, 0);n1 = n(ix, iy, iz, 1);n3 = n(ix, iy, iz, 3);
    n2 = n(ix, iy, iz, 2);n4 = n(ix, iy, iz, 4);n5 = n(ix, iy, iz, 5);
    n6 = n(ix, iy, iz, 6);
    fnew[n0] = 0;fnew[n1] = 0;fnew[n2] = 0;
    fnew[n3] = 0;fnew[n4] = 0;
    fnew[n5] = 0;fnew[n6] = 0;
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
                if (Pared(ix,iy,iz) || Techo(ix,iy,iz)){
                    Rebote(ix,iy,iz,D_paredes);
                }
                else if(Silla(ix,iy,iz) || Mesa(ix,iy,iz)){
                    Rebote(ix,iy,iz,D_acero);
                }
                else if (ix==Lx-1 || ix==0 || iy==Ly-1 || iy==0 || iz==Lz-1 || iz==0){
                    Frontera(ix,iy,iz);
                }
                else{
                    for(i=0;i<Q;i++){ //En cada direccion
                        n0=n(ix,iy,iz,i);
                        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);
                    }
                }
            
            }
}
void LatticeBoltzman::ImponerCampos(int t){
    int i, ix, iy, iz, n0;
    double lambda, omega, rho0, Jx0, Jy0, Jz0; lambda=7; omega=2*M_PI/lambda*C;
    //Una fuente oscilante en el medio
    ix=Lx/2; iy=Ly/2; iz=4;
    rho0=10*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
    for(i=0;i<Q;i++){
        n0=n(ix,iy,iz,i);
        fnew[n0]=feq(rho0,Jx0,Jy0,Jz0,i);
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
void LatticeBoltzman::Print(const char * NameFile,int z, bool UseNew){
    std::ofstream MyFile(NameFile); double rho0; int ix, iz, iy;
    //true----->plano XZ
    //false----->plano XY
    if(UseNew){
        iy = z;
        for(ix=0;ix<Lx;ix++){
            for(iz=0;iz<Lz;iz++){
                rho0=rho(ix,iy,iz,true);
                MyFile<<(float)ix*deltax<<" "<<(float)iz*deltax<<" "<<rho0<<std::endl;
            }
            MyFile<<std::endl;
        }
        MyFile.close();
    }
    else{
        iz = z;
        for(ix=0;ix<Lx;ix++){
            for(iy=0;iy<Ly;iy++){
                rho0=rho(ix,iy,iz,true);
                MyFile<<(float)ix*deltax<<" "<<(float)iy*deltax<<" "<<rho0<<std::endl;
            }
            MyFile<<std::endl;
        }
        MyFile.close();
    }
}
//-------------------------------------------------Clase Fuentes-------------------------------------------------
class Fuentes {
private:
    std::string archivotxt;
    int ix, iy, iz;
    std::vector<double> sonido;
    LatticeBoltzman &LB;
public:
    Fuentes(std::string nombreArchivo, LatticeBoltzman& LBn, int Ix, int Iy, int Iz, int tmax)
        : archivotxt(nombreArchivo), ix(Ix), iy(Iy), iz(Iz), LB(LBn), sonido(int(tmax/deltaT),0) {
        std::ifstream archivo(nombreArchivo);
        if (!archivo.is_open()) {
            std::cerr << "Error al abrir el archivo: " << nombreArchivo << std::endl;
            return;
        }
        double valor;
        int i = 0;
        while (archivo >> valor && i < tmax) {
            sonido[i] = 20e-6*std::pow(10,valor/20.0);
            i++;
        }
        archivo.close();
        // std::cout << "Archivo " << nombreArchivo << " leído con éxito. Valores almacenados en el vector 'sonido'." << std::endl;
        // std::cout << "Contenido del vector 'sonido':" << std::endl;
        // for (int i = 0; i < tmax; i++) {
        //     std::cout << "sonido[" << i << "] = " << sonido[i] << std::endl;
        // }

    }
    void ImponerFuente(int t);
    friend class LatticeBoltzman;
};



void Fuentes::ImponerFuente(int t) {
    double rho0 =sonido[t];//amplitud * sin(omega * t);
    double Jx0 = LB.Jx(ix, iy, iz, false);
    double Jy0 = LB.Jy(ix, iy, iz, false);
    double Jz0 = LB.Jz(ix, iy, iz, false);
    for (int i = 0; i < Q; i++) {
        int n0 = LB.n(ix, iy, iz, i);
        LB.fnew[n0] = LB.feq(rho0, Jx0, Jy0, Jz0, i);
    }
}

//-------------------------------------------------Función principal-------------------------------------------------

int main(void){


    // Establecer el número de hilos de forma explícita
    int num_threads = 2;
    omp_set_num_threads(num_threads);
    LatticeBoltzman Ondas;
    int t;
    double tmax=20.0;//segundos 
    double rho0=0.0, Jx0=0, Jy0=0, Jz0=0;
    bool plano=false;

    // Imprimir la cantidad de celdas en cada eje y los valores de deltax y deltaT
    std::cout << "Cantidad de celdas en el eje X (Lx): " << Lx << std::endl;
    std::cout << "Cantidad de celdas en el eje Y (Ly): " << Ly << std::endl;
    std::cout << "Cantidad de celdas en el eje Z (Lz): " << Lz << std::endl;
    std::cout << "Valor de deltax: " << deltax << " metros por celda" << std::endl;
    std::cout << "Valor de deltaT: " << deltaT << " segundos por click" << std::endl;
    std::cout << "Valor de tiempo total: " << tmax << " segundos "<<tmax/deltaT<<" click"<< std::endl;


    //INICIE

    Ondas.Inicie(rho0,Jx0,Jy0,Jz0);
    
    //Corre
    // ...

    auto start = std::chrono::high_resolution_clock::now(); // Start timer
    // Fuentes
    Crandom ran64(23);
    const int Numero_fuentes=0;
    Fuentes* fuentes[Numero_fuentes];
    int random_number_x;
    int random_number_y;
    int txt_number;
    
    for(int r=0;r<Numero_fuentes;r++){// Get a random integer between 0 and 120
        random_number_x= ran64.intRange(Lx/4,Lx*3/4);
        random_number_y= ran64.intRange(Lx/4,Ly*3/4);
        txt_number=ran64.intRange(1,5);
        fuentes[r] = new Fuentes("Fuentes/fuente_" + std::to_string(txt_number) + ".txt", Ondas, random_number_x, random_number_y, 4,tmax); 
    }
    
    for(t=0; t<int(tmax/deltaT); t++){
        Ondas.Colision();
        Ondas.ImponerCampos(t);
        for(int r=0; r<Numero_fuentes; r++) {
            fuentes[r]->ImponerFuente(t);
        }
        Ondas.Adveccion();
        if(t % int(1/deltaT) == 0){
            std::cout << "Imprimendo: " << t << " click "<<(double)t*deltaT<<" segundos"<< std::endl;
            // Crear la carpeta D3/z si no existe
            char directory[30];
            sprintf(directory, "D3/%d", int(3));

            // Crear el archivo y guardar los datos
            char filename[50];
            sprintf(filename, "D3/%d/Ondas_%d.txt", int(3), int(1000*t*deltaT));
            Ondas.Print(filename, int(3/deltax),plano);
        }
        std::clog << t*deltaT << "     " << Ondas.rho(Lz/2, Ly/2, Lz/2, true) << std::endl; //Lugar de medicion
    }


    auto end = std::chrono::high_resolution_clock::now(); // End timer
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start); // Calculate duration in milliseconds
    std::cout << "Total time for all iterations: " << duration.count() << " seconds" << std::endl;
    //Print
    
    return 0;
}