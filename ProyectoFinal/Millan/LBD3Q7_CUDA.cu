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


const double deltax=0.1;//metro por celda
const double Lx_real=19.7;
const double Ly_real=26.5;
const double Lz_real=8;
const int Lx=Lx_real/deltax+2;
const int Ly=Ly_real/deltax+2;
const int Lz=Lz_real/deltax+2;
const double deltaT=0.50000*deltax/300.0;//segundo por click 



const int Q=7;
const double W0=1.0/4;

const double C=0.5;
const double C2=C*C;
const double AUX0=1-4*C2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau= 1-Utau;

const double D = 0.987;




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
__global__ void ColisionKernel(double *f, double *fnew, int Lx, int Ly, int Lz, double Utau, double UmUtau, double D) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix < Lx && iy < Ly && iz < Lz) {
        int index = (ix * Ly * Lz + iy * Lz + iz) * 7;
        double rho0 = 0.0, Jx0 = 0.0, Jy0 = 0.0, Jz0 = 0.0;
        for (int i = 0; i < 7; i++) {
            rho0 += f[index + i];  // Calcula la densidad
        }
        // Aplicar las condiciones y actualizaciones (ejemplo)
        if (ix == Lx - 2 || ix == 1 || iy == Ly - 2 || iy == 1 || iz == Lz - 2 || iz == 1) {
            fnew[index + 0] = D * f[index + 0];
            fnew[index + 1] = D * f[index + 2];
            // Completar el resto de las asignaciones
        } else if (ix == Lx - 1 || ix == 0 || iy == Ly - 1 || iy == 0 || iz == Lz - 1 || iz == 0) {
            for (int i = 0; i < 7; i++) {
                fnew[index + i] = 0.0;
            }
        } else {
            for (int i = 0; i < 7; i++) {
                fnew[index + i] = UmUtau * f[index + i] + Utau * /* feq(...) */;
            }
        }
    }
}
void LatticeBoltzman::Colision(void) {
    double *d_f, *d_fnew;

    // Tamaño de los arreglos
    int size = Lx * Ly * Lz * Q * sizeof(double);

    // Reservar memoria en la GPU
    cudaMalloc((void**)&d_f, size);
    cudaMalloc((void**)&d_fnew, size);

    // Copiar datos desde la CPU a la GPU
    cudaMemcpy(d_f, f, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_fnew, fnew, size, cudaMemcpyHostToDevice);

    // Configuración de la cuadrícula y bloques
    dim3 threadsPerBlock(8, 8, 8);  // 512 hilos por bloque
    dim3 numBlocks((Lx + 7) / 8, (Ly + 7) / 8, (Lz + 7) / 8);  // Dividir el trabajo en bloques

    // Llamar al kernel
    ColisionKernel<<<numBlocks, threadsPerBlock>>>(d_f, d_fnew, Lx, Ly, Lz, Utau, UmUtau, D);

    // Sincronizar para asegurar que el kernel ha terminado
    cudaDeviceSynchronize();

    // Copiar los resultados de vuelta a la CPU
    cudaMemcpy(fnew, d_fnew, size, cudaMemcpyDeviceToHost);

    // Liberar la memoria de la GPU
    cudaFree(d_f);
    cudaFree(d_fnew);
}
__global__ void AdveccionKernel(double *f, double *fnew, int *Vx, int *Vy, int *Vz, int Lx, int Ly, int Lz) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix < Lx && iy < Ly && iz < Lz) {
        for (int i = 0; i < 7; i++) {
            int ixnext = (ix + Vx[i] + Lx) % Lx;
            int iynext = (iy + Vy[i] + Ly) % Ly;
            int iznext = (iz + Vz[i] + Lz) % Lz;
            int n0 = (ix * Ly * Lz + iy * Lz + iz) * 7 + i;
            int n0next = (ixnext * Ly * Lz + iynext * Lz + iznext) * 7 + i;
            f[n0next] = fnew[n0];
        }
    }
}


void LatticeBoltzman::ImponerCampos(int t){}
void LatticeBoltzman::Adveccion(void) {
    double *d_f, *d_fnew;
    int *d_Vx, *d_Vy, *d_Vz;

    // Tamaño de los arreglos
    int size = Lx * Ly * Lz * Q * sizeof(double);
    int Vsize = Q * sizeof(int);

    // Reservar memoria en la GPU
    cudaMalloc((void**)&d_f, size);
    cudaMalloc((void**)&d_fnew, size);
    cudaMalloc((void**)&d_Vx, Vsize);
    cudaMalloc((void**)&d_Vy, Vsize);
    cudaMalloc((void**)&d_Vz, Vsize);

    // Copiar los datos desde la CPU a la GPU
    cudaMemcpy(d_f, f, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_fnew, fnew, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vx, Vx, Vsize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vy, Vy, Vsize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vz, Vz, Vsize, cudaMemcpyHostToDevice);

    // Configurar la cuadrícula y los bloques
    dim3 threadsPerBlock(8, 8, 8);
    dim3 numBlocks((Lx + 7) / 8, (Ly + 7) / 8, (Lz + 7) / 8);

    // Llamar al kernel
    AdveccionKernel<<<numBlocks, threadsPerBlock>>>(d_f, d_fnew, d_Vx, d_Vy, d_Vz, Lx, Ly, Lz);

    // Sincronizar para asegurarse de que el kernel ha terminado
    cudaDeviceSynchronize();

    // Copiar los resultados de vuelta a la CPU
    cudaMemcpy(f, d_f, size, cudaMemcpyDeviceToHost);

    // Liberar memoria de la GPU
    cudaFree(d_f);
    cudaFree(d_fnew);
    cudaFree(d_Vx);
    cudaFree(d_Vy);
    cudaFree(d_Vz);
}

void LatticeBoltzman::Print(const char * NameFile,int z){
    std::ofstream MyFile(NameFile); double rho0; int ix, iy;
    int iz = z;
    for(ix=0;ix<Lx;ix++){
        for(iy=0;iy<Ly;iy++){
            rho0=rho(ix,iy,iz,true);
            MyFile<<(float)ix*deltax<<" "<<(float)iy*deltax<<" "<<rho0<<std::endl;
        }
        MyFile<<std::endl;
    }
    MyFile.close();
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
    int num_threads = 4;
    omp_set_num_threads(num_threads);
    LatticeBoltzman Ondas;
    int t;
    double tmax=20.0;//segundos 
    double rho0=0.0, Jx0=0, Jy0=0, Jz0=0;

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
    const int Numero_fuentes=5;
    Fuentes* fuentes[Numero_fuentes];
    int random_number_x;
    int random_number_y;
    int txt_number;
    
    for(int r=0;r<Numero_fuentes;r++){// Get a random integer between 0 and 120
        random_number_x= ran64.intRange(Lx/4,Lx*3/4);
        random_number_y= ran64.intRange(Lx/4,Ly*3/4);
        txt_number=ran64.intRange(1,5);
        fuentes[r] = new Fuentes("Fuentes/fuente_" + std::to_string(txt_number) + ".txt", Ondas, Lx/2, Ly/2, Lz/2,tmax); 
    }
    
    for(t=0; t<int(tmax/deltaT); t++){
        Ondas.Colision();
        Ondas.ImponerCampos(t);
        for(int r=0; r<Numero_fuentes; r++) {
            fuentes[r]->ImponerFuente(t);
        }
        Ondas.Adveccion();
        if(t % int(0.4/deltaT) == 0){
            std::cout << "Imprimendo: " << t << " click "<<(double)t*deltaT<<" segundos"<< std::endl;
            #pragma omp parallel
            for(int z=Lz/4; z<Lz; z=z+Lz/4){
                // Crear la carpeta D3/z si no existe
                char directory[30];
                sprintf(directory, "D3/%d", z);

                // Verificar si el directorio existe, si no, crearlo
                struct stat st = {0};
                if (stat(directory, &st) == -1) {
                    mkdir(directory, 0700);  // Crear el directorio con permisos de lectura/escritura
                }

                // Crear el archivo y guardar los datos
                char filename[50];
                sprintf(filename, "D3/%d/Ondas_%d.txt", z, int(1000*t*deltaT));
                Ondas.Print(filename, z);
            }
        }
        std::clog << t*deltaT << "     " << Ondas.rho(Lz/2, Ly/2, Lz/2, true) << std::endl;
    }


    auto end = std::chrono::high_resolution_clock::now(); // End timer
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start); // Calculate duration in milliseconds
    std::cout << "Total time for all iterations: " << duration.count() << " seconds" << std::endl;
    //Print
    
    return 0;
}