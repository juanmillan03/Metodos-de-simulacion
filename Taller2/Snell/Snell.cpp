#include<fstream>
#include  <cmath>

const int Lx=400;
const int Ly=200;

const int Q=5;//
const double W0=1.0/3;// peso del central

const double theta=45.0;

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;



class LatticeBoltzmann{
    private:
      double C=0.5; // Velocidad de la onda SONIDO C<1 criterio courant c<0.707 click
      double w[Q];      //Weights 
      int Vx[Q],Vy[Q];  //Velocity vectors
      double *f, *fnew; //Distribution Functions
    public:
        LatticeBoltzmann(void);
        ~LatticeBoltzmann(void);
        int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
        //---- campos macroscopicos
        double rho(int ix,int iy, bool UseNew);
        double Jx(int ix,int iy, bool UseNew);
        double Jy(int ix,int iy, bool UseNew);
        double Ccelda(int ix,int iy);
        // funciones de equilibrio
        double feq(double rho0,double Jx0,double Jy0,int i);
        //------------- evolucion temporal
        void Start(double rho0,double Jx0,double Jy0);
        void Collision(void);
        void ImposeFields(int t);
        void Advection(void);
        void Print(const char * NameFile);
}; 

LatticeBoltzmann::LatticeBoltzmann(void){
    //Set the weights
    w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
    //Set the velocity vectors
    Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
    Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
    //Create the dynamic arrays
    int ArraySize=Lx*Ly*Q;
    f=new double [ArraySize]; fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[]f;delete []fnew;
}
double LatticeBoltzmann::rho(int ix,int iy, bool UseNew){
    double sum; int i,n0;
    for(sum=0,i=0;i<Q;i++){
        n0=n(ix,iy,i);
        if(UseNew)sum+=fnew[n0];else sum+=f[n0];
    }
    return sum;
}
double LatticeBoltzmann::Jx(int ix,int iy, bool UseNew){
    double sum; int i,n0;
    for(sum=0,i=0;i<Q;i++){
        n0=n(ix,iy,i);
        if(UseNew)sum+=Vx[i]*fnew[n0];else sum+=Vx[i]*f[n0];
    }
    return sum;
}
double LatticeBoltzmann::Jy(int ix,int iy, bool UseNew){
    double sum; int i,n0;
    for(sum=0,i=0;i<Q;i++){
        n0=n(ix,iy,i);
        if(UseNew) sum+=Vy[i]*fnew[n0];else sum+=Vy[i]*f[n0];
    }
    return sum;
}
double  LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i){
  if(i>0)
    return 3*w[i] *(C*C*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return (1.0 - 3*C*C*(1 - W0))*rho0;
}
void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      C=Ccelda(ix,iy);
      for(i=0;i<Q;i++){ //on each direction
        n0=n(ix,iy,i);
        f[n0]=feq(rho0,Jx0,Jy0,i); 
      }
}
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      C=Ccelda(ix,iy);
      for(i=0;i<Q;i++){ //for each velocity vector
        n0=n(ix,iy,i);  
        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i);
      }
    }  
}
double LatticeBoltzmann::Ccelda(int ix, int iy){
  double Cx;
  int x_center = 100;
  int ix0 = tan(theta*M_PI/180)*(iy - x_center) + x_center+theta;
  Cx = -tanh(ix-ix0)/8 + 0.375;
  return Cx;
}
void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0;
  double lambda,omega,rho0,Jx0,Jy0; lambda=10;
  //an oscillating source in the middle

  ix=1;
  for(iy=0;iy<200;iy++){
    C=Ccelda(ix,iy);
    omega=2*M_PI/lambda*C;
    rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
    
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(rho0,Jx0,Jy0,i);
    }
  }
} 
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	    n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	    f[n0next]=fnew[n0]; //periodic boundaries
      }
} 
void LatticeBoltzmann::Print(const char * NameFile){
  std::ofstream MyFile(NameFile); double rho0; int ix,iy;
  for(ix=0;ix<200;ix++){
    for(iy=0;iy<200;iy++){
      rho0=rho(ix,iy,true); 
      //MyFile<<ix<<" "<<iy<<" "<<0.5/Ccelda(ix,iy)<<std::endl;
      MyFile<<ix<<" "<<iy<<" "<<rho0<<std::endl;
    }
    MyFile<<std::endl;
  }
  MyFile.close();
}


int main(void){
    LatticeBoltzmann Ondas;
    int t,tmax=400;
    double rho0=0,Jx0=0,Jy0=0;
    Ondas.Start(rho0,Jx0,Jy0);
    //run
    for(t=0;t<tmax;t++){
        Ondas.Collision();
        Ondas.ImposeFields(t);
        Ondas.Advection();
    }
    Ondas.Print("Ondas.dat");
    return 0;
}