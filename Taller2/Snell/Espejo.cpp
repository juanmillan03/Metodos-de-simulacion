#include<fstream>
#include  <cmath>

const int Lx=600;
const int Ly=200;

const int Q=5;//
const double W0=1.0/3;// peso del central

const double theta=45.0;

double C=0.5;
double C2=C*C;

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;



class LatticeBoltzmann{
    private:
      // Velocidad de la onda SONIDO C<1 criterio courant c<0.707 click
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
    return 3*w[i] *(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return (1.0 - 3*C2*(1 - W0))*rho0;
}
void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
        n0=n(ix,iy,i);
        f[n0]=feq(rho0,Jx0,Jy0,i); 
      }
}
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0,n1,n2,n3,n4; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
        //compute the macroscopic fields on the cell
        rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
        if (ix>=50 && ix<=200 && std::pow(ix-50,2)+std::pow(iy-100,2)>=std::pow(100,2)){
            n0 = n(ix, iy, 0);
            n1 = n(ix, iy, 1);
            n3 = n(ix, iy, 3);
            n2 = n(ix, iy, 2);
            n4 = n(ix, iy, 4);
            //fnew[n0] = f[n0];
            // fnew[n1] = f[n3];
            // fnew[n2] = f[n4];
            // fnew[n3] = f[n1];
            // fnew[n4] = f[n2];
            fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,0);
            fnew[n3] = UmUtau*f[n1]+Utau*feq(rho0,Jx0,Jy0,1);   
            fnew[n4] = UmUtau*f[n2]+Utau*feq(rho0,Jx0,Jy0,2); 
            fnew[n1] = UmUtau*f[n3]+Utau*feq(rho0,Jx0,Jy0,3); 
            fnew[n2] = UmUtau*f[n4]+Utau*feq(rho0,Jx0,Jy0,4); 
        }  
        else{
            for(i=0;i<Q;i++){ //En cada direccion
                n0=n(ix,iy,i);
                fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i);
            }
        }
    }  
}

void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0;
  double lambda,omega,rho0,Jx0,Jy0; lambda=10;
  //an oscillating source in the middle
  ix=0;
  for(iy=0;iy<200;iy++){
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
    int t,tmax=500;
    double rho0=0,Jx0=0,Jy0=0;
    Ondas.Start(rho0,Jx0,Jy0);
    //run
    for(t=0;t<tmax;t++){
        Ondas.Collision();
        Ondas.ImposeFields(t);
        Ondas.Advection();
    }
    Ondas.Print("Espejo.dat");
    return 0;
}