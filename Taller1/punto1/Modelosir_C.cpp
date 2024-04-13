#include<iostream>
#include <cmath>
#include<fstream>

const double i0=5.326666666666667e-05;
const double s0= 1.0-i0;
const double a= 1.035900845530171;
const double ERR=1.0e-6;
double funcion_gamma(double gamma){return max_modelosir(gamma)-162;}

double s_t(double t, double i, double s,double beta);
double i_t(double t, double i, double s,double beta, double gamma);
void  stepRungekutta4_dos(double & t,double & s,double & i,double  dt,double beta, double gamma);
double max_modelosir(double gamma);
// void  stepRungekutta4(double & t,double & r,double i,double & dt);

int main(){
    std::ofstream outfile;
    outfile.open("datos_C.dat");
// --------------------------------------------
    double gammai{0.001},gammaf{1.5},m,fa,fm;
    fa=funcion_gamma(gammai);
    while ((gammaf-gammai)>ERR)
    {
        m=(gammai+gammaf)/2.0;fm=funcion_gamma(m);
        if (fa*fm<0){gammaf=m;}
        else{gammai=m; fa=fm;}   
    }
    std::cout<<"-----------------------"<<"\n"<<std::endl;
    std::cout<<"El cero se encuentra en gamma="<<(gammai+gammaf)/2.0<<std::endl;
    std::cout<<"\n"<<"\n"<<"-----------------------"<<std::endl;

// --------------------------------------------------------
    double gamma=(gammai+gammaf)/2.0;
    double beta=gamma+std::log(a);
    double t,i,s,r,sum=0;
    double dt=4;
    int T=300;
    for (t=0,s=s0,i=i0,r=0; t<=T;)
    {
        outfile<<t<<" "<<s<<" "<<i<<" "<<r<<std::endl;
        stepRungekutta4_dos(t,i,s,dt,beta,gamma);
        r=1-i-s;
    }
    outfile.close();
    return 0;


}
double s_t(double t, double i, double s,double beta){
    return (-1.0)*beta*s*i;
}
double i_t(double t, double i, double s,double beta, double gamma){
    return beta*s*i-(gamma*i);
}
void stepRungekutta4_dos(double & t,double & i,double & s,
    double dt,double beta, double gamma){

    double di1,di2,di3,di4;                     double ds1,ds2,ds3,ds4;
    di1=dt*i_t(t,i,s,beta,gamma);                           ds1=dt*s_t(t,i,s,beta);
    di2=dt*i_t(t+dt/2,i+di1/2,s+ds1/2,beta,gamma);          ds2=dt*s_t(t+dt/2,i+di1/2,s+ds1/2,beta);
    di3=dt*i_t(t+dt/2,i+di2/2,s+ds2/2,beta,gamma);          ds3=dt*s_t(t+dt/2,i+di2/2,s+ds2/2,beta);
    di4=dt*i_t(t*dt,i+di3,s+ds3,beta,gamma);                ds4=dt*s_t(t*dt,i+di3,s+ds3,beta);
    i+=(di1+2*(di3+di2)+di4)/6;                 s+=(ds1+2*(ds3+ds2)+ds4)/6;
    t+=dt;
}
double max_modelosir(double gamma){
    double t,r,imax=0;
    double s=s0,i=i0;
    double beta=gamma+std::log(a);
    double dt=1,tmax=0;int T=300;
   for (t=0;t<T;){
    if(imax<i){
        imax=i;
        tmax=t;
    }
    stepRungekutta4_dos(t,i,s,dt,beta,gamma);
    r=1-i-s;
    }
    return tmax;
}