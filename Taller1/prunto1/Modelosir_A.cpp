#include<iostream>
#include <cmath>
#include<fstream>

double s_t(double t, double i, double s);
double i_t(double t, double i, double s);
const double beta=0.35;
const double Gamma=0.08;
void  stepRungekutta4_dos(double & t,double & s,double & i,double  dt);
// void  stepRungekutta4(double & t,double & r,double i,double & dt);
int main(){
    double t,i,s,r,sum=0;
    double dt=2;
    int T=70;
    
    for (t=0,s=0.999,i=0.001,r=0; t<=T;)
    {
        std::cout<<t<<" "<<s<<" "<<i<<" "<<r<<std::endl;
        stepRungekutta4_dos(t,i,s,dt);
        r=1-i-s;
        // stepRungekutta4(t,r,i,dt);
    }
    return 0;
}
double s_t(double t, double i, double s){
    return (-1.0)*beta*s*i;
}
double i_t(double t, double i, double s){
    return beta*s*i-(Gamma*i);
}
// double r_t(double t, double i){return gamma*i;}

void stepRungekutta4_dos(double & t,double & i,double & s,double dt){
    double di1,di2,di3,di4;                     double ds1,ds2,ds3,ds4;
    di1=dt*i_t(t,i,s);                           ds1=dt*s_t(t,i,s);
    di2=dt*i_t(t+dt/2,i+di1/2,s+ds1/2);          ds2=dt*s_t(t+dt/2,i+di1/2,s+ds1/2);
    di3=dt*i_t(t+dt/2,i+di2/2,s+ds2/2);          ds3=dt*s_t(t+dt/2,i+di2/2,s+ds2/2);
    di4=dt*i_t(t*dt,i+di3,s+ds3);                ds4=dt*s_t(t*dt,i+di3,s+ds3);
    i+=(di1+2*(di3+di2)+di4)/6;                 s+=(ds1+2*(ds3+ds2)+ds4)/6;
    t+=dt;
}
