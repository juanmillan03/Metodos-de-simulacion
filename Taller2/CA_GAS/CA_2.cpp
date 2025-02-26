//CA de Difusion 1D en C++
#include  <iostream>
#include <vector>
#include  <cmath>
#include <numeric>
using namespace std;

const int Lx=256;
const int L2=Lx*Lx;
const double p0=0.25;


class LatticeGas{
private:
  double f[L2],fnew[L2]; // n[ix][i]
public:
  void Borrese(void);
  int Lista(int x,int y);
  void Inicie(int N, double mu,double sigma);
  std::pair<double, double>  GetSigma2(void);
  int Frontera(int x, int y);
  double rho(int ix);
  void Colisione(double p);
  void Adveccione(void);
  void Show(void);
};

void LatticeGas::Borrese(void){
  for(int ix=0;ix<L2;ix++)
    {f[ix]=0;fnew[ix]=0;}
}
int LatticeGas::Lista(int x,int y){
     return y*Lx+x;
}
void LatticeGas::Inicie(int N,double mu,double sigma){
  int ix,i,x,y;
  for(y = 0; y < Lx; y++) 
    for(x = 0; x < Lx; x++){ 
      f[Lista(x,y)]=fnew[Lista(x,y)]=0.5/(2*M_PI*std::pow(sigma,2))*std::exp(-0.5*(std::pow((x-mu)/sigma,2.0)+std::pow((y-mu)/sigma,2.0)));

  }
}

int LatticeGas::Frontera(int x, int y){
  if (x >= Lx) x -= Lx;
  else if (x < 0) x += Lx;
  if (y >= Lx) y -= Lx;
  else if (y < 0) y += Lx;
  return Lista(x, y);
}

void LatticeGas::Colisione(double p){
  int x, y;
  for(y = 0; y < Lx; y++) 
    for(x = 0; x < Lx; x++){
      fnew[Lista(x,y)]=p0*f[Frontera(x,y+1)]+p*f[Frontera(x+1,y)]+p*f[Frontera(x-1,y)]+(1-2*p-p0)*f[Frontera(x,y-1)];
  } 
}

double LatticeGas::rho(int ix){
  return f[ix];
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<L2;ix++){
    f[ix]=fnew[ix];
  }
}
void LatticeGas::Show(void){
  for(int ix=0;ix<L2;ix++){
    if(ix%Lx==0 && ix!=0 )cout<<std::endl;
    std::cout<<f[ix]<<" ";
    
  }
  cout<<endl;
}
//------------------- FUNCIONES GLOBALES -------

std::pair<double, double>  LatticeGas::GetSigma2(void){
  double N=0;
  for(int celda=0;celda<L2;celda++) N +=rho(celda);
  //Calcular la posición promedio  
  double xprom=0,yprom=0;
  for(int y=0;y<Lx;y++){
    for(int x=0;x<Lx;x++){
      xprom+=x*rho(Lista(x,y));
      yprom+=y*rho(Lista(x,y));
    }
  }
  xprom/=N;
  yprom/=N;
  double Sigma2x=0,Sigma2y=0;
  for(int y=0;y<Lx;y++){
    for(int x=0;x<Lx;x++){
      Sigma2x+=pow(x-xprom,2.0)*rho(Lista(x,y));
      Sigma2y+=pow(y-yprom,2.0)*rho(Lista(x,y));
    }
  }
  Sigma2x/=N;
  Sigma2y/=N;
  return std::make_pair(Sigma2x, Sigma2y);
}

pair<double, double> linearRegression(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    double sumX = std::accumulate(x.begin(), x.end(), 0.0);
    double sumY = std::accumulate(y.begin(), y.end(), 0.0);
    double sumXY = 0.0;
    double sumXX = 0.0;

    for (int i = 0; i < n; i++) {
        sumXY += x[i] * y[i];
        sumXX += x[i] * x[i];
    }

    double m = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
    double b = (sumY - m * sumX) / n;

    return make_pair(m, b);
  }
int main(int argc, char **argv){
  double p=std::atof(argv[1]);
  double tmax = 350;
  LatticeGas Difusion;
  double mu=Lx/2, sigma=16;
  int t;
  int N=2400;
  vector<double> T(tmax,0.0);
  vector<double> S2x(tmax,0.0);
  vector<double> S2y(tmax,0.0);
  
  Difusion.Borrese();
  Difusion.Inicie(N,mu,sigma);
  for(t=0;t<tmax;t++){
    std::clog<<t<<" "<<Difusion.GetSigma2().first<<" "<<Difusion.GetSigma2().second<<std::endl;
    T[t]=t;S2x[t]=Difusion.GetSigma2().first;S2y[t]=Difusion.GetSigma2().second;

    Difusion.Colisione(p);
    Difusion.Adveccione();
  }
  pair<double, double> result1= linearRegression(T, S2x);
  pair<double, double> result2= linearRegression(T, S2y);
  std::cout<<p<<" "<<(p+p0)/(2*(1-(p+p0)))<<" "<<result1.first<<" "<<result2.first<<std::endl;
  //Difusion.Show();
  
  
  return 0;
}

