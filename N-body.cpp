#include<iostream>
#include<cmath>
#include<random>
//structura
struct Particulas{
double x,y=0;
double Vx,Vy=0;
double Fx=0;
double M=1;
};

void fuerza(double & F); // Fuerza gravitacional

int int main(int argc, char const *argv[]) {
  int seed=1;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(-1.0, 1.0);
  std::vector<Particulas>masas;
  masas.resize(N);

  for (int ii = 0; ii < N; ++ii) {
    double masas[N].x = dis(gen);
    double masas[N].y = dis(gen);

}

  for (int ii = 0; ii < N; ++ii){
      std::cout << masa[ii].x << "\t"<< masa[ii].y<< std::endl;

  }



//fuerza Total sobre la particula



  for (int ii = 0; ii < N; ++ii) {
    for(int jj = 0; jj < N; ++ii){
      sum=0;
      // para una particula suma las fuerzas debida a las demas  excepto ella misma
        if(jj=!ii){
          double M1T=masa[ii].M;
          double M2T=masa[jj].M;
          double distancia= std::pow(masa[ii].x-std::pow(masa[jj].x,2),2)-std::pow(masa[ii].x-std::pow(masa[jj].x,2),2);
          sum=sum + fuerza(M1T,M2T,distancia)
        }
        //fuerza neta sobre la particula
        masa[ii].F=sum;
        }
      }
    }

  return 0;
}






double fuerza(double m1,double m2,double d){
  double G= 30;// constante Gravitacional
  double F=(G*m1*m2)/std::pow(d,2);// fuerza gravitacional


return double;
}
