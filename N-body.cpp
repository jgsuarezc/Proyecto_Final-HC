#include<iostream>
#include<cmath>
#include<random>
#include<vector>

//estructura que define las propiedad de la particula
struct Particulas{
  double x,y;// posición del la particula
  double Vx,Vy;//velocidad
  double F;//fuerza sobre ella
  double M=1;// masa
};

double fuerza(double m1,double m2,double d); // Fuerza gravitacional

int main(int argc, char *argv[]) {
  //Numero de particulas
  int N=10;
  int seed=3;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(0, 1.0);
  //vector con N particulas
  std::vector<Particulas>planeta;
  planeta.resize(N);
  //bucle for que me genera posiciones aleatorias para cada particula del vector
  for (int ii = 0; ii < N; ++ii) {
    planeta[N].x = dis(gen);
    planeta[N].y = dis(gen);
}

  for (int ii = 0; ii< N; ++ii){
      std::cout << planeta[ii].x << "\t"<< planeta[ii].y<< std::endl;

  }


  //fuerza Total sobre la particula

  for (int ii = 0; ii < N; ii+1) {
    double sum=0;
    for(int jj = 0; jj < N; jj+1){

      // para una particula suma las fuerzas debida a las demas  excepto ella misma
        if(jj=!ii){
          double M1T=planeta[ii].M;
          double M2T=planeta[jj].M;
          double distancia= std::pow(planeta[ii].x-std::pow(planeta[jj].x,2),2)-std::pow(planeta[ii].x-std::pow(planeta[jj].x,2),2); //(xi-xj)^2+(xj,yj)^2
          sum=sum + fuerza(M1T,M2T,3);//variable va sumando la fuerza sobre la particula Mi debida a las Mj con jdistinto que i 
        }
        //fuerza neta sobre la particula
        }
   planeta[ii].F=sum;
      }


  return 0;
}
