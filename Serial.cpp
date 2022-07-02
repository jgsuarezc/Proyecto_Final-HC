#include<iostream>
#include <cmath>
#include <vector>
#include<random>
#include "D.hpp"

//serial

//fuerza gravitacional
double fuerza(double m1,double m2,double d){

  double G= 1;// constante de proporcionalidad de la fuerza
  double F=(G*m1*m2)/(d*d);// fuerza gravitacional

return F;
}
//genera una posicion x e y aleatoria para las particulas
void posicion (std::vector<Particulas> &planeta,int N,double seed){

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(0, 1.0);

  //bucle for que me genera posiciones aleatorias para cada particula del vector
  for (int ii = 0; ii < N; ++ii) {
    planeta[ii].x = dis(gen);
    planeta[ii].y = dis(gen);

  }

}

void FuerzaT (std::vector<Particulas> &planeta,int N){
  double G=1;
  double K=9*10E9;//carga electrica
  //fuerza Total sobre la particula


  for (int ii = 0; ii < N; ii++) {
    double sumx=0;
    double sumy=0;
    for(int jj = 0; jj < N; jj++){

      // para una particula suma las fuerzas debida a las demas  excepto ella misma
        if(jj==ii){continue;}
          double M1T=planeta[ii].M;
          double M2T=planeta[jj].M;
          double dx=planeta[ii].x-planeta[jj].x;//distancia xi-xj
          double dy=planeta[ii].y-planeta[jj].y;//distancia yi-yj
          double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
          double d3=d*d*d;
          sumx=(G*(M1T*M2T)/(d3))*dx+sumx;//sumando la fuerza en x debida a las otras particulas
          sumy=(G*(M1T*M2T)/(d3))*dy+sumy;//sumando la fuerza en y debida a las otras particulas
          ;
    }

    planeta[ii].Fx=sumx;//fuerza en x
    planeta[ii].Fy=sumy;//fuerza en y
  }

}
void imprimir(const std::vector<Particulas> &planeta,int N){
  for (int ii = 0; ii< N; ++ii){
    // imprime la posicion x e y junto a la fuerza
    std::cout << planeta[ii].x << "\t"<< planeta[ii].y<<"\t"<< planeta[ii].Fx<<"\t"<< planeta[ii].Fy<<std::endl;
  }
}
