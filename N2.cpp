#include <iostream>
#include <cstdlib>
#include <cmath>
#include<vector>
#include"D.hpp"
#include"Serial.h"



int main(int argc, char **argv)
{
  //MPI_Init(&argc, &argv); /* Mandatory */
  //recibe por consola
    int N=atoi(argv[1]);//Numero de particulas
    double seed=atoi(argv[2]);//semilla genera la posicion aleatoria
    std::vector<Particulas>planeta;
    planeta.resize(N);
    posicion(planeta,N,seed);//llena aleatoriamente
    FuerzaT(planeta,N);//calcula la component x e y de la fuerza
    imprimir(planeta,N);//imprime en pantalla
    int M=2;// M columnas (xi,Yi)
    double posiciones[N*M];
    for(int ii=0;ii<N;++ii){
      for(int jj=0;jj<M;++jj){
        if(jj==0){
          posiciones[M*ii+jj]=planeta[ii].x;
        }
        else{
          posiciones[M*ii+jj]=planeta[ii].y;
        }
      }
    }

FLocal(posiciones,2*N);



  return 0;
}





/*
//calcula la fuerza de las particulas propias del proceso con las nuevas particulas en buffer
void FBufer(double P[],double Buf[],int Tamaño){
double G=6.67E-11;

  for(int ii=0;ii<Tamaño/2;ii++){
    double sumax=0;
    double sumay=0;
    for(int jj=0;jj<Tamaño/2;jj++){

      double dx=P[2*ii]-B[2*jj];
      double dy=P[2*ii+1]-B[2*jj+1];
      double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
      double d3=d*d*d;
      sumax=(-G/d3)*dx+sumax;//sumando la fuerza en x debida a las otras particulas
      sumay=(-G/d3)*dy+sumay;//sumando la fuerza en y debida a las otras particulas
        }
    std::cout<<sumax<<"\t"<<sumay<<"\n";

  }
}
*/
