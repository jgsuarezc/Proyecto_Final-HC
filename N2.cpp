#include <iostream>
#include <cstdlib>
#include <cmath>
#include<vector>
#include"D.hpp"
#include"Serial.h"
#include "mpi.h"


void ring(double P[],int N,int pid, int np);
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

FuerzaTA(posiciones,2*N);
MPI_Init(&argc, &argv); /* Mandatory */

int pid; /* rank of process */
int np;  /* number of processes */

MPI_Comm_size(MPI_COMM_WORLD, &np);
MPI_Comm_rank(MPI_COMM_WORLD, &pid);

ring(pid, np);

MPI_Finalize(); /* Mandatory */






  return 0;
}






//calcula la fuerza de las particulas propias del proceso con las nuevas particulas en buffer
void FBufer(double P[],double Buf[],int Tamaño){
double G=6.67E-11;

  for(int ii=0;ii<Tamaño/2;ii++){
    double sumax=0;
    double sumay=0;
    for(int jj=0;jj<Tamaño/2;jj++){

      double dx=P[2*ii]-Buf[2*jj];
      double dy=P[2*ii+1]-Buf[2*jj+1];
      double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
      double d3=d*d*d;
      sumax=(-G/d3)*dx+sumax;//sumando la fuerza en x debida a las otras particulas
      sumay=(-G/d3)*dy+sumay;//sumando la fuerza en y debida a las otras particulas
        }
    std::cout<<sumax<<"\t"<<sumay<<"\n";

  }
}



void ring(double P[],int N,int pid, int np)
{
  int Tam=2*N/np;//numero de particulas que le corresponderia a cada proceso por dos coordenadas
  int tag =0;
  double val[Tam];
  double buf[Tam];
      for(int ii=0;ii<Tam;ii++){
       double tmp=0;
       val[ii]=P[Tam*pid+ii];

   }

  int next = (pid + 1) % np;
  int prev = (pid - 1 + np) % np;
  if (pid == 0)
  {
    MPI_Send(val, Tam, MPI_INT, next, tag, MPI_COMM_WORLD);

    MPI_Recv(buf, Tam, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else
  { // pid != 0
    MPI_Recv(buf, Tam, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Send(val, Tam, MPI_INT, next, tag, MPI_COMM_WORLD);
  }

}
