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

ring(posiciones,N,pid, np);

MPI_Finalize(); /* Mandatory */






  return 0;
}




void ring(double P[],int N,V &tmp,int pid, int np)
{
  int Tam=2*N/np;//numero de particulas que le corresponderia a cada proceso por dos coordenadas
  int tag =0;
  double val[Tam];
  double buf[Tam];
      for(int ii=0;ii<Tam;ii++){
       double tmp=0;
       val[ii]=P[Tam*pid+ii];

   }
  V tmpF;
  tmpF.resize(Tam);
  int next = (pid + 1) % np;
  int prev = (pid - 1 + np) % np;
  if (pid == 0)
  {
    FuerzaTA(val,Tam);
    MPI_Send(val, Tam, MPI_INT, next, tag, MPI_COMM_WORLD);
    MPI_Recv(buf, Tam, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    FBufer(val,buf,Tam);
  }
  else
  { // pid != 0
    FuerzaTA(val,Tam);
    MPI_Recv(buf, Tam, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    FBufer(val,buf,Tam);
    MPI_Send(val, Tam, MPI_INT, next, tag, MPI_COMM_WORLD);
  }

}
