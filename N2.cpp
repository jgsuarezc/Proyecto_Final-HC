#include<iostream>
#include<cmath>
#include<vector>
#include"D.hpp"
#include"Serial.h"
#include "mpi.h"

void FParallelo(int N,int np,int pid);

int main(int argc, char *argv[]){
  int N=atoi(argv[1]);//Numero de particulas
  double seed=atoi(argv[2]);//semilla genera la posicion aleatoria


    std::vector<Particulas>planeta;
    planeta.resize(N);
    posicion(planeta,N,seed);//llena aleatoriamente
    int M=2;// M columnas (xi,Yi)
    double posiciones[N*M];
    imprimir(planeta,N);//imprime en pantalla


    for(int ii=0;ii<N;++ii){
      for(int jj=0;jj<2;++jj){
          if(jj==0){
          posiciones[M*ii+jj]=planeta[ii].x;
          }
          else{
            posiciones[M*ii+jj]=planeta[ii].y;
          }
      }
    }

    MPI_Init(&argc, &argv);
    int np, pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    //FParallelo(planeta,N,pid,np);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
      FParallelo(N,np,pid);

    MPI_Finalize();



    return 0;
  }


















  void FParallelo(int N,int np,int pid){
    //calcula la furza en paralelo por medio de un anillo
    int Tam=N/np;//numero de particulas que le corresponderia a cada proceso
    int tag =0;
    double Local[2*Tam];//dos datos que voy a pasar x e y por lo tanto  2*(N/np)
    double buffer[2*Tam];
    int siguiente = (pid+1)%np;// proceso siguiente
    int anterior = (pid-1+np)%np;// proceso anterior
    if(0==pid){
        MPI_Send(Local,2*Tam,MPI_INT, siguiente, tag, MPI_COMM_WORLD);

        MPI_Recv(buffer,2*Tam,MPI_INT, anterior, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      }
      else {
        MPI_Recv(buffer,2*Tam,MPI_INT, anterior, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(Local,2*Tam,MPI_INT,siguiente, tag, MPI_COMM_WORLD);
      }

  }
