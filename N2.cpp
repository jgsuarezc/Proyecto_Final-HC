#include <iostream>
#include <cstdlib>
#include <cmath>
#include "mpi.h"
#include<vector>
#include"D.hpp"
#include"Serial.h"



void ring(double P[],int N,int pid, int np);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv); /* Mandatory */
  //recibe por consola
    int N=atoi(argv[1]);//Numero de particulas
    double seed=atoi(argv[2]);//semilla genera la posicion aleatoria
    std::vector<Particulas>planeta;
    planeta.resize(N);
    posicion(planeta,N,seed);//llena aleatoriamente
    int M=2;// M columnas (xi,Yi)
    double posiciones[N*M];
    //imprimir(planeta,N);//imprime en pantalla


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






  int pid; /* rank of process */
  int np;  /* number of processes */

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  ring(posiciones,N,pid, np);

  MPI_Finalize(); /* Mandatory */

  return 0;
}

void ring(double P[],int N,int pid, int np)
{
  int Tam=2*N/np;//numero de particulas que le corresponderia a cada proceso por dos coordenadas
  int tag =0;
  double totaltime = 0;
  int val[Tam];
  int buf[Tam];
    for(int ii=0;ii<Tam;ii++){
       val[ii]=P[Tam*pid+ii];

   }
for(int ii=0;ii<N;ii++) {


	std::cout<<P[ii]<<"\n";}

 for(int ii=0;ii<Tam;ii++){
	 std::cout<<val[ii]<<"\n";
 }


  int next = (pid + 1) % np;
  int prev = (pid - 1 + np) % np;
  double start = MPI_Wtime();
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
  double end = MPI_Wtime();
  totaltime += end - start;

  if (pid == 0)
  {
    std::cout << "totaltime: " << totaltime << std::endl;
    std::cout << "0 buf: " << buf << std::endl;
    for (int i = 0; i < np; i++)
    {
      std::cout << buf[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "0 val1: " << val << std::endl;
    for (int i = 0; i < np; i++)
    {
      std::cout << val[i] << "\t";
    }
    std::cout << std::endl;
  }
    if (pid == 1)
  {
    std::cout << "totaltime: " << totaltime << std::endl;
    std::cout << "1 buf: " << buf << std::endl;
    for (int i = 0; i < np; i++)
    {
      std::cout << buf[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "1 val1: " << val << std::endl;
    for (int i = 0; i < np; i++)
    {
      std::cout << val[i] << "\t";
    }
    std::cout << std::endl;
  }
      if (pid == 2)
  {
    std::cout << "totaltime: " << totaltime << std::endl;
    std::cout << "2 buf: " << buf << std::endl;
    for (int i = 0; i < np; i++)
    {
      std::cout << buf[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "2 val1: " << val << std::endl;
    for (int i = 0; i < np; i++)
    {
      std::cout << val[i] << "\t";
    }
    std::cout << std::endl;
  }
}
