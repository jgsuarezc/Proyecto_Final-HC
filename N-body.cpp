        #include<iostream>
        #include<cmath>
        #include<vector>
        #include"D.hpp"
        #include"Serial.h"
        #include "mpi.h"

        void XPropios(const std::vector<Particulas> PR,int N,int pid ,int np);
        void XPrestados(const std::vector<Particulas> PR,int N,const std::vector<Particulas> Bu,int pid ,int np);

        void FParallelo(const std::vector<Particulas> planeta,int N,int pid,int np);

        int main(int argc, char *argv[]) {

        //recibe por consola
          int N=atoi(argv[1]);//Numero de particulas
          double seed=atoi(argv[2]);//semilla genera la posicion aleatoria
          //vector con N particulas
          std::vector<Particulas>planeta;
          planeta.resize(N);
          std::cout.precision(5);// precision de 5 cifras
          posicion(planeta,N,seed);
          FuerzaT(planeta,N);
          imprimir(planeta,N);
          /* MPI setup */
          MPI_Init(&argc, &argv);
          MPI_Comm_rank(MPI_COMM_WORLD, &pid);
          MPI_Comm_size(MPI_COMM_WORLD, &np);


          MPI_Finalize();




          return 0;
        }





        // paralelizado
        // Suma lafuerza en paralelo
        //Funcion calcula las fuerzas de las particulas perteneciente al proceso
        void XPropios  (const std::vector<Particulas> PR,int N,int pid ,int np){
          double sum;
        //fuerza debida a las n  particulas propias del proceso
        for (int ii = N*pid/np; ii < N*(pid+1)/np; ii++) {
          double sum=0;
          for(int jj = N*pid/np; jj <N*(pid+1)/np; jj++){
              // para una particula suma las fuerzas debida a las demas  excepto ella misma
              if(jj==ii){continue;}
                double M1T=PR[ii].M;
                double M2T=PR[jj].M;
                double dx=PR[ii].x-PR[jj].x;//distancia xi-xj
                double dy=PR[ii].y-PR[jj].y;//distancia yi-yj
                double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
                sum=sum + fuerza(M1T,M2T,d);//variable va sumando la fuerza sobre la particula Mi debida a las Mj con j distinto que i

              }

          }

        //calcula la fuerza debida a los posiciones recibidas que  guardo  en buffer
        }
        void XPrestados(const std::vector<Particulas> PR,int N,const std::vector<Particulas> Bu,int pid ,int np){
          double sum;
          //fuerza debida a las n  particulas recibidas
          for (int ii = 0; ii <= N*(pid+1)/np; ii++) {
            double sum=0;
            for(int jj = N/np; jj <=(N+1)/np; jj++){
                // para una particula suma las fuerzas debida a las demas  excepto ella misma
                  double M1T=PR[ii].M;
                  double M2T=Bu[jj].M;
                  double dx=PR[ii].x-Bu[jj].x;//distancia xi-xj
                  double dy=PR[ii].y-Bu[jj].y;//distancia yi-yj
                  double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
                  sum=sum + fuerza(M1T,M2T,d);//variable va sumando la fuerza sobre la particula Mi debida a las Mj con j distinto que i

                }
            }

        }

        void FParallelo(const std::vector<Particulas> planeta,int N,int pid,int np)
        {
          //Nuevo Tipo de dato en MPI por medio de la estructura:
          const int nitems=6;
          int len[6] = {1,2,3,4,5,6};
          MPI_Datatype types[6] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
          MPI_Datatype mpi_struct_type;
          MPI_Aint offset[6];
          offset[0] = offsetof(Particulas, x);
          offset[1] = offsetof(Particulas, y);
          offset[2] = offsetof(Particulas, Vx);
          offset[3] = offsetof(Particulas, Vy);
          offset[4] = offsetof(Particulas, F);
          offset[5] = offsetof(Particulas, M);
          MPI_Type_create_struct(nitems, len, offset, types, &mpi_struct_type);
          MPI_Type_commit(&mpi_struct_type);

          int tag =0;
          //N/np particulas para cada proceso
          std::vector<Particulas>P;
          P.resize(N/np);

          for(int ii=0;ii<N/np;ii++){
          P[ii]=planeta[(N*pid)+ii];
        }

          std::vector<Particulas>buffer;
          buffer.resize(N/np);
          //ring
          int next = (pid+1)%np;//siguiente proceso
          int prev = (pid-1+np)%np;//anterior proceso
          if (pid == 0){
              //XPropios(P,N,pid,np);//calcula la fuerza de los N/np primeras particulas
              MPI_Send(&P,N/np,mpi_struct_type, next, tag, MPI_COMM_WORLD);

              MPI_Recv(&buffer,N/np,mpi_struct_type, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            }
            else {
              //XPropios(P,N,pid,np);
              MPI_Recv(&buffer,N/np,mpi_struct_type, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              //XPrestados(P,N,buffer,pid,np);
              MPI_Send(&P,N/np,mpi_struct_type ,next, tag, MPI_COMM_WORLD);
            }

        }
