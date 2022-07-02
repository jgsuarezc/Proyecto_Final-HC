        #include<iostream>
        #include<cmath>
        #include<vector>
        #include"D.hpp"
        #include"Serial.h"
        #include "mpi.h"

        //void XPropios(const std::vector<Particulas> PR,int N);
        //void XPrestados(const std::vector<Particulas> PR,const std::vector<Particulas> Bu,int N);

        void FParallelo(int N,int np,int pid);

        int main(int argc, char *argv[]) {

        //recibe por consola
          int N=atoi(argv[1]);//Numero de particulas
          double seed=atoi(argv[2]);//semilla genera la posicion aleatoria
          //vector con N particulas
          std::vector<Particulas>planeta;
          planeta.resize(N);
          posicion(planeta,N,seed);//llena aleatoriamente
          //FuerzaT(planeta,N);//calcula la component x e y de la fuerza
          //imprimir(planeta,N);//imprime en pantalla

          MPI_Init(&argc, &argv);
          int np, pid;
          MPI_Comm_rank(MPI_COMM_WORLD, &pid);
          FParallelo(N,np,pid);
          MPI_Comm_size(MPI_COMM_WORLD, &np);


          MPI_Finalize();




          return 0;
        }


/*


        // paralelizado  Suma lafuerza en paralelo


        //Funcion calcula las fuerzas de las particulas perteneciente al proceso propio
        void XPropios  (const std::vector<Particulas> PR,int N){
          double G=1;
          double K=9*10E9;//carga electrica
          //fuerza Total sobre la particula


          for (int ii = 0; ii < N; ii++) {
            double sumx=0;
            double sumy=0;
            for(int jj = 0; jj < N; jj++){

              // para una particula suma las fuerzas debida a las demas  excepto ella misma
                if(jj==ii){continue;}
                  double M1T=PR[ii].M;
                  double M2T=PR[jj].M;
                  double dx=PR[ii].x-PR[jj].x;//distancia xi-xj
                  double dy=PR[ii].y-PR[jj].y;//distancia yi-yj
                  double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
                  double d3=d*d*d;
                  sumx=-(G*(M1T*M2T)/(d3))*dx+sumx;//sumando la fuerza en x debida a las otras particulas
                  sumy=-(G*(M1T*M2T)/(d3))*dy+sumy;//sumando la fuerza en y debida a las otras particulas

            }

            PR[ii].Fx=sumx;//fuerza en x
            PR[ii].Fy=sumy;//fuerza en y
          }



        }

        //calcula la fuerzas de las particulas nuevas
        void XPrestados(const std::vector<Particulas> PR,const std::vector<Particulas> Bu,int tamaño){
          double G=10;
          double K=9*10E9;//carga electrica
          //fuerza Total sobre la particula


          for (int ii = 0; ii < tamaño; ii++) {
            double sumx=0;
            double sumy=0;
            for(int jj = 0; jj < tamaño; jj++){

              // para una particula suma las fuerzas debida a las demas  excepto ella misma
                if(jj==ii){continue;}
                  double M1T=PR[ii].M;
                  double M2T=Bu[jj].M;
                  double dx=PR[ii].x-Bu[jj].x;//distancia xi-xj
                  double dy=PR[ii].y-Bu[jj].y;//distancia yi-yj
                  double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
                  double d3=d*d*d;
                  sumx=-(G*(M1T*M2T)/(d3))*dx+sumx;//sumando la fuerza en x debida a las otras particulas
                  sumy=-(G*(M1T*M2T)/(d3))*dy+sumy;//sumando la fuerza en y debida a las otras particulas

            }

            PR[ii].Fx=sumx;//fuerza en x
            PR[ii].Fy=sumy;//fuerza en y
        }
*/

  void FParallelo(int N,int np,int pid)
        {

          //Nuevo Tipo de dato en MPI por medio de la estructura:
          const int nitems=2;
          int len[2] = {1,1};
          MPI_Datatype types[2] = {MPI_DOUBLE,MPI_DOUBLE};
          MPI_Datatype mpi_SistemaP;
          MPI_Aint offset[2];
          offset[0] = offsetof(Particulas, x);
          offset[1] = offsetof(Particulas, y);
          MPI_Type_create_struct(nitems, len, offset, types, &mpi_SistemaP);
          MPI_Type_commit(&mpi_SistemaP);
          int tag = 0;
          double totaltime = 0;
          std::vector<Particulas> val;
          val.resize(np);

          std::vector<Particulas> buf;
          buf.resize(np);

          int next = (pid + 1) % np;
          int prev = (pid - 1 + np) % np;
          double start = MPI_Wtime();
          if (pid == 0)
          {
            MPI_Send(val, np, mpi_SistemaP, next, tag, MPI_COMM_WORLD);
            MPI_Recv(buf, np, mpi_SistemaP, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }
          else
          { // pid != 0
            MPI_Recv(buf, np, mpi_SistemaP, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < np; i++)
            {
              val[i].x = buf[i].y + 1;
            }
            MPI_Send(val, np, mpi_SistemaP, next, tag, MPI_COMM_WORLD);
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



FLocal(double A[array];int tamaño;jj)
