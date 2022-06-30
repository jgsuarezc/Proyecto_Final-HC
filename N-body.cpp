#include<iostream>
#include<cmath>
#include<random>
#include<vector>
#include "mpi.h"


//estructura que define las propiedad de la particula
struct Particulas{
  double x,y;// posición del la particula
  double Vx,Vy;//velocidad
  double F;//fuerza sobre ella
  double M=1;//masa  1(provicional)
};


double fuerza(double m1,double m2,double d); // Fuerza gravitacional
void posicion (std::vector<Particulas> &planeta,int N,double seed);//posición (modifica el vector por referencia)
//serie
void FuerzaT (std::vector<Particulas> &planeta,int N);//Fuerza Total
void imprimir(const std::vector<Particulas> &planeta,int N);// imprime en pantalla la posicion y la fuerza
//void FuerzaT2(double % M,const std::vector<Particulas> planeta, int N);//fuerza total array
// imprimir2(const double M,int N);

void XManda(std::vector<Particulas> PR,int pid ,int np);
void XRecibe(std::vector<Particulas> PR,std::vector<Particulas> &Bu,int pid ,int np);

void FParallelo(std::vector<Particulas> &planeta,int N,int pid,int np);




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



  MPI_Init(&argc, &argv);
  int pid, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  MPI_Finalize();



  return 0;
}


//serial

//fuerza gravitacional
double fuerza(double m1,double m2,double d){

  double G= 30;// constante de proporcionalidad de la fuerza
  double F=(G*m1*m2)/std::pow(d,2);// fuerza gravitacional

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
  //fuerza Total sobre la particula
  for (int ii = 0; ii < N; ii++) {
    double sum=0;
    for(int jj = 0; jj < N; jj++){

      // para una particula suma las fuerzas debida a las demas  excepto ella misma
        if(jj==ii){continue;}
          double M1T=planeta[ii].M;
          double M2T=planeta[jj].M;
          double dx=planeta[ii].x-planeta[jj].x;//distancia xi-xj
          double dy=planeta[ii].y-planeta[jj].y;//distancia yi-yj
          double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
          sum=sum + fuerza(M1T,M2T,d);//variable va sumando la fuerza sobre la particula Mi debida a las Mj con j distinto que i
    }//fuerza x e y
    planeta[ii].F=sum;
  }

}
void imprimir(const std::vector<Particulas> &planeta,int N){
  for (int ii = 0; ii< N; ++ii){
    // imprime la posicion x e y junto a la fuerza
    std::cout << planeta[ii].x << "\t"<< planeta[ii].y<< "\t"<< planeta[ii].F<<std::endl;
  }
}

/*
void FuerzaT2(double & M,const std::vector<Particulas> planeta,int N){ //fuerza total array

  for (int ii = 0; ii < N; ii++) {
    double sum=0;
    for(int jj = 0; jj < N; jj++){

      // para una particula suma las fuerzas debida a las demas  excepto ella misma
        if(jj==ii){M[ii+N*jj]=0}
          double M1T=planeta[ii].M;
          double M2T=planeta[jj].M;
          double dx=planeta[ii].x-planeta[jj].x;//distancia xi-xj
          double dy=planeta[ii].y-planeta[jj].y;//distancia yi-yj
          double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
          double d3=sqrt(d*d)
          M[ii+N*jj]= fuerza(M1T,M2T,d);// guarda en un arreglo un la fuerza sobre i debido a j
    }
  }
}


void imprimir2(const double M,int N){
  for (int ii = 0; ii< N; ++ii){
    for(int jj=0;jj<N;++jj)
    // imprime la posicion x e y junto a la fuerza
    std::cout << M[ii+N*jj] <<std::endl;
  }
}

*/

// paralelizado
// Suma lafuerza en paralelo
void FParallelo(std::vector<Particulas> &planeta,int N,int pid,int np)
{

  const int nitems=6;
  int len[6] = {1,2,3,4,5,6};
  MPI_Datatype types[6] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  MPI_Datatype mpi_struct_type;//Nuevo Tipo de dato en MPI
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
  std::vector<Particulas>P;
  P.resize(N/np);

  std::vector<Particulas>buffer;
  buffer.resize(N/np);

  std::vector<double> Fuerzas;
  Fuerzass.resize(N/np);

  int next = (pid+1)%np;
  int prev = (pid-1+np)%np;
  if (pid == 0){
      XManda(P,pid,np);
      MPI_Send(&P, N/np,mpi_struct_type, next, tag, MPI_COMM_WORLD);
      MPI_Recv(&buffer, N/np,mpi_struct_type, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
     else {
      XManda(P,pid,np);
      MPI_Recv(&P,N/np,mpi_struct_type, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      XRecibe(P,buffer,pid,np);
      MPI_Send(&buffer, N/np,mpi_struct_type , next, tag, MPI_COMM_WORLD);
    }

}

//Funcion calcula las fuerzas de las particulas perteneciente al proceso
void XManda(std::vector<Particulas> PR,int pid ,int np){
  double sum;
//fuerza debida a las n  particulas propias del proceso
for (int ii = N*pid/np; ii < N*(pid+1)/np; ii++) {
  double sum=0;
  for(int jj = N*pid/np; jj <N*(pid+1)/np; jj++){
      // para una particula suma las fuerzas debida a las demas  excepto ella misma
      if(jj==ii){continue;}
        double M1T=PR[ii].M;
        double M2T=PR[jj].M;
        double dx=PR[ii].x-planeta[jj].x;//distancia xi-xj
        double dy=planeta[ii].y-planeta[jj].y;//distancia yi-yj
        double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
        sum=sum + fuerza(M1T,M2T,d);//variable va sumando la fuerza sobre la particula Mi debida a las Mj con j distinto que i

      }

  }

//calcula la fuerza debida a los datos recibido del proceso anterior guardos  en buffer
}
void XRecibe(std::vector<Particulas> PR,std::vector<Particulas> Bu,int pid ,int np){

  //fuerza debida a las n  particulas propias del proceso
  for (int ii = 0; ii <= N*(pid+1)/np; ii++) {
    double sum=0;
    for(int jj = N/np; jj <=(N+1)/Np; jj++){
        // para una particula suma las fuerzas debida a las demas  excepto ella misma
          double M1T=Bu[ii].M;
          double M2T=PR[jj].M;
          double dx=Bu[ii].x-planeta[jj].x;//distancia xi-xj
          double dy=PR[ii].y-planeta[jj].y;//distancia yi-yj
          double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
          sum=sum + fuerza(M1T,M2T,d);//variable va sumando la fuerza sobre la particula Mi debida a las Mj con j distinto que i

        }
    }

}
