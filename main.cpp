#include<iostream>
#include<cmath>
#include<vector>
#include"D.hpp"
#include"Serial.h"
//metodo de Euler
void tiempo(Sistema & F,int N);//actualiza la posicion  a; transcurrir el tiempo
double DT=0.01;
int main(int argc, char *argv[]) {

double N =  std::atof(argv[1]);// Numero de Particulas
double NP = std::atof(argv[2]);// Numero de Particulas
double seed=std::atof(argv[3]);// Semilla genera posiciones aleatorias
Sistema S;
S.resize(N);
//condiciones Iniciales
posicion(S,N,seed);//llena aleatoriamente
for (int tstep = 0; tstep < NP; ++tstep){
    FuerzaT(S,N);
    tiempo(S,N);
      for(int ii=0;ii<N;ii++){
        std::cout << tstep*DT << "\t"<< S[ii].x << "\t"<< S[ii].y<< std::endl;
      }
    }




return 0;
}



//metodo de Euler
void tiempo(Sistema & F,int N){
  for(int ii=0;ii<N;ii++){
      F[ii].y = F[ii].y + DT*F[ii].Vy;
      F[ii].x = F[ii].x + DT*F[ii].Vx;

      F[ii].Vy = F[ii].Vy + DT*F[ii].Fy/F[ii].M;
      F[ii].Vx = F[ii].Vx + DT*F[ii].Fx/F[ii].M;
  }
}
