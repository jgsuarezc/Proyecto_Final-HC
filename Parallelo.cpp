#include<iostream>
#include <cmath>

//calcula la fuerza de las particulas propias del proceso
void FuerzaTA(double P[],int Tamaño){
double G=6.67E-11;

  for(int ii=0;ii<Tamaño/2;ii++){
    double sumax=0;
    double sumay=0;
    for(int jj=0;jj<Tamaño/2;jj++){
      if(ii==jj){continue;}
      double dx=P[2*ii]-P[2*jj];
      double dy=P[2*ii+1]-P[2*jj+1];
      double d= sqrt(dx*dx+dy*dy); // distancia al cuadrado
      double d3=d*d*d;
      sumax=(-G/d3)*dx+sumax;//sumando la fuerza en x debida a las otras particulas
      sumay=(-G/d3)*dy+sumay;//sumando la fuerza en y debida a las otras particulas
        }
    //std::cout<<sumax<<"\t"<<sumay<<"\n";

  }
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
    // std::cout<<sumax<<"\t"<<sumay<<"\n";

  }
}
