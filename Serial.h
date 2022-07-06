#pragma once
//estructura que define las propiedad de la particula

double fuerza(double m1,double m2,double d); // Fuerza gravitacional
void posicion (std::vector<Particulas> &planeta,int N,double seed);//posición (modifica el vector por referencia)
//serie
void FuerzaT (std::vector<Particulas> &planeta,int N);//Fuerza Total
void imprimir(const std::vector<Particulas> &planeta,int N);// imprime en pantalla la posicion y la fuerza
void FuerzaTA(double P[],int Tamaño);//calcula Fuerza pero usando array
