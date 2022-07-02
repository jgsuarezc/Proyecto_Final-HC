#include <iostream>
#include <cmath>
#include <random>
#include <vector>

// estructura que define las propiedad de la particula
struct Particulas
{
  double x, y;   // posición del la particula
  double Vx, Vy; // velocidad
  double F;      // fuerza sobre ella
  double M = 1;  // masa 1 (provicional)
};

double fuerza(double m1, double m2, double d);                       // Fuerza gravitacional
void posicion(std::vector<Particulas> &planeta, int N, double seed); // posición (modifica el vector por referencia)
// serie
void FuerzaT(std::vector<Particulas> &planeta, int N);        // Fuerza Total
void imprimir(const std::vector<Particulas> &planeta, int N); // imprime en pantalla la posicion y la fuerza
// paralelo
void FuerzaT(std::vector<Particulas> &planeta, int N, int pid, int np);        // Fuerza Total
void imprimir(const std::vector<Particulas> &planeta, int N, int pid, int np); // imprime en pantalla la posicion y la fuerza

int main(int argc, char *argv[])
{

  // recibe por consola
  int N = atoi(argv[1]);       // Numero de particulas
  double seed = atoi(argv[2]); // semilla genera la posicion aleatoria
  // vector con N particulas
  std::vector<Particulas> planeta;
  planeta.resize(N);
  std::cout.precision(5); // precision de 5 cifras
  posicion(planeta, N, seed);
  FuerzaT(planeta, N);
  imprimir(planeta, N);

  return 0;
}

// serial

// fuerza gravitacional
double fuerza(double m1, double m2, double d)
{

  double G = 30;                             // constante de proporcionalidad de la fuerza
  double F = (G * m1 * m2) / std::pow(d, 2); // fuerza gravitacional

  return F;
}
// genera una posicion x e y aleatoria para las particulas
void posicion(std::vector<Particulas> &planeta, int N, double seed)
{

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(0, 1.0);

  // bucle for que me genera posiciones aleatorias para cada particula del vector
  for (int ii = 0; ii < N; ++ii)
  {
    planeta[ii].x = dis(gen);
    planeta[ii].y = dis(gen);
  }
}

void FuerzaT(std::vector<Particulas> &planeta, int N)
{
  // fuerza Total sobre la particula
  for (int ii = 0; ii < N; ii++)
  {
    double sum = 0;
    for (int jj = 0; jj < N; jj++)
    {

      // para una particula suma las fuerzas debida a las demas  excepto ella misma
      if (jj == ii)
      {
        continue;
      }
      double M1T = planeta[ii].M;
      double M2T = planeta[jj].M;
      double dx = planeta[ii].x - planeta[jj].x; // distancia xi-xj
      double dy = planeta[ii].y - planeta[jj].y; // distancia yi-yj
      double d = sqrt(dx * dx + dy * dy);        // distancia al cuadrado
      sum = sum + fuerza(M1T, M2T, d);           // variable va sumando la fuerza sobre la particula Mi debida a las Mj con j distinto que i
    }
    planeta[ii].F = sum;
  }
}
void imprimir(const std::vector<Particulas> &planeta, int N)
{
  for (int ii = 0; ii < N; ++ii)
  {
    // imprime la posicion x e y junto a la fuerza
    std::cout << planeta[ii].x << "\t" << planeta[ii].y << "\t" << planeta[ii].F << std::endl;
  }
}
