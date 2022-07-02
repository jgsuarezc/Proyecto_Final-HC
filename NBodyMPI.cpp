#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "mpi.h"

// estructura que define las propiedad de la particula
struct Particulas
{
    double x, y;   // posición del la particula
    double Vx, Vy; // velocidad
    double F;      // fuerza sobre ella
    double M = 1;  // masa 1 (provicional)
};

/*
MPI_Datatype mpi_struct_p;//Nuevo Tipo de dato en MPI
MPI_Datatype types[2] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  int len[6] = {1,2,3,4,5,6};
  MPI_Aint offset[6] = {(MPI_Aint)offsetof(Star, m),(MPI_Aint)offsetof(Star, x)};
  MPI_Type_create_struct(2, len, offset, types, &mpi_star);
  MPI_Type_commit(&mpi_star);
*/

double fuerza(double m1, double m2, double d);                       // Fuerza gravitacional
void posicion(std::vector<Particulas> &planeta, int N, double seed); // posición (modifica el vector por referencia)
// paralelo
void FuerzaT(std::vector<Particulas> &planeta, int N, int pid, int np);        // Fuerza Total
void imprimir(const std::vector<Particulas> &planeta, int N, int pid, int np); // imprime en pantalla la posicion y la fuerza
void ring(int pid, int np);                                                    // Implementación de MPIring.

int main(int argc, char *argv[])
{

    // recibe por consola
    int N = atoi(argv[1]);       // Numero de particulas
    double seed = atoi(argv[2]); // semilla genera la posicion aleatoria
    // vector con N particulas SERIALIZADO
    std::vector<Particulas> planeta;
    planeta.resize(N);
    std::cout.precision(5); // precision de 5 cifras
    posicion(planeta, N, seed);
    FuerzaT(planeta, N);
    imprimir(planeta, N);

    // PARALELIZACIÓN.
    MPI_Init(&argc, &argv); /* Mandatory */

    int pid; /* rank of process */
    int np;  /* number of processes */

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // Construcción del ring.
    ring(pid, np);

    MPI_Finalize();

    return 0;
}

/**
 * @brief Usando estructura de Ring MPI para paralelización.
 * Se permite utilizar MPIring para hallar las fuerzas en proceso
 * cíclico, para cada N Cuerpo y suma las fuerzas.
 * Implementación para N = 2.
 *
 * @param pid
 * @param np
 */
void ringFuerzas(int pid, int np)
{
    int tag = 0;
    // 
    std::vector<Double> val = {planeta[0].x, planeta[0].y, planeta[1].x, planeta[1].y, 0.0, 0.0, 0.0};
    std::vector<Double> buf;
    std::vector<Particulas> FPlaneta;
    FPlanetas.resize(2);
    buf.resize(7);
    int next = (pid + 1) % np;
    int prev = (pid - 1 + np) % np;

    double start = MPI_Wtime();
    if (pid == 0)
    {
        MPI_Send(&pid, 1, MPI_INT, next, tag, MPI_COMM_WORLD);
        Particulas p1, p2;
        p1.x = val[0];
        p1.y = val[1];
        p2.x = val[2];
        p2.y = val[3];
        FPlaneta.push_back(p1);
        FPlaneta.push_back(p2);
        double dx = p1.x - p2.x; // distancia xi-xj
        double dy = p1.y - p2.y; // distancia yi-yj
        double d = sqrt(dx * dx + dy * dy);
        double F = fuerza(1, 1, d);
        val 
        MPI_Recv(&buf, 1, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    { // pid != 0
        MPI_Recv(&buf, 1, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        FuerzaT();
        MPI_Send(&val, 1, MPI_INT, next, tag, MPI_COMM_WORLD);
    }
    // Cálculo del tiempo de cómputo total.
    double end = MPI_Wtime();
    totaltime = end - start;

    if (pid == 0)
    {
        std::cout << "Fuerza aplicada a cada cuerpo:" << std::endl;
        std::cout << "Fuerza inicial: "
                  << "\t" <<
    }
    else
    {
        std::cout << "Tiempo total: " << totaltime << std::endl;
        // std::cout << "totaltime/nsteps: " << totaltime / nsteps << std::endl;
    }
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
