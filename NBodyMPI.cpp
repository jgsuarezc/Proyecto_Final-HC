#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "mpi.h"

/**
 * Advertencias:
 * Código se espera 2 particulas para cada proceso, de forma momentanea hasta que se
 * pueda implementar generalidad. Por tanto, según N particulas, se deben poner N/2 = np
 * número de procesos
 */

// estructura que define las propiedad de la particula
struct Particulas
{
    double x, y;   // posición del la particula
    double Vx, Vy; // velocidad
    double F;
    double Fx; // fuerza sobre ella
    double Fy;
    double M = 1; // masa 1 (provicional)
};

/*
MPI_Datatype mpi_struct_p;//Nuevo Tipo de dato en MPI
MPI_Datatype types[2] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  int len[6] = {1,2,3,4,5,6};
  MPI_Aint offset[6] = {(MPI_Aint)offsetof(Star, m),(MPI_Aint)offsetof(Star, x)};
  MPI_Type_create_struct(2, len, offset, types, &mpi_star);
  MPI_Type_commit(&mpi_star);
*/

std::vector<Particulas> planeta;
// double[3] FuerzaN (double r1x, double r1y, double r2x, double r2y)
double fuerza(double m1, double m2, double d);                       // Fuerza gravitacional
void posicion(std::vector<Particulas> &planeta, int N, double seed); // posición (modifica el vector por referencia)
// paralelo
void FuerzaT(std::vector<Particulas> &planeta, int N);        // Fuerza Total
void imprimir(const std::vector<Particulas> &planeta, int N); // imprime en pantalla la posicion y la fuerza
void ringFuerzas(int pid, int np, int N);                     // Implementación de MPIring.

int main(int argc, char *argv[])
{

    // recibe por consola
    int N = atoi(argv[1]);       // Numero de particulas
    double seed = atoi(argv[2]); // semilla genera la posicion aleatoria
    // vector con N particulas SERIALIZADO
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
    ringFuerzas(pid, np, N);

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
void ringFuerzas(int pid, int np, int N)
{
    int tag = 0;
    int size = 4 + 2 * N;
    int totaltime = 0;
    // Orden (r: Posición, F: Fuerza): r1_x, r1_y, r2_x, r2_y, F1_x, F1_y, F2_x, F2_y... Fn_
    double buf[size];
    double val[size];
    int next = (pid + 1) % np;
    int prev = (pid - 1 + np) % np;

    double start = MPI_Wtime();
    if (pid == 0)
    {
        // Para el proceso 1: Enviar a val, buf tiene que traerme el dato después de la vuelta.
        val[4] = 0;
        val[5] = 0;
        for (int jj = 0; jj < (N / np); jj++)
        {
            val[0] = planeta[jj].x;
            val[1] = planeta[jj].y;
            for (int ii = 0; ii < N; ii++)
            {
                if (ii == jj)
                {
                    continue;
                }
                val[2] = planeta[ii].x;
                val[3] = planeta[ii].y;
                double dx = val[2] - val[0]; // distancia xi-xj
                double dy = val[3] - val[1]; // distancia yi-yj
                double d = sqrt(dx * dx + dy * dy);
                double F = fuerza(1, 1, d);
                val[2 * jj + 4] += -F * dx / d;
                val[2 * jj + 5] += -F * dy / d;
            }
        }

        val[0] = planeta[2].x;
        val[1] = planeta[2].y;
        val[2] = planeta[3].x;
        val[3] = planeta[3].y;
        MPI_Send(val, size, MPI_DOUBLE, next, tag, MPI_COMM_WORLD);
        MPI_Recv(buf, size, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Recv(buf, size, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int jj = 0; jj < (N / np); jj++)
        {
            val[0] = buf[2 * jj];
            val[1] = buf[2 * jj + 1];
            std::cout << val[0] << "\t" << val[1] << std::endl;
            for (int ii = 0; ii < N; ii++)
            {
                val[2] = planeta[ii].x;
                val[3] = planeta[ii].y;
                double dx = val[2] - val[0]; // distancia xi-xj
                double dy = val[3] - val[1]; // distancia yi-yj
                if (dx == 0 || dy == 0)
                {
                    continue;
                }

                double d = sqrt(dx * dx + dy * dy);
                double F = fuerza(1, 1, d);
                val[2 * (jj + pid) + 4] += -F * dx / d;
                val[2 * (jj + pid) + 5] += -F * dy / d;
                std::cout << val[2 * (jj + pid) + 4] << std::endl;
            }
        }
        val[0] = planeta[2 * pid + 2].x;
        val[1] = planeta[2 * pid + 2].y;
        val[2] = planeta[2 * pid + 3].x;
        val[3] = planeta[2 * pid + 3].y;
        // for (int ii = 0; ii < sizeof(val); ii++)
        // {
        //     std::cout << val[ii] << " ";
        // }
        // std::cout << std::endl;

        MPI_Send(val, size, MPI_DOUBLE, next, tag, MPI_COMM_WORLD);
    }
    // Cálculo del tiempo de cómputo total.
    double end = MPI_Wtime();
    totaltime += end - start;

    if (pid == 0)
    {
        std::cout << "Fuerza aplicada a cada cuerpo: (pid == 0, cuando dió la vuelta) " << std::endl;
        std::cout << "F -   Fx1 -   Fy1 -   Fx2 -   Fy2" << std::endl;
        for (int i = 0; i < size; i++)
        {
            std::cout << buf[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Tiempo total: " << totaltime << std::endl;
        std::cout << "Tiempo respecto cada proceso : " << totaltime / np << std::endl;
    }
    else
    {
        std::cout << "Fuerza aplicada a cada cuerpo: " << pid << ", cuando dió la vuelta) " << std::endl;
        std::cout << "F -   Fx1 -   Fy1 -   Fx2 -   Fy2" << std::endl;
        for (int i = 0; i < size; i++)
        {
            std::cout << buf[i] << " ";
        }
        std::cout << std::endl;
    }
}

// Funcion calcula las fuerzas de las particulas perteneciente al proceso propio
// void FPropios(const std::vector<Particulas> PR, int N)
// {
//     double G = 1;
//     double K = 9 * 10E9; // carga electrica
//     // fuerza Total sobre la particula

//     for (int ii = 0; ii < N; ii++)
//     {
//         double sumx = 0;
//         double sumy = 0;
//         for (int jj = 0; jj < N; jj++)
//         {

//             // para una particula suma las fuerzas debida a las demas  excepto ella misma
//             if (jj == ii)
//             {
//                 continue;
//             }
//             double M1T = 1;
//             double M2T = 1;
//             double dx = PR[ii].x - PR[jj].x;    // distancia xi-xj
//             double dy = PR[ii].y - PR[jj].y;    // distancia yi-yj
//             double d = sqrt(dx * dx + dy * dy); // distancia al cuadrado
//             double d3 = d * d * d;
//             sumx = -(G * (M1T * M2T) / (d3)) * dx + sumx; // sumando la fuerza en x debida a las otras particulas
//             sumy = -(G * (M1T * M2T) / (d3)) * dy + sumy; // sumando la fuerza en y debida a las otras particulas
//         }

//         PR[ii].Fx = sumx; // fuerza en x
//         PR[ii].Fy = sumy; // fuerza en y
//     }

// SERIAL
// fuerza gravitacional
double fuerza(double m1, double m2, double d)
{

    double G = 6.67384 * pow(10, -11);         // constante de proporcionalidad de la fuerza
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

// Fuerza para paralelización.
// double[3] FuerzaN (double r1x, double r1y, double r2x, double r2y){

//     double fuerzas[3] = {F, Fx, Fy};
//     return
// }

void FuerzaT(std::vector<Particulas> &planeta, int N)
{
    // fuerza Total sobre la particula
    for (int ii = 0; ii < N; ii++)
    {
        double sum = 0;
        double sumx = 0;
        double sumy = 0;
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
            sumx = sumx - fuerza(M1T, M2T, d) * dx / d;
            sumy = sumy - fuerza(M1T, M2T, d) * dy / d;
        }
        planeta[ii].F = sum;
        planeta[ii].Fx = sumx;
        planeta[ii].Fy = sumy;
    }
}
void imprimir(const std::vector<Particulas> &planeta, int N)
{
    for (int ii = 0; ii < N; ++ii)
    {
        // imprime la posicion x e y junto a la fuerza
        std::cout << planeta[ii].x << "\t" << planeta[ii].y << "\t" << planeta[ii].Fx << "\t" << planeta[ii].Fy << std::endl;
    }
}
