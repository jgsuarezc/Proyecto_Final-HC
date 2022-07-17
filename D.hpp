struct Particulas{
  double x,y;// posición del la particula
  double Vx,Vy=0;//velocidad
  double Fx,Fy;
  double F;//modulo de la fuerza
  double M=1;//masa  1(provicional)
};
typedef std::vector<double> V;

typedef std::vector<Particulas> Sistema;
