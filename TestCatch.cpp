#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "D.hpp"
#include "Serial.h"
#include <vector>



TEST_CASE( "Verifica la fuerza Gravitacional", "Gravitacional" ) {
  REQUIRE( fuerza(1,1,1) == 1);
  REQUIRE( std::abs(fuerza(1,1,4)-0.0625)<10E-10);
  REQUIRE( std::abs(fuerza(1,1,4)-0.0625)<10E-10);
  REQUIRE( std::abs(fuerza(120,124,49)-6.19)<0.3);
}



TEST_CASE( "Fuerza entre dos particulas", "Finterna" ) {
  //test verifica si  la fuerza entre dos o tres particulas es correcta
  std::vector<Particulas> prueba;
  prueba.resize(2);
  Particulas M1;
  M1.x=111341;
  M1.y=214343;
  Particulas M2;
  M1.x=141234;
  M1.y=133;
  prueba[0]=M1;
  prueba[1]=M2;
  FuerzaT(prueba,2);
  

  REQUIRE(std::abs(prueba[0].Fx==prueba[1].Fx)<0.001);
  REQUIRE(std::abs(prueba[0].Fx==prueba[1].Fx)<0.001);

}
