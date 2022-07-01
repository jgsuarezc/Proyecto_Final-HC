#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "D.hpp"
#include "Serial.h"


TEST_CASE( "Verifica la fuerza Gravitacional", "Gravitacional" ) {
  REQUIRE( fuerza(1,1,1) == 1);
  REQUIRE( fuerza(1,1,4)==1/2);

}
