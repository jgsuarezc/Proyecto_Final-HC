CXX = g++
CXXFLAGS = -fsanitize=address -fsanitize=leak -fsanitize=undefined N-body.cpp


data.txt : ejecutable.x
  ./$^ > $@

ejecutable.x : N-body.cpp
	$(CXX) $(CXXFLAGS) -o $@
