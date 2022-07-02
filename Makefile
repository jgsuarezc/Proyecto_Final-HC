SHELL:=/bin/zsh
CXX=g++
CXXFLAGS =-fsanitize=address -fsanitize=leak -fsanitize=undefined


ejecutable.x : NBody.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
clean :
	rm datos.txt ejecutable.x
