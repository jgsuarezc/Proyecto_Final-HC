CXX=g++
CXXFLAGS =-fsanitize=address -fsanitize=leak -fsanitize=undefined


ejecutable.x : N-body.cpp
			$(CXX) $(CXXFLAGS) $< -o $@
clean :
			rm datos.txt ejecutable.x
