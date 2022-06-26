CXX=g++
CXXFLAGS =-fsanitize=address -fsanitize=leak -fsanitize=undefined

datos.txt : ejecutable.x
			./$^ > $@

ejecutable.x : N-body.cpp
			$(CXX) $(CXXFLAGS) $< -o $@
clean :
			rm -f datos.txt ejecutable.x
