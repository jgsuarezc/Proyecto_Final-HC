CXX=g++
CXXFLAGS =-fsanitize=address -fsanitize=leak -fsanitize=undefined


ejecutable.x : D.hpp Serial.cpp N-body.cpp
			$(CXX) $(CXXFLAGS) $< -o $@


#%.x: %.o .hpp.o
#	source $$HOME/repos/spack/share/spack/setup-env.sh; \spack load catch2; \
		g++ $$(pkg-config --cflags catch2) $^ -o $@

clean:
	rm -f *.o *~ *.x

#g++ $(pkg-config --cflags catch2) TestCatch.cpp
