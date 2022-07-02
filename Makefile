CXX =mpic++
@CXXFLAGS =-fsanitize=address -fsanitize=leak -fsanitize=undefined
A= Serial.cpp NBodyMPI.cpp

ejecutable.x : $(A)
	$(CXX) $(CXXFLAGS ) $(A) -o ejecutable.x
