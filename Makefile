CXX =g++
CXXFLAGS =-fsanitize=address -fsanitize=leak -fsanitize=undefined
A= Serial.cpp N-body.cpp

ejecutable.x : $(A)
	$(CXX) $(CXXFLAGS ) $(A) -o ejecutable.x
