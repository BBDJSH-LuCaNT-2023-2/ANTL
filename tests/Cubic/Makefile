main: cubicTestMain.o generalFunctions.o GlobalCubicField.o CubicNumberField.o\
					RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderNF.o \
					CubicElement.o CubicElementNF.o CubicIdeal.o MultiplyStrategyWilliams.o \
					xgcd_plain.o CubicIdealNF.o #QQ.o common.o debug.o

	g++ -g -o testFile -O2 -I/c/include generalFunctions.o cubicTestMain.o GlobalCubicField.o CubicNumberField.o RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderNF.o CubicElement.o CubicElementNF.o CubicIdeal.o MultiplyStrategyWilliams.o xgcd_plain.o CubicIdealNF.o -L/c/lib -lntl -lgmp -lm -lpthread

cubicTestMain.o: cubicTestMain.cpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/ComplexCubicNumberField.hpp ../../include/ANTL/Cubic/RealCubicNumberField.hpp
	g++ -g -c -O2 -I/c/include cubicTestMain.cpp -L/c/lib -lntl -lgmp -lm -lpthread

generalFunctions.o: ../../src/Cubic/generalFunctions.cpp ../../include/ANTL/Cubic/generalFunctions.hpp
	g++ -g -c -O2 -I/c/include ../../src/Cubic/generalFunctions.cpp -L/c/lib -lntl -lgmp -lm -lpthread

GlobalCubicField.o: ../../src/Cubic/GlobalCubicField.cpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
	g++ -g -c -O2 -I/c/include ../../src/Cubic/GlobalCubicField.cpp -L/c/lib -lntl -lgmp -lm -lpthread

CubicNumberField.o: ../../src/Cubic/CubicNumberField.cpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/CubicNumberField.cpp -L/c/lib -lntl -lgmp -lm -lpthread

RealCubicNumberField.o: ../../src/Cubic/RealCubicNumberField.cpp ../../include/ANTL/Cubic/RealCubicNumberField.hpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/RealCubicNumberField.cpp -L/c/lib -lntl -lgmp -lm -lpthread

ComplexCubicNumberField.o: ../../src/Cubic/ComplexCubicNumberField.cpp ../../include/ANTL/Cubic/ComplexCubicNumberField.hpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/ComplexCubicNumberField.cpp -L/c/lib -lntl -lgmp -lm -lpthread

CubicOrder.o: ../../src/Cubic/CubicOrder.cpp ../../include/ANTL/Cubic/CubicOrder.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/CubicOrder.cpp -L/c/lib -lntl -lgmp -lm -lpthread

CubicOrderNF.o: ../../src/Cubic/CubicOrderNF.cpp ../../include/ANTL/Cubic/CubicOrderNF.hpp ../../include/ANTL/Cubic/CubicOrder.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/CubicOrderNF.cpp -L/c/lib -lntl -lgmp -lm -lpthread

CubicElement.o: ../../src/Cubic/CubicElement.cpp ../../include/ANTL/Cubic/CubicElement.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/CubicElement.cpp -L/c/lib -lntl -lgmp -lm -lpthread

CubicElementNF.o: ../../src/Cubic/CubicElementNF.cpp ../../include/ANTL/Cubic/CubicElementNF.hpp ../../include/ANTL/Cubic/CubicElement.hpp #include/ANTL/Arithmetic/QQ.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/CubicElementNF.cpp -L/c/lib -lntl -lgmp -lm -lpthread

CubicIdeal.o: ../../src/Cubic/CubicIdeal.cpp ../../include/ANTL/Cubic/CubicIdeal.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/CubicIdeal.cpp -L/c/lib -lntl -lgmp -lm -lpthread

CubicIdealNF.o: ../../src/Cubic/CubicIdealNF.cpp ../../include/ANTL/Cubic/CubicIdealNF.hpp ../../include/ANTL/Cubic/CubicIdeal.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/CubicIdealNF.cpp -L/c/lib -lntl -lgmp -lm -lpthread

MultiplyStrategyWilliams.o: ../../src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp ../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp ../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp
		g++ -g -c -O2 -I/c/include ../../src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp -L/c/lib -lntl -lgmp -lm -lpthread

xgcd_plain.o: ../../src/XGCD/xgcd_plain.cpp ../../include/ANTL/XGCD/xgcd_plain.hpp
	g++ -g -c -O2 -I/c/include ../../src/XGCD/xgcd_plain.cpp -L/c/lib -lntl -lgmp -lm -lpthread



objects = cubicTestMain.o generalFunctions.o GlobalCubicField.o CubicNumberField.o\
					RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderNF.o \
					CubicElement.o CubicElementNF.o CubicIdeal.o MultiplyStrategyWilliams.o \
					CubicIdealNF.o
clean:
	rm $(objects)
