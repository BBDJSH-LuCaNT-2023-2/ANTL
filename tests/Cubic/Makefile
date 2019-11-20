CC=g++
CFLAGS=-g -O2
LFLAGS= -lntl -lgmp -lmpfi -lmpfr -lm -lpthread
SRC_FILES := $(wildcard ../../src/Cubic/*.cpp)
OBJ_FILES := $(notdir $(SRC_FILES:.cpp=.o))

main: cubicTestMain.o generalFunctions.o GlobalCubicField.o CubicNumberField.o\
					RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderReal.o\
					CubicElement.o CubicIdeal.o IdealMultiplicationStrategy.o MultiplyStrategyWilliams.o \
					VoronoiMethods.o VoronoiReal.o VoronoiComplex.o BasicVoronoi.o BSGSVoronoi.o
					#QQ.o common.o debug.o xgcd.o

	$(CC) $(CFLAGS) -o testFile -I/c/include generalFunctions.o cubicTestMain.o GlobalCubicField.o CubicNumberField.o RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderReal.o CubicElement.o CubicIdeal.o IdealMultiplicationStrategy.o MultiplyStrategyWilliams.o VoronoiMethods.o VoronoiReal.o VoronoiComplex.o BasicVoronoi.o BSGSVoronoi.o -L/c/lib $(LFLAGS)

cubicTestMain.o: cubicTestMain.cpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/ComplexCubicNumberField.hpp ../../include/ANTL/Cubic/RealCubicNumberField.hpp
	$(CC) $(CFLAGS) -c -I/c/include cubicTestMain.cpp -L/c/lib $(LFLAGS)

generalFunctions.o: ../../src/Cubic/generalFunctions.cpp ../../include/ANTL/Cubic/generalFunctions.hpp
	$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/generalFunctions.cpp -L/c/lib $(LFLAGS)

GlobalCubicField.o: ../../src/Cubic/GlobalCubicField.cpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
	$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/GlobalCubicField.cpp -L/c/lib $(LFLAGS)

CubicNumberField.o: ../../src/Cubic/CubicNumberField.cpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/CubicNumberField.cpp -L/c/lib $(LFLAGS)

RealCubicNumberField.o: ../../src/Cubic/RealCubicNumberField.cpp ../../include/ANTL/Cubic/RealCubicNumberField.hpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/RealCubicNumberField.cpp -L/c/lib $(LFLAGS)

ComplexCubicNumberField.o: ../../src/Cubic/ComplexCubicNumberField.cpp ../../include/ANTL/Cubic/ComplexCubicNumberField.hpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/ComplexCubicNumberField.cpp -L/c/lib $(LFLAGS)

CubicOrder.o: ../../src/Cubic/CubicOrder.cpp ../../include/ANTL/Cubic/CubicOrder.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/CubicOrder.cpp -L/c/lib $(LFLAGS)

CubicOrderReal.o: ../../src/Cubic/CubicOrderReal.cpp ../../include/ANTL/Cubic/CubicOrderReal.hpp ../../include/ANTL/Cubic/CubicOrder.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/CubicOrderReal.cpp -L/c/lib $(LFLAGS)

CubicElement.o: ../../src/Cubic/CubicElement.cpp ../../include/ANTL/Cubic/CubicElement.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/CubicElement.cpp -L/c/lib $(LFLAGS)

CubicElementNF.o: ../../src/Cubic/CubicElementNF.cpp ../../include/ANTL/Cubic/CubicElementNF.hpp ../../include/ANTL/Cubic/CubicElement.hpp #include/ANTL/Arithmetic/QQ.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/CubicElementNF.cpp -L/c/lib $(LFLAGS)

CubicIdeal.o: ../../src/Cubic/CubicIdeal.cpp ../../include/ANTL/Cubic/CubicIdeal.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/CubicIdeal.cpp -L/c/lib $(LFLAGS)

CubicIdealNF.o: ../../src/Cubic/CubicIdealNF.cpp ../../include/ANTL/Cubic/CubicIdealNF.hpp ../../include/ANTL/Cubic/CubicIdeal.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/CubicIdealNF.cpp -L/c/lib $(LFLAGS)

IdealMultiplicationStrategy.o: ../../src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp ../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp -L/c/lib $(LFLAGS)

MultiplyStrategyWilliams.o: ../../src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp ../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp ../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp -L/c/lib $(LFLAGS)

VoronoiMethods.o: ../../src/Cubic/VoronoiMethods.cpp ../../include/ANTL/Cubic/VoronoiMethods.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/VoronoiMethods.cpp -L/c/lib $(LFLAGS)

VoronoiReal.o: ../../src/Cubic/VoronoiReal.cpp ../../include/ANTL/Cubic/VoronoiReal.hpp ../../include/ANTL/Cubic/VoronoiMethods.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/VoronoiReal.cpp -L/c/lib $(LFLAGS)

VoronoiComplex.o: ../../src/Cubic/VoronoiComplex.cpp ../../include/ANTL/Cubic/VoronoiComplex.hpp ../../include/ANTL/Cubic/VoronoiMethods.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/VoronoiComplex.cpp -L/c/lib $(LFLAGS)

BasicVoronoi.o: ../../src/Cubic/fund_unit/BasicVoronoi.cpp ../../include/ANTL/Cubic/fund_unit/BasicVoronoi.hpp ../../include/ANTL/Cubic/fund_unit/FundUnitStrategy.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/fund_unit/BasicVoronoi.cpp -L/c/lib $(LFLAGS)

BSGSVoronoi.o: ../../src/Cubic/fund_unit/BSGSVoronoi.cpp ../../include/ANTL/Cubic/fund_unit/BSGSVoronoi.hpp ../../include/ANTL/Cubic/fund_unit/FundUnitStrategy.hpp
		$(CC) $(CFLAGS) -c -I/c/include ../../src/Cubic/fund_unit/BSGSVoronoi.cpp -L/c/lib $(LFLAGS)

xgcd.o: ../../src/XGCD/xgcd.cpp ../../include/ANTL/XGCD/xgcd.hpp
	$(CC) $(CFLAGS) -c -I/c/include ../../src/XGCD/xgcd.cpp -L/c/lib $(LFLAGS)



objects = cubicTestMain.o generalFunctions.o GlobalCubicField.o CubicNumberField.o\
					RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o \
					CubicElement.o  CubicIdeal.o  IdealMultiplicationStrategy.o \
					MultiplyStrategyWilliams.o VoronoiMethods.o VoronoiReal.o VoronoiComplex.o\
					CubicOrderReal.o BasicVoronoi.o BSGSVoronoi.o

mul_objects = cubicTestMain.o IdealMultiplicationStrategy.o MultiplyStrategyWilliams.o cubicTestMain.o
voronoi_objects = VoronoiMethods.o VoronoiReal.o VoronoiComplex.o CubicIdeal.o cubicTestMain.o BasicVoronoi.o BSGSVoronoi.o
core_obj = CubicIdeal.o CubicOrder.o CubicElement.o cubicTestMain.o
unit_obj = BasicVoronoi.o CubicOrderReal.o CubicOrder.o
clean:
	rm $(objects)

clean_voronoi:
	rm $(voronoi_objects)
clean_core:
	rm $(core_obj)
clean_unit:
	rm $(unit_obj)
