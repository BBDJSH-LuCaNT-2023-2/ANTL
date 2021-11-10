CC=g++
CFLAGS= -g -O2 -Werror
LFLAGS= -lntl -lgmp -lmpfi -lm -lpthread
IFLAGS= -I/c/include -I/home/cs-gb2khv2/ANTL/include
SRC_FILES := $(wildcard ../../src/Cubic/*.cpp)
OBJ_FILES := $(notdir $(SRC_FILES:.cpp=.o))
#For GCC/G++, you can pass compiler flags -std=c++11, -std=c++14, -std=c++17, or -std=c++2a to enable C++11/14/17/2a support respectively.
main: cubicTestMain.o GlobalCubicField.o CubicNumberField.o\
					RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderReal.o\
					CubicElement.o CubicIdeal.o IdealMultiplicationStrategy.o MultiplyStrategyWilliams.o \
					VoronoiMethods.o VoronoiReal.o VoronoiComplex.o BasicVoronoi.o BSGSVoronoi.o
					#QQ.o common.o debug.o xgcd.o

	$(CC) $(CFLAGS) -o bin/testFile $(IFLAGS) cubicTestMain.o GlobalCubicField.o CubicNumberField.o RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderReal.o CubicElement.o CubicIdeal.o IdealMultiplicationStrategy.o MultiplyStrategyWilliams.o VoronoiMethods.o VoronoiReal.o VoronoiComplex.o BasicVoronoi.o BSGSVoronoi.o -L/c/lib $(LFLAGS)

main2: belabas-extract.o
	$(CC) $(CFLAGS) -o bin/belabasrun $(IFLAGS) belabas-extract.o GlobalCubicField.o CubicNumberField.o RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o CubicOrderReal.o CubicElement.o CubicIdeal.o IdealMultiplicationStrategy.o MultiplyStrategyWilliams.o VoronoiMethods.o VoronoiReal.o VoronoiComplex.o BasicVoronoi.o BSGSVoronoi.o -L/c/lib $(LFLAGS)

belabas-extract.o:
	 $(CC) $(CFLAGS) -c $(IFLAGS) belabas-extract.cpp -L/c/lib $(LFLAGS)

cubicTestMain.o: cubicTestMain.cpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/ComplexCubicNumberField.hpp ../../include/ANTL/Cubic/RealCubicNumberField.hpp
	$(CC) $(CFLAGS) -c $(IFLAGS) cubicTestMain.cpp -L/c/lib $(LFLAGS)

GlobalCubicField.o: ../../src/Cubic/GlobalCubicField.cpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
	$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/GlobalCubicField.cpp -L/c/lib $(LFLAGS)

CubicNumberField.o: ../../src/Cubic/CubicNumberField.cpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/CubicNumberField.cpp -L/c/lib $(LFLAGS)

RealCubicNumberField.o: ../../src/Cubic/RealCubicNumberField.cpp ../../include/ANTL/Cubic/RealCubicNumberField.hpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/RealCubicNumberField.cpp -L/c/lib $(LFLAGS)

ComplexCubicNumberField.o: ../../src/Cubic/ComplexCubicNumberField.cpp ../../include/ANTL/Cubic/ComplexCubicNumberField.hpp ../../include/ANTL/Cubic/CubicNumberField.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/ComplexCubicNumberField.cpp -L/c/lib $(LFLAGS)

CubicOrder.o: ../../src/Cubic/CubicOrder.cpp ../../include/ANTL/Cubic/CubicOrder.hpp ../../include/ANTL/Cubic/GlobalCubicField.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/CubicOrder.cpp -L/c/lib $(LFLAGS)

CubicOrderReal.o: ../../src/Cubic/CubicOrderReal.cpp ../../include/ANTL/Cubic/CubicOrderReal.hpp ../../include/ANTL/Cubic/CubicOrder.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/CubicOrderReal.cpp -L/c/lib $(LFLAGS)

CubicElement.o: ../../src/Cubic/CubicElement.cpp ../../include/ANTL/Cubic/CubicElement.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/CubicElement.cpp -L/c/lib $(LFLAGS)

CubicElementNF.o: ../../src/Cubic/CubicElementNF.cpp ../../include/ANTL/Cubic/CubicElementNF.hpp ../../include/ANTL/Cubic/CubicElement.hpp #include/ANTL/Arithmetic/QQ.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/CubicElementNF.cpp -L/c/lib $(LFLAGS)

CubicIdeal.o: ../../src/Cubic/CubicIdeal.cpp ../../include/ANTL/Cubic/CubicIdeal.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/CubicIdeal.cpp -L/c/lib $(LFLAGS)

CubicIdealNF.o: ../../src/Cubic/CubicIdealNF.cpp ../../include/ANTL/Cubic/CubicIdealNF.hpp ../../include/ANTL/Cubic/CubicIdeal.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/CubicIdealNF.cpp -L/c/lib $(LFLAGS)

IdealMultiplicationStrategy.o: ../../src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp ../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp -L/c/lib $(LFLAGS)

MultiplyStrategyWilliams.o: ../../src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp ../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp ../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp -L/c/lib $(LFLAGS)

VoronoiMethods.o: ../../src/Cubic/VoronoiMethods.cpp ../../include/ANTL/Cubic/VoronoiMethods.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/VoronoiMethods.cpp -L/c/lib $(LFLAGS)

VoronoiReal.o: ../../src/Cubic/VoronoiReal.cpp ../../include/ANTL/Cubic/VoronoiReal.hpp ../../include/ANTL/Cubic/VoronoiMethods.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/VoronoiReal.cpp -L/c/lib $(LFLAGS)

VoronoiComplex.o: ../../src/Cubic/VoronoiComplex.cpp ../../include/ANTL/Cubic/VoronoiComplex.hpp ../../include/ANTL/Cubic/VoronoiMethods.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/VoronoiComplex.cpp -L/c/lib $(LFLAGS)

BasicVoronoi.o: ../../src/Cubic/FundamentalUnits/BasicVoronoi.cpp ../../include/ANTL/Cubic/FundamentalUnits/BasicVoronoi.hpp ../../include/ANTL/Cubic/FundamentalUnits/FundUnitStrategy.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/FundamentalUnits/BasicVoronoi.cpp -L/c/lib $(LFLAGS)

BSGSVoronoi.o: ../../src/Cubic/FundamentalUnits/BSGSVoronoi.cpp ../../include/ANTL/Cubic/FundamentalUnits/BSGSVoronoi.hpp ../../include/ANTL/Cubic/FundamentalUnits/FundUnitStrategy.hpp
		$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/Cubic/FundamentalUnits/BSGSVoronoi.cpp -L/c/lib $(LFLAGS)

xgcd.o: ../../src/XGCD/xgcd.cpp ../../include/ANTL/XGCD/xgcd.hpp
	$(CC) $(CFLAGS) -c $(IFLAGS) ../../src/XGCD/xgcd.cpp -L/c/lib $(LFLAGS)



objects = cubicTestMain.o GlobalCubicField.o CubicNumberField.o\
					RealCubicNumberField.o ComplexCubicNumberField.o CubicOrder.o \
					CubicElement.o  CubicIdeal.o  IdealMultiplicationStrategy.o \
					MultiplyStrategyWilliams.o VoronoiMethods.o VoronoiReal.o VoronoiComplex.o\
					CubicOrderReal.o BasicVoronoi.o BSGSVoronoi.o

mul_objects = cubicTestMain.o IdealMultiplicationStrategy.o MultiplyStrategyWilliams.o cubicTestMain.o
voronoi_objects = VoronoiMethods.o VoronoiReal.o VoronoiComplex.o CubicIdeal.o cubicTestMain.o BasicVoronoi.o BSGSVoronoi.o
core_obj = CubicIdeal.o CubicOrder.o CubicElement.o cubicTestMain.o CubicOrderReal.o
unit_obj = BasicVoronoi.o CubicOrderReal.o CubicOrder.o
belabas_obj = belabas-extract.o VoronoiReal.o CubicOrderReal.o CubicElement.o
clean:
	rm $(objects)

clean_voronoi:
	rm $(voronoi_objects)
clean_core:
	rm $(core_obj)
clean_unit:
	rm $(unit_obj)
clean_belabas:
	rm $(belabas_obj)
