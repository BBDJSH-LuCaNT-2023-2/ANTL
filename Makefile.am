AUTOMAKE_OPTIONS = subdir-objects


AM_LIBS = @PTHREAD_LIBS@ @LIBS@
AM_CXXFLAGS = @CXXFLAGS@ @PTHREAD_CFLAGS@
AM_CC = @PTHREAD_CC@
AM_CPPFLAGS = @CPPFLAGS@ -I$(top_srcdir)/include

bin_PROGRAMS = test main

test_SOURCES= tests/IndexCalculus/IndCalc_test.cpp 
main_SOURCES= tests/HeaderTest.cpp

bin_PROGRAMS2 = main main2 testelem cubic

cubic_SOURCES=tests/Cubic/cubicTestMain.cpp belabas-extract.cpp src/Cubic/generalFunctions.cpp src/Cubic/GlobalCubicField.cpp src/Cubic/CubicNumberField.cpp src/Cubic/RealCubicNumberField.cpp src/Cubic/ComplexCubicNumberField.cpp src/Cubic/CubicOrder.cpp src/Cubic/CubicOrderReal.cpp src/Cubic/CubicElement.cpp src/Cubic/CubicIdeal.cpp src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp src/Cubic/VoronoiMethods.cpp src/Cubic/VoronoiReal.cpp src/Cubic/VoronoiComplex.cpp src/Cubic/fund_unit/BasicVoronoi.cpp src/Cubic/fund_unit/BSGSVoronoi.cpp 

testelem_SOURCES=tests/Cubic/catch-tests/test-element.cpp tests/Cubic/catch-tests/testingMain.cpp 

