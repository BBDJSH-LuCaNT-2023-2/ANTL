AUTOMAKE_OPTIONS = subdir-objects


AM_LIBS = @PTHREAD_LIBS@ @LIBS@
AM_CXXFLAGS = @CXXFLAGS@ @PTHREAD_CFLAGS@
AM_CC = @PTHREAD_CC@
AM_CPPFLAGS = @CPPFLAGS@ -I$(top_srcdir)/include


bin_PROGRAMS = test main cubic testorder


test_SOURCES= tests/IndexCalculus/IndCalc_test.cpp
main_SOURCES= tests/HeaderTest.cpp src/Quadratic/QuadraticOrder_ZZ.cpp src/Quadratic/QuadraticOrder_long.cpp

cubic_SOURCES=tests/Cubic/cubicTestMain.cpp src/Cubic/generalFunctions.cpp src/Cubic/GlobalCubicField.cpp src/Cubic/CubicNumberField.cpp src/Cubic/RealCubicNumberField.cpp src/Cubic/ComplexCubicNumberField.cpp src/Cubic/CubicOrder.cpp src/Cubic/CubicOrderReal.cpp src/Cubic/CubicElement.cpp src/Cubic/CubicIdeal.cpp src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp src/Cubic/VoronoiMethods.cpp src/Cubic/VoronoiReal.cpp src/Cubic/VoronoiComplex.cpp src/Cubic/fund_unit/BasicVoronoi.cpp src/Cubic/fund_unit/BSGSVoronoi.cpp

testorder_SOURCES=tests/Cubic/catch-tests/test-order.cpp tests/Cubic/catch-tests/testingMain.cpp src/Cubic/generalFunctions.cpp src/Cubic/CubicOrder.cpp

