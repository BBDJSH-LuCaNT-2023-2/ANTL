AUTOMAKE_OPTIONS = subdir-objects

AM_LIBS = @PTHREAD_LIBS@ @LIBS@
AM_CXXFLAGS = @CXXFLAGS@ @PTHREAD_CFLAGS@ -g
AM_CC = @PTHREAD_CC@
AM_CPPFLAGS = @CPPFLAGS@ -I$(top_srcdir)/include

if EXTRADEPS
bin_PROGRAMS = test main cubic
else
bin_PROGRAMS = test main
endif

test_SOURCES= tests/UnitTests.cpp tests/IndexCalculus/IndCalc_Tests.cpp tests/Quadratic/QuadraticOrder_ZZ_Tests.cpp tests/Quadratic/QuadraticOrder_long_Tests.cpp src/Quadratic/QuadraticOrder_ZZ.cpp src/Quadratic/QuadraticOrder_long.cpp

main_SOURCES= tests/HeaderTest.cpp src/Quadratic/QuadraticOrder_ZZ.cpp src/Quadratic/QuadraticOrder_long.cpp

cubic_SOURCES=tests/Cubic/cubicTestMain.cpp src/Cubic/GlobalCubicField.cpp src/Cubic/CubicNumberField.cpp src/Cubic/RealCubicNumberField.cpp src/Cubic/ComplexCubicNumberField.cpp src/Cubic/CubicOrder.cpp src/Cubic/CubicOrderReal.cpp src/Cubic/CubicElement.cpp src/Cubic/CubicIdeal.cpp src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp src/Cubic/VoronoiMethods.cpp src/Cubic/VoronoiReal.cpp src/Cubic/VoronoiComplex.cpp src/Cubic/FundamentalUnits/BasicVoronoi.cpp src/Cubic/FundamentalUnits/BSGSVoronoi.cpp

