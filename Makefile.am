AUTOMAKE_OPTIONS = subdir-objects

AM_LIBS = @PTHREAD_LIBS@ @LIBS@
AM_CXXFLAGS = @CXXFLAGS@ @PTHREAD_CFLAGS@ -g
AM_CC = @PTHREAD_CC@
AM_CPPFLAGS = @CPPFLAGS@ -I$(top_srcdir)/include

if EXTRADEPS
bin_PROGRAMS = regulator
else
bin_PROGRAMS = regulator
endif

#main_SOURCES = tests/HeaderTest.cpp                  \
               src/Quadratic/QuadraticOrder_long.cpp \
               src/Quadratic/QuadraticOrder_ZZ.cpp

regulator_SOURCES = tests/UnitTests.cpp                                 \
                    tests/Quadratic/Regulator/RegulatorLenstra_Test.cpp \
                    src/common.cpp                                      \
                    src/L_function/Character.cpp                        \
                    src/L_function/L_function_util.cpp                  \
                    src/L_function/L_function_ZZ.cpp                    \
                    src/Quadratic/QuadraticIdealBase_ZZ.cpp             \
                    src/Quadratic/QuadraticInfElement_ZZ.cpp            \
                    src/Quadratic/QuadraticNumber_ZZ.cpp                \
                    src/Quadratic/QuadraticOrder_ZZ.cpp                 \
                    src/Quadratic/Cube/CubePlain_ZZ.cpp                 \
                    src/Quadratic/Cube/CubeNucube_ZZ.cpp                \
                    src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp        \
                    src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp         \
                    src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp         \
                    src/Quadratic/Square/SquareNudupl_ZZ.cpp            \
                    src/Quadratic/Square/SquarePlain_ZZ.cpp             \
                    src/thresholds.cpp                                  \
                    src/XGCD/hxgcd.cpp                                  \
                    src/XGCD/xgcd.cpp                                   \
                    src/XGCD/xgcd_iter.cpp                              \
                    src/XGCD/xgcd_plain.cpp

#quadinfelm_SOURCES = tests/UnitTests.cpp                                \
                    tests/Quadratic/QuadraticInfElement_ZZ_Tests.cpp    \
                    src/common.cpp                                      \
                    src/L_function/Character.cpp                        \
                    src/L_function/L_function_util.cpp                  \
                    src/L_function/L_function_ZZ.cpp                    \
                    src/Quadratic/QuadraticIdealBase_ZZ.cpp             \
                    src/Quadratic/QuadraticInfElement_ZZ.cpp            \
                    src/Quadratic/QuadraticNumber_ZZ.cpp                \
                    src/Quadratic/QuadraticOrder_ZZ.cpp                 \
                    src/Quadratic/Cube/CubePlain_ZZ.cpp                 \
                    src/Quadratic/Cube/CubeNucube_ZZ.cpp                \
                    src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp        \
                    src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp         \
                    src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp         \
                    src/Quadratic/Square/SquareNudupl_ZZ.cpp            \
                    src/Quadratic/Square/SquarePlain_ZZ.cpp             \
                    src/thresholds.cpp                                  \
                    src/XGCD/hxgcd.cpp                                  \
                    src/XGCD/xgcd.cpp                                   \
                    src/XGCD/xgcd_iter.cpp                              \                                 \
                    src/XGCD/xgcd_plain.cpp

#cubic_SOURCES=tests/Cubic/cubicTestMain.cpp src/Cubic/generalFunctions.cpp src/Cubic/GlobalCubicField.cpp src/Cubic/CubicNumberField.cpp src/Cubic/RealCubicNumberField.cpp src/Cubic/ComplexCubicNumberField.cpp src/Cubic/CubicOrder.cpp src/Cubic/CubicOrderReal.cpp src/Cubic/CubicElement.cpp src/Cubic/CubicIdeal.cpp src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp src/Cubic/VoronoiMethods.cpp src/Cubic/VoronoiReal.cpp src/Cubic/VoronoiComplex.cpp src/Cubic/FundamentalUnits/BasicVoronoi.cpp src/Cubic/FundamentalUnits/BSGSVoronoi.cpp

#testorder_SOURCES=tests/Cubic/catch-tests/test-order.cpp tests/Cubic/catch-tests/testingMain.cpp src/Cubic/generalFunctions.cpp src/Cubic/CubicOrder.cpp

#test_SOURCES = tests/UnitTests.cpp                                   \
               tests/IndexCalculus/IndCalc_Tests.cpp                  \
               tests/common_Tests.cpp                                 \
               tests/LFunction/LFunction_ZZ_Tests.cpp                 \
               tests/Quadratic/QuadraticIdealBase_long_Tests.cpp      \
               tests/Quadratic/QuadraticIdealBase_ZZ_Tests.cpp        \
               tests/Quadratic/QuadraticInfElement_ZZ_Tests.cpp       \
               tests/Quadratic/QuadraticNumber_ZZ_Tests.cpp           \
               tests/Quadratic/QuadraticOrder_ZZ_Tests.cpp            \
               tests/Quadratic/QuadraticOrder_long_Tests.cpp          \
               tests/Quadratic/Cube/CubeNucube_long_Tests.cpp         \
               tests/Quadratic/Cube/CubeNucube_ZZ_Tests.cpp           \
               tests/Quadratic/Cube/CubePlain_long_Tests.cpp          \
               tests/Quadratic/Multiply/MultiplyNucomp_long_Tests.cpp \
               tests/Quadratic/Multiply/MultiplyNucomp_ZZ_Tests.cpp   \
               tests/Quadratic/Multiply/MultiplyPlain_long_Tests.cpp  \
               tests/Quadratic/Multiply/MultiplyPlain_ZZ_Tests.cpp    \
               tests/Quadratic/Reduce/ReducePlainReal_long_Tests.cpp  \
               tests/Quadratic/Reduce/ReducePlainReal_ZZ_Tests.cpp    \
               tests/Quadratic/Square/SquareNudupl_long_Tests.cpp     \
               tests/Quadratic/Square/SquareNudupl_ZZ_Tests.cpp       \
               tests/Quadratic/Square/SquarePlain_ZZ_Tests.cpp        \
               tests/Quadratic/Square/SquarePlain_long_Tests.cpp      \
               src/common.cpp                                         \
               src/L_function/Character.cpp                           \
               src/L_function/L_function_util.cpp                     \
               src/L_function/L_function_ZZ.cpp                       \
               src/Quadratic/QuadraticIdealBase_long.cpp              \
               src/Quadratic/QuadraticIdealBase_ZZ.cpp                \
               src/Quadratic/QuadraticInfElement_ZZ.cpp               \
               src/Quadratic/QuadraticNumber_ZZ.cpp                   \
               src/Quadratic/QuadraticOrder_long.cpp                  \
               src/Quadratic/QuadraticOrder_ZZ.cpp                    \
               src/Quadratic/Cube/CubePlain_long.cpp                  \
               src/Quadratic/Cube/CubePlain_ZZ.cpp                    \
               src/Quadratic/Cube/CubeNucube_long.cpp                 \
               src/Quadratic/Cube/CubeNucube_ZZ.cpp                   \
               src/Quadratic/Multiply/MultiplyNucomp_long.cpp         \
               src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp           \
               src/Quadratic/Multiply/MultiplyPlain_long.cpp          \
               src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp            \
               src/Quadratic/Reduce/ReducePlainReal_long.cpp          \
               src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp            \
               src/Quadratic/Square/SquareNudupl_long.cpp             \
               src/Quadratic/Square/SquareNudupl_ZZ.cpp               \
               src/Quadratic/Square/SquarePlain_long.cpp              \
               src/Quadratic/Square/SquarePlain_ZZ.cpp                \
               src/thresholds.cpp                                     \
               src/XGCD/hxgcd.cpp                                     \
               src/XGCD/xgcd.cpp                                      \
               src/XGCD/xgcd_iter.cpp                                 \
               src/XGCD/xgcd_plain.cpp
