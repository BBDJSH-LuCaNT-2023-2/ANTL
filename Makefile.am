AUTOMAKE_OPTIONS = subdir-objects

AM_LIBS = @PTHREAD_LIBS@ @LIBS@ -L/home/austin.docherty/ntl/lib
AM_CXXFLAGS = @CXXFLAGS@ -I/home/austin.docherty/ntl/include -L/home/austin.docherty/ntl/lib
AM_CC =
AM_CPPFLAGS = @CPPFLAGS@ -I$(top_srcdir)/include

#extra libraries when building with non-standard lib locations
interval_LDADD = -lntl -lpthread

if EXTRADEPS
bin_PROGRAMS = arc regulator classgroupBSGS classgroupBS
else
bin_PROGRAMS = arc regulator classgroupBSGS classgroupBS
endif

arc_SOURCES = tests/UnitTests.cpp \
              tests/ARC_test.cpp

#regulator classgroupBSGS classgroupBS interval tabulation harvest

#main_SOURCES = tests/HeaderTest.cpp                  \
               src/Quadratic/QuadraticOrder_long.cpp \
               src/Quadratic/QuadraticOrder_ZZ.cpp

#tabulation_SOURCES = tests/Quadratic/Applications/Tabulation/InvariantsTabulation_Tests.cpp     \
                     tests/Quadratic/Applications/Tabulation/timer.cpp                          \
                     src/common.cpp                                                             \
                     src/LinearAlgebra/Smith.cpp                                                \
                     src/L_function/Character.cpp                                               \
                     src/L_function/L_function_util.cpp                                         \
                     src/L_function/L_function_ZZ.cpp                                           \
                     src/Quadratic/QuadraticIdealBase_ZZ.cpp                                    \
                     src/Quadratic/QuadraticInfElement_ZZ.cpp                                   \
                     src/Quadratic/QuadraticNumber_ZZ.cpp                                       \
                     src/Quadratic/QuadraticOrder_ZZ.cpp                                        \
                     src/Quadratic/Cube/CubePlain_ZZ.cpp                                        \
                     src/Quadratic/Cube/CubeNucube_ZZ.cpp                                       \
                     src/Quadratic/ClassGroup/ClassGroupBSGSReal.cpp                            \
                     src/Quadratic/Multiply/MultiplyComp_ZZ.cpp                                 \
                     src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp                               \
                     src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp                                \
                     src/Quadratic/Regulator/RegulatorLenstra_ZZ.cpp                            \
                     src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp                                \
                     src/Quadratic/Square/SquareNudupl_ZZ.cpp                                   \
                     src/Quadratic/Square/SquarePlain_ZZ.cpp                                    \
                     src/thresholds.cpp                                                         \
                     src/XGCD/hxgcd.cpp                                                         \
                     src/XGCD/xgcd.cpp                                                          \
                     src/XGCD/xgcd_iter.cpp                                                     \
                     src/XGCD/xgcd_plain.cpp

#interval_SOURCES = tests/Quadratic/Applications/Tabulation/QCL_Tabulation_Interval.cpp        \
                     tests/Quadratic/Applications/Tabulation/timer.cpp                          \
                     src/common.cpp                                                             \
                     src/LinearAlgebra/Smith.cpp                                                \
                     src/L_function/Character.cpp                                               \
                     src/L_function/L_function_util.cpp                                         \
                     src/L_function/L_function_ZZ.cpp                                           \
                     src/Quadratic/QuadraticIdealBase_ZZ.cpp                                    \
                     src/Quadratic/QuadraticInfElement_ZZ.cpp                                   \
                     src/Quadratic/QuadraticNumber_ZZ.cpp                                       \
                     src/Quadratic/QuadraticOrder_ZZ.cpp                                        \
                     src/Quadratic/Cube/CubePlain_ZZ.cpp                                        \
                     src/Quadratic/Cube/CubeNucube_ZZ.cpp                                       \
                     src/Quadratic/ClassGroup/ClassGroupBSGSReal.cpp                            \
                     src/Quadratic/Multiply/MultiplyComp_ZZ.cpp                                 \
                     src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp                               \
                     src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp                                \
                     src/Quadratic/Regulator/RegulatorLenstra_ZZ.cpp                            \
                     src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp                                \
                     src/Quadratic/Square/SquareNudupl_ZZ.cpp                                   \
                     src/Quadratic/Square/SquarePlain_ZZ.cpp                                    \
                     src/thresholds.cpp                                                         \
                     src/XGCD/hxgcd.cpp                                                         \
                     src/XGCD/xgcd.cpp                                                          \
                     src/XGCD/xgcd_iter.cpp                                                     \
                     src/XGCD/xgcd_plain.cpp

#harvest_SOURCES =    tests/Quadratic/Applications/Tabulation/QCL_Tabulation_Harvest.cpp         \
                     tests/Quadratic/Applications/Tabulation/timer.cpp                          \
                     tests/Quadratic/Applications/Tabulation/iq-data.cpp                        \
                     src/common.cpp                                                             \
                     src/LinearAlgebra/Smith.cpp                                                \
                     src/L_function/Character.cpp                                               \
                     src/L_function/L_function_util.cpp                                         \
                     src/L_function/L_function_ZZ.cpp                                           \
                     src/Quadratic/QuadraticIdealBase_ZZ.cpp                                    \
                     src/Quadratic/QuadraticInfElement_ZZ.cpp                                   \
                     src/Quadratic/QuadraticNumber_ZZ.cpp                                       \
                     src/Quadratic/QuadraticOrder_ZZ.cpp                                        \
                     src/Quadratic/Cube/CubePlain_ZZ.cpp                                        \
                     src/Quadratic/Cube/CubeNucube_ZZ.cpp                                       \
                     src/Quadratic/ClassGroup/ClassGroupBSGSReal.cpp                            \
                     src/Quadratic/Multiply/MultiplyComp_ZZ.cpp                                 \
                     src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp                               \
                     src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp                                \
                     src/Quadratic/Regulator/RegulatorLenstra_ZZ.cpp                            \
                     src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp                                \
                     src/Quadratic/Square/SquareNudupl_ZZ.cpp                                   \
                     src/Quadratic/Square/SquarePlain_ZZ.cpp                                    \
                     src/thresholds.cpp                                                         \
                     src/XGCD/hxgcd.cpp                                                         \
                     src/XGCD/xgcd.cpp                                                          \
                     src/XGCD/xgcd_iter.cpp                                                     \
                     src/XGCD/xgcd_plain.cpp


regulator_SOURCES = tests/UnitTests.cpp                                     \
                    tests/Quadratic/Regulator/RegulatorLenstra_Test.cpp     \
                    tests/Quadratic/TestData/TestData.cpp                   \
                    src/common.cpp                                          \
                    src/LinearAlgebra/Smith.cpp                             \
                    src/L_function/Character.cpp                            \
                    src/L_function/L_function_util.cpp                      \
                    src/L_function/L_function_ZZ.cpp                        \
                    src/L_function/L_function_long.cpp                      \
                    src/Quadratic/QuadraticIdealBase_ZZ.cpp                 \
                    src/Quadratic/QuadraticIdealBase_long.cpp               \
                    src/Quadratic/QuadraticInfElement_ZZ.cpp                \
                    src/Quadratic/QuadraticInfElement_long.cpp              \
                    src/Quadratic/QuadraticNumber_ZZ.cpp                    \
                    src/Quadratic/QuadraticNumber_long.cpp                  \
                    src/Quadratic/QuadraticOrder_ZZ.cpp                     \
                    src/Quadratic/QuadraticOrder_long.cpp                   \
                    src/Quadratic/Cube/CubePlain_ZZ.cpp                     \
                    src/Quadratic/Cube/CubePlain_long.cpp                   \
                    src/Quadratic/Cube/CubeNucube_ZZ.cpp                    \
                    src/Quadratic/Cube/CubeNucube_long.cpp                  \
                    src/Quadratic/ClassGroup/ClassGroupBSGSReal.cpp         \
                    src/Quadratic/Multiply/MultiplyComp_ZZ.cpp              \
                    src/Quadratic/Multiply/MultiplyComp_long.cpp            \
                    src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp            \
                    src/Quadratic/Multiply/MultiplyNucomp_long.cpp          \
                    src/Quadratic/Multiply/MultiplyNucomp_Opt_ZZ.cpp        \
                    src/Quadratic/Multiply/MultiplyNucomp_Opt_long.cpp      \
                    src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp             \
                    src/Quadratic/Multiply/MultiplyPlain_long.cpp           \
                    src/Quadratic/Regulator/RegulatorLenstra_ZZ.cpp         \
                    src/Quadratic/Regulator/RegulatorLenstra_long.cpp       \
                    src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp             \
                    src/Quadratic/Reduce/ReducePlainReal_long.cpp           \
                    src/Quadratic/Reduce/ReducePlainReal_Opt_ZZ.cpp         \
                    src/Quadratic/Reduce/ReducePlainReal_Opt_long.cpp       \
                    src/Quadratic/Square/SquareNudupl_ZZ.cpp                \
                    src/Quadratic/Square/SquareNudupl_long.cpp              \
                    src/Quadratic/Square/SquareNudupl_Opt_ZZ.cpp            \
                    src/Quadratic/Square/SquareNudupl_Opt_long.cpp          \
                    src/Quadratic/Square/SquarePlain_ZZ.cpp                 \
                    src/Quadratic/Square/SquarePlain_long.cpp               \
                    src/thresholds.cpp                                      \
                    src/XGCD/hxgcd.cpp                                      \
                    src/XGCD/xgcd.cpp                                       \
                    src/XGCD/xgcd_iter.cpp                                  \
                    src/XGCD/xgcd_plain.cpp

classgroupBSGS_SOURCES = tests/UnitTests.cpp                                \
                    tests/Quadratic/ClassGroup/ClassGroupBSGSReal_Tests.cpp \
                    tests/Quadratic/TestData/TestData.cpp                   \
                    src/common.cpp                                          \
                    src/LinearAlgebra/Smith.cpp                             \
                    src/L_function/Character.cpp                            \
                    src/L_function/L_function_util.cpp                      \
                    src/L_function/L_function_ZZ.cpp                        \
                    src/L_function/L_function_long.cpp                      \
                    src/Quadratic/QuadraticIdealBase_ZZ.cpp                 \
                    src/Quadratic/QuadraticIdealBase_long.cpp               \
                    src/Quadratic/QuadraticInfElement_ZZ.cpp                \
                    src/Quadratic/QuadraticInfElement_long.cpp              \
                    src/Quadratic/QuadraticNumber_ZZ.cpp                    \
                    src/Quadratic/QuadraticNumber_long.cpp                  \
                    src/Quadratic/QuadraticOrder_ZZ.cpp                     \
                    src/Quadratic/QuadraticOrder_long.cpp                   \
                    src/Quadratic/Cube/CubePlain_ZZ.cpp                     \
                    src/Quadratic/Cube/CubePlain_long.cpp                   \
                    src/Quadratic/Cube/CubeNucube_ZZ.cpp                    \
                    src/Quadratic/Cube/CubeNucube_long.cpp                  \
                    src/Quadratic/ClassGroup/ClassGroupBSGSReal.cpp         \
                    src/Quadratic/Multiply/MultiplyComp_ZZ.cpp              \
                    src/Quadratic/Multiply/MultiplyComp_long.cpp            \
                    src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp            \
                    src/Quadratic/Multiply/MultiplyNucomp_long.cpp          \
                    src/Quadratic/Multiply/MultiplyNucomp_Opt_ZZ.cpp        \
                    src/Quadratic/Multiply/MultiplyNucomp_Opt_long.cpp      \
                    src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp             \
                    src/Quadratic/Multiply/MultiplyPlain_long.cpp           \
                    src/Quadratic/Regulator/RegulatorLenstra_ZZ.cpp         \
                    src/Quadratic/Regulator/RegulatorLenstra_long.cpp       \
                    src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp             \
                    src/Quadratic/Reduce/ReducePlainReal_long.cpp           \
                    src/Quadratic/Reduce/ReducePlainReal_Opt_ZZ.cpp         \
                    src/Quadratic/Reduce/ReducePlainReal_Opt_long.cpp       \
                    src/Quadratic/Square/SquareNudupl_ZZ.cpp                \
                    src/Quadratic/Square/SquareNudupl_long.cpp              \
                    src/Quadratic/Square/SquareNudupl_Opt_ZZ.cpp            \
                    src/Quadratic/Square/SquareNudupl_Opt_long.cpp          \
                    src/Quadratic/Square/SquarePlain_ZZ.cpp                 \
                    src/Quadratic/Square/SquarePlain_long.cpp               \
                    src/thresholds.cpp                                      \
                    src/XGCD/hxgcd.cpp                                      \
                    src/XGCD/xgcd.cpp                                       \
                    src/XGCD/xgcd_iter.cpp                                  \
                    src/XGCD/xgcd_plain.cpp

classgroupBS_SOURCES = tests/UnitTests.cpp                                 \
                    tests/Quadratic/TestData/TestData.cpp                    \
                    tests/Quadratic/ClassGroup/ClassGroupBSReal_Tests.cpp  \
                    src/common.cpp                                          \
                    src/LinearAlgebra/Smith.cpp                             \
                    src/L_function/Character.cpp                            \
                    src/L_function/L_function_util.cpp                      \
                    src/L_function/L_function_ZZ.cpp                        \
                    src/L_function/L_function_long.cpp                      \
                    src/Quadratic/QuadraticIdealBase_ZZ.cpp                 \
                    src/Quadratic/QuadraticIdealBase_long.cpp               \
                    src/Quadratic/QuadraticInfElement_ZZ.cpp                \
                    src/Quadratic/QuadraticInfElement_long.cpp              \
                    src/Quadratic/QuadraticNumber_ZZ.cpp                    \
                    src/Quadratic/QuadraticNumber_long.cpp                  \
                    src/Quadratic/QuadraticOrder_ZZ.cpp                     \
                    src/Quadratic/QuadraticOrder_long.cpp                   \
                    src/Quadratic/Cube/CubePlain_ZZ.cpp                     \
                    src/Quadratic/Cube/CubePlain_long.cpp                   \
                    src/Quadratic/Cube/CubeNucube_ZZ.cpp                    \
                    src/Quadratic/Cube/CubeNucube_long.cpp                  \
                    src/Quadratic/ClassGroup/ClassGroupBSGSReal.cpp         \
                    src/Quadratic/Multiply/MultiplyComp_ZZ.cpp              \
                    src/Quadratic/Multiply/MultiplyComp_long.cpp            \
                    src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp            \
                    src/Quadratic/Multiply/MultiplyNucomp_long.cpp          \
                    src/Quadratic/Multiply/MultiplyNucomp_Opt_ZZ.cpp        \
                    src/Quadratic/Multiply/MultiplyNucomp_Opt_long.cpp      \
                    src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp             \
                    src/Quadratic/Multiply/MultiplyPlain_long.cpp           \
                    src/Quadratic/Regulator/RegulatorLenstra_ZZ.cpp         \
                    src/Quadratic/Regulator/RegulatorLenstra_long.cpp       \
                    src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp             \
                    src/Quadratic/Reduce/ReducePlainReal_long.cpp           \
                    src/Quadratic/Reduce/ReducePlainReal_Opt_ZZ.cpp         \
                    src/Quadratic/Reduce/ReducePlainReal_Opt_long.cpp       \
                    src/Quadratic/Square/SquareNudupl_ZZ.cpp                \
                    src/Quadratic/Square/SquareNudupl_long.cpp              \
                    src/Quadratic/Square/SquareNudupl_Opt_ZZ.cpp            \
                    src/Quadratic/Square/SquareNudupl_Opt_long.cpp          \
                    src/Quadratic/Square/SquarePlain_ZZ.cpp                 \
                    src/Quadratic/Square/SquarePlain_long.cpp               \
                    src/thresholds.cpp                                      \
                    src/XGCD/hxgcd.cpp                                      \
                    src/XGCD/xgcd.cpp                                       \
                    src/XGCD/xgcd_iter.cpp                                  \
                    src/XGCD/xgcd_plain.cpp

#quadinfelm_SOURCES = tests/UnitTests.cpp                                \
                    tests/Quadratic/TestData/TestData.cpp               \
                    tests/Quadratic/QuadraticInfElement_ZZ_Tests.cpp    \
                    tests/Quadratic/Multiply/MultiplyPlain_ZZ_Tests.cpp \
                    tests/Quadratic/Reduce/ReducePlainReal_ZZ_Tests.cpp \
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
                    src/Quadratic/Multiply/MultiplyComp_ZZ.cpp          \
                    src/Quadratic/Multiply/MultiplyNucomp_ZZ.cpp        \
                    src/Quadratic/Multiply/MultiplyNucomp_Opt_ZZ.cpp    \
                    src/Quadratic/Multiply/MultiplyPlain_ZZ.cpp         \
                    src/Quadratic/Reduce/ReducePlainReal_ZZ.cpp         \
                    src/Quadratic/Reduce/ReducePlainReal_Opt_ZZ.cpp     \
                    src/Quadratic/Square/SquareNudupl_ZZ.cpp            \
                    src/Quadratic/Square/SquareNudupl_Opt_ZZ.cpp        \
                    src/Quadratic/Square/SquarePlain_ZZ.cpp             \
                    src/thresholds.cpp                                  \
                    src/XGCD/hxgcd.cpp                                  \
                    src/XGCD/xgcd.cpp                                   \
                    src/XGCD/xgcd_iter.cpp                              \
                    src/XGCD/xgcd_plain.cpp

#cubic_SOURCES=tests/Cubic/cubicTestMain.cpp src/Cubic/generalFunctions.cpp src/Cubic/GlobalCubicField.cpp src/Cubic/CubicNumberField.cpp src/Cubic/RealCubicNumberField.cpp src/Cubic/ComplexCubicNumberField.cpp src/Cubic/CubicOrder.cpp src/Cubic/CubicOrderReal.cpp src/Cubic/CubicElement.cpp src/Cubic/CubicIdeal.cpp src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp src/Cubic/VoronoiMethods.cpp src/Cubic/VoronoiReal.cpp src/Cubic/VoronoiComplex.cpp src/Cubic/FundamentalUnits/BasicVoronoi.cpp src/Cubic/FundamentalUnits/BSGSVoronoi.cpp

#testorder_SOURCES=tests/Cubic/catch-tests/test-order.cpp tests/Cubic/catch-tests/testingMain.cpp src/Cubic/generalFunctions.cpp src/Cubic/CubicOrder.cppquadinfelm

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
