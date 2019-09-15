main:  xgcd_test.o xgcd_plain.o   # xgcd.o hxgcd.o xgcd_iter.o xgcd_pseudo.o
	g++ -g -o testFile -O2 -I/c/include xgcd_test.o xgcd_plain.o  -L/c/lib -lntl -lgmp -lm -lpthread
#xgcd.o hxgcd.o xgcd_iter.o xgcd_pseudo.o

xgcd_test.o: ../../include/ANTL/XGCD/xgcd_plain.hpp
	g++ -g -c -O2 -I/c/include xgcd_test.cpp -L/c/lib -lntl -lgmp -lm -lpthread

xgcd_plain.o: ../../src/XGCD/xgcd_plain.cpp ../../include/ANTL/XGCD/xgcd_plain.hpp
	g++ -g -c -O2 -I/c/include ../../src/XGCD/xgcd_plain.cpp -L/c/lib -lntl -lgmp -lm -lpthread


#xgcd.o: ../../src/XGCD/xgcd.cpp ../../include/ANTL/XGCD/xgcd.hpp
#	g++ -g -c -O2 -I/c/include ../../src/XGCD/xgcd.cpp -L/c/lib -lntl -lgmp -lm -lpthread

#hxgcd.o: ../../src/XGCD/hxgcd.cpp ../../include/ANTL/XGCD/hxgcd.hpp
#	g++ -g -c -O2 -I/c/include ../../src/XGCD/hxgcd.cpp -L/c/lib -lntl -lgmp -lm -lpthread

#xgcd_iter.o: ../../src/XGCD/xgcd_iter.cpp ../../include/ANTL/XGCD/xgcd_iter.hpp
#	g++ -g -c -O2 -I/c/include ../../src/XGCD/xgcd_iter.cpp -L/c/lib -lntl -lgmp -lm -lpthread

#xgcd_pseudo.o: ../../src/XGCD/xgcd_pseudo.cpp ../../include/ANTL/XGCD/xgcd_pseudo.hpp
#	g++ -g -c -O2 -I/c/include ../../src/XGCD/xgcd_pseudo.cpp -L/c/lib -lntl -lgmp -lm -lpthread
