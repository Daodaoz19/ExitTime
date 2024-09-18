OBJS = exittime

all: ${OBJS}

exittime: exittime.o
	 g++ -o exittime exittime.o

exittime.o: exittime.cpp
	gcc -c exittime.cpp -std=c++0x

clean : 
	rm -f ${OBJS} *.o core*
