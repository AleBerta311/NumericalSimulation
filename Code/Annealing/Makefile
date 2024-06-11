CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: anneal.exe
anneal.exe : anneal.o random.o
	$(CC) random.o anneal.o -o anneal.exe $(AFLAGS)
anneal.o: anneal.cpp random.h
	$(CC) -c anneal.cpp -o anneal.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulator.exe
remove:
	rm *.dat
