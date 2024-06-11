CC = g++
CFLAGS = -O0 --std=c++11 -g
AFLAGS = -larmadillo
all: remove clean genetic.exe
gen: genetic.exe
genetic.exe : genetic.o random.o city.o path.o population.o
	$(CC) random.o path.o population.o city.o genetic.o -o genetic.exe $(AFLAGS)
genetic.o: genetic.cpp random.h city.h Path.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)
path.o : Path.cpp Path.h
	$(CC) -c Path.cpp -o path.o $(CFLAGS)
population.o : population.cpp Population.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
clean :
	rm -f *.o
remove:
	rm -f *.dat
