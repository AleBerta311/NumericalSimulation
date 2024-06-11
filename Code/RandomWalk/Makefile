CC = g++
CFLAGS = -Wall -O3 --std=c++11

randomwalk.exe : randomwalk.o random.o
	$(CC) random.o randomwalk.o -o randomwalk.exe
randomwalk.o : random_walk.cpp
	$(CC) -c random_walk.cpp -o randomwalk.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt 
