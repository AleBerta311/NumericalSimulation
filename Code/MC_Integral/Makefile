CC = g++
CFLAGS = -Wall -O3 --std=c++11

mcintegral.exe : mcintegral.o random.o
	$(CC) random.o mcintegral.o -o mcintegral.exe
mcintegral.o : MC_Integral.cpp
	$(CC) -c MC_Integral.cpp -o mcintegral.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt 
