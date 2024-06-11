CC = g++
CFLAGS = -Wall -O3 --std=c++11

distributions.exe : distributions.o random.o
	$(CC) random.o distributions.o -o distributions.exe
distributions.o : distributions.cpp
	$(CC) -c distributions.cpp -o distributions.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o test.exe *.txt chi2.exe
