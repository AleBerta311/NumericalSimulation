CC = g++
CFLAGS = -Wall -O3 --std=c++11

buffon.exe : buffon.o random.o
	$(CC) random.o buffon.o -o buffon.exe
buffon.o : Buffon.cpp
	$(CC) -c Buffon.cpp -o buffon.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt 
