CC = g++
CFLAGS = -Wall -O3 --std=c++11

BlackScholes.exe : BlackScholes.o random.o dudu.o
	$(CC) random.o dudu.o BlackScholes.o -o BlackScholes.exe
BlackScholes.o : BlackScholes.cpp
	$(CC) -c BlackScholes.cpp -o BlackScholes.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
dudu.o : DuduMath.cpp DuduMath.h
	$(CC) -c DuduMath.cpp -o dudu.o $(CFLAGS)

BlackScholes_discrete.exe : BlackScholes_discrete.o random.o dudu.o
	$(CC) random.o dudu.o BlackScholes_discrete.o -o BlackScholes_discrete.exe
BlackScholes_discrete.o: BlackScholes_discrete_steps.cpp
	$(CC) -c BlackScholes_discrete_steps.cpp -o BlackScholes_discrete.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt 
