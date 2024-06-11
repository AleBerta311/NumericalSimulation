CC = g++
CFLAGS = -Wall -O3 --std=c++11

test.exe : RandomNumberTest.o random.o
	$(CC) random.o RandomNumberTest.o -o test.exe
test.o : RandomNumberTest.cpp
	$(CC) -c RandomNumberTest.cpp -o RandomNumberTest.o $(CFLAGS)
chi2.exe : chi2.o random.o
	$(CC) random.o chi2.o -o chi2.exe
chi2.o : Chi2Test.cpp
	$(CC) -c Chi2Test.cpp -o chi2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o test.exe *.txt chi2.exe
