CC = g++
CFLAG = -O3 -Wall -std=c++0x

hirgc: hirgc.cpp
	$(CC) hirgc.cpp -o hirgc $(CFLAG)

de_hirgc: de_hirgc.cpp
	$(CC) de_hirgc.cpp -o de_hirgc $(CFLAG)

clean:
	rm *.o