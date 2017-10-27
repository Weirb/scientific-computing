all: makeme

makeme: functions.o main.o
	g++ -o out functions.o main.o

functions.o: functions.cpp
	g++ -c functions.cpp

main.o: main.cpp
	g++ -c main.cpp

.PHONY: clean
clean:
	rm *.o out