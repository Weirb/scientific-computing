CC=g++
CFLAGS_DEBUG=-c -Wall
CFLAGS_RELEASE=-c -O2 -ffast-math -march=native
CFLAGS=$(CFLAGS_DEBUG)
OBJ_PATH=obj
SRC_PATH=src
EXECUTABLE=cg.exe
FILES=functions.cpp main.cpp
OBJECTS = $(addprefix $(OBJ_PATH)/, $(FILES:.cpp=.o))

all: $(OBJECTS)
	$(CC) -o $(EXECUTABLE) $(OBJECTS)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.cpp
	$(CC) $(CFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm $(OBJECTS) $(EXECUTABLE)
