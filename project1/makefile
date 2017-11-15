CC=g++
CFLAGS_DEBUG=-c -Wall
CFLAGS_RELEASE=-c -O2 -ffast-math -march=native
CFLAGS=$(CFLAGS_RELEASE)
OBJ_PATH=obj
SRC_PATH=src
FIGURE_PATH=figures
EXECUTABLE=cg.exe
FILES=functions.cpp main.cpp
OBJECTS = $(addprefix $(OBJ_PATH)/, $(FILES:.cpp=.o))

all: $(OBJECTS)
	$(CC) -o $(EXECUTABLE) $(OBJECTS)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.cpp
	mkdir -p $(OBJ_PATH)
	$(CC) $(CFLAGS) -o $@ $<

figures:
	# $(EXECUTABLE)
	mkdir -p $(FIGURE_PATH);
	matlab -r 'gen_plots(); exit;';

.PHONY: clean
clean:
	rm -rf $(OBJ_PATH) $(FIGURE_PATH) $(EXECUTABLE)
