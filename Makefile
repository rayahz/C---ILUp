MPI=mpichcc
CC=gcc
RM=rm
AR=ar
LDFLAGS=
CFLAGS=-Iinclude/
LFLAGS=-lm
OBJ_FILES=$(OBJ_DIR)main.o $(OBJ_DIR)projet.o
TARGET=projet
EXEC=proj
SRC_DIR=./src/
OBJ_DIR=./obj/
BIN_DIR=./bin/

all: $(TARGET)

$(TARGET): $(OBJ_FILES)
	$(CC) $(LDFLAGS) -o $(BIN_DIR)$(TARGET) $(OBJ_FILES) $(LFLAGS) -Wall -fopenmp

$(OBJ_DIR)%.o: $(SRC_DIR)%.c
	$(CC) $(CFLAGS) -c -o $@ $^ 

$(EXEC): 
	$(BIN_DIR)$(TARGET) 

clean:
	$(RM) $(BIN_DIR)$(TARGET) $(OBJ_FILES)
	
