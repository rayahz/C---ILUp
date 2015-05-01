MPI=mpichcc
MPI_RUN=mpichexec -n $(NBR)
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
	$(MPI) $(LDFLAGS) -o $(BIN_DIR)$(TARGET) $(OBJ_FILES) $(LFLAGS) -Wall 

$(OBJ_DIR)%.o: $(SRC_DIR)%.c
	$(MPI) $(CFLAGS) -c -o $@ $^ 

$(EXEC): 
	$(MPI_RUN) $(BIN_DIR)$(TARGET) 

clean:
	$(RM) $(BIN_DIR)$(TARGET) $(OBJ_FILES)
	
