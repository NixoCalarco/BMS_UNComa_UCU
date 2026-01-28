# compilador
CC = gcc
SHELL := /bin/bash

BUILD_DIR = build

# archivitos
SRC_DIR = .
SRC = $(addprefix $(SRC_DIR)/, main.c matrix.c emf.c )

# ejecutable principal
BIN = $(BUILD_DIR)/bin

FLAGS_LINKING = -lm
FLAGS_COMMON  = -Werror -Wall -Wextra -Wpedantic -std=c23 
FLAGS_RELEASE = -O3 -march=native
FLAGS_TESTING = -g -O0
FLAGS_SANITIZER = -fsanitize=address,undefined 

_release: FLAGS = $(FLAGS_COMMON) $(FLAGS_RELEASE) $(FLAGS_LINKING)
_release: $(BIN)

_testing: FLAGS = $(FLAGS_COMMON) $(FLAGS_TESTING) $(FLAGS_LINKING)
_testing: $(BIN)

_sanitizer: FLAGS = $(FLAGS_COMMON) $(FLAGS_TESTING) $(FLAGS_SANITIZER) $(FLAGS_LINKING)
_sanitizer: $(BIN)


$(BIN): $(SRC)
	$(CC) $^ -o $@ $(FLAGS)

clean:
	-rm $(BIN)

test: $(BIN)
	valgrind --leak-check=full ./$(BIN)

run:
	./$(BIN)

release:
	$(MAKE) clean
	$(MAKE) _release

testing:
	$(MAKE) clean
	$(MAKE) _testing

sanitizer:
	$(MAKE) clean
	$(MAKE) _sanitizer

.PHONY: plot
plot:
	./$(BIN) > ./plot/output.py
	cd plot; uv run plt_c_sim.py

