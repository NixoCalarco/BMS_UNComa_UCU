#!/bin/bash
gcc main.c matrix.c emf.c -Wall -Werror -Wextra -Wpedantic -std=c23 -Werror=vla -O3 -o build/bin
