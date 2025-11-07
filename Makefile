# Makefile para Sistema de Refrigeración por Absorción

# Compilador y flags
CC = gcc
CFLAGS = -Wall -Wno-deprecated-declarations
LIBS = -framework OpenGL -framework GLUT -lm

# Targets principales
all: refrigerador

# Ejecutable principal
refrigerador: main.o
	$(CC) -o $@ $^ $(LIBS)

# Objeto principal
main.o: main.c
	$(CC) $(CFLAGS) -c main.c

# Limpieza
clean:
	rm -f *.o refrigerador

# Ejecutar
run: refrigerador
	./refrigerador

.PHONY: all clean run