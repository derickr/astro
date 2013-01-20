CFLAGS=-DHAVE_CBRT -ggdb3 -O0 -Wall
LDFLAGS=-lm

all: main

main: astro.o

clean:
	rm *.o main
