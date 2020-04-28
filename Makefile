CFLAGS=-DHAVE_CBRT -ggdb3 -O0 -Wall
LDFLAGS=-lm
OBJECTS=astro.o main.o

all: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o main $(LDFLAGS)

clean:
	rm *.o main
