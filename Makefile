CC = clang
CFLAGS = -O1 -std=gnu11 -Wall -Werror
LDFLAGS = -pthread

.PHONY: all clean

all: matrix

matrix: main.c matrix.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

profile: main.c matrix.c
	gcc -pg -g -pthread -std=gnu11 main.c matrix.c -o matrix

clean:
	-rm -f *.o
	-rm -f matrix
