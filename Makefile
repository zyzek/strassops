CC = clang
CFLAGS = -O1 -std=gnu11 -Wall -Werror
LDFLAGS = -pthread

.PHONY: all clean

all: matrix

matrix: main.c matrix.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

clean:
	-rm -f *.o
	-rm -f matrix
