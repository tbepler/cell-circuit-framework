CC = gcc
CFLAGS = -Wall -Ofast -lm
OBJS = c_utils/array.o c_utils/file.o c_utils/optimize.o

all: RelaxationOscillator.o
	$(CC) $(OBJS) RelaxationOscillator.o -o simulate.out $(CFLAGS)

%.o: %.c
	$(CC) -c $<


