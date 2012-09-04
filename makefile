#	Makefile -- Sat Apr  5 14:50:06 MET DST 1997
#	Copyright (c) 1992, 1997 Axel T. Schreiner
include /people/isie002/include/make/c.h

D =	Fourier-Transformation

T = test$x grf$x grf.o

all:	$T

run::	test$x;  ./test

test$x:	test.c;
	$(CC) $(CFLAGS) -o $@  test.c $(GSL)

grf$x:	grf.c;
	$(CC) $(CFLAGS) -DWHITE -DTEST -o $@  grf.c $(GSL)

grf.o:	grf.c;
	$(CC) $(CFLAGS) -c -o $@  grf.c

grftest$x:	grf.o grftest.c;
	$(CC) $(CFLAGS) -o $@  grf.o grftest.c $(GSL)
