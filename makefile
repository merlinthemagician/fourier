#	Makefile -- Sat Apr  5 14:50:06 MET DST 1997
#	Copyright (c) 1992, 1997 Axel T. Schreiner

# include /people/isie002/include/make/c.h
# include /Users/merlin/c/include/make/c.h
include ${HOME}/c/include/make/c.h

D =	Fourier-Transformation

T = test$x grf$x grf.o

INC=/sw/include/

all:	$T

run::	test$x;  ./test

test$x:	test.c;
	$(CC) $(CFLAGS) -o $@  test.c $(GSL) -I$(INC)

grf$x:	grf.c;
	$(CC) $(CFLAGS) -DWHITE -DTEST -o $@  grf.c -I$(INC) $(GSL)

grf.o:	grf.c;
	$(CC) $(CFLAGS) -c -o $@  grf.c -I$(INC)

grftest$x:	grf.o grftest.c;
	$(CC) $(CFLAGS) -o $@  grf.o grftest.c -I$(INC) $(GSL)
