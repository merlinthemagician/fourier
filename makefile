##	Makefile
#

CC      = gcc#
CFLAGS	=-Wall -pedantic#

rm	= rm -f#				delete file

all:	_
clean::		;@ $(MAKE) T='$T' _clean
_clean:	_	;  $(rm) *.o $T a.out core *.tmp *.ps *.bak
run::	_
_:		;@ echo -------------------- $D --------------------

D =	Generate coloured noise

T = test$x grf$x grf.o testSpatTemp$x grf_field.o

INC=/sw/include/
GSL = `gsl-config --cflags` `gsl-config --libs` -lm

all:	$T

run::	test$x;  ./test

test$x:	test.c;
	$(CC) $(CFLAGS) -o $@  test.c $(GSL) -I$(INC)

grf$x:	grf.c;
	$(CC) $(CFLAGS) -DWHITE -DTEST -o $@  grf.c -I$(INC) $(GSL)

grf.o:	grf.c;
	$(CC) $(CFLAGS) -c -o $@  grf.c -I$(INC)

grf_field.o:	grf_field.c;
		$(CC) $(CFLAGS) -c -o $@  grf_field.c -I$(INC)


grftest$x:	grf.o grftest.c;
	$(CC) $(CFLAGS) -o $@  grf.o grftest.c -I$(INC) $(GSL)

testSpatTemp$x:	grf.o testSpatTemp.c;
	$(CC) $(CFLAGS) -o $@ grf.o testSpatTemp.c $(GSL) -I$(INC)
