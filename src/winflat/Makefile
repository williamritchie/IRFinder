all: winflat test 

winflat: winflat_with_beta.c
	$(CC) -o winflat winflat_with_beta.c -lm 

test:	winflat
	sh runtest.sh 2> /dev/null
