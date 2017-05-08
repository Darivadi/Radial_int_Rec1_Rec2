CC = gcc
CFLAGSDEBUG = -g -Wall -c -I/home/$(USER)/local/include/ -I/usr/include/ 
CFLAGS = -c -O3 -I$(HOME)/local/include/ -I/usr/include/  -DCOMPLETERAY -DCOMPLETERAY_REC1 
CFLAGSUNTIL = -c -O3 -I$(HOME)/local/include/ -I/usr/include/ -DUNTIL90MPC -DCOMPLETERAY_REC1 #-DUNTIL90MPC_REC1
CFLAGSFROM = -c -O3 -I$(HOME)/local/include/ -I/usr/include/ -DFROM90MPC -DCOMPLETERAY_REC1 #-DFROM90MPC_REC1
LFLAGS = -lm -L$(HOME)/local/lib -Wl,"-R /export/$(USER)/local/lib"

PROGRAM = main_radial_ISW_integral


$(PROGRAM):
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -lgsl -lgslcblas -lm -o $@
	mv main_radial_ISW_integral main_radial_ISW.x


UNTIL90MPC:
	$(CC) $(CFLAGSUNTIL) $(PROGRAM).c -o $(PROGRAM).o
	$(CC) $(PROGRAM).o $(LFLAGS) -lgsl -lgslcblas -lm -o $(PROGRAM)
	mv main_radial_ISW_integral main_radial_ISW.x


FROM90MPC:
	$(CC) $(CFLAGSFROM) $(PROGRAM).c -o $(PROGRAM).o
	$(CC) $(PROGRAM).o $(LFLAGS) -lgsl -lgslcblas -lm -o $(PROGRAM)
	mv main_radial_ISW_integral main_radial_ISW.x

debug:
	echo Compiling for debug $(PROGRAM).c
	$(CC) $(CFLAGSDEBUG) $(PROGRAM).c -o $(PROGRAM).o
	$(CC) $(PROGRAM).o $(LFLAGS) -lgsl -lgslcblas -lm -o $(PROGRAM).x


clean:
	rm -rf $(PROGRAM)
	rm -rf *~
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a      
	rm -rf *.so
	rm *.x
