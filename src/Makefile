CC = gcc
CC_FLAGS =  -O2  

#Parallel version
#CC =mpicc 
#CC_FLAGS =  -O2  -DPARALLEL

LD_FLAGS =-lm 


OBJECTS= main.o initialize.o createLattice.o ngbrs.o cellMethod.o fileio.o  potentialWww.o mcSubs.o conjgrad.o nrMinSubs.o  parallel.o randFunc.o  bondList.o  miscSubs.o rings.o writeMat.o wwwMc.o tabulatedValues.o

Final:	$(OBJECTS)
	$(CC) $(CC_FLAGS) *.o -o bomc $(LD_FLAGS)

%.o: %.c 
	$(CC) $(CC_FLAGS) -c   $< 

clean:
	rm -f *.o 



