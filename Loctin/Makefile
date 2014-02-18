############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Compilation & chargement de ex7
##
## Auteur          : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Créé le         : Mon Nov 30 16:43:26 1998
## Dern. mod. par  : Jalel Chergui (CNRS/IDRIS) <Jalel.Chergui@idris.fr>
## Dern. mod. le   : Wed Jan  4 11:48:52 2006
###############################################################################

ROOT    = ex7
CMD     = $(ROOT).x
SRC     = $(ROOT).c
OBJ     = $(ROOT).o

FC	= mpicc
FFLAGS	= -O3
LDFLAGS	=
NP	= 4

.SUFFIXES : .c .o

all     : $(CMD)

clean   :; rm -f *.o $(CMD) core

$(OBJ) : $(SRC)

.c.o    :; $(FC) $(FFLAGS) -c $<

$(CMD)  : $(OBJ)
	$(FC) $(LDFLAGS) -o $(@) $(OBJ)

exe: $(CMD)
	mpirun --mca btl_tcp_if_include eth0 -np $(NP) --hostfile ../hostfile $(CMD)
