ALLFLAGS= -m64
CFLAGS	= 
FFLAGS 	= 

ALL	= sac2helm window bin2sac sac2bin fromHelm mkHelm putmech putmt decon sac2columns sac2columns_mt
PROGLIST= get_resp sac2helm window bin2sac sac2bin fromHelm mkHelm putmech putmtsac2columns sac2columns_mt 

all:	$(ALL)

sac2helm: sac2helm.c getpar.c
	gcc $(ALLFLAGS) -o sac2helm sac2helm.c getpar.c
	\mv sac2helm ../BIN

sac2bin: sac2bin.c getpar.c
	gcc $(ALLFLAGS) -o sac2bin sac2bin.c getpar.c
	\mv sac2bin ../BIN

putmech:  putmech_iso.c readhelm.c getpar.c
	gcc $(ALLFLAGS) -g -o putmech putmech_iso.c readhelm.c getpar.c -lm
	\mv putmech ../BIN

putmt:  putmt.c readhelm.c getpar.c
	gcc $(ALLFLAGS) -g -o putmt putmt.c readhelm.c getpar.c -lm
	\mv putmt ../BIN

window: window.c getpar.c
	gcc $(ALLFLAGS) -o window window.c getpar.c 
	\mv window ../BIN


fromHelm: fromHelm_new.c
	gcc $(ALLFLAGS) -o fromHelm fromHelm_new.c
	\mv fromHelm ../BIN

mkHelm: mkHelm2.c getpar.c
	gcc $(ALLFLAGS) -o mkHelm mkHelm2.c getpar.c
	\mv mkHelm ../BIN

decon: deconnew.c deconnew_subs2.c getpar.c
	gcc $(ALLFLAGS) -o decon deconnew.c deconnew_subs2.c getpar.c -lm
	\mv decon ../BIN

bin2sac: bin2sac.c getpar.c
	gcc $(ALLFLAGS) -o bin2sac bin2sac.c getpar.c
	\mv bin2sac ../BIN

sac2columns: sac2columns.c
	gcc $(ALLFLAGS) -o sac2columns sac2columns.c
	\mv sac2columns ../BIN

sac2columns_mt: sac2columns_mt.c
	gcc $(ALLFLAGS) -o sac2columns_mt sac2columns_mt.c
	\mv sac2columns_mt ../BIN

##Cleanup routine

clean:
	\rm *.o
