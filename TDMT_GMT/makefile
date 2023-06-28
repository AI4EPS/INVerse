# Standalone Makefile for Solaris and Sun Studio 12 compiler suite

cc	= gcc
CFLAGS  = -m64
FFLAGS  = 
FF      = gfortran

########################################################################
# If you want plotting within tdmt_invc using the stdplt library, 
# define the following; otherwise comment them out:
#:: STDPLTDIR = ../stdplt/plotsubs
#:: PLTLIB    = -L$(STDPLTDIR) -lstdplt
#:: TDMT_PLOT_SOURCE = mt_plot6iso2_linux.c
#:: TDMT_PLOT_SWITCH = -DTDMT_PLOT

# To turn off all exceptt errors being written to stderr, 
# define the following; otherwise comment it out
#TDMT_STDERR_SWITCH = -DONLY_ERRORS_TO_STDERR
########################################################################

TDMT_OPT_FLAGS = $(TDMT_PLOT_SWITCH) $(TDMT_STDERR_SWITCH)
LIBS 	= $(PLTLIB) -lm

TARGET1 = sac2helm 
CSRCS1 = sac2helm.c getpar.c 
TARGET2 = tdmt_invc_iso
CSRCS2 = tdmt_invc_iso.c minvdbl.c correl2b.c readhelm.c fitcheck2_iso.c \
     mt_plot_dat.c ${TDMT_PLOT_SOURCE}
CSRCS3 = tdmt_invc_iso_fwd.c minvdbl.c correl2b.c readhelm.c fitcheck2_iso.c \
     mt_plot_dat.c ${TDMT_PLOT_SOURCE}
CSRCS4 = tdmt_invc_iso_out.c minvdbl.c correl2b.c readhelm.c fitcheck2_iso.c \
     mt_plot_dat.c ${TDMT_PLOT_SOURCE}
CSRCS5 = tdmt_invc_iso_zcor.c minvdbl.c correl2b.c readhelm.c fitcheck2_iso.c \
     mt_plot_dat.c ${TDMT_PLOT_SOURCE}
CSRCS6 = tdmt_invc_iso_outputmtline.c minvdbl.c correl2b.c readhelm.c fitcheck2_iso.c \
     mt_plot_dat.c ${TDMT_PLOT_SOURCE}
TARGET5 = tdmt_invc_iso_fwd
TARGET6 = tdmt_invc_iso_out
TARGET7 = tdmt_invc_iso_zcor
TARGET8 = tdmt_invc_iso_outputmtline
TARGET3	= fmap
FSRCS3	= fmap_new.f
TARGET4	= mtmanip
FSRCS4	= mtmanip_new.f
FSRCSUB = mtmanip_module.f
FSRCSUB2 = fmap_subs_new.f

ALL	= $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8)
OBJS1	= $(CSRCS1:%.c=%.o)
OBJS2	= $(CSRCS2:%.c=%.o)
OBJS3	= $(FSRCS3:%.f=%.o)
OBJS4	= $(FSRCS4:%.f=%.o)
OBJS5	= $(CSRCS3:%.c=%.o)
OBJS6	= $(CSRCS4:%.c=%.o)
OBJS7	= $(CSRCS5:%.c=%.o)
OBJS8	= $(CSRCS6:%.c=%.o)
FOBJSUB	= $(FSRCSUB:%.f=%.o)
FOBJSUB2= $(FSRCSUB2:%.f=%.o)

all:	$(ALL)

$(TARGET1): $(OBJS1)
	$(cc) $(CFLAGS) -o $@ $(OBJS1)

$(TARGET2): $(FOBJSUB) $(OBJS2) $(FOBJSUB2)
	$(FF) $(CFLAGS) -o $@ $(OBJS2) $(FOBJSUB2) $(FOBJSUB) $(LIBS)

$(TARGET3): $(FOBJSUB) $(OBJS3) $(FOBJSUB2)
	$(FF) $(CFLAGS) -o $@ $(OBJS3) $(FOBJSUB2) $(FOBJSUB)

$(TARGET4): $(FOBJSUB2) $(OBJS4) $(FOBJSUB)
	$(FF) $(CFLAGS) -o $@ $(OBJS4) $(FOBJSUB) $(FOBJSUB2)
	
$(TARGET5): $(FOBJSUB) $(OBJS5) $(FOBJSUB2)
	$(FF) $(CFLAGS) -o $@ $(OBJS5) $(FOBJSUB2) $(FOBJSUB) $(LIBS)

$(TARGET6): $(FOBJSUB) $(OBJS6) $(FOBJSUB2)
	$(FF) $(CFLAGS) -o $@ $(OBJS6) $(FOBJSUB2) $(FOBJSUB) $(LIBS)

$(TARGET7): $(FOBJSUB) $(OBJS7) $(FOBJSUB2)
	$(FF) $(CFLAGS) -o $@ $(OBJS7) $(FOBJSUB2) $(FOBJSUB) $(LIBS)

$(TARGET8): $(FOBJSUB) $(OBJS8) $(FOBJSUB2)
	$(FF) $(CFLAGS) -o $@ $(OBJS8) $(FOBJSUB2) $(FOBJSUB) $(LIBS)

.c.o:
	$(cc) $< -c $(CFLAGS) $(TDMT_OPT_FLAGS)

.f.o:
	$(FF) $< -c $(CFLAGS)

clean:
	\rm *.o tdmt_invc_iso tdmt_invc_iso_fwd mtmanip sac2helm fmap

install:	$(ALL) $(SCRIPTS)
	\mv tdmt_invc_iso mtmanip tdmt_invc_iso_fwd tdmt_invc_iso_out tdmt_invc_iso_outputmtline tdmt_invc_iso_zcor ../BIN

