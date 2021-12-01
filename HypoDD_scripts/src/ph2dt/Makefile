CMD	= ph2dt
CC	= gcc
FC	= gfortran
#FC	= f77
#FC	= gfortran
SRCS	= $(CMD).f \
	  cluster.f datetime.f delaz.f ifindi.f \
	  indexx.f indexxi.f sorti.f trimlen.f
CSRCS	= atoangle_.c atoangle.c rpad_.c sscanf3_.c
OBJS	= $(SRCS:%.f=%.o) $(CSRCS:%.c=%.o)
INCLDIR	= ../../include
CFLAGS	= -I$(INCLDIR) -O

# Flags for GNU g77 compiler
#FFLAGS	= -O -I$(INCLDIR) -fno-silent -Wall -implicit

# Flags for SUN f77 compiler
FFLAGS	= -O -I$(INCLDIR)

# Following line needed on HP-UX (hasn't been tested, though).
#LDFLAGS	= +U77

all: $(CMD)

$(CMD): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@

clean:
	-rm -f $(CMD) *.o core a.out junk

atoangle_.o : ../../include/compat.h
atoangle_.o : ../../include/f77types.h
ph2dt.o : ../../include/ph2dt.inc
