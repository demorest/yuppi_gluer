# define some Makefile variables for the compiler and compiler flags
# to use Makefile variables later in the Makefile: $()
#
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#
# for C++ define  CC = g++
CC = gcc
CFLAGS  = -g -Wall -O3 #-I.
LIBS = -lm -L/usr/lib -lcfitsio
OBJS = main.o chkio.o multifiles.o sigproc_fb2.o main_cmd.o vectors2.o \
	write_psrfits.o read_psrfits.o polyco.o downsample.o fill_psrfits_struct.o


# typing 'make' will invoke the first target entry in the file 
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
default: fil2fits

fil2fits:  $(OBJS)
	$(CC) $(CFLAGS) -o fil2fits $(OBJS) $(LIBS)

main.o: main.c sigproc_fb2.h presto2.h main_cmd.h vectors2.h psrfits.h
	$(CC) $(CFLAGS) -c main.c

main_cmd.o: main_cmd.c main_cmd.h
	$(CC) $(CFLAGS) -c main_cmd.c

chkio.o:  chkio.c chkio.h
	$(CC) $(CFLAGS) -c chkio.c

multifiles.o: multifiles.c multifiles.h
	$(CC) $(CFLAGS) -c multifiles.c

sigproc_fb2.o: sigproc_fb2.c sigproc_fb2.h presto2.h
	$(CC) $(CFLAGS) -c sigproc_fb2.c

vectors2.o: vectors2.c vectors2.h
	$(CC) $(CFLAGS) -c vectors2.c

write_psrfits.o: write_psrfits.c psrfits.h polyco.h downsample.h
	$(CC) $(CFLAGS) -c write_psrfits.c

read_psrfits.o: read_psrfits.c psrfits.h downsample.h
	$(CC) $(CFLAGS) -c read_psrfits.c

polyco.o: polyco.c polyco.h polyco_struct.h psrfits.h
	$(CC) $(CFLAGS) -c polyco.c

downsample.o: downsample.c psrfits.h downsample.h
	$(CC) $(CFLAGS) -c downsample.c

fill_psrfits_struct.o: fill_psrfits_struct.c chkio.h vectors2.h psrfits.h
	$(CC) $(CFLAGS) -c fill_psrfits_struct.c

.PHONY: clean

clean: 
	$(RM) fil2fits *.o *~
