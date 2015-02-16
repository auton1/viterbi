# Make file for viterbi
# Author: Adam Auton
# ($Revision: 230 $)

# Compiler
CC = gcc
CPP = g++
# Output executable
EXECUTABLE = viterbi
# Compiler flags
CFLAGS = -O3 -m64 -D_FILE_OFFSET_BITS=64
CPPFLAGS = -D_FILE_OFFSET_BITS=64 -O3 -std=c++11 -pthread
#CPPFLAGS = -O2 -Wall -g -D_FILE_OFFSET_BITS=64
# Included libraries (zlib)
LIB = -lz #-lrt
#LIB = -lz -I/opt/local/include/ -L/opt/local/lib/
#RLIB = -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/Library/Frameworks/R.framework/Resources/lib -lRlapack  /Library/Frameworks/R.framework/Versions/3.0/Resources/library/RInside/lib/libRInside.a -framework R
#RINC = -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/R.framework/Versions/3.0/Resources/library/Rcpp/include -I/Library/Frameworks/R.framework/Versions/3.0/Resources/library/RInside/include -I/usr/local/include -F/Library/Frameworks/R.framework/.. 

OBJS = viterbi.o bcf_file.o vcf_file.o variant_file.o \
	header.o bcf_entry.o vcf_entry.o entry.o entry_getters.o entry_setters.o \
	vcf_entry_setters.o bcf_entry_setters.o entry_filters.o variant_file_filters.o parameters.o \
	output_log.o 

OBJS2 = bgzf.oc

all: viterbi

viterbi: $(OBJS) $(OBJS2)
	$(CPP) $(CPPFLAGS) $(OBJS) $(OBJS2) -o viterbi $(LIB)

bgzf: bgzf.c
	$(CC) -c $(CFLAGS) bgzf.c $(LIB) -o bgzf.oc

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $*.cpp -o $*.o $(RINC)

%.oc: %.c
	$(CC) $(CFLAGS) -c $*.c -o $*.oc $(RINC)

# remove compilation products
clean:
	@rm -f viterbi *.o *.d 
