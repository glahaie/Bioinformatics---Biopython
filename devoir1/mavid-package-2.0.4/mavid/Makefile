
.SUFFIXES += .exe .tst
CPP = g++
BASE_CPP_OPTIONS = ${IDIRS}
DEBUG_CPP_OPTIONS = -g -fno-inline -DDEBUG -Wall -Wmissing-prototypes -Wpointer-arith -Wcast-align
OPT_CPP_OPTIONS = -O4 -ffast-math -fexpensive-optimizations -finline-functions
IDIRS = -I../bioc
BASE_LINK_OPTIONS = -lstdc++ -L../bioc -lbioc -lm

objects = mavid.o output.o

# comment out the following line if you don't want to clean automatically
#clean_option = clean

# Debug version
debug:
	make mavid "CPP_OPTIONS = ${BASE_CPP_OPTIONS} ${DEBUG_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Optimized version
opt:
	make mavid "CPP_OPTIONS = ${OPT_CPP_OPTIONS} ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Profiling version
prof:
	make mavid "CPP_OPTIONS = -g -pg ${OPT_CPP_OPTIONS} ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Memory checking version
mem:
	make mavid "CPP_OPTIONS = -g ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = -g -lgc_dbg ${BASE_LINK_OPTIONS}"


# Optimized version
mavid: $(clean_option) $(objects)
	$(CPP) $(IDIRS) $(CPP_OPTIONS) -o mavid $(objects) $(LINK_OPTIONS)

mavid.o:         mavid.cc output.cc output.h ../bioc/libbioc.a
output.o:         output.cc output.h


%.o: %.cc
	$(CPP) $(IDIRS) $(CPP_OPTIONS) -c $< 

.PHONY: clean
clean:
	rm -f *.o mavid
 
# DO NOT DELETE THIS LINE -- make depend depends on it.



