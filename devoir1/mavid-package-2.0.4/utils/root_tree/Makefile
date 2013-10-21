
.SUFFIXES += .exe .tst
CPP = g++
PROJNAME = root_tree
IDIRS = -I../../bioc
BASE_CPP_OPTIONS =
DEBUG_CPP_OPTIONS = -Wall -Wmissing-prototypes -Wtraditional -Wpointer-arith -Wcast-align
BASE_LINK_OPTIONS = -lstdc++ -lbioc -L../../bioc


objects = root_tree.o

# comment out the following line if you don't want to have it clean
# automatically
#clean_option = clean

# Debug version
debug:
	make main "CPP_OPTIONS = -g -DDEBUG ${DEBUG_CPP_OPTIONS} ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Optimized version
opt:
	make main "CPP_OPTIONS = -O4 -ffast-math -fexpensive-optimizations -finline-functions ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Profiling version
prof:
	make main "CPP_OPTIONS = -O4 -g -ffast-math -fexpensive-optimizations -finline-functions -pg ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Memory checking version
mem:
	make main "CPP_OPTIONS = -g ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = -g -lgc_dbg ${BASE_LINK_OPTIONS}"


# actually compile it
main: $(clean_option) $(objects)
	c++ $(IDIRS) -o $(PROJNAME) $(objects) $(LINK_OPTIONS)

%.o: %.cc
	$(CPP) $(IDIRS) $(CPP_OPTIONS) -c $< 

.PHONY: clean
clean:
	rm -f *.o ${PROJNAME}
 
# DO NOT DELETE THIS LINE -- make depend depends on it.



