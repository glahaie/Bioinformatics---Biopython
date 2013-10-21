
.SUFFIXES += .exe .tst
CPP = g++
BASE_CPP_OPTIONS = 
DEBUG_CPP_OPTIONS = -g -fno-inline -DDEBUG -Wall -Wmissing-prototypes -Wpointer-arith -Wcast-align
OPT_CPP_OPTIONS = -O4 -ffast-math -fexpensive-optimizations -finline-functions
BASE_LINK_OPTIONS = -lstdc++ -lm

objects = fasta.o file.o alignments.o seq.o tree.o matrices.o axt.o common.o match.o align.o sufftree.o wobble.o anchor.o malign.o constraints.o refine.o prune.o

# comment out the following line if you don't want to have it clean
# automatically
#clean_option = clean

# Debug version
debug:
	make bioc "CPP_OPTIONS = ${BASE_CPP_OPTIONS} ${DEBUG_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Optimized version
opt:
	make bioc "CPP_OPTIONS = ${BASE_CPP_OPTIONS} ${OPT_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Profiling version
prof:
	make bioc "CPP_OPTIONS = -g -pg ${BASE_CPP_OPTIONS} ${OPT_CPP_OPTIONS}" "LINK_OPTIONS = ${BASE_LINK_OPTIONS}"

# Memory checking version
mem:
	make bioc "CPP_OPTIONS = -g ${BASE_CPP_OPTIONS}" "LINK_OPTIONS = -g -lgc_dbg ${BASE_LINK_OPTIONS}"


# actually compile it
bioc: $(clean_option) $(objects)
	ar r libbioc.a $(objects)

fasta.o: 	fasta.cc fasta.h
file.o: 	file.cc file.h
alignments.o: 	alignments.cc alignments.h tree.h
seq.o:		seq.cc seq.h
tree.o:		tree.cc tree.h
matrices.o:	matrices.cc matrices.h
axt.o:		axt.cc axt.h
common.o:	common.cc common.h
match.o:	match.cc match.h sufftree.cc sufftree.h
align.o:	align.cc align.h match.cc match.h
sufftree.o:	sufftree.cc sufftree.h
wobble.o:	wobble.cc wobble.h
anchor.o:	anchor.cc anchor.h align.cc
malign.o:	malign.cc malign.h
constraints.o:	constraints.cc constraints.h
refine.o:	refine.h malign.cc tree.cc

%.o: %.cc
	$(CPP) $(IDIRS) $(CPP_OPTIONS) -c $<

.PHONY: clean
clean:
	rm -f *.o libbioc.a
 
# DO NOT DELETE THIS LINE -- make depend depends on it.



