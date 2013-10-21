DIRS = checkfasta cut_alignment extract_seq extract_tree fasta2phylip phylip2fasta project_alignment randtree root_tree tree_dists

opt: 
	for i in $(DIRS); do cd $$i; make opt; cd ..; done;

debug: 
	for i in $(DIRS); do cd $$i; make; cd ..; done;

clean:
	for i in $(DIRS); do cd $$i; make clean; cd ..; done;