DIRS = bioc mavid utils

opt:
	for i in $(DIRS); do cd $$i; make opt; cd ..; done;

debug:
	for i in $(DIRS); do cd $$i; make; cd ..; done;

clean:
	for i in $(DIRS); do cd $$i; make clean; cd ..; done;