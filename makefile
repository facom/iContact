CC=gcc
CFLAGS=-I. -Iinclude
CSPICE_LIBS=lib/csupport.a lib/cspice.a -lgsl -lgslcblas

clean:
	@rm -rf *.o *.out
	@find . -name "*~" -exec rm -rf {} \;

%.out:%.o
	$(CC) $^ $(CSPICE_LIBS) -lm -o $(@:.o=.out)

%.o:%.c
	$(CC) $(CFLAGS) -c -g $^ -o $@ 

commit:
	@echo "Commiting..."
	@git commit -am "Update"
	@git push origin master

pull:
	@git pull