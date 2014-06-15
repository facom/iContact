INCLUDEDIR=util/include
LIBDIR=util/lib
CC=gcc
CFLAGS=-I. -I$INCLUDEDIR/include
CLIBS=$LIBDIR/lib/csupport.a $LIBDIR/lib/cspice.a -L$LIBDIR -lgsl -lgslcblas
USER=root
GROUP=www-data

clean:
	@rm -rf *.o *.out *.dat *.png
	@find . -name "*~" -exec rm -rf {} \;

cleanall:clean
	@rm -rf util/lib/* util/include/* util/bin/* util/share/*

%.out:%.o
	$(CC) $^ $(CLIBS) -lm -o $(@:.o=.out)

%.o:%.c
	$(CC) $(CFLAGS) -c -g $^ -o $@ 

commit:
	@echo "Commiting..."
	@git commit -am "Update"
	@git push origin master

pull:
	@git pull reset --hard HEAD

perms:
	@chown -R root.www-data .
	chmod -R g+w .
