INCLUDEDIR=util/include
LIBDIR=util/lib
CC=gcc
CFLAGS=-I. -I$(INCLUDEDIR)
CLIBS=$(LIBDIR)/csupport.a $(LIBDIR)/cspice.a -L$(LIBDIR) -lgsl -lgslcblas
USER=www-data
GROUP=www-data

clean:
	@echo "Basic cleaning..."
	@rm -rf *.o *.out *.dat *.png
	@find . -name "*~" -exec rm -rf {} \;

cleandata:clean
	@echo "Cleaning data..."
	@rm -rf data/*

cleangsl:
	@echo "Cleaning gsl..."
	@rm -rf $(LIBDIR)/*gsl*
	@rm -rf $(INCLUDEDIR)/*gsl*

cleancspice:
	@echo "Cleaning cspice..."
	@rm -rf $(LIBDIR)/cspice.a $(LIBDIR)/csupport.a
	@rm -rf $(INCLUDEDIR)/*.h

cleanall:clean cleandata cleangsl cleancspice
	@echo "Cleaning all..."
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
	@git reset --hard HEAD
	@git pull 

perms:
	@chown -R $(USER).$(GROUP) .
	@chmod -R g+w .
