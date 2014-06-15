INCLUDEDIR=util/include
LIBDIR=util/lib
CC=gcc
CFLAGS=-I. -I$(INCLUDEDIR)
CLIBS=$(LIBDIR)/csupport.a $(LIBDIR)/cspice.a -L$(LIBDIR) -lgsl -lgslcblas
USER=www-data
GROUP=www-data

clean:
	@rm -rf *.o *.out *.dat *.png
	@find . -name "*~" -exec rm -rf {} \;

cleanall:cleangsl cleancspice
	@rm -rf util/lib/* util/include/* util/bin/* util/share/*

cleangsl:
	@rm -rf $(LIBDIR)/*gsl*
	@rm -rf $(INCLUDEDIR)/*gsl*

cleancspice:
	@rm -rf $(LIBDIR)/cspice.a $(LIBDIR)/csupport.a
	@rm -rf $(INCLUDEDIR)/*.h

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
