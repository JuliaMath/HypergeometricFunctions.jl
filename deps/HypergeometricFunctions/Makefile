include Make.inc

all:
	make lib

lib:
	$(CC) $(CFLAGS) $(LIBFLAGS) $(SRC) $(LDFLAGS) $(LDLIBS) -o lib$(LIB).$(SLIB)

coverage:
	$(COV) $(SRC) -o $(LIBDIR)

clean:
	rm -f lib$(LIB).*
	rm -f test_*
	find . -type f -name '*.s' -delete
	find . -type f -name '*.gcda' -delete
	find . -type f -name '*.gcno' -delete
	find . -type f -name '*.gcov' -delete

.PHONY: all lib tests runtests coverage clean
