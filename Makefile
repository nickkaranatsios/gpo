
.PHONY: all clean test

all clean:
	$(MAKE) -C build $@

test: all
	build/test_gpo.exe
