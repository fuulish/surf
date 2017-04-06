default:
	@echo 'Choose one of the following targets'
	@echo '* xdrfile downloads and compiles xdrfile 1.1.4'
	@echo '* einspline downloads and compiles EINSPLINE 0.9.2'
	@echo '* nlopt downloads and compiles nlopt 2.4.2'

.PHONY: nlopt xdrfile einspline

nlopt: lib/libnlopt.a

gsl: lib/libgsl.a

xdrfile: lib/libxdrfile.a

einspline: lib/libeinspline.a

all: xdrfile einspline nlopt

# Build einspline from scratch

lib/libeinspline.a: einspline-0.9.2
	cd einspline-0.9.2 && CONFIG_SHELL=/bin/bash ./configure --prefix=$(PWD)
	cd einspline-0.9.2 && $(MAKE) install

einspline-0.9.2: einspline-0.9.2.tar.gz
	tar -xzvf einspline-0.9.2.tar.gz

einspline-0.9.2.tar.gz:
	wget 'http://sourceforge.net/projects/einspline/files/einspline-0.9.2.tar.gz/download' -O einspline-0.9.2.tar.gz

# Build xdrfile from scratch

lib/libxdrfile.a: xdrfile-1.1.4
	cd xdrfile-1.1.4 && CONFIG_SHELL=/bin/bash ./configure --prefix=$(PWD) --enable-shared
	cd xdrfile-1.1.4 && $(MAKE) install

xdrfile-1.1.4: xdrfile-1.1.4.tar.gz
	tar -xzvf xdrfile-1.1.4.tar.gz

xdrfile-1.1.4.tar.gz:
	wget 'ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz' -O xdrfile-1.1.4.tar.gz

lib/libnlopt.a: nlopt-2.4.2
	cd nlopt-2.4.2&& CONFIG_SHELL=/bin/bash ./configure --prefix=$(PWD)
	cd nlopt-2.4.2&& $(MAKE) install

nlopt-2.4.2: nlopt-2.4.2.tar.gz
	tar -xzvf nlopt-2.4.2.tar.gz

nlopt-2.4.2.tar.gz:
	wget 'http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz' -O nlopt-2.4.2.tar.gz

lib/libgsl.a: gsl-2.3
	cd gsl-2.3 && CONFIG_SHELL=/bin/bash ./configure --prefix=$(PWD)
	cd gsl-2.3 && $(MAKE) install

gsl-2.3: gsl-2.3.tar.gz
	tar -xzvf gsl-2.3.tar.gz

gsl-2.3.tar.gz:
	wget 'ftp://ftp.gnu.org/gnu/gsl/gsl-2.3.tar.gz' -O gsl-2.3.tar.gz

clean:
	-rm -rf einspline-0.9.2* xdrfile-1.1.4* nlopt-2.4.2* bin include lib
