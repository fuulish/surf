default:
	@echo 'Choose one of the following targets'
	@echo '* xdrfile downloads and compiles xdrfile 1.1.4'
	@echo '* einspline downloads and compiles EINSPLINE 0.9.2'

xdrfile: lib/libxdrfile.a

einspline: lib/libeinspline.a

all: xdrfile einspline

# Build einspline from scratch

lib/libeinspline.a: einspline-0.9.2
	cd einspline-0.9.2 && CONFIG_SHELL=/bin/sh ./configure --prefix=$(PWD)
	cd einspline-0.9.2 && $(MAKE) install

einspline-0.9.2: einspline-0.9.2.tar.gz
	tar -xzvf einspline-0.9.2.tar.gz

einspline-0.9.2.tar.gz:
	wget 'http://sourceforge.net/projects/einspline/files/einspline-0.9.2.tar.gz/download' -O einspline-0.9.2.tar.gz

# Build xdrfile from scratch

lib/libxdrfile.a: xdrfile-1.1.4
	cd xdrfile-1.1.4 && CONFIG_SHELL=/bin/sh ./configure --prefix=$(PWD)
	cd xdrfile-1.1.4 && $(MAKE) install

xdrfile-1.1.4: xdrfile-1.1.4.tar.gz
	tar -xzvf xdrfile-1.1.4.tar.gz

xdrfile-1.1.4.tar.gz:
	wget 'ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz' -O xdrfile-1.1.4.tar.gz

clean:
	rm -rf einspline-0.9.2* xdrfile-1.1.4* bin include lib
