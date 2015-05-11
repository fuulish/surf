default:
	@echo 'Choose one of the following targets'
	@echo '* xdrfile downloads and compiles xdrfile 1.1.4'
	@echo '* einspline downloads and compiles EINSPLINE 0.9.2'

xdrfile: xdrfile-1.1.4/src/.libs/libxdrfile.a

einspline: einspline-0.9.2/src/.libs/libeinspline.a

all: xdrfile einspline

# Build xdrfile from scratch

einspline-0.9.2/src/.libs/libeinspline.a: einspline-0.9.2
	cd einspline-0.9.2 && CONFIG_SHELL=/bin/sh ./configure
	cd einspline-0.9.2 && $(MAKE)

einspline-0.9.2: einspline-0.9.2.tar.gz
	tar -xzvf einspline-0.9.2.tar.gz

einspline-0.9.2.tar.gz:
	wget 'http://sourceforge.net/projects/einspline/files/einspline-0.9.2.tar.gz/download' -O einspline-0.9.2.tar.gz

# Build libxc from scratch

xdrfile-1.1.4.tar.gz:
	wget 'ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz' -O xdrfile-1.1.4.tar.gz

xdrfile-1.1.4: xdrfile-1.1.4.tar.gz
	tar -xzvf xdrfile-1.1.4.tar.gz

xdrfile-1.1.4/src/.libs/libxdrfile.a: xdrfile-1.1.4
	cd xdrfile-1.1.4 && CONFIG_SHELL=/bin/sh ./configure
	cd xdrfile-1.1.4 && $(MAKE)


