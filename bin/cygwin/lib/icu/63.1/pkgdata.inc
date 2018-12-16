GENCCODE_ASSEMBLY_TYPE=-a gcc-mingw64
SO=dll
SOBJ=dll
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=63.dll
COMPILE=gcc -D_REENTRANT  -DU_HAVE_ELF_H=1 -DU_HAVE_STRTOD_L=1 -DU_HAVE_XLOCALE_H=1  -DU_ATTRIBUTE_DEPRECATED= -ggdb -O2 -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fstack-protector-strong --param=ssp-buffer-size=4 -fdebug-prefix-map=/home/kbrown/src/cygicu/icu-63.1-1.x86_64/build=/usr/src/debug/icu-63.1-1 -fdebug-prefix-map=/home/kbrown/src/cygicu/icu-63.1-1.x86_64/src/icu/source=/usr/src/debug/icu-63.1-1 -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -c
LIBFLAGS=-I/usr/include -DPIC 
GENLIB=gcc -ggdb -O2 -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fstack-protector-strong --param=ssp-buffer-size=4 -fdebug-prefix-map=/home/kbrown/src/cygicu/icu-63.1-1.x86_64/build=/usr/src/debug/icu-63.1-1 -fdebug-prefix-map=/home/kbrown/src/cygicu/icu-63.1-1.x86_64/src/icu/source=/usr/src/debug/icu-63.1-1 -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings    -shared -Wl,-Bsymbolic -Wl,--enable-auto-import -Wl,--out-implib=all.a
LDICUDTFLAGS=
LD_SONAME=
RPATH_FLAGS=
BIR_LDFLAGS=-Wl,-Bsymbolic
AR=ar
ARFLAGS=r
RANLIB=ranlib
INSTALL_CMD=/usr/bin/install -c
