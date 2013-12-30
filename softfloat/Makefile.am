## Process this file with automake to produce Makefile.in
#
# Makefile for building the SoftFloat package
# for use with Hercules S/370, ESA/390 and z/Architecture emulator
#
# This file was added by the Hercules project.
# It is not part of the original SoftFloat distribution.
#

lns=@LN_S@

LDADD = @LIBS@ 

AM_CPPFLAGS = -I$(top_srcdir)

softfloat_SRC = softfloat.c

if BUILD_SHARED
   XSTATIC =
else
   XSTATIC = -static
endif

if OPTION_DYNAMIC_LOAD
  LTDL  	 = ../ltdl.c
  LIB_LD_FLAGS    = -export-dynamic  \
		   $(XSTATIC)       \
		   -no-undefined    \
		   -avoid-version
else
  LTDL  	 =
  LIB_LD_FLAGS    = $(XSTATIC)    \
		   -no-undefined \
		   -avoid-version
endif

  HERCLIBS =
  HERCLIBS2 = libsoftfloat.la
  noinst_LTLIBRARIES = $(HERCLIBS)
  lib_LTLIBRARIES = $(HERCLIBS2)

  libsoftfloat_la_SOURCES = $(softfloat_SRC)
  libsoftfloat_la_LDFLAGS = $(LIB_LD_FLAGS)
  libsoftfloat_la_LIBADD  = $(LDADD)

noinst_HEADERS = milieu.h          \
		 processor.h       \
		 softfloat.h       \
		 softfloat-macros  \
		 softfloat-specialize

EXTRA_DIST = README.txt            \
	     SoftFloat.txt         \
	     SoftFloat-history.txt \
	     SoftFloat-source.txt

%.s: %.c
	$(COMPILE) -S $<
