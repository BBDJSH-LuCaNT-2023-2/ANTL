AUTOMAKE_OPTIONS = subdir-objects

AM_LIBS = @PTHREAD_LIBS@ @LIBS@
AM_CXXFLAGS = @CXXFLAGS@ @PTHREAD_CFLAGS@
AM_CC = @PTHREAD_CC@
AM_CPPFLAGS = @CPPFLAGS@ -I$(top_srcdir)/include

bin_PROGRAMS = test main

test_SOURCES= tests/IndexCalculus/IndCalc_test.cpp 
main_SOURCES= tests/HeaderTest.cpp
