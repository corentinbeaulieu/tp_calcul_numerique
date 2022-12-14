##############################################
# archlinux.mk
# Default options for my archlinux computer
##############################################

CC=gcc
LIBSLOCAL=-L/usr/lib -llapack -lcblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native
