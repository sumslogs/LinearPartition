################################
# Makefile
#
# author: He Zhang
# edited by: 03/2019
################################

# NOTE:
# -ggdb for debugging
# -O3 to not have names mangled during debugging
CC=g++
DEPS=src/bpp.cpp src/LinearPartition.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h
CFLAGS=-std=c++11 -ggdb
# -O3
.PHONY : all clean linearpartition
objects=bin/linearpartition_v bin/linearpartition_c

all: linearpartition

linearpartition: src/LinearPartition.cpp $(DEPS) 
		chmod +x linearpartition draw_bpp_plot draw_heatmap
		mkdir -p bin
		$(CC) src/LinearPartition.cpp $(CFLAGS) -Dlpv -o bin/linearpartition_v 
		$(CC) src/LinearPartition.cpp $(CFLAGS) -o bin/linearpartition_c

clean:
	-rm $(objects)