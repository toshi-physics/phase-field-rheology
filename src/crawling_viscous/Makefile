
# CC := gcc
# CC := gcc-10
CC := $(CC)
CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11 -fopenmp -flto -g
#CFLAGS := -Wall -Wshadow -Winline -O3 -std=c11

# CXX := g++
# CXX := g++-10
CXX := $(CXX)
CXXFLAGS := -Wall -Wshadow -Winline -O3 -std=c++11 -flto -g

SHELL := /bin/bash
BIN := ../../bin/crawling_viscous/
LDIR := ../../lib/
CONDA_DIR := /home/s1309877/.local/bin/miniconda3/envs/myUtil/
INC :=  -I$(CONDA_DIR)/include/ -I$(LDIR)/
LFLAGS := -L$(CONDA_DIR)/lib/ -lm -lgomp -flto -O3 -llapack -lblas -g

PRG_SUFFIX_FLAG := 0
PRG_SUFFIX =

# 0 = no armadillo and 1 = with armadillo
HAS_ARMA := 1

ODIR := $(BIN)/obj/
EDIR := $(BIN)/exe/
LODIR := $(BIN)/lib/

SRC_C := cell.c model.c array.c image.c overlap.c dump.c dump_cm.c dump_bulk_cm.c dump_field.c dump_cell_field.c dump_gyration.c dump_gyration_field.c dump_velocity.c dump_velocity_field.c dump_deform.c dump_deform_field.c dump_neighbour.c dump_energy.c dump_overlap.c dump_overlap_field.c dump_index_field.c random_mtwister.c util.c linear_lapack.c
SRC_CXX := neighbour.cpp
ifeq ($(HAS_ARMA),1)
	SRC_C := $(SRC_C) shape.c dump_shape.c
	SRC_CXX := $(SRC_CXX) sgolay_filter.cpp
	LFLAGS := $(LFLAGS) -larmadillo
endif
OBJ_C := $(patsubst %.c, $(ODIR)/%.o, $(SRC_C))
OBJ_CXX := $(patsubst %.cpp, $(ODIR)/%.o, $(SRC_CXX))
PRG_SRC := run_model.c
PRG_OBJ := $(patsubst %.c, $(ODIR)/%.o, $(PRG_SRC))
LSRC := arralloc.c mtwister.c
LOBJ := $(patsubst %.c, $(LODIR)/%.o, $(LSRC))
PRG := $(patsubst %.c, $(EDIR)/%$(PRG_SUFFIX), $(PRG_SRC))

all : check_arma make_bin $(LOBJ) $(OBJ) $(PRG_OBJ) $(PRG)

%.o : %.c %.h
	$(CC) -o $@ -c $< $(CFLAGS) $(INC)

$(LODIR)/%.o : $(LDIR)/%.c
	$(CC) -o $@ -c $< $(CFLAGS) $(INC)

$(ODIR)/%.o : %.c %.h
	$(CC) -o $@ -c $< $(CFLAGS) $(INC) 

$(ODIR)/%.o : %.c
	$(CC) -o $@ -c $< $(CFLAGS) $(INC) 

$(ODIR)/%.o : %.cpp %.hpp
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(INC)

$(ODIR)/%.o : %.cpp %.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(INC)

$(ODIR)/%.o : %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(INC)

$(PRG) : $(LOBJ) $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ)
	$(CXX) $^ -o $@ $(LFLAGS)

.PHONY : clean check_arma make_bin

clean :
	rm -rf $(OBJ_C) $(OBJ_CXX) $(PRG_OBJ) $(PRG) $(LOBJ) $(LODIR) $(ODIR) $(EDIR)

check_arma :
ifeq ($(HAS_ARMA),1)
	if [ $$(grep 'PF_HAS_ARMA 1' dump.h | wc -l ) != 1 ]; \
	then perl -pi -e 's/PF_HAS_ARMA [0,1]/PF_HAS_ARMA 1/g' dump.h; \
	fi;
else
	if [ $$(grep 'PF_HAS_ARMA 0' dump.h | wc -l ) != 1 ]; \
	then perl -pi -e 's/PF_HAS_ARMA [0,1]/PF_HAS_ARMA 0/g' dump.h; \
	fi;
endif

make_bin : $(ODIR) $(EDIR) $(LODIR)

$(ODIR) : 
	mkdir -p $(ODIR)

$(EDIR) :
	mkdir -p $(EDIR)

$(LODIR) :
	mkdir -p $(LODIR)
