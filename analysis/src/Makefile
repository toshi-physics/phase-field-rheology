
#CC := gcc
CC := $(CC)
CFLAGS := -O3 -Wall -Wshadow -std=c11 -fopenmp

#CXX := g++
CXX := $(CXX)
CXXFLAGS := -O3 -Wall -Wshadow -std=c++11 -fopenmp

PRG_SUFFIX_FLAG := 0

LDIR := ../lib
LODIR = ../bin/lib
ODIR := ../bin/obj
EDIR := ../bin/exe

CONDA_DIR := /home/s1309877/.local/bin/miniconda3/envs/myUtil/
INC = -I$(CONDA_DIR)/include/ -I$(LDIR)
LDFLAGS := -L$(CONDA_DIR)/lib/ -lm -lgomp -larmadillo

LSRCC := $(wildcard $(LDIR)/*.c)
LOBJC := $(patsubst %.c, $(LODIR)/%.o, $(notdir $(LSRCC)))
LSRCCXX := $(wildcard $(LDIR)/*.cpp)
LOBJCXX := $(patsubst %.cpp, $(LODIR)/%.o, $(notdir $(LSRCCXX)))
SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp, $(ODIR)/%.o, $(SRC))
PRG := $(patsubst %.cpp,$(EDIR)/%,$(SRC))
PRG_SUFFIX=
BIN := $(patsubst %,%$(PRG_SUFFIX),$(PRG))

all : make_bin $(LOBJC) $(LOBJCXX) $(OBJ) $(BIN)

#% : %.cpp
#	make $(EDIR)/$@

$(ODIR)/%.o : %.cpp 
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

$(EDIR)/% : $(ODIR)/%.o $(LOBJC) $(LOBJCXX)
	$(CXX) $(LDFLAGS) $^ -o $@

$(LODIR)/%.o : $(LDIR)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(LODIR)/%.o : $(LDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY : clean make_bin
clean : 
	rm -rf $(EDIR) $(ODIR) $(LODIR)

make_bin : $(LODIR) $(ODIR) $(EDIR)

$(LODIR) : 
	mkdir -p $(LODIR)

$(ODIR) :
	mkdir -p $(ODIR)

$(EDIR) :
	mkdir -p $(EDIR)
