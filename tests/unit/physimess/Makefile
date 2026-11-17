PROGRAM_NAME := physimess_unit_tests

CC := g++

# Check for environment definitions of compiler
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

# The tests only need the header PhysiMeSS_geometry.h: no PhysiCell object,
# no OpenMP. -ffp-contract=off keeps the results independent of the machine.
CFLAGS := -O2 -ffp-contract=off -std=c++11 -Wall -Wextra

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): test_geometry.cpp ../../../addons/PhysiMeSS/PhysiMeSS_geometry.h
	$(CC) $(CFLAGS) -o $(PROGRAM_NAME) test_geometry.cpp

run: $(PROGRAM_NAME)
	./$(PROGRAM_NAME)

clean:
	rm -f $(PROGRAM_NAME) $(PROGRAM_NAME).exe

.PHONY: all run clean
