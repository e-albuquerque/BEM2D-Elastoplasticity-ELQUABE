# ==============================================================================
# Makefile for Modernized ELQUABE (2D Elastoplastic BEM)
# ==============================================================================

FC      = gfortran
FFLAGS  = -fopenmp -O3 -fdefault-real-8 -Wall -Wextra
LIBS    = -llapack -lblas
TARGET  = elquabe_modern.exe

ifeq ($(OS),Windows_NT)
    LDFLAGS = -static
else
    LDFLAGS =
endif

OBJS    = data_mod.o io_mod.o solver_mod.o assembly_mod.o main.o

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

data_mod.o: data_mod.f90
	$(FC) $(FFLAGS) -c data_mod.f90

io_mod.o: io_mod.f90 data_mod.o
	$(FC) $(FFLAGS) -c io_mod.f90

solver_mod.o: solver_mod.f90
	$(FC) $(FFLAGS) -c solver_mod.f90

assembly_mod.o: assembly_mod.f90 data_mod.o solver_mod.o
	$(FC) $(FFLAGS) -c assembly_mod.f90

main.o: main.f90 data_mod.o io_mod.o solver_mod.o assembly_mod.o
	$(FC) $(FFLAGS) -c main.f90

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o *.mod $(TARGET)
