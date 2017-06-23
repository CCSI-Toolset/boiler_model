#!/bin/sh
#please select the compilers
#using g++ and gfortran
CXX=g++
F77=gfortran
CFLAGS=-Wall -I.
#using intel's icc and ifort
#CXX=g++
#F77=gfortran
#CFLAGS=-I.
#end of compiler selection

OBJS=main.o \
	BoundaryFaceProperty.o \
	CoalKinetics.o \
	DiscreteOrdinateEquation.o \
	HexCell.o \
	IdealGas.o \
	IdealGasMixture.o \
	LiquidFuel.o \
	MaterialStream.o \
	Mesh.o \
	ParticleRadiationProperty.o \
	PSD.o \
	QuadFace.o \
	RadiantFurnace.o \
	SolidFuel.o \
	SolveBoilerModel.o \
	SuperHeater.o \
	UtilityFunctions.o \
	Vertex.o \
	GasEmissivity.o

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

.f.o:
	$(F77) -c $<

boilermodel: $(OBJS)
	$(CXX) $(CFLAGS) -o boilermodel $(OBJS)

$(OBJS): Constants.h \
	BoundaryFaceProperty.h \
	CoalKinetics.h \
	DiscreteOrdinateEquation.h \
	HexCell.h \
	IdealGas.h \
	IdealGasMixture.h \
	LiquidFuel.h \
	MaterialStream.h \
	Mesh.h \
	ParticleRadiationProperty.h \
	PSD.h \
	QuadFace.h \
	RadiantFurnace.h \
	SolidFuel.h \
	SolveBoilerModel.h \
	SuperHeater.h \
	UtilityFunctions.h \
	Vertex.h

clean:
	@rm -rf *.o boilermodel
