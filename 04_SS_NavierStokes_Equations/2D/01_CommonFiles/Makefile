CC=gcc
IDIR= 
CFLAGS= -DSSNavierStokesEquations2D -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DSSNavierStokesEquations2D -lm -lrt -Wall ../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc64/libmc64.a \
-lrt -lm -lgfortran -Wall -lblas -llapack -lquadmath -Ofast -march=native
LDFLAGS= 
SOURCES= main.c \
Preprocess.c \
Process.c \
Fill_ID.c \
Fill_LM.c \
setProblem.c \
setSolver.c \
setStabilizationForm.c \
setPreconditioner.c \
setMatrixVectorProductType.c \
csr_Initialization.c \
ebe_assembly.c \
csr_assembly.c \
Build_K_F_SUPG_PSPG.c \
eta_newton.c\
Postprocess.c \
Paraview_Output.c \
../CAVITY/CAVITY_v1presc.c \
../CAVITY/CAVITY_v2presc.c \
../CAVITY/CAVITY_ppresc.c \
../CHANNEL/CHANNEL_v1presc.c \
../CHANNEL/CHANNEL_v2presc.c \
../CHANNEL/CHANNEL_ppresc.c \
../EXATA/EXATA_v1presc.c \
../EXATA/EXATA_v2presc.c \
../EXATA/EXATA_ppresc.c \
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/BLAS_Operations/izero.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/BlockDiagDOF3_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/BlockDiagDOF3_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemvNDOF3.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/Reordering/reordering.c
#../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemvMK.c \
#../../../00_CommonFiles/MatrixVector_Operations/EDE/edemv.c \
#.../../../00_CommonFiles/MatrixVector_Operations/EDE/edemvMK.c \
#.../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
#../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmvMK.c \
#ede_assembly.c \
#ede_end_assembly.c \
#csr_assembly.c \
#csr_end_assembly.c \
#Postprocess.c \
#Predictor_New.c \
#Predictor_Old.c \
#setzeros.c \
#Build_M_K_F_SUPG.c \
#Build_M_K_F_DD_Transiente.c \
#Build_M_K_F_DD_Estatico.c \
#Build_M_F_DD_Transiente.c \
#Delta_CAU.c \
#Delta_DD.c \
#Delta_YZBeta.c \
#Paraview_Output_3D.c \
#Paraview_Output_DeltaT.c \
#Paraview_Output_3D_DeltaT.c \
#calculate_DaB.c \
#../OBLIQUO/OBLIQUO_InitialSolution.c \
#../OBLIQUO/OBLIQUO_gamma.c \
#../OBLIQUO/OBLIQUO_cv.c \
#../OBLIQUO/OBLIQUO_rhopresc.c \
#../OBLIQUO/OBLIQUO_v1presc.c \
#../OBLIQUO/OBLIQUO_v2presc.c \
#../OBLIQUO/OBLIQUO_epresc.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=SSNavierStokesEquations2D 

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DSSNavierStokesEquations2D -g -Ddebug -c -lm -lrt -Wall
debug: GCFLAGS= -DSSNavierStokesEquations2D -g -Ddebug -lm -lrt -Wall ../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-lrt -lm -lgfortran -Wall -lblas -llapack -lquadmath 
debug: all

fast: CFLAGS= -DSSNavierStokesEquations2D -c -lm -lrt -Wall -Ofast -march=native
fast: GCFLAGS= -DSSNavierStokesEquations2D -lm -lrt -Wall ../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-lrt -lm -lgfortran -Wall -lblas -llapack -lquadmath -Ofast -march=native
fast: all

gprof: CFLAGS+=-pg  
gprof: LDFLAGS+=-pg
gprof: all
 
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)



