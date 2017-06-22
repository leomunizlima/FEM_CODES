CC=gcc
#CC=tau_cc.sh
IDIR= 
CFLAGS= -DSSTranspEquation2D -Ddebug -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DSSTranspEquation2D -Ddebug -lm -lrt -Wall ../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-lrt -lm -lgfortran -Wall -lblas -llapack -lquadmath -Ofast -march=native
LDFLAGS=
SOURCES= main.c \
Preprocess.c \
Process.c \
Postprocess.c \
Fill_LM.c \
ede_Initialization.c \
csr_Initialization.c \
eval_U.c \
F_assembly.c \
Build_K_F_SUPG.c \
Build_K_F_DD.c \
ebe_assembly.c \
ede_assembly.c \
csr_assembly.c \
Paraview_Output.c \
CAU_ShockCapture.c \
CAU_DD_ShockCapture.c \
YZBeta_ShockCapture.c \
YZBeta_ShockCapture2.c \
h_shock_2sqrtArea.c \
h_shock_Option1.c \
h_shock_Option2.c \
set_h_Shock.c \
setStabilizationForm.c \
setSolver.c \
setPreconditioner.c \
setProblem.c \
setMatrixVectorProductType.c \
setzeros.c \
../CARTOLA/CARTOLA_Condutivity.c \
../CARTOLA/CARTOLA_Font.c \
../CARTOLA/CARTOLA_Reaction.c \
../CARTOLA/CARTOLA_Velocity.c \
../CARTOLA/CARTOLA_upresc.c \
../HEMKER/HEMKER_Condutivity.c \
../HEMKER/HEMKER_Font.c \
../HEMKER/HEMKER_Reaction.c \
../HEMKER/HEMKER_Velocity.c \
../HEMKER/HEMKER_upresc.c \
../PUDIM/PUDIM_Condutivity.c \
../PUDIM/PUDIM_Font.c \
../PUDIM/PUDIM_Reaction.c \
../PUDIM/PUDIM_Velocity.c \
../PUDIM/PUDIM_upresc.c \
../TESTE/TESTE_Condutivity.c \
../TESTE/TESTE_Font.c \
../TESTE/TESTE_Reaction.c \
../TESTE/TESTE_Velocity.c \
../TESTE/TESTE_upresc.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemv.c \
../../../00_CommonFiles/MatrixVector_Operations/EDE/edemv.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pcg.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond_CSR_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/Jacobi_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/Jacobi_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGS_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SGS_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSOR_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SSOR_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/Jacobi_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/Jacobi_precond_EDE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGS_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SGS_precond_EDE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSOR_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SSOR_precond_EDE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/ILUp_precond.c \
../../../00_CommonFiles/Reordering/reordering.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=SSTranspEquation 

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DSSTranspEquation2D -Ddebug -c -lm -lrt -Wall -g
debug: GCFLAGS= -DSSTransEquation2D -Ddebug -lm -lrt -Wall \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-lgfortran -lblas -llapack -lquadmath -g
debug: all

fast: CFLAGS= -DSSTranspEquation2D -c -lm -lrt -Wall -Ofast -march=native 
fast: GCFLAGS= -DSSTransEquation2D -lm -lrt -Wall -Ofast -march=native \
../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o \
../../../00_CommonFiles/Reordering/hsl_mc73/libhsl_mc73.a ../../../00_CommonFiles/Reordering/hsl_mc60/libmc60.a \
-lgfortran -lblas -llapack -lquadmath 
fast: all

gprof: CFLAGS+=-pg -Ddebug 
gprof: LDFLAGS+=-pg -Ddebug
gprof: all
 
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)



