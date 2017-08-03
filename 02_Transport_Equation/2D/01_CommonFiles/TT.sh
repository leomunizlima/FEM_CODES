#!/bin/bash
make clean && make fast
echo CSR DD NOT
./TranspEquation 0_Parameters/CSR/DD/NOT.dat >> CSR_DD_NOT.txt
echo CSR DD ILU
./TranspEquation 0_Parameters/CSR/DD/ILU.dat >> CSR_DD_ILU.txt
echo CSR DD DIAG
./TranspEquation 0_Parameters/CSR/DD/DIAG.dat >> CSR_DD_DIAG.txt
################################################################################
echo CSR SUPG NOT
./TranspEquation 0_Parameters/CSR/SUPG/NOT.dat >> CSR_SUPG_NOT.txt
echo CSR SUPG ILU
./TranspEquation 0_Parameters/CSR/SUPG/ILU.dat >> CSR_SUPG_ILU.txt
echo CSR SUPG DIAG
./TranspEquation 0_Parameters/CSR/SUPG/DIAG.dat >> CSR_SUPG_DIAG.txt
################################################################################
echo EBE DD NOT
./TranspEquation 0_Parameters/EBE/DD/NOT.dat >> EBE_DD_NOT.txt
echo EBE DD DIAG
./TranspEquation 0_Parameters/EBE/DD/DIAG.dat >> EBE_DD_DIAG.txt
echo EBE DD JACOBI
./TranspEquation 0_Parameters/EBE/DD/JACOBI.dat >> EBE_DD_JACOBI.txt
echo EBE DD SOR05
./TranspEquation 0_Parameters/EBE/DD/SOR05.dat >> EBE_DD_SOR05.txt
echo EBE DD SOR10
./TranspEquation 0_Parameters/EBE/DD/SOR10.dat >> EBE_DD_SOR10.txt
echo EBE DD SOR15
./TranspEquation 0_Parameters/EBE/DD/SOR15.dat >> EBE_DD_SOR15.txt
################################################################################
echo EBE SUPG NOT
./TranspEquation 0_Parameters/EBE/SUPG/NOT.dat >> EBE_SUPG_NOT.txt
echo EBE SUPG DIAG
./TranspEquation 0_Parameters/EBE/SUPG/DIAG.dat >> EBE_SUPG_DIAG.txt
echo EBE SUPG JACOBI
./TranspEquation 0_Parameters/EBE/SUPG/JACOBI.dat >> EBE_SUPG_JACOBI.txt
echo EBE SUPG SOR05
./TranspEquation 0_Parameters/EBE/SUPG/SOR05.dat >> EBE_SUPG_SOR05.txt
echo EBE SUPG SOR10
./TranspEquation 0_Parameters/EBE/SUPG/SOR10.dat >> EBE_SUPG_SOR10.txt
echo EBE SUPG SOR15
./TranspEquation 0_Parameters/EBE/SUPG/SOR15.dat >> EBE_SUPG_SOR15.txt
