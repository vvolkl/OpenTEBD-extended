#Makefile for the Heisenberg Case Study Wrapper
#----------------------------------------------
OFLAGS = -O3 -ffree-line-length-none #Set the optimization flags
LIBFLAGS = # -llapack -lblas #Link the LAPACK and BLAS libraries
lgfortran=gfortran  ${F95ROOT}/lib/intel64/libmkl_blas95_lp64.a -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm  -I${F95ROOT}/include/intel64/lp64 -m64 -I${MKLROOT}/include
MC = $(lgfortran) #Specify the compiler

#-----------------------------------
#List objects
#-----------------------------------
OBJS = system_parameters.o TEBDtools_module.o io_module.o  Hamiltonian_tools_module.o bose_hubbard_module.o fermi_hubbard_module.o heisenberg_module.o spinS_module.o rotation_module.o local_operations_module.o observables_module.o  propagation_module.o  

#-----------------------------------
#Object compilation rules
#-----------------------------------


#Modules	
system_parameters.o : system_parameters.f90
		$(MC) $(OFLAGS) -c system_parameters.f90

TEBDtools_module.o : ../../VidalCore_v_2.0/TEBDtools_module.f90
		$(MC) $(OFLAGS) -c ../../VidalCore_v_2.0/TEBDtools_module.f90

io_module.o : ../../ToolsCore_v_2.0/io_module.f90
		$(MC) $(OFLAGS) -c ../../ToolsCore_v_2.0/io_module.f90
		
Hamiltonian_tools_module.o : ../../HamiltonianCore_v_2.0/Hamiltonian_tools_module.f90
		$(MC) $(OFLAGS) -c ../../HamiltonianCore_v_2.0/Hamiltonian_tools_module.f90

bose_hubbard_module.o : ../../HamiltonianCore_v_2.0/bose_hubbard_module.f90
		$(MC) $(OFLAGS) -c ../../HamiltonianCore_v_2.0/bose_hubbard_module.f90

fermi_hubbard_module.o : ../../HamiltonianCore_v_2.0/fermi_hubbard_module.f90
		$(MC) $(OFLAGS) -c ../../HamiltonianCore_v_2.0/fermi_hubbard_module.f90

heisenberg_module.o : ../../HamiltonianCore_v_2.0/heisenberg_module.f90
		$(MC) $(OFLAGS) -c ../../HamiltonianCore_v_2.0/heisenberg_module.f90

spinS_module.o : ../../HamiltonianCore_v_2.0/spinS_module.f90
		$(MC) $(OFLAGS) -c ../../HamiltonianCore_v_2.0/spinS_module.f90

rotation_module.o : ../../HamiltonianCore_v_2.0/rotation_module.f90
		$(MC) $(OFLAGS) -c ../../HamiltonianCore_v_2.0/rotation_module.f90

local_operations_module.o : ../../VidalCore_v_2.0/local_operations_module.f90
		$(MC) $(OFLAGS) -c ../../VidalCore_v_2.0/local_operations_module.f90

observables_module.o : ../../VidalCore_v_2.0/observables_module.f90
		$(MC) $(OFLAGS) -c ../../VidalCore_v_2.0/observables_module.f90

propagation_module.o : ../../VidalCore_v_2.0/propagation_module.f90
		$(MC) $(OFLAGS) -c ../../VidalCore_v_2.0/propagation_module.f90

#Main Programs	
XXDynamics.o : XXDynamics.f90
		$(MC) $(OFLAGS) -c XXDynamics.f90

TEBDSolve.o : TEBDSolve.f90
		$(MC) $(OFLAGS) -c TEBDSolve.f90
XXZ_ITP_block.o : XXZ_ITP_block.f90
		$(MC) $(OFLAGS) -c XXZ_ITP_block.f90

#-----------------------------------
#Primary rules
#-----------------------------------

XX :  $(OBJS) XXDynamics.o
	$(MC) $(OFLAGS)  -o Execute_XX XXDynamics.o $(OBJS) $(LIBFLAGS)

Solver :  $(OBJS) TEBDSolve.o
	$(MC) $(OFLAGS)  -o Execute_TEBDSolve TEBDSolve.o $(OBJS) $(LIBFLAGS)
XXZblock :  $(OBJS) XXZ_ITP_block.o
	$(MC) $(OFLAGS)  -o Execute_XXZ_block XXZ_ITP_block.o $(OBJS) $(LIBFLAGS)

clean: #removes all old *.o and *.mod files
	rm *.o
	rm *.mod
	rm Execute*

touch : #go through all of the *.f90 source files in the Core folders, touching them (updating the timestamp, forcing 'make' to recompile them)
	touch *.f90
	cd ../../VidalCore_v_2.0/
	touch *.f90
	cd ../../ToolsCore_v_2.0/
	touch *.f90
	cd ../../HamiltonianCore_v_2.0/
	touch *.f90


	
