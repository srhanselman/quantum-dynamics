FC = gfortran


qd: sym_tree_ify.o linalg_algorithms.o qd_procedures.o main.o
	$(FC) sym_tree_ify.o linalg_algorithms.o qd_procedures.o main.o -o qd


main.o:
	$(FC) -c src/main.f95

qd_procedures.o:
	$(FC) -c src/qd_procedures.f95

sym_tree_ify.o:
	$(FC) -c src/sym_tree_ify.f90

linalg_algorithms.o:
	$(FC) -c src/linalg_algorithms.f90



debug:
	$(FC) -c -g src/sym_tree_ify.f90
	$(FC) -c -g src/linalg_algorithms.f90
	$(FC) -c -g src/qd_procedures.f95
	$(FC) -c -g src/main.f95
	$(FC) sym_tree_ify.o linalg_algorithms.o qd_procedures.o main.o -o qd-debug
       
