MAIN = main.f90
MDIR = ./modules/
MODULES_names := physicals.f90 numericals.f90 lele_filter.f90 derivative.f90 IO.f90
MODULES := $(addprefix $(MDIR), $(MODULES_names))
COMPILE_BASE := gfortran -J $(MDIR) $(MODULES)

help:
	@echo "Available tasks:"
	@echo "compile  -> Compile fortran main and modules"
	@echo "run      -> Run the compiled codes"
	@echo "all      -> Call compile, run and plot results"
	@echo "clean    -> Remove *.mod *.so *.c *.o *.html build"

all: compile run plot

compile:
	$(COMPILE_BASE) $(MAIN) -o run

.PHONY: run
run:
	./run

plot:
	python plots/main.py

# Phony targets for cleanup and similar uses
.PHONY: clean
clean:
	rm -rf *~ *.mod $(MDIR)/*.mod *.out *.o *.html build __pycache__

clean-data:
	rm -rf data/*

test-penta: 
	$(COMPILE_BASE) tests/test_penta.f90 -llapack -lblas -o tests/test_penta.out
	clear
	./tests/test_penta.out
	
