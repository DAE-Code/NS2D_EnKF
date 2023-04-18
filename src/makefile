FC        = gfortran
RM        = rm -rf
LINKER    = $(FC)
LIBS      = 

ifeq ($(MAKECMDGOALS),debug) 
  COPT    = -c -O0 -cpp -fno-range-check -g -Wall -fbounds-check 
  CFLAG   =
  LOPT    =
else
  COPT    = -c -O3 -cpp -fno-range-check 
  CFLAG   =
  LOPT    = 
endif

EXECDIR   = ./

VPATH     = src

APNAME    = enkf 


TARGET    = mod_variables.o m_random3.o m_ranmean3.o m_ranvar3.o    \
            sub_kfilter.o sub_measure.o sub_utils.o sub_bc_outer.o  \
            sub_bc_wall.o sub_hsmac.o sub_initial.o sub_plot3d.o    \
            sub_rhs3rd.o NS2D_EnKF.o svd.o pythag.o


CLNTARGET = *.o *.mod *.d *.pc *.L 


$(APNAME) : $(TARGET)
	$(LINKER) $(CFLAG) -o $(EXECDIR)$(APNAME) $(TARGET) $(LIBS)


%.o : %.F90
	$(FC) $(COPT) -c -o $@ $<
%.o : %.f90
	$(FC) $(COPT) -c -o $@ $<
%.o : %.f
	$(FC) $(COPT) -c -o $@ $<
svd.f :
	test -f $@ || wget -O svd.f --no-check-certificate https://netlib.org/eispack/svd.f
pythag.f :
	test -f $@ || wget -O pythag.f --no-check-certificate https://netlib.org/eispack/3090vf/double/pythag.f


debug : $(TARGET)
	$(LINKER) $(CFLAG) -o $(EXECDIR)$(APNAME) $(TARGET) $(LIBS)

clean :
	$(RM) $(CLNTARGET)

realclean :
	$(RM) $(CLNTARGET) $(APNAME) *~
