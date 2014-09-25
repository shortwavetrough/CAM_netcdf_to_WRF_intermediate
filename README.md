CAM_netcdf_to_WRF_intermediate
==============================
Fortran program for reading CAM/CLM/POP/CIC NetCDF output and writing in the WPS intermediate format.

To compile and run the program on the ERDC Garnet machine:

(1) load the Intel module: >module load intel/14.0.2.144,

(2) compile the program: >ifort -c -CB -CU -ftrapuv -par_report0 -vec_report0 -heap-arrays -convert big_endian -I/opt/cray/netcdf/4.3.0/INTEL/130/include/ CAM_netcdf_to_WRF_intermediate.f90 ; ifort CAM_netcdf_to_WRF_intermediate.o -L/opt/cray/netcdf/4.3.0/INTEL/130/lib -lnetcdf -lnetcdff

Note* - the compile flag -heap-arrays is needed, else compile errors.

Note** - be sure the netcdf include and library locations are up to date.

Note*** - debug compile flags include the following: -O0 -check all -traceback -fstack-protector -assume protect_parens -implicitnone -debug -gen-interfaces -check arg_temp_created -ftrapuv -g -traceback 

(3) run the program: >a.out
