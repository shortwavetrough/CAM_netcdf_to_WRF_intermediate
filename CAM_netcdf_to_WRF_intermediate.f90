program CAM_netcdf_to_WRF_intermediate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fortran 90 program for reading CAM netcdf output and writing in
  ! the WPS intermediate format.
  !
  ! REQUIRED STEPS: (1) Run CAM and produce netcdf output files with
  ! the variables specified in the cam-header.txt and clm-header.txt
  ! files in the /Doc subdirectory here.  (2) Run the
  ! Input/create_read_file.sh to create an input file for this FORTRAN
  ! program.  (3) Run this FORTRAN program.  (4) Run metgrid.exe.  (5)
  ! Run real.exe.  You're done.
  !
  ! To compile and run this program on our harvard cluster, use the
  ! following shell command (use crest rather than swell): Note that
  ! this uses the ifort compiler.  If you use another one, make sure
  ! you write in big endian format using either compiler switches or
  ! open statement below:
  !
  !! cd ~/WRF_Mauricio/; ifort -c -CB -CU -ftrapuv -par_report0 -vec_report0 -heap-arrays -I/opt/netcdf-3.6.0-p1/include/ CAM_netcdf_to_WRF_intermediate.f90; ifort CAM_netcdf_to_WRF_intermediate.o -L/opt/netcdf-3.6.0-p1/lib/ -lnetcdff; ./a.out
  !
  ! BJF - compile
  ! ifort -c -CB -CU -ftrapuv -par_report0 -vec_report0 -heap-arrays -convert big_endian
  ! -I/opt/cray/netcdf/4.3.0/INTEL/130/include/ CAM_netcdf_to_WRF_intermediate.f90 ;
  ! ifort CAM_netcdf_to_WRF_intermediate.o -L/opt/cray/netcdf/4.3.0/INTEL/130/lib -lnetcdf -lnetcdff
  !
  ! BJF - add for debug
  ! -O0 -check all -traceback -fstack-protector -assume protect_parens -implicitnone
  ! -debug -gen-interfaces -check arg_temp_created -ftrapuv -g -traceback
  !
  ! cd ~/WRF_Mauricio/Output/; ./metgrid.exe
  ! ln -s met_em*.nc real/
  ! cd ~/WRF_Mauricio/Output/real/; ./real.exe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use netcdf
  implicit none

  ! Declarations:
  integer, parameter :: outfile_intermediate = 10
  integer, parameter :: outfile_intermediate_SST = 11
  integer, parameter :: outfile_diagnostics = 16
  integer, parameter :: infile_CAM_files_and_dates = 15
  character(len=24) :: HDATE

  ! dimensions:
  integer, parameter :: nx_CAM=288,ny_CAM=192,nz_CAM=26 &
       ,nfields=5,nfields2d=9,nfields2d_to_read=5 &
       ,nz_soil=4,nz_CLM=1,nfields_soil=2
  integer, parameter :: nz_WRF=38
!  integer, parameter :: nz_WRF=27 ! (26 levels plus surface value)
  character(len=128) :: netcdf_cam_filename,netcdf_clm_filename,netcdf_pop_filename
  character(len=128) :: netcdf_ice_filename
  integer :: iEOF
  logical :: EOF
  logical :: EOD

  ! (1) Things read from the netcdf file need to be real*8 if declared
  ! there as DOUBLE and real if declared FLOAT; those written to
  ! intermediate format are real.  (2) Also, note that the order of
  ! dimension in the fortran program (nx_CAM,ny_CAM,nz_CAM,1)
  ! needs to be the reverse of that in the netcdf header
  ! (time,lev,lat,lon)!!
  real(8) :: lon(nx_CAM), lat(ny_CAM)
  real :: PS(nx_CAM,ny_CAM)
  ! pressure variables:
  real(8) :: P(nx_CAM,ny_CAM,nz_CAM), log_P(nx_CAM,ny_CAM,nz_CAM)
  ! interpolation variables:
  real(8) :: P_int(nz_WRF),log_P_int(nz_WRF)
  real :: field_data_int(nfields,nx_CAM,ny_CAM,nz_WRF,1)
  real, dimension(nfields,nx_CAM,ny_CAM,nz_CAM,1) :: field_data
  real, dimension(nfields2d,nx_CAM,ny_CAM,1) :: field2d_data
  real, dimension(nfields_soil,nx_CAM,ny_CAM,nz_soil,1) :: field_soil_data
  ! presure levels to which data are interpolated, these levels are
  ! taken from a standard case of WRF/WPS:
  data P_int /200100,100000,97500,95000,92500,90000,87500,85000,82500,80000 &
       ,77500,75000,70000,65000,60000,55000,50000,45000,40000,35000,30000 &
       ,25000,22500,20000,17500,15000,12500,10000,7000,5000,3000,2000,1000 &
       ,700,500,300,200,100/
!  data P_int /200100,100000,97500,95000,92500,90000,85000,80000 &
!       ,75000,70000,65000,60000,55000,50000,45000,40000,35000 &
!       ,30000,25000,20000,15000,10000,7000,5000,3000,2000,1000/
  character(len=25) :: field_units(nfields),field2d_units(nfields2d) &
       ,field_soil_units(nfields_soil)
  character(len=46) :: field_DESC(nfields),field2d_DESC(nfields2d) &
       ,field_soil_DESC(nfields_soil)
  character(len=9) :: field2d_name_to_output(nfields2d) &
       ,field_name_to_output(nfields),field_soil_name_to_output(nfields,nz_soil)

  ! specify the order of 2d and 3d variables in the data, units and title
  ! arrayes by specifying the integer pointers here:
  integer, parameter :: i2d_PSFC=1,i2d_PMSL=2,i2d_landsea=3,i2d_SKINTEMP=4 &
       ,i2d_TT=5,i2d_RH=6,i2d_UU=7,i2d_VV=8,i2d_SEAICE=9 ! XX ,i2d_SOILHGT=9
  integer, parameter :: i3d_TT=1,i3d_RH=2,i3d_UU=3,i3d_VV=4,i3d_GHT=5
  integer :: ifield,k

  ! 2d names expected by metgridl
  field2d_name_to_output(i2d_PSFC)    =   'PSFC     '  !
  field2d_name_to_output(i2d_PMSL)    =   'PMSL     '  !
  field2d_name_to_output(i2d_landsea) =   'LANDSEA  '  !
  field2d_name_to_output(i2d_SKINTEMP)=   'SKINTEMP '  !
  field2d_name_to_output(i2d_TT)      =   'TT       '  !  ! at 2m
  field2d_name_to_output(i2d_RH)      =   'RH       '  !  ! at 2m
  field2d_name_to_output(i2d_UU)      =   'UU       '  !  ! at 10M
  field2d_name_to_output(i2d_VV)      =   'VV       '  !  ! at 10M
  field2d_name_to_output(i2d_SEAICE)  =   'SEAICE   '  !
  !XX  field2d_name_to_output(i2d_SOILHGT) =   'SOILHGT  '  !

  ! names of 3d variables expected by metgrid
  ! (temp, relative humidity, u, v, geopotential height):
  field_name_to_output(i3d_TT) ='TT       '
  field_name_to_output(i3d_RH) ='RH       '
  field_name_to_output(i3d_UU) ='UU       '
  field_name_to_output(i3d_VV) ='VV       '
  field_name_to_output(i3d_GHT)='GHT      '

  ! names of soil variables required by Noah LSM (land surface model):
  ifield=1; k=1; field_soil_name_to_output(ifield,k)='SM000010'
  ifield=1; k=2; field_soil_name_to_output(ifield,k)='SM010040'
  ifield=1; k=3; field_soil_name_to_output(ifield,k)='SM040100'
  ifield=1; k=4; field_soil_name_to_output(ifield,k)='SM100200'
  ifield=2; k=1; field_soil_name_to_output(ifield,k)='ST000010'
  ifield=2; k=2; field_soil_name_to_output(ifield,k)='ST010040'
  ifield=2; k=3; field_soil_name_to_output(ifield,k)='ST040100'
  ifield=2; k=4; field_soil_name_to_output(ifield,k)='ST100200'

  ! the required variables by WRF are given at
  ! http://www.mmm.ucar.edu/wrf/OnLineTutorial/Basics/UNGRIB/ungrib_req_fields.htm

  ! open outpuf log file:
  open(outfile_diagnostics,form='formatted',file="Output/CCSM2WRF.log")

  ! read the first date and netcdf file name from the input file:
  open(infile_CAM_files_and_dates,form='formatted',file="Input/CCSM2WRF.input")
  read(infile_CAM_files_and_dates,*,iostat=iEOF) netcdf_cam_filename,netcdf_clm_filename,&
                        netcdf_pop_filename,netcdf_ice_filename,hdate
  if (iEOF<0) then;
     EOF=.true.;
  else;
     EOF=.false.;
  end if

  ! Loop over all CAM netcdf files and dates
  ! specified in unit infile_CAM_files_and_dates:
  ! =============================================
  do while (.not.EOF)

     write(*,*) "processing date=",hdate
     write(outfile_diagnostics,*) "processing CAM/CLM/POP/CIC files=" &
          ,netcdf_cam_filename,netcdf_clm_filename,netcdf_pop_filename,netcdf_ice_filename,"; date=",hdate

EOD=.TRUE.
     call read_netcdf_files(nz_WRF,hdate,outfile_diagnostics,outfile_intermediate &
                    ,outfile_intermediate_SST,netcdf_cam_filename &
                    ,netcdf_clm_filename,netcdf_pop_filename &
                    ,netcdf_ice_filename,nx_CAM,ny_CAM,nz_CAM,nfields,nfields2d &
                    ,i3d_TT,i3d_RH,i3d_UU,i3d_VV,i3d_GHT,field_units,field_DESC &
                    ,i2d_PSFC,i2d_PMSL,i2d_landsea,i2d_SKINTEMP,i2d_TT,i2d_RH,i2d_UU &
                    ,i2d_VV,i2d_SEAICE,field2d_units,field2d_DESC,field_soil_units &
                    ,field_soil_DESC,nfields_soil,field_data,field_name_to_output &
                    ,field2d_data,PS,field2d_name_to_output,P,log_P,log_P_int,P_int &
                    ,field_soil_data,nz_soil,field_soil_name_to_output)
     call interpolate_to_pressure(ifield,nfields,nx_CAM,ny_CAM,nz_CAM,log_P &
                     ,field_data,nz_WRF,P_int,PS,log_P_int &
                     ,outfile_diagnostics,field_data_int &
                     ,field_name_to_output)
     call write_intermediate(nz_WRF,outfile_diagnostics,nfields &
                             ,field_name_to_output,outfile_intermediate &
                             ,field_units,field_DESC,P_int,nx_CAM &
                             ,ny_CAM,HDATE,lat,lon,field_data_int &
                             ,nfields2d,field2d_name_to_output &
                             ,field2d_units,field2d_DESC,field2d_data &
                             ,nz_soil,nfields_soil,field_soil_name_to_output &
                             ,field_soil_units,field_soil_DESC,field_soil_data &
                             ,outfile_intermediate_SST,i2d_SKINTEMP)
 if (EOD) then
   STOP 'End of dummy'
 end if
     ! read next CAM netcdf filename and date:
     read(infile_CAM_files_and_dates,*,iostat=iEOF) netcdf_cam_filename,netcdf_clm_filename &
                   ,netcdf_pop_filename,netcdf_ice_filename,hdate
     if (iEOF<0) then;
        EOF=.true.;
        write(outfile_diagnostics,*) ,"reached EOF for unit infile_CAM_files_and_dates."
     else
        EOF=.false.;
     end if

  end do ! while loop over reading of CAM filenames and dates
  write(outfile_diagnostics,'(/,"End of read loop.  Program finished.")')

  stop
end program CAM_netcdf_to_WRF_intermediate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE spline(X,Y,nz_CAM,yp1,ypn,Y2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  INTEGER :: nz_CAM
  INTEGER, parameter :: NMAX=500
  REAL(8) :: yp1,ypn,X(nz_CAM),Y(nz_CAM),Y2(nz_CAM)
  INTEGER :: i,k
  REAL(8) :: p,qn,sig,un,u(NMAX)
  if (yp1.gt..99e30) then
     Y2(1)=0.
     u(1)=0.
  else
     Y2(1)=-0.5
     u(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-yp1)
  end if
  do i=2,nz_CAM-1
     sig=(X(i)-X(i-1))/(X(i+1)-X(i-1))
     p=sig*Y2(i-1)+2.
     Y2(i)=(sig-1.)/p
     u(i)=(6.*((Y(i+1)-Y(i))/(X(i+1)-X(i))-(Y(i)-Y(i-1)) &
          /(X(i)-X(i-1)))/(X(i+1)-X(i-1))-sig*u(i-1))/p
  end do

  if (ypn.gt..99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(X(nz_CAM)-X(nz_CAM-1)))*(ypn-(Y(nz_CAM)-Y(nz_CAM-1))/(X(nz_CAM)-X(nz_CAM-1)))
  end if
  Y2(nz_CAM)=(un-qn*u(nz_CAM-1))/(qn*Y2(nz_CAM-1)+1.)
  do k=nz_CAM-1,1,-1
     Y2(k)=Y2(k)*Y2(k+1)+u(k)
  end do
  return
end SUBROUTINE spline


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE splint(X,Y,Y2,nz_CAM,XINT,YINT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  INTEGER :: nz_CAM
  REAL(8) :: XINT,YINT,X(nz_CAM),Y2(nz_CAM),Y(nz_CAM)
  INTEGER :: k,khi,klo
  REAL(8) :: a,b,h

  ! Eli: avoid actual extrapolation by using the end values:
  if (XINT<X(1)) then
     YINT=Y(1);
  elseif (XINT>X(nz_CAM)) then
     YINT=Y(nz_CAM);
  else
     ! Eli: end of my addition here.
     klo=1
     khi=nz_CAM
1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(X(k).gt.XINT)then
           khi=k
        else
           klo=k
        end if
        GO TO 1
     end if
     h=X(khi)-X(klo)
     if (h.eq.0.) then
        write(*,'("bad X input in splint, type [enter] to continue")')
        read(*,*)
     end if
!     if (h.eq.0.) pause 'bad X input in splint'
     a=(X(khi)-XINT)/h
     b=(XINT-X(klo))/h
     YINT=a*Y(klo)+b*Y(khi)+ &
          ((a**3-a)*Y2(klo)+(b**3-b)*Y2(khi))*(h**2)/6.
  end if

  return
end SUBROUTINE splint



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE HANDLE_ERR(STATUS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use netcdf
  implicit none
  INTEGER STATUS
  IF (STATUS .NE. NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     STOP 'Stopped'
  ENDIF
END SUBROUTINE HANDLE_ERR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine read_netcdf_files &
     (nz_WRF,hdate,outfile_diagnostics,outfile_intermediate &
     ,outfile_intermediate_SST,netcdf_cam_filename &
     ,netcdf_clm_filename,netcdf_pop_filename &
     ,netcdf_ice_filename,nx_CAM,ny_CAM,nz_CAM,nfields,nfields2d &
     ,i3d_TT,i3d_RH,i3d_UU,i3d_VV,i3d_GHT,field_units,field_DESC &
     ,i2d_PSFC,i2d_PMSL,i2d_landsea,i2d_SKINTEMP,i2d_TT,i2d_RH,i2d_UU &
     ,i2d_VV,i2d_SEAICE,field2d_units,field2d_DESC,field_soil_units &
     ,field_soil_DESC,nfields_soil,field_data,field_name_to_output &
     ,field2d_data,PS,field2d_name_to_output,P,log_P,log_P_int,P_int &
     ,field_soil_data,nz_soil,field_soil_name_to_output)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !DEC$ ATTRIBUTES REFERENCE :: HDATE
  use netcdf
  implicit none
  integer :: nz_WRF
  integer :: nx_CAM,ny_CAM,nz_CAM,nfields,nfields2d,nfields_soil,nz_soil
  character(len=128) :: filename
  character(len=24) :: HDATE
  integer :: outfile_diagnostics,outfile_intermediate,outfile_intermediate_SST
  integer :: STATUS, NCID, NCID_clm, NCID_pop, NCID_ice
  character(len=128) :: netcdf_cam_filename, netcdf_clm_filename, netcdf_pop_filename
  character(len=128) :: netcdf_ice_filename
  integer :: lat_var_id,lon_var_id,lev_var_id,time_var_id,nlon_var_id &
            ! 3d fields
            ,T_var_id,RH_var_id,U_var_id,V_var_id,GEOP_var_id &
            ! 2d fields
            ,PS_var_id,PSL_var_id,LANDFRAC_var_id,TS_var_id,SEAICE_var_id &
            ! pressure variables
            ,P0_var_id,hyam_var_id,hybm_var_id &
            ! soil variables
            ,SM_var_id,ST_var_id &
            ! surface geopotential
            ,TOPO_var_id
  integer :: field_var_id(nfields),field2d_var_id(nfields2d) &
            ,lon_lat_netcdf_units_length,lon_lat_netcdf_title_length &
            ,ifield,i,j,k
  character(len=12) :: lon_netcdf_units
  character(len=9) :: lon_netcdf_title
  character(len=13) :: lat_netcdf_units
  character(len=8) :: lat_netcdf_title
  integer :: i3d_TT,i3d_RH,i3d_UU,i3d_VV,i3d_GHT
  integer :: i2d_PSFC,i2d_PMSL,i2d_landsea,i2d_SKINTEMP,i2d_TT,i2d_RH &
            ,i2d_UU,i2d_VV,i2d_SEAICE
  character(len=25) :: field_units(nfields),field2d_units(nfields2d) &
                      ,field_soil_units(nfields_soil)
  character(len=46) :: field_DESC(nfields),field2d_DESC(nfields2d) &
                      ,field_soil_DESC(nfields_soil)
  character(len=9) :: field_name_to_output(nfields),field2d_name_to_output(nfields2d) &
                     ,field_soil_name_to_output(nfields,nz_soil)
  real(8) :: lon(nx_CAM),lat(ny_CAM)
  real, dimension(nx_CAM,ny_CAM,nz_CAM,1) :: field3d ! T(time,lev,lat,lon)
  real, dimension(nx_CAM,ny_CAM,1) :: field2d ! T(time,lat,lon)
  real, dimension(nx_CAM,ny_CAM,1) :: field_soil ! TSOI(time,lev,lat,lon)
  real, dimension(nfields,nx_CAM,ny_CAM,nz_CAM,1) :: field_data
  real, dimension(nfields2d,nx_CAM,ny_CAM,1) :: field2d_data
  real, dimension(nfields_soil,nx_CAM,ny_CAM,nz_soil,1) :: field_soil_data
  real, dimension(nx_CAM,ny_CAM,nz_CAM,1) :: spechumd
  real, dimension(nx_CAM,ny_CAM) :: topo,t2m,tsfc
  real, dimension(nx_CAM,ny_CAM) :: tbot,ubot,vbot,z3bot,esat,pressure,qsat
  real :: PS(nx_CAM,ny_CAM)
  real :: field_max,field_min
  real(8) :: P0,hyam(nz_CAM),hybm(nz_CAM),P(nx_CAM,ny_CAM,nz_CAM),P_avg(nz_CAM) &
            ,log_P(nx_CAM,ny_CAM,nz_CAM),log_P_avg(nz_CAM)
  real(8) :: P_int(nz_WRF),log_P_int(nz_WRF)

  ! open output files for metgrid in WRF/WPS intermediate format:
  write(filename,'("Output/FILE:",A13)') hdate(1:13)
  write(outfile_diagnostics,*) "output intermediate file filename=",filename
  open(outfile_intermediate,form='unformatted',file=filename)

  write(filename,'("Output/SST:",A13)') hdate(1:13)
  write(outfile_diagnostics,*) "output intermediate SST file filename=",filename
  open(outfile_intermediate_SST,form='unformatted',file=filename)

  ! CAM
  STATUS = NF90_OPEN(netcdf_cam_filename, 0, NCID)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  ! CLM
  STATUS = NF90_OPEN(netcdf_clm_filename, 0, NCID_clm)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  ! POP
  STATUS = NF90_OPEN(netcdf_pop_filename, 0, NCID_pop)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  ! ICE
  STATUS = NF90_OPEN(netcdf_ice_filename, 0, NCID_ice)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  write(outfile_diagnostics,*) "done nf90_open"

  ! read netcdf data for all levels:
  ! ================================

  ! get dimension IDs: (only needed if we don't know the dimension and
  ! want to read them from the netcdf file)
  ! STATUS = NF_INQ_DIMID(NCID, 'T', dim_id)
  ! IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)

  ! get variable IDs
  ! 3d:
  STATUS = NF90_INQ_VARID(NCID, 'lat', lat_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'lon', lon_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'lev', lev_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'time', time_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'nlon', nlon_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'T', T_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'Q', RH_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'U', U_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'V', V_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'Z3', GEOP_var_id)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

  ! 2d:
  STATUS = NF90_INQ_VARID(NCID, 'PS', PS_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'PSL', PSL_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'LANDFRAC', LANDFRAC_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'PHIS', TOPO_var_id)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

  ! pressure variables:
  STATUS = NF90_INQ_VARID(NCID, 'P0', P0_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'hyam', hyam_var_id)
  STATUS = NF90_INQ_VARID(NCID, 'hybm', hybm_var_id)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

  ! soil temp and moisture:
  STATUS = NF90_INQ_VARID(NCID_clm, 'SOILWATER_10CM', SM_var_id)
  STATUS = NF90_INQ_VARID(NCID_clm, 'TSOI_10CM', ST_var_id)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

  ! Ying.Liu sst:
  STATUS = NF90_INQ_VARID(NCID_pop, 'tos', TS_var_id)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

  ! Ying.Liu ice:
  STATUS = NF90_INQ_VARID(NCID_ice, 'aice_d', SEAICE_var_id)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

  ! 3d:
  field_var_id(i3d_TT) = T_var_id
  field_var_id(i3d_RH) = RH_var_id
  field_var_id(i3d_UU) = U_var_id
  field_var_id(i3d_VV) = V_var_id
  field_var_id(i3d_GHT)= GEOP_var_id

  ! 2d:
  field2d_var_id(i2d_PSFC)    = PS_var_id
  field2d_var_id(i2d_PMSL)    = PSL_var_id
  field2d_var_id(i2d_landsea) = LANDFRAC_var_id
  field2d_var_id(i2d_SKINTEMP)= TS_var_id
  field2d_var_id(i2d_TT)      = T_var_id
  field2d_var_id(i2d_SEAICE)  = SEAICE_var_id

  write(outfile_diagnostics,*) "done NF90_INQ_VARID, field_var_id=",field_var_id
  write(outfile_diagnostics,*) "field2d_var_id=",field2d_var_id

  ! get attribute values: title and units
  ! =====================================
  ! get units and titles for lon:
  STATUS = NF90_INQUIRE_ATTRIBUTE(NCID, lon_var_id, 'units', lon_lat_netcdf_units_length)
  STATUS = NF90_INQUIRE_ATTRIBUTE(NCID, lon_var_id, 'long_name', lon_lat_netcdf_title_length)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  write(outfile_diagnostics,*) " lon: netcdf units length=",lon_lat_netcdf_units_length &
       ,"; netcdf title length=",lon_lat_netcdf_title_length
  STATUS = NF90_GET_ATT(NCID, lon_var_id, 'units', lon_netcdf_units)
  STATUS = NF90_GET_ATT(NCID, lon_var_id, 'long_name', lon_netcdf_title)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  write(outfile_diagnostics,*) "netcdf title=",lon_netcdf_title,"; units=",lon_netcdf_units
  ! get units and titles for lat:
  STATUS = NF90_INQUIRE_ATTRIBUTE(NCID, lat_var_id, 'units', lon_lat_netcdf_units_length)
  STATUS = NF90_INQUIRE_ATTRIBUTE(NCID, lat_var_id, 'long_name', lon_lat_netcdf_title_length)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  write(outfile_diagnostics,*) " lat: netcdf units length=",lon_lat_netcdf_units_length &
       ,"; netcdf title length=",lon_lat_netcdf_title_length
  STATUS = NF90_GET_ATT(NCID, lat_var_id, 'units', lat_netcdf_units)
  STATUS = NF90_GET_ATT(NCID, lat_var_id, 'long_name', lat_netcdf_title)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  write(outfile_diagnostics,*) "netcdf title=",lat_netcdf_title,"; units=",lat_netcdf_units

  !Ying.Liu set units(len=25) and descriptions(len=46) for 3d fields
  field_units(i3d_TT)="K                        "  ! TT
  field_DESC(i3d_TT)="CAM temperature                               "
  field_units(i3d_RH)="%                        "  ! RH
  field_DESC(i3d_RH)="RH computed from CAM specific humidity        "
  field_units(i3d_UU)="m/s                      "  ! UU
  field_DESC(i3d_UU)="CAM U                                         "
  field_units(i3d_VV)="m/s                      "  ! VV
  field_DESC(i3d_VV)="CAM V                                         "
  field_units(i3d_GHT)="m                        "  ! GHT
  field_DESC(i3d_GHT)="Geopotential Height (above sea level)         "

  !Ying.Liu set units(len=25) and descriptions(len=46) for 2d fields
  field2d_units(i2d_PSFC)     ="Pa                       "  !
  field2d_DESC(i2d_PSFC)      ="Surface pressure                              "
  field2d_units(i2d_PMSL)     ="Pa                       "  !
  field2d_DESC(i2d_PMSL)      ="Sea level pressure                            "
  field2d_units(i2d_landsea)  ="percent                  "
  field2d_DESC(i2d_landsea)   ="land mask from CLM                            "
  field2d_units(i2d_SKINTEMP) ="K                        "  ! SKINTEMP
  field2d_DESC(i2d_SKINTEMP)  ="regridded SST from POP to CAM                 "
  field2d_units(i2d_TT)       ="K                        "  ! T2m
  field2d_DESC(i2d_TT)        ="VerticalExtrapolate to 2m from CAM T at lev=26"
  field2d_units(i2d_RH)       ="%                        "  ! RH
  field2d_DESC(i2d_RH)        ="computed RH from CAM SPECHUMD at lev=26       "
  field2d_units(i2d_UU)       ="m/s                      "  ! UU
  field2d_DESC(i2d_UU)        ="VerticalExtrapolate to 10m from CAM U at lev26"
  field2d_units(i2d_VV)       ="m/s                      "  ! VV
  field2d_DESC(i2d_VV)        ="VerticalExtrapolate to 10m from CAM V at lev26"
  field2d_units(i2d_SEAICE)   ="percent                  "  !
  field2d_DESC(i2d_SEAICE)    ="CICE ice area (aggregate)                     "

  ! Ying.Liu set units and titles for soil variables:
  write(outfile_diagnostics,*) " getting units and title for soil fields"
  ifield=1
  field_soil_units(ifield)="mm3/mm3                  "  ! moisture
  field_soil_DESC(ifield)="volumetric soil water                         "
  ifield=2
  field_soil_units(ifield)="K                        "  ! temperature
  field_soil_DESC(ifield)="soil temperature                              "

  ! get values of the needed variables from netcdf file:
  ! ====================================================

  ! first, get lon/ lat:
  STATUS = NF90_GET_VAR(NCID, lon_var_ID, lon)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  write(outfile_diagnostics,*) "lon=",lon
  STATUS = NF90_GET_VAR(NCID, lat_var_ID, lat)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  write(outfile_diagnostics,*) "lat=",lat

  ! get 3d T,RH,U,V:
  do ifield=1,nfields
     STATUS = NF90_GET_VAR(NCID, field_var_ID(ifield), field3d)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
     if(ifield.ne.i3d_RH) then
       field_data(ifield,:,:,:,1)=field3d(:,:,:,1)
     else    !specific humidity to Relative humidity
       spechumd(:,:,:,1)=field3d(:,:,:,1)
     endif
! Uncomment to check values, otherwise comment out for optimal IO performance
!     if(ifield.ne.i3d_RH) then
!       write(outfile_diagnostics,*) ,"ifield=",ifield,"; name=",field_name_to_output(ifield) &
!            ,"; field3d(64,32,:,1)=",field3d(64,32,:,1) &
!            ,"; field3d(64,:,nz_CAM,1)=",field3d(64,:,nz_CAM,1)
!     endif
  end do
  write(outfile_diagnostics,*) "done getting T,RH,U,V from netcdf"

  ! get 2d fields:
  ! first those that are in the netcdf file as 2d fields:
  do ifield=1,2   ! 1=i2d_PSFC  2=i2d_PMSL
     STATUS = NF90_GET_VAR(NCID, field2d_var_id(ifield), field2d)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
     field2d_data(ifield,:,:,1)=field2d(:,:,1)
  end do
  PS=field2d_data(i2d_PSFC,:,:,1)

  STATUS = NF90_GET_VAR(NCID, field2d_var_id(i2d_landsea), field2d)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  field2d_data(i2d_landsea,:,:,1)=field2d(:,:,1)

  ! Ying.Liu fix land mask: WRF cannot deal with mask values that are not 0 or 1:
!  write(outfile_diagnostics,*) ,"land mask before fixing: i2d_landsea(64,:,1)=" &
!       ,field2d_data(i2d_landsea,64,:,1)
  do i=1,nx_CAM; do j=1,ny_CAM;
     if (field2d_data(i2d_landsea,i,j,1)<0.5) then
        field2d_data(i2d_landsea,i,j,1)=0; else; field2d_data(i2d_landsea,i,j,1)=1;
     end if
  end do; end do
!  write(outfile_diagnostics,*) ,"land mask after  fixing: i2d_landsea(64,:,1)=" &
!       ,field2d_data(i2d_landsea,64,:,1)
  write(outfile_diagnostics,*) ,"done fixing land mask to be only 0 or 1."

  STATUS = NF90_GET_VAR(NCID_ice, field2d_var_id(i2d_SEAICE), field2d)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  field2d_data(i2d_SEAICE,:,:,1)=0.01*field2d(:,:,1)

  STATUS = NF90_GET_VAR(NCID_pop, field2d_var_id(i2d_SKINTEMP), field2d)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  field2d_data(i2d_SKINTEMP,:,:,1)=field2d(:,:,1)+273.15

  STATUS = NF90_INQ_VARID(NCID, 'PHIS', TOPO_var_id)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  STATUS = NF90_GET_VAR(NCID, TOPO_var_id, field2d)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  topo=field2d(:,:,1)/9.81    ! "PHIS" units=m^2/s^2

  ! get fields needed to calculate 3d pressure:
  STATUS = NF90_GET_VAR(NCID, hyam_var_id, hyam)
  STATUS = NF90_GET_VAR(NCID, hybm_var_id, hybm)
  STATUS = NF90_GET_VAR(NCID, P0_var_id, P0)
  write(outfile_diagnostics,*) "P0=",P0
  write(outfile_diagnostics,*) "hyam=",hyam
  write(outfile_diagnostics,*) "hybm=",hybm
  write(outfile_diagnostics,*) "done getting pressure fields from netcdf"

  ! Ying.Liu next, those 2d fields that need to be obtained by extrapolating
  ! some 3d Field to the surface:
   tbot = field_data(i3d_TT,:,:,nz_CAM,1) ! T
   ubot = field_data(i3d_UU,:,:,nz_CAM,1) ! UU
   vbot = field_data(i3d_VV,:,:,nz_CAM,1) ! VV
  z3bot = field_data(i3d_GHT,:,:,nz_CAM,1) ! Geop Height

  field2d_data(i2d_UU,:,:,1)= 10.*ubot/(z3bot-topo) ! U10m
  field2d_data(i2d_VV,:,:,1)= 10.*vbot/(z3bot-topo) ! V10m
  t2m=0.0065*(z3bot-(topo+2))+tbot ! T2m
  tsfc=0.0065*(z3bot-topo)+tbot ! surface temperature, use it on land as skin temperature
  field2d_data(i2d_TT,:,:,1)= t2m ! T2m

  do k=1,nz_CAM
    esat=611.2*exp(17.67*(field_data(i3d_TT,:,:,k,1)-273.15)/(field_data(i3d_TT,:,:,k,1)-273.15+243.5)) !saturated vapor pressure
    pressure=hyam(k)*P0+hybm(k)*PS !pressure
    qsat=0.622*esat/(pressure-0.378*esat)
    field_data(i3d_RH,:,:,k,1)=(spechumd(:,:,k,1)/qsat)*100. !relative humidity
  end do
  field2d_data(i2d_RH,:,:,1)=field_data(i3d_RH,:,:,nz_CAM,1)  !RH2m
!  write(outfile_diagnostics,*) ,"ifield=",i3d_RH,"; name=",field_name_to_output(i3d_RH) &
!          ,"; field3d(64,32,:,1)=",field_data(i3d_RH,64,32,:,1) &
!          ,"; field3d(64,:,nz_CAM,1)=",field_data(i3d_RH,64,:,nz_CAM,1)

  ! Ying.Liu substitude skin temperature on land with Tsfc: use landsea
  do i=1,nx_CAM; do j=1,ny_CAM;
     if (field2d_data(i2d_landsea,i,j,1).eq.1) then
        field2d_data(i2d_SKINTEMP,i,j,1)=tsfc(i,j)
     end if
  end do; end do

  ! print some diagnostics about 2d fields:
  do ifield=1,nfields2d
     field_max=-1.e10
     do i=1,nx_CAM;do j=1,ny_CAM;
        field_max=amax1(field_max,field2d_data(ifield,i,j,1))
     end do; end do;
     field_min=1.e10
     do i=1,nx_CAM;do j=1,ny_CAM;
        field_min=amin1(field_min,field2d_data(ifield,i,j,1))
     end do; end do;
     write(outfile_diagnostics,*) ,"2d fields: ifield=",ifield &
          ,"; field name=",field2d_name_to_output(ifield) &
          ,"; max=",field_max,"; min=",field_min
!     write(outfile_diagnostics,*) ,"ifield=",ifield,"; 2d field name to output=" &
!          ,field2d_name_to_output(ifield) &
!          ,"; field2d_data(ifield,1, 1,1)=",field2d_data(ifield,1,1, 1) &
!          ,"; field2d_data(ifield,64,:,1)=",field2d_data(ifield,64,:,1)
  end do
  write(outfile_diagnostics,*) ,"done getting 2d fields from netcdf"

  ! calculate 3d pressure field from sigma coordinates:
  ! P(k,j,i)=hyam(k)*P0+hybm(k)*PS(j,i)
  do k=1,nz_CAM; do i=1,nx_CAM; do j=1,ny_CAM;
     P(i,j,k)=hyam(k)*P0+hybm(k)*PS(i,j)
  end do; end do; end do
!  write(outfile_diagnostics,*) ,"P(64,32,:)=",P(64,32,:)

  ! calculate average pressure profile:
  do k=1,nz_CAM;
     P_avg(k)=0
     do i=1,nx_CAM; do j=1,ny_CAM;
        P_avg(k)=P_avg(k)+P(i,j,k)
     end do; end do;
     P_avg(k)=P_avg(k)/(nx_CAM*ny_CAM)
  end do

  ! Calculate log pressure and log average pressure, in preparation
  ! for interpolation from sigma to pressure coordinates:
  log_P=log10(P);
  log_P_avg=log10(P_avg);
  log_P_int=log10(P_int);
  write(outfile_diagnostics,*) ,"P_avg=",P_avg
  write(outfile_diagnostics,*) ,"log_P_avg=",log_P_avg
  write(outfile_diagnostics,*) ,"P_int=",P_int
  write(outfile_diagnostics,*) ,"log_P_int=",log_P_int
  write(outfile_diagnostics,*) ,"done calculating pressure."

  ! get soil variables from CLM file:
  ! first, set soil moisture:
  STATUS = NF90_GET_VAR(NCID_clm, SM_var_id, field_soil)
  write(outfile_diagnostics,*) ,"read in soil moisture from CLM"
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  ifield=1
  ! XX set WRF soil moisture values at all levels to be equal to
  ! the first CAM soil moisture level, note 1000 conversion factor
  ! from CLM's H2OSOI(mm3/mm3) field to WRF's SMXXXXXX (kg/m3):
  ! PMA calculates soil moisture at 4 levels using a stupid averaging method
  ! PMA says no need to mutiply by 1000. WRF Vtable.GFS is wrong by using Kg/m3.
  ! All soil moisture should have the unit as "fraction"
  ! CCSM soil moisture is in kg water/ m2

  field_soil_data(ifield,:,:,1,1)=field_soil(:,:,1)/100.0
  field_soil_data(ifield,:,:,2,1)=field_soil(:,:,1)/100.0
  field_soil_data(ifield,:,:,3,1)=field_soil(:,:,1)/100.0
  field_soil_data(ifield,:,:,4,1)=field_soil(:,:,1)/100.0
  ! second, set soil temperature:
  STATUS = NF90_GET_VAR(NCID_clm, ST_var_id, field_soil)
  IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
  ifield=2
  ! XX set WRF soil values at all levels to be equal to the first CAM soil level:
  ! PMA calculates soil temperature at 4 levels using a stupid averaging method
!  print *, field_soil(121,21,:)
  field_soil_data(ifield,:,:,1,1)=field_soil(:,:,1)
  field_soil_data(ifield,:,:,2,1)=field_soil(:,:,1)
  field_soil_data(ifield,:,:,3,1)=field_soil(:,:,1)
  field_soil_data(ifield,:,:,4,1)=field_soil(:,:,1)
!  print *, field_soil_data(2,121,21,:,1)
!  print *, field_soil_data(1,121,21,:,1)

  ! deal with missing soil values (1.e36):
  ! for moisture, set missing values to 0:
  ifield=1
  do i=1,nx_CAM; do j=1,ny_CAM; do k=1,nz_soil
     if (field_soil_data(ifield,i,j,k,1)>1.0e30) then
        field_soil_data(ifield,i,j,k,1)= 0.0;
     elseif (field_soil_data(ifield,i,j,k,1)>1.0) then
        field_soil_data(ifield,i,j,k,1)= 1.0;
     end if
  end do; end do; end do;
  ! temperature: set missing vaues to surface temperature:
  ifield=2
  do i=1,nx_CAM; do j=1,ny_CAM; do k=1,nz_soil
     if (field_soil_data(ifield,i,j,k,1)>1.0e30) then
        field_soil_data(ifield,i,j,k,1)=field2d_data(i2d_SKINTEMP,i,j,1)
     end if
  end do; end do; end do;
!  print*,field_soil_data(2,121,21,:,1)
!  print*,field_soil_data(1,121,21,:,1)

  ! print some diagnostics about soil fields:
  do ifield=1,nfields_soil
     do k=1,nz_soil
        write(outfile_diagnostics,*) ,"soil diagnostics for ifield=",ifield,"; k=",k &
             ,"; variable=",field_soil_name_to_output(ifield,k)
        field_max=-1.e10
        do i=1,nx_CAM; do j=1,ny_CAM;
           field_max=amax1(field_max,field_soil_data(ifield,i,j,k,1))
        end do; end do;
        field_min=1.e10
        do i=1,nx_CAM; do j=1,ny_CAM;
           field_min=amin1(field_min,field_soil_data(ifield,i,j,k,1))
        end do; end do;
        write(outfile_diagnostics,*) ,"max=",field_max,"; min=",field_min
!        write(outfile_diagnostics,*) ,"; field_soil_data(ifield,64,32,:,1)=" &
!             ,field_soil_data(ifield,64,32,:,1) &
!             ,"; field_soil_data(ifield,64,:,nz_soil,1)=" &
!             ,field_soil_data(ifield,64,:,nz_soil,1)
     end do
  end do

  ! close netCDF dataset
  status=NF90_CLOSE(NCID)
  status=NF90_CLOSE(NCID_clm)
  status=NF90_CLOSE(NCID_pop)
  status=NF90_CLOSE(NCID_ice)

! print *, "Leaving read subroutine, going to MAIN"
  return
end Subroutine read_netcdf_files


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine interpolate_to_pressure(ifield,nfields,nx_CAM,ny_CAM,nz_CAM,log_P &
                                  ,field_data,nz_WRF,P_int,PS,log_P_int &
                                  ,outfile_diagnostics,field_data_int &
                                  ,field_name_to_output)
  use netcdf
  implicit none
  integer :: nz_WRF,outfile_diagnostics
  integer :: ifield,nfields,i,j,k,nx_CAM,ny_CAM,nz_CAM
  character(len=9) :: field_name_to_output(nfields)
  real, dimension(nfields,nx_CAM,ny_CAM,nz_CAM,1) :: field_data
  real(8) :: yp1,ypn,X(nz_CAM),Y(nz_CAM),XINT,YINT,Y2(nz_CAM)
  real(8) :: log_P(nx_CAM,ny_CAM,nz_CAM)
  real(8) :: P_int(nz_WRF),log_P_int(nz_WRF)
  real :: PS(nx_CAM,ny_CAM)
  real :: field_data_int(nfields,nx_CAM,ny_CAM,nz_WRF,1)
  real :: field_max,field_min

  ! interpolate from sigma to pressure coordinates:
  ! "natural" b.c.:
  yp1=2.d30; ypn=2.d30
  do ifield=1,nfields
     do j=1,ny_CAM; do i=1,nx_CAM;
        X=log_P(i,j,:)
        Y=field_data(ifield,i,j,:,1)
        call spline(X,Y,nz_CAM,yp1,ypn,Y2)
        do k=1,nz_WRF
           if (abs(P_int(k)-200100).le.0.01) then
              ! surface level:
              XINT=log(PS(i,j))
           else
              XINT=log_P_int(k)
           end if
           call splint(X,Y,Y2,nz_CAM,XINT,YINT)
           ! make sure RH is not negative:
           if (ifield==2 .and. YINT<0) then
              write(outfile_diagnostics,*) ,"*** RH<0: RH=",YINT,"; i,j,k=",i,j,k
              YINT=0;
           end if
           field_data_int(ifield,i,j,k,1)=YINT
        end do
     end do; end do;

     ! print some diagnostics about interpolated 3d fields:
     field_max=-1.e10
     do i=1,nx_CAM;do j=1,ny_CAM;do k=1,nz_WRF;
        field_max=amax1(field_max,field_data_int(ifield,i,j,k,1))
     end do; end do; end do;
     field_min=1.e10
     do i=1,nx_CAM;do j=1,ny_CAM;do k=1,nz_WRF;
        field_min=amin1(field_min,field_data_int(ifield,i,j,k,1))
     end do; end do; end do;
     write(outfile_diagnostics,*) ,"interpolated fields: ifield=",ifield &
          ,"; field name=",field_name_to_output(ifield) &
          ,"; max=",field_max,"; min=",field_min
!     write(outfile_diagnostics,*) ,"; field_data_int(ifield,64,32,:,1)=" &
!          ,field_data_int(ifield,64,32,:,1) &
!          ,"; field_data_int(ifield,64,:,5,1)=" &
!          ,field_data_int(ifield,64,:,5,1)

  end do

  write(outfile_diagnostics,*) ,"done interpolating to pressure coordinates"

! print *, "Leaving interpolate subroutine, going to MAIN"
  return
end Subroutine interpolate_to_pressure


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine write_intermediate(nz_WRF,outfile_diagnostics,nfields &
                             ,field_name_to_output,outfile_intermediate &
                             ,field_units,field_DESC,P_int,nx_CAM &
                             ,ny_CAM,HDATE,lat,lon,field_data_int &
                             ,nfields2d,field2d_name_to_output &
                             ,field2d_units,field2d_DESC,field2d_data &
                             ,nz_soil,nfields_soil,field_soil_name_to_output &
                             ,field_soil_units,field_soil_DESC,field_soil_data &
                             ,outfile_intermediate_SST,i2d_SKINTEMP)

  ! the intermediate format is described in:
  ! http://www.mmm.ucar.edu/wrf/OnLineTutorial/Basics/IM_files/IM_wps.htm

  !====================================================================================!
  ! READ in your data from the original source - you need to add the reading code here !
  !                                                                                    !
  ! You need to allocate SLAB (this is a 2D array) and place each 2D slab here before  !
  ! you can write it out to into the intermadiate file format                          !
  !                                                                                    !
  ! Other information you need to know about your data:                                !
  !    Time at which data is valid                                                     !
  !    Forecast time of the data                                                       !
  !    Source of data - you can make something up, it is never used                    !
  !    Field name - NOTE THEY NEED TO MATCH THOSE EXPECTED BY METGRID                  !
  !    Units of field                                                                  !
  !    Description of data                                                             !
  !    Level of data - Pa, 200100 Pa is used for surface, and 201300 Pa is used        !
  !          for sea-level pressure                                                    !
  !    X dimension                                                                     !
  !    Y dimension                                                                     !
  !    Data projection - only recognize                                                !
  !         0:  Cylindrical Equidistant (Lat/lon) projection.                          !
  !         1:  Mercator projection.                                                   !
  !         3:  Lambert-conformal projection.                                          !
  !         4:  Gaussian projection.                                                   !
  !         5:  Polar-stereographic projection.                                        !
  !    Start location of data - "CENTER", "SWCORNER". "SWCORNER" is typical            !
  !    Start lat & long of data                                                        !
  !    Lat/Lon increment                                                               !
  !    Number of latitudes north of equator (for Gaussian grids)                       !
  !    Grid-spacing in x/y                                                             !
  !    Center long                                                                     !
  !    truelat1/2                                                                      !
  !    Has the winds been rotated                                                      !
  !====================================================================================!

  use netcdf
  implicit none
  integer :: nz_WRF,k,outfile_diagnostics,ifield,nfields,outfile_intermediate &
            ,outfile_intermediate_SST,nx_CAM,ny_CAM,i,j,nfields2d,nz_soil &
            ,nfields_soil,IERR,i2d_SKINTEMP
  character(len=9) :: field_name_to_output(nfields),field2d_name_to_output(nfields2d) &
                     ,field_soil_name_to_output(nfields,nz_soil)
  character(len=25) :: field_units(nfields),field2d_units(nfields2d) &
                      ,field_soil_units(nfields_soil)
  character(len=46) :: field_DESC(nfields),field2d_DESC(nfields2d) &
                      ,field_soil_DESC(nfields_soil)
  character(len=24) :: HDATE
  real, dimension(nfields2d,nx_CAM,ny_CAM,1) :: field2d_data
  real, dimension(nfields_soil,nx_CAM,ny_CAM,nz_soil,1) :: field_soil_data
  real(8) :: P_int(nz_WRF)
  real(8) :: lon(nx_CAM),lat(ny_CAM)
  real :: field_data_int(nfields,nx_CAM,ny_CAM,nz_WRF,1)

  character(len=32) :: MAP_SOURCE
  character(len=9) :: FIELD
  character(len=25) :: UNITS
  character(len=46) :: DESC
  character(len=8) :: STARTLOC
  integer :: NX
  integer :: NY
  integer :: IPROJ
  integer :: IFV=5
  real, dimension(nx_CAM,ny_CAM) :: SLAB
  real :: XFCST
  real :: XLVL
  real :: STARTLAT
  real :: STARTLON
  real :: NLATS
  real :: DELTALON
  real :: EARTH_RADIUS = 6367470. * .001
  logical :: IS_WIND_EARTH_REL = .FALSE.

  ! loop over all levels and output 3d fields to the intermediate format
  ! ====================================================================

  do k=1,nz_WRF

     write(outfile_diagnostics,*) ," outputing 3d fields for k=",k

     ! go over all 3d fields and output them to the intermediate format
     ! =================================================================
     do ifield=1,nfields

        write(outfile_diagnostics,*) ,"writing 3d fields in intermediate format: ifield=" &
             ,ifield, "; field=",field_name_to_output(ifield)

        write (outfile_intermediate, IOSTAT=IERR) IFV
        write(outfile_diagnostics,*) ,"done writing IFV"

        ! WRITE the second record, common to all projections:
        ! HDATE=  this is read from input file.
        XFCST=0.0
        MAP_SOURCE="CCSMrun_b40.20th.track1.1deg.012"
        FIELD=field_name_to_output(ifield)
        UNITS=field_units(ifield)
        DESC=field_DESC(ifield)
        XLVL=P_int(k)
        NX=nx_CAM
        NY=ny_CAM
        IPROJ=4 ! CAM data are on a gaussian grid  XX right?
        write (outfile_intermediate) HDATE,XFCST,MAP_SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(outfile_diagnostics,*) ,"done writing second record, date:" &
             , HDATE//"  ", FIELD,"; xlvl=",xlvl

        ! WRITE the third record, which depends on the projection:

        if (IPROJ == 4) then
           ! Gaussian projection
           STARTLOC="SWCORNER"
!           STARTLAT=lat(1)
           STARTLAT=89.284
           STARTLON=lon(1)
!           NLATS=180./192. ! number of latitudes north of equator
           NLATS=ny_CAM/2 ! number of latitudes north of equator
           DELTALON=lon(2)-lon(1)
           WRITE (outfile_intermediate) STARTLOC,STARTLAT,STARTLON,NLATS,DELTALON,EARTH_RADIUS
        else
           write(outfile_diagnostics,*) ," *** error: wrong projection"
           stop
        end if
        write(outfile_diagnostics,*) ,"done writing third record"

        WRITE (outfile_intermediate) IS_WIND_EARTH_REL
        write(outfile_diagnostics,*) ,"done writing IS_WIND_EARTH_REL=",IS_WIND_EARTH_REL

        do i=1,nx_CAM; do j=1,ny_CAM
           SLAB(i,j)=field_data_int(ifield,i,j,k,1)
        end do; end do
        WRITE (outfile_intermediate) SLAB
!        write(outfile_diagnostics,*) ,"done writing slab; slab(1,1)=",slab(1,1)

     end do ! do loop over 3d fields

  end do ! do loop over levels

  ! go over all 2d fields and output them to the intermediate format
  ! =================================================================
  do ifield=1,nfields2d

     write(outfile_diagnostics,*) ,"writing 2d fields in intermediate format: ifield=" &
          ,ifield, "; field=",field2d_name_to_output(ifield)

     write (outfile_intermediate, IOSTAT=IERR) IFV
     write(outfile_diagnostics,*) ,"done writing IFV"

     ! WRITE the second record, common to all projections:
     ! HDATE=  this is read from input file.
     XFCST=0.0
     MAP_SOURCE="CCSMrun_b40.20th.track1.1deg.012"
     FIELD=field2d_name_to_output(ifield)
     UNITS=field2d_units(ifield)
     DESC=field2d_DESC(ifield)
     if (ifield==1) then
        XLVL=200100.0
     elseif (ifield==2) then
        XLVL=201300.0
     elseif (ifield==3) then
        XLVL=200100.0
     elseif (ifield==4) then
        XLVL=200100.00
     elseif (ifield==5) then
        XLVL=200100.0
     elseif (ifield==6) then
        XLVL=200100.0
     elseif (ifield==7) then
        XLVL=200100.0
     elseif (ifield==8) then
        XLVL=200100.0
     elseif (ifield==9) then
        XLVL=200100.0
     end if
     NX=nx_CAM
     NY=ny_CAM
     IPROJ=4 ! CAM data are on a gaussian grid  XX right?
     write (outfile_intermediate) HDATE,XFCST,MAP_SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
     write(outfile_diagnostics,*) ,"done writing second record, date:" &
          , HDATE//"  ", FIELD,"; xlvl=",XLVL

     ! WRITE the third record, which depends on the projection:

     if (IPROJ == 4) then

        ! Gaussian projection
        STARTLOC="SWCORNER"
!        STARTLAT=lat(1)
        STARTLAT=89.284
        STARTLON=lon(1)
!        NLATS=180./192. ! number of latitudes north of equator
        NLATS=ny_CAM/2 ! number of latitudes north of equator
        DELTALON=lon(2)-lon(1)
        WRITE (outfile_intermediate) STARTLOC,STARTLAT,STARTLON,NLATS,DELTALON,EARTH_RADIUS
     else
        write(outfile_diagnostics,*) ," *** error: wrong projection"
        stop
     endif
     write(outfile_diagnostics,*) ,"done writing third record"

     WRITE (outfile_intermediate) IS_WIND_EARTH_REL
     write(outfile_diagnostics,*) ,"done writing IS_WIND_EARTH_REL=",IS_WIND_EARTH_REL

     do i=1,nx_CAM; do j=1,ny_CAM
        SLAB(i,j)=field2d_data(ifield,i,j,1)
     end do; end do
     WRITE(outfile_intermediate) SLAB
     write(outfile_diagnostics,*) ,"done writing slab"

  end do ! do loop over 2d fields

  ! loop over all soil levels and output 3d soil fields to the intermediate format
  ! ==============================================================================
  do k=1,nz_soil

     write(outfile_diagnostics,*) ," outputing soil variables for k=",k

     ! go over all soil fields and output them to the intermediate format
     ! =================================================================
     do ifield=1,nfields_soil

        write(outfile_diagnostics,*) ,"writing soil fields in intermediate format: ifield=" &
             ,ifield, "; field=",field_soil_name_to_output(ifield,k)

        write (outfile_intermediate, IOSTAT=IERR) IFV
        write(outfile_diagnostics,*) ,"done writing IFV"

        ! WRITE the second record, common to all projections:
        ! HDATE=  this is read from input file.
        XFCST=0.0
        MAP_SOURCE="CCSMrun_b40.20th.track1.1deg.012"
        FIELD=field_soil_name_to_output(ifield,k)
        UNITS=field_soil_units(ifield)
        DESC=field_soil_DESC(ifield)
        XLVL=200100
        NX=nx_CAM
        NY=ny_CAM
        IPROJ=4 ! CAM data are on a gaussian grid  XX right?
        write (outfile_intermediate) HDATE,XFCST,MAP_SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(outfile_diagnostics,*) ,"done writing second record, date:" &
             , HDATE//"  ", FIELD,"; xlvl=",XLVL

        ! WRITE the third record, which depends on the projection:

        if (IPROJ == 4) then

           ! Gaussian projection
           STARTLOC="SWCORNER"
!           STARTLAT=lat(1)
           STARTLAT=89.284
           STARTLON=lon(1)
!           NLATS=180./192. ! number of latitudes north of equator
           NLATS=ny_CAM/2 ! number of latitudes north of equator
           DELTALON=lon(2)-lon(1)
           WRITE (outfile_intermediate) STARTLOC,STARTLAT,STARTLON,NLATS,DELTALON,EARTH_RADIUS
        else
           write(outfile_diagnostics,*) ," *** error: wrong projection"
           stop
        end if
        write(outfile_diagnostics,*) ,"done writing third record"

        WRITE (outfile_intermediate) IS_WIND_EARTH_REL
        write(outfile_diagnostics,*) ,"done writing IS_WIND_EARTH_REL=",IS_WIND_EARTH_REL

        do i=1,nx_CAM; do j=1,ny_CAM
           SLAB(i,j)=field_soil_data(ifield,i,j,k,1)
        end do; end do
        WRITE (outfile_intermediate) SLAB
        write(outfile_diagnostics,*) ,"done writing slab"

     end do ! do loop over soil fields

  end do ! do loop over soil levels

  ! write SST to the intermediate format:
  ! =====================================

  write(outfile_diagnostics,*) ,"writing SST in intermediate format, date=" &
          ,HDATE
  write (outfile_intermediate_SST, IOSTAT=IERR) IFV
  write(outfile_diagnostics,*) ,"done writing IFV"

  ! WRITE the second record, common to all projections:
  ! HDATE=  this is read from input file.
  XFCST=0.0
  MAP_SOURCE="CCSMrun_b40.20th.track1.1deg.012"
  FIELD=field2d_name_to_output(i2d_SKINTEMP)
  UNITS=field2d_units(i2d_SKINTEMP)
  DESC=field2d_DESC(i2d_SKINTEMP)
  XLVL=200100
  NX=nx_CAM
  NY=ny_CAM
  IPROJ=4 ! CAM data are on a gaussian grid  XX right?
  write (outfile_intermediate_SST) HDATE,XFCST,MAP_SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
  write(outfile_diagnostics,*) ,"done writing second record, date:" &
          , HDATE//"  ", FIELD,"; xlvl=",xlvl

  ! WRITE the third record, which depends on the projection:

  if (IPROJ == 4) then
     ! Gaussian projection
     STARTLOC="SWCORNER"
!     STARTLAT=lat(1)
     STARTLAT=89.284
     STARTLON=lon(1)
     NLATS=ny_CAM/2 ! number of latitudes north of equator
     DELTALON=lon(2)-lon(1)
     WRITE (outfile_intermediate_SST) STARTLOC,STARTLAT,STARTLON,NLATS,DELTALON,EARTH_RADIUS
  else
     write(outfile_diagnostics,*) " *** error: wrong projection"
     stop
  endif
  write(outfile_diagnostics,*) ,"done writing third record"

  WRITE (outfile_intermediate_SST) IS_WIND_EARTH_REL
  write(outfile_diagnostics,*) ,"done writing IS_WIND_EARTH_REL=",IS_WIND_EARTH_REL

  SLAB=field2d_data(i2d_SKINTEMP,:,:,1)
  WRITE (outfile_intermediate_SST) SLAB
  write(outfile_diagnostics,*) ,"done writing slab; slab(1,1)=",slab(1,1)



  ! close output intermediate format files:
  close(outfile_intermediate)
  close(outfile_intermediate_SST)

! print *, "Leaving write subroutine, going to MAIN"
  return
end Subroutine write_intermediate
