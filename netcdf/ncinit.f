C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.

      subroutine ncinit(FILE_NAME,NLONS,NLATS,lats,lons,A,f,g,H0,hb)
      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      integer ncid,i,j

C     We are writing 3D data, a M x N lat-lon grid, with
C     unlimited timesteps of data.
      integer NDIMS, NRECS
      parameter (NDIMS = 3)
      integer NLATS, NLONS
C      parameter (NLATS = M, NLONS = N)
      integer lon_dimid, lat_dimid, rec_dimid

C     The start and count arrays will tell the netCDF library where to
C     write our data.
      integer start(NDIMS-1), count(NDIMS-1)
      double precision g,H0,A,f
      double precision hb(NLONS,NLATS)

C constants
      integer g_varid,A_varid,f_varid,H0_varid

C     These program variables hold the latitudes and longitudes.
      double precision lats(NLATS), lons(NLONS)
      integer lon_varid, lat_varid, rec_varid

C     We will create six netCDF variables, one each for h, hb, d, z, u and v
      integer h_varid
      integer hb_varid
      integer u_varid
      integer v_varid
      integer z_varid
      integer d_varid

      integer dimids(NDIMS)
      integer dimids2(NDIMS-1)

C     Error handling.
      integer retval

C     Create the file. 
      retval = nf_create(FILE_NAME, nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the dimensions. The record dimension is defined to have
C     unlimited length - it can grow as needed. In this example it is
C     the time dimension.
      retval = nf_def_dim(ncid, 'LAT', NLATS, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, 'LON', NLONS, lon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, 'TIME', NF_UNLIMITED, rec_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the coordinate variables. We will only define coordinate
C     variables for lat and lon.  Ordinarily we would need to provide
C     an array of dimension IDs for each variable's dimensions, but
C     since coordinate variables only have one dimension, we can
C     simply provide the address of that dimension ID (lat_dimid) and
C     similarly for (lon_dimid).
      retval = nf_def_var(ncid, 'LAT', NF_DOUBLE, 1, lat_dimid, 
     +     lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'LON', NF_DOUBLE, 1, lon_dimid, 
     +     lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'TIME', NF_DOUBLE, 1, rec_dimid, 
     +     rec_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C constants
      retval = nf_def_var(ncid, 'g', NF_DOUBLE, 0, 0, 
     +     g_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'A', NF_DOUBLE, 0, 0, 
     +     a_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'f', NF_DOUBLE, 0, 0, 
     +     f_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'H0', NF_DOUBLE, 0, 0, 
     +     H0_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to coordinate variables.
      retval = nf_put_att_text(ncid, lat_varid, 'units', 7, 
     +     'degrees')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, lon_varid, 'units', 7, 
     +     'degrees')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, rec_varid, 'units', 1, 
     +     's')
      if (retval .ne. nf_noerr) call handle_err(retval)

C constants
      retval = nf_put_att_text(ncid, g_varid, 'units', 4, 
     +     'm/s2')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, A_varid, 'units', 1, 
     +     'm')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, f_varid, 'units', 3, 
     +     '1/s')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, H0_varid, 'units', 1, 
     +     'm')
      if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the dimids of the dimensions of
C     the netCDF variables. Both of the netCDF variables we are creating
C     share the same four dimensions. In Fortran, the unlimited
C     dimension must come last on the list of dimids.
      dimids(1) = lon_dimid
      dimids(2) = lat_dimid
      dimids(3) = rec_dimid
      dimids2(1) = lon_dimid
      dimids2(2) = lat_dimid

C     Define the netCDF variables for the data.
      retval = nf_def_var(ncid, 'HB', NF_DOUBLE, NDIMS-1, dimids2, 
     +     hb_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'H', NF_DOUBLE, NDIMS, dimids, 
     +     h_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'U', NF_DOUBLE, NDIMS, dimids, 
     +     u_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'V', NF_DOUBLE, NDIMS, dimids, 
     +     v_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'ROT', NF_DOUBLE, NDIMS, dimids, 
     +     z_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, 'DIV', NF_DOUBLE, NDIMS, dimids, 
     +     d_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to the netCDF variables.
      retval = nf_put_att_text(ncid, hb_varid, 'units', 1, 
     +     'm')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, h_varid, 'units', 1, 
     +     'm')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, u_varid, 'units', 3, 
     +     'm/s')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, v_varid, 'units', 3, 
     +     'm/s')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, z_varid, 'units', 3, 
     +     '1/s')
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, d_varid, 'units', 3, 
     +     '1/s')
      if (retval .ne. nf_noerr) call handle_err(retval)

C     End define mode.
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Write the coordinate variable data. This will put the latitudes
C     and longitudes of our data grid into the netCDF file.
      retval = nf_put_var_double(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)
C constants
      retval = nf_put_var_double(ncid, g_varid, g)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid, f_varid, f)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid, A_varid, A)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid, H0_varid, H0)
      if (retval .ne. nf_noerr) call handle_err(retval)

      count(1) = NLONS
      count(2) = NLATS
      start(1) = 1
      start(2) = 1
      retval = nf_inq_varid(ncid, 'HB', hb_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid, hb_varid, hb)
      if (retval .ne. nf_noerr) call handle_err(retval)

C      retval = nf_put_var1_double(ncid, rec_varid, 0.d0)
C      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      write(*,*)
      print *,'*** SUCCESS CREATING FILE ', FILE_NAME, '!'
      write(*,*)

      end

      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end
