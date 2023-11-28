C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.

      subroutine ncreadhist(H,U,V,Z,D,t,timestep,N,M,FILE_NAME,recno)
      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      integer ncid

      integer NDIMS
      parameter (NDIMS = 3)
      integer M, N

C     The start and count arrays will tell the netCDF library where to
C     write our data.
      integer start(NDIMS), count(NDIMS)

      integer h_varid
      integer u_varid
      integer v_varid
      integer z_varid
      integer d_varid
      integer t_varid

C     Program variables to hold the data we will write out. We will only
C     need enough space to hold one timestep of data; one record.
      double precision H(N, M)
      double precision U(N, M)
      double precision V(N, M)
      double precision Z(N, M)
      double precision D(N, M)
      double precision t

C     
      integer recno,timestep

C     Error handling.
      integer retval


C     open the file. 
      retval = nf_open(FILE_NAME, nf_write, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     These settings tell netcdf to write one timestep of data. (The
C     setting of start(3) from the time step input tells netCDF which
C     timestep to write.)
      count(1) = N
      count(2) = M
      count(3) = 1
      start(1) = 1
      start(2) = 1
      start(3) = recno

         retval = nf_inq_varid(ncid, 'H', h_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_vara_double(ncid, h_varid, start, count, 
     +        H)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_inq_varid(ncid, 'U', u_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_vara_double(ncid, u_varid, start, count, 
     +        U)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_inq_varid(ncid, 'V', v_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_vara_double(ncid, v_varid, start, count, 
     +        V)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_inq_varid(ncid, 'ROT', z_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_vara_double(ncid, z_varid, start, count, 
     +        Z)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_inq_varid(ncid, 'DIV', d_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_vara_double(ncid, d_varid, start, count, 
     +        D)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_inq_varid(ncid, 'TIME', t_varid)
         if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_get_var1_double(ncid, t_varid, recno,
     +        t)
         if (retval .ne. nf_noerr) call handle_err(retval)


C     Close the file. This causes netCDF to flush all buffers and make
C     sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      write(*,*)   
      print *,'*** INITIAL DATA READ FROM    '
     +				, FILE_NAME
      write(*,*)   

      end

