C     Local grid writer for NCAR STSWM: time in hours, winds in m/s,
C     height in metres and Gaussian latitude in radians. Float64 fields
C     follow two int32 dimensions; longitude varies fastest.
      SUBROUTINE GRDOUT(U,V,H,TAU)
      INCLUDE 'params.i'
      REAL U(NLON+2,NLAT), V(NLON+2,NLAT), H(NLON+2,NLAT), TAU
      REAL GLAT
      EXTERNAL GLAT
      INTEGER I, J
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
      IF (FIRST) THEN
         OPEN(77,FILE='grid_out.bin',ACCESS='STREAM',FORM='UNFORMATTED',
     $        STATUS='REPLACE')
         WRITE(77) NLON, NLAT
         WRITE(77) (GLAT(J), J=1,NLAT)
         FIRST = .FALSE.
      ENDIF
      WRITE(77) TAU
      WRITE(77) ((U(I,J), I=1,NLON), J=1,NLAT)
      WRITE(77) ((V(I,J), I=1,NLON), J=1,NLAT)
      WRITE(77) ((H(I,J), I=1,NLON), J=1,NLAT)
      FLUSH(77)
      RETURN
      END
