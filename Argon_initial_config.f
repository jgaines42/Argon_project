      FUNCTION  random(iseed)
      ! When first call, iseed must be a large positive integer.
      ! iseed will be changed when exit and be used for next calling.
      ! The range of the generated random number is between 1 and -1
      implicit none
      integer, intent(inout) :: iseed
      real :: random
      iseed = mod(8121*iseed+28411, 134456) ! 0 =< iseed < 134456
      random = real(iseed)/134456. ! 0 < random < 1
      end FUNCTION random

      PROGRAM main
      IMPLICIT NONE

      INTEGER I,J,K,L,M,N

      ! Initiate constants
      REAL*8 NA
      PARAMETER (NA=6.02214076D23)
      REAL*8 PI
      PARAMETER (PI=4.0*atan(1.0))
      REAL*8 BOLTZ
      PARAMETER (BOLTZ=1.38064852D-23) ! in J/K

      ! Set initial parameters for our system
      INTEGER natom ! the number of atoms in the simulation
      PARAMETER (natom=864) ! From rahman64

      REAL*8 massAr ! mass of the particle
      PARAMETER (massAr=39.95/1000*1.6747E-24) ! in kg from kg/mol rahman64

      REAL*8 Tref !the reference temprerature for the thermostat
      PARAMETER (Tref=94.400000) ! in K from Rahman64

      REAL*8 Length ! Length of the box
      PARAMETER (Length=3.478) ! in nm from Raham64

      REAL vrms ! target average final_velocity
      PARAMETER (vrms = sqrt(3.0*BOLTZ*Tref/massAr)/1000) ! velocity in nm/ps

      REAL*8 sigma ! sigma value for calculating LJ potential
      PARAMETER (sigma=0.34) ! sigma in nm, from rahman64

      REAL*8 dist_part ! distance between particles in unit cell
      PARAMETER (dist_part=Length/12.0)!sigma/sqrt(2.0))

      REAL*8 dist_unit ! distance between unit cells
      PARAMETER (dist_unit=Length/6.0)

      INTEGER repeats
      PARAMETER (repeats = 6)

      ! Create unit cell
      REAL*8 unit_cell(4,3) ! sum of the total force on each particle


      REAL phi,theta,costheta ! for generating random vectors


      REAL scaled ! the value that the random vectors are scaled by

      REAL vx,vy,vz !velocity vectors in the x,y,z directions
      REAL rand_mag ! the magnitude of the x,y,z velocity vectors

      CHARACTER*5 resname,atomname
      PARAMETER (resname='   Ar',atomname='   Ar')

      INTEGER icount
      REAL x,y,z
      REAL dx,dy,dz
      REAL random
      INTEGER iseed

      ! Create unit cell
      DO I=1,4
         DO J=1,3
           unit_cell(I,J) = 0.0
         END DO
      END DO
      unit_cell(1,2) = dist_part
      unit_cell(1,3) = dist_part
      unit_cell(3,1) = dist_part
      unit_cell(3,3) = dist_part
      unit_cell(4,1) = dist_part
      unit_cell(4,2) = dist_part

      print*,massAr
      print*,Length
      print*,(massAr*natom/(Length*Length*Length))
      iseed=36

      OPEN(91,FILE="Argon_864_initial_coordinates.gro") !open output file

      WRITE(91,*)'A box of liquid Ar'
      WRITE(91,*)natom

      OPEN(92,FILE="Argon_864_initial_velocity.gro") !open output file
      WRITE(92,*)'A box of liquid Ar'
      WRITE(92,*)natom

      icount=1      ! stores which atom we are on
      ! Loop over all repeats of the unit cell
      DO I=1,repeats
         dx=REAL(I-1)*dist_unit
         DO J=1,repeats
            dy=REAL(J-1)*dist_unit
            DO K=1,repeats
              dz=REAL(K-1)*dist_unit

              ! Now that we know the displacements of the unit cell, move each atom
              DO M=1,4
                x=dx+unit_cell(M,1)
                y=dy+unit_cell(M,2)
                z=dz+unit_cell(M,3)

                ! Get random velocity
                phi = random(iseed)*PI*2.0
                costheta = -1.0 + random(iseed)*2.0

                theta = acos( costheta )
                vx = sin( theta) * cos( phi )
                vy = sin( theta) * sin( phi )
                vz = cos( theta )

                ! calculate the magnitude of the random x,y,z vectors
                rand_mag = sqrt((vx**2)+(vy**2)+(vz**2))

                ! calculate how much the random vector needs to be scaled
                scaled=vrms/rand_mag

                ! Scale the random vector
                vx=vx*vrms ! scale the x vector
                vy=vy*vrms ! scale the y vector
                vz=vz*vrms ! scale the z vector


                WRITE(91,31)icount,resname,atomname,icount,x,y,z,
     :vx,vy,vz
                WRITE(92,31)icount,resname,atomname,icount,10.0*vx,
     :10.0*vy,10.0*vz,vx,vy,vz
                icount=icount+1
C                IF(icount.GT.natom) THEN
C                  GO TO 51
C                END IF
              END DO
          END DO
        END DO
      END DO

  51  WRITE(91,*)Length,Length,Length
      WRITE(92,*)Length,Length,Length

  31  FORMAT(i5,2a5,i5,3f8.4,3f8.4)

      CLOSE(91)
      CLOSE(92)

      END
