      !-------------------------------------------------------------
      ! function kinetic_energy: Calculate kinetic energy based on velocities
      !
      ! Input:
      !   natoms: number of atoms in the system
      !   vel1: (natoms,3) velocities of all atoms
      !   massAr: Mass of 1 atom
      !
      ! Output:
      !   ke: Total kinetic energy in the system
      !-------------------------------------------------------------
      function kinetic_energy(natoms, vel1, mass) result(ke)
        INTEGER, intent(in):: natoms          ! number of atoms
        REAL*8, intent(in) :: vel1(natoms,3)  ! velocity array
        REAL*8, intent(in) :: mass            ! mass of 1  atom
        REAL*8             :: ke              ! output
        INTEGER I,K
        ke=0.0

        ! Loop over all atoms
        DO I=1,natoms
          ! Loop over x,y,z components of velocity
           DO K=1,3
              ke=ke+vel1(I,K)**2    ! sum vel^2
           END DO
        END DO
        ke = 0.5*ke*mass            ! calculate kinetic energy
      end function


      !-------------------------------------------------------------
      ! function shift_vel_com
      !-------------------------------------------------------------
      function shift_vel_com(natom, vel) result(amount_moved)  
        INTEGER, intent(in):: natom           ! number of atoms
        REAL*8             :: vel(natom,3)    ! velocity array
        REAL*8 vel_sum(3)                     ! value to shift by
        REAL*8 amount_moved                   ! vel_sum(1) 
        INTEGER I,K
        vel_sum(1) = 0.0
        vel_sum(2) = 0.0
        vel_sum(3) = 0.0

        DO I=1,natom
           DO K=1,3                           ! Loop over x,y,z components
              vel_sum(K)=vel_sum(K)+vel(I,K)  ! Calculate total velocity of system
           END DO
        END DO

        ! Calculate how much system needs to shift to keep from flying ice cube
        vel_sum(1) = vel_sum(1)/natom
        vel_sum(2) = vel_sum(2)/natom
        vel_sum(3) = vel_sum(3)/natom
   
        ! shift velocity to keep whole system from moving in a direction
        DO I=1,natom                          ! Loop over all atoms
           DO K=1,3                           ! Loop over x,y,z components
             vel(I,K) = vel(I,K)-vel_sum(K) ! shift velocity
           END DO
        END DO
        amount_moved=vel_sum(1)
      end function

      !-------------------------------------------------------------
      ! Argon_simulation
      !
      ! Simulate a set of argon atoms
      ! Based on code from Stephanie Kearing
      !
      ! Input:
      !   Argon_864_initial_coordinates.gro : file with atomic coordinates and
      !       velocities for to initiate simulation
      !
      ! Ouput:
      !   A bunch of files
      !-------------------------------------------------------------

      PROGRAM Argon_simulation

      ! Based on code from Stephanie Kearing
      IMPLICIT NONE

      ! Declare functions
      REAL*8 kinetic_energy
      REAL*8 shift_vel_com

      !-------------------------------------------------------------
      ! Initiate constants
      !-------------------------------------------------------------
      REAL*8 NA
      PARAMETER (NA=6.02214076D23)

      REAL*8 BOLTZ
      PARAMETER (BOLTZ=1.38064852D-23)  ! in J/K

      INTEGER natom                     ! the number of atoms in the simulation
      PARAMETER (natom=864)

      REAL*8 eps                        ! epsilon value in J
      PARAMETER (eps=120.0*BOLTZ)

      REAL*8 sigma                      ! sum of radii, in m
      PARAMETER (sigma=0.34000000E-9)
      
      REAL*8 sigma2                     ! sigma squared, in m^2
      PARAMETER (sigma2=sigma**2)

      REAL*8 massAr                     ! mass of on Ar particle, in kg
      PARAMETER (massAr=39.95/1000*1.6747E-24)

      INTEGER nstep                     ! number of steps in the simulation
      PARAMETER (nstep=100000)

      INTEGER nsave                     ! frequency to save data
      PARAMETER (nsave=10)

      INTEGER comShiftTime              ! frequency to shift com of velocity
      PARAMETER (comShiftTime=150)

      REAL*8 DT                         ! time step, in seconds
      PARAMETER (DT=10.0E-15) 

      REAL*8 cutoff,cutoffSQ            ! cutoffs for energy equations
      PARAMETER (cutoff=2.25*sigma,cutoffSQ=cutoff**2)

      REAL*8 Tref                       ! the reference temprerature, in K
      PARAMETER (Tref=94.400000)

      REAL*8 eps4                       ! epsilon*4 for potential energy
      PARAMETER (eps4 = 4.0*eps)

      REAL*8 eps24                      ! epsilon*24 for force
      PARAMETER (eps24=24.0*eps)

      INTEGER nDOF                      ! DOF of each atom
      PARAMETER (nDOF = 3)

      REAL*8 KE_Temp                    ! conversion between KE and Temp
      PARAMETER (KE_Temp = 2.0/(REAL(nDOF)*BOLTZ*natom))

      !-------------------------------------------------------------
      ! Declare variables
      !-------------------------------------------------------------
      INTEGER I,J,K,L,M,N,time_loop

      REAL*8 Length                     ! Length of the box

      ! Create arrays for coordinates and velocities
      REAL*8 pos(natom,3),vel(natom,3) ! 3 for x,y,z directions
      REAL*8 vel_old(natom,3),vel_half(natom,3)

      ! Set to NVE or NVT, 0 if NVE, 1 for NVT
      INTEGER is_NVT

      ! System variables
      REAL*8 time                      ! the step multiplied by DT to give the time in sec
      REAL*8 temp                      ! the temperature of the system to be calc from KE
      REAL*8 KE                        ! Kinetic energy in J
      REAL*8 V_tot                     ! The total potential energy in the system

      ! preventing the flying ice cube
      REAL*8 vel_sum(3),vel_scale(3)

      ! calculating VLJ and FLJ in x,y,z directions
      REAL*8 rij(3)                    ! the vector between atoms I and J (I-J)
      REAL*8 dist_ij                   ! total distance between atoms I and J
      REAL*8 V_ij                      ! The LJ potential between atoms I and J
      REAL*8 F_ij                      ! the magnitude of force between atoms I and J
      REAL*8 sr_12                     ! (sigma/r)^12
      REAL*8 sr_6                      ! (sigma/r)^6
      REAL*8 force(natom,3)            ! the force on each particles

      ! Variables for leapfrog algorithm 
      REAL*8 this_force(3)            ! force on this particle
      REAL*8 accel(natom,3)           ! acceleration
      REAL*8 is_NVT_scale             ! value for scaling velocity in NVT


      ! reading in/values from input file
      INTEGER idum                    ! dummy integer
      CHARACTER*5 resname,atomname    ! names in input file

      ! File for storing energy data
      OPEN(92,FILE="Energies_Ar_NVT.txt")
      OPEN(91,FILE="Ar_positions_NVT.gro") ! open output traj file
      OPEN(94,FILE="Temp_NVT.txt") ! open output temperature
      OPEN(95,FILE="Velocity_NVT.txt")
      OPEN(96,FILE="Velocity_vector_NVT.txt")
      OPEN(97,FILE="G_r_all_distances_NVT.txt")
      OPEN(98,FILE="Argon_Equ.gro")

      !-------------------------------------------------------------
      ! Set initial system variables
      !-------------------------------------------------------------

      is_NVT=1
      KE = 0.0
      V_tot = 0.0
      time = 0.0

      !-------------------------------------------------------------
      ! Read in initial file
      !-------------------------------------------------------------
      OPEN(11,FILE="Argon_864_initial_coordinates.gro") ! open input file
      READ(11,*)
      READ(11,*)
      DO I=1,natom
         READ(11,31)idum,resname,atomname,idum,
     :(pos(I,K),K=1,3),(vel(I,K),K=1,3)
      END DO
      READ(11,*)Length,Length,Length
      CLOSE(11)

      !--------------------------------------------
      ! unit conversion for pos/vel/box read in from the gro file; to SI
      ! Also calculate average velocity and initialize some variables
      !--------------------------------------------
      Length=Length*1.0E-9 ! from nm to m

      ! Loop over all atoms
      DO I=1,natom
         DO K=1,3                           ! Loop over x,y,z components
            pos(I,K)=pos(I,K)*1.0E-9        ! from nm to m
            vel(I,K)=vel(I,K)*1.0E2         ! from 10*nm/ps to m/s
            accel(I,K) = 0.0                ! initialize acceleration to 0
            force(I,K) = 0.0                ! initialize forces to 0
         END DO
      END DO

      KE = kinetic_energy(natom, vel,massAr)
      temp = (KE*KE_Temp)
      print*,temp
      ! Prevent flying ice cube
      vel_sum(1) = shift_vel_com(natom, vel)
      print*,vel_sum(1)

      !--------------------------------------------
      ! Progress system in time
      !--------------------------------------------
      DO time_loop=1,nstep
        time = time+dt                       ! Calculate current time

        !--------------------------------------------
        ! Calculate forces and potential energy
        !--------------------------------------------

        V_tot = 0.0                     ! Initialize total velocity to 0

        ! loop over all pairs of atoms
        DO I=1,natom-1
          DO J=i+1, natom

            dist_ij = 0.0             ! reset dist_ij to 0

            ! Loop over x,y,z components
            DO K=1,3
              ! get distance vector between atoms
              rij(K)= pos(I,K)-pos(J,K)
              rij(K)= rij(K)-length*ANINT(rij(K)/Length)
              dist_ij = dist_ij + rij(K)*rij(K)
            END DO ! end K: calculate dist_ij

            ! If within a certain distance, calculate force
            if (dist_ij <= cutoffSQ) then
              sr_6 = (sigma2/dist_ij)**3
              sr_12 = sr_6**2

              ! Calculate potential energy and add to total of system
              V_ij=eps4*(sr_12-sr_6)
              V_tot=V_tot+V_ij

              ! Calculate magnitude of force
              F_ij=-eps24*(-2.0*sr_12+sr_6)/dist_ij

              ! Add force vector to i and j
              DO K=1,3
                 this_force(K) = F_ij*rij(K)
                 force(i,K) = force(i,K) + this_force(K)
                 force(j,K) = force(j,K) - this_force(K)
              END DO ! end K: force vectors

            END if  ! end if dist_ij < cutoffSQ
          END do    ! end J
        END DO      ! end I


        !---------------------------------------------------------------------
        ! Use the force to update the positions
        ! This is done using leap frog algorithm
        ! Remember that vel is really vel at t-0.5dt
        !----------------------------------------------------------------------------------------------

        DO I=1,natom
           DO K=1,3

             accel(I,K) =force(I,K)/massAr                  ! Acceleration based on force
             vel_old(I,K) = vel(I,K)                        ! Store old velocity
             vel(I,K)=vel(I,K) + DT*accel(I,K)              ! Get new velocity (at t+0.tdt)
             pos(I,K) = pos(I,K) + DT*vel(I,K)              ! Get new position (at t+dt)
             vel_half(I,K) = (vel(I,K) + vel_old(I,K))*0.5  ! Get velocity at (t+dt) for KE
             force(I,K) = 0.0                               ! reset force array for next time loop
           END DO

        END DO

        ! Calculate kinetic energy and temperature
        KE = kinetic_energy(natom, vel_half,massAr)
        temp = (KE*KE_Temp)
        WRITE(94,*)temp

        !---------------------------------------------------------------------
        ! If we're in NVT, scale velocity to keep set temperature
        !---------------------------------------------------------------------
        IF (is_NVT == 1) THEN
            KE = kinetic_energy(natom, vel,massAr)
            temp = (KE*KE_Temp)
            is_NVT_scale = sqrt(Tref/temp)  ! Scale for velocity

            ! Loop over all atoms and scale velocity
            DO I=1,natom
              DO K=1,3
                vel(I,K)=vel(I,K)*is_NVT_scale
              end do
            end do
        end if

        !---------------------------------------------------------------------
        ! Every 100 time points, shift com of veloctiy to be 0
        !---------------------------------------------------------------------
        IF(MOD(time_loop,comShiftTime).EQ.0) THEN
          vel_sum = shift_vel_com(natom, vel)
        END IF


        !---------------------------------------------------------------------
        ! Store system data if we're in a save time
        !---------------------------------------------------------------------
        IF(MOD(time_loop,nsave).EQ.0) THEN
           IF(is_NVT == 0) THEN
             WRITE(91,*)'After step ',time_loop
             WRITE(91,*)natom
             DO I=1,natom
               WRITE(91,31)I,resname,atomname,I,
     :(pos(I,K)*1.0E9,K=1,3),
     :(vel(I,K)*1.0E-3,K=1,3)
                WRITE(96,*)time_loop,vel(I,1),vel(I,2),vel(I,3)
              END DO
              WRITE(91,*)Length*1.0E9,Length*1.0E9,Length*1.0E9
            END IF
        END IF
        WRITE(92,*)time,V_tot,KE,V_tot+KE,temp
      END DO ! End of time loop

      ! Write final coordinates
      WRITE(98,*)'A box of liquid Ar'
      WRITE(98,*)natom
      DO I=1,natom
        WRITE(98,31)I,resname,atomname,I,
     :(pos(I,K)*1.0E9,K=1,3),
     :(vel(I,K)*1.0E-2,K=1,3)
      END DO
      WRITE(98,*)Length*1.0E9,Length*1.0E9,Length*1.0E9


  31  FORMAT(i5,2a5,i5,3f8.4,3f8.4) ! format of input file

      CLOSE(92)
      CLOSE(91)
      CLOSE(94)
      CLOSE(95)
      CLOSE(96)
      CLOSE(97)
      CLOSE(98)

      END
