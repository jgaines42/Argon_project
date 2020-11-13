      !-------------------------------------------------------------
      ! function kinetic_energy: Calculate kinetic energy based on velocities
      !
      ! Input:
      !   natoms: number of atoms in the system
      !   vel1: (natoms,3) velocities of all atoms
      !   mass: Mass of 1 atom
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
      ! function shift_vel_com: Shift all velocites to remove global shift
      !
      ! Input:
      !   natom: number of atoms in the sytem
      !   vel: (natom,3) velocity of each atom
      !
      ! Modifies:
      !   vel: updated with new velocities
      !   
      ! Output:
      !   amount_moved: the value of the change in v_x 
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

        ! Loop over all atoms and find center of mass of velocties
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

        amount_moved = vel_sum(1)
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
      PARAMETER (massAr=39.95/1000*1.6747E-24) ! Taken directly from the paper

      INTEGER nstep                     ! number of steps in the simulation
      PARAMETER (nstep=300000)

      INTEGER nsave                     ! frequency to save data
      PARAMETER (nsave=10)

      INTEGER comShiftTime              ! frequency to shift com of velocity
      PARAMETER (comShiftTime=101)

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

      INTEGER nDOF                      ! DOF per atom
      PARAMETER (nDOF=3)

      REAL*8 KE_Temp                    ! conversion between KE and Temp
      PARAMETER (KE_Temp = 2.0/(REAL(nDOF)*BOLTZ*natom))

      !-------------------------------------------------------------
      ! Declare variables
      !-------------------------------------------------------------
      INTEGER I,J,K,L,M,N,time_loop

      REAL*8 Length                     ! Length of the box

      ! Create arrays for coordinates and velocities
      REAL*8 pos(natom,3),vel(natom,3)
      REAL*8 vel_old(natom,3),vel_half(natom,3)

      ! Set to NVE or NVT, 0 if NVE, 1 for NVT
      INTEGER is_NVT

      ! System variables
      REAL*8 time                      ! the step multiplied by DT to give the time in sec
      REAL*8 temp                      ! the temperature of the system to be calc from KE
      REAL*8 KE                        ! Kinetic energy in J
      REAL*8 V_tot                     ! The total potential energy in the system
      INTEGER start_NVE                ! How many loops before starting NVE

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

      ! Variables for g(r)
      INTEGER nGrBins,ibin
      PARAMETER (nGrBins=500)             ! number of g(r) bins
      REAL*8 bin_ends(nGrBins)            ! the edge value of each bin for g(r)
      REAL*8 g_of_r(nGrBins)              ! g(r) data
      REAL*8 gr_bin_W,dmax                ! width of each bin (A)
      PARAMETER (gr_bin_W=0.05E-10,dmax=REAL(nGrBins)*gr_bin_W)

      ! Variables for CVV
      INTEGER CVV_size                    ! Number of delta T time points
      PARAMETER (CVV_size=300)
      REAL*8 vel_store(CVV_size,natom,3)  ! Store last N velocities
      REAL*8 CVV_data(CVV_size)           ! CVV values
      INTEGER CVV_I1, CVV_I2,CVV_loop     ! Looping variables
      INTEGER CVV_count                   ! Keep track of how many CVV

      ! Variables for velocity distribution
      REAL*8 vel_min                      ! left most bin
      PARAMETER (vel_min=0.0)
      REAL*8 v_bin_width                  ! width of velocity bins
      PARAMETER (v_bin_width = 1.0)
      REAL*8 this_vel                     ! speed of this atom
      INTEGER vnbin                       ! number of velocity bins
      PARAMETER (vnbin=600)
      INTEGER vel_bin                     ! Loop variable
      INTEGER vel_hist(vnbin)             ! histogram data

      ! Variables for mean squared discplcement
      INTEGER num_mds_steps               ! Number of time steps to use
      PARAMETER (num_mds_steps=400)
      INTEGER MSD_I1, MSD_I2, MSD_loop    ! index and loop variables
      REAL*8 msd(num_mds_steps)           ! msd data
      REAL*8 p(num_mds_steps,natom,3)     ! store last N positions
      Integer MSD_Count                   ! Keep track of how many MSD

      ! reading in/values from input file
      INTEGER idum ! loop/array variables / dummy integer
      CHARACTER*5 resname,atomname

      ! File for storing energy data
      OPEN(91,FILE="Ar_positions.gro")         ! output traj file
      OPEN(92,FILE="Energies_Ar.txt")          ! Energy 
      OPEN(94,FILE="Temp_Ar.txt")              ! Temperature
      OPEN(95,FILE="Velocity_Ar.txt")          ! CVV data
      OPEN(96,FILE="Mean_squared_displacement_Ar.txt")
      OPEN(97,FILE="G_r_Ar.txt")               ! g(r) data
      OPEN(98,FILE="Vel_distribution_Ar.txt")  ! Speed histogram


      !-------------------------------------------------------------
      ! Set initial system variables
      !-------------------------------------------------------------
      is_NVT=1
      KE = 0.0
      V_tot = 0.0
      start_NVE = 20000                       ! Start in NVT to check equilibration
      time = 0.0
      CVV_count = 0
      MSD_Count = 0

      !-------------------------------------------------------------
      ! Read in initial file
      !-------------------------------------------------------------
      OPEN(11,FILE="Argon_Equ.gro") ! open input file
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
            vel(I,K)=vel(I,K)*1.0E2         ! to m/s
            accel(I,K) = 0.0                ! initialize acceleration to 0
            force(I,K) = 0.0                ! initialize forces to 0
         END DO
      END DO

      KE = kinetic_energy(natom, vel,massAr)
      temp = (KE*KE_Temp)
      print*,temp
      print*,massAr
      print*,Length

      ! Initialize g_of_r data
      DO I=1,nGrBins
         g_of_r(i) = 0.0
         bin_ends(i) = (gr_bin_W*I)
       END DO
       
       ! Initialize CVV data
       DO I=1,CVV_size
          CVV_data(I)=0.0
       end DO

       ! Initialize msd data
       DO I=1,num_mds_steps
          msd(I)=0.0
       end DO

      !--------------------------------------------
      ! Progress system in time
      !--------------------------------------------
      DO time_loop=1,nstep
        time = time+dt                       ! Calculate current time

        ! Change to NVE after start_NVE
        if (time_loop > start_NVE) then
          is_NVT=0
        end if

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
              rij(K) = pos(I,K) - pos(J,K)
              rij(K) = rij(K) - Length*ANINT(rij(K)/Length)
              dist_ij = dist_ij + rij(K)*rij(K)
            END DO ! end K: calculate dist_ij

            ! If within a certain distance, calculate force
            if (dist_ij <= cutoffSQ) then
              sr_6 = (sigma2/dist_ij)**3
              sr_12 = sr_6**2

              ! Calculate potential energy and add to total of system
              V_ij = eps4*(sr_12 - sr_6)
              V_tot = V_tot + V_ij

              ! Calculate magnitude of force
              F_ij = -eps24*(-2.0*sr_12 + sr_6)/dist_ij

              ! Add force vector to i and j
              DO K=1,3
                 this_force(K) = F_ij*rij(K)
                 force(i,K) = force(i,K) + this_force(K)
                 force(j,K) = force(j,K) - this_force(K)
              END DO ! end K: force vectors
            END IF  ! end if dist_ij < cutoffSQ

            IF (is_NVT == 0) THEN

              ! Add to g(r) data if in save time
              IF(MOD(time_loop,nsave).EQ.0) THEN
                dist_ij = sqrt(dist_ij)

                ! Figure out what g(r) bin we are in
                ibin=FLOOR((dist_ij)/gr_bin_W)+1
                IF (ibin.LE.nGrBins) THEN
                  g_of_r(ibin)=g_of_r(ibin)+2 ! add 2 (i-j, j-i)
                END IF
              end IF ! if save_loop
            END IF ! end if NVE

          END do    ! end J
        END DO      ! end I


        !---------------------------------------------------------------------
        ! Use the force to update the positions
        ! This is done using leap frog algorithm
        ! Remember that vel is really vel at t-0.5dt
        !----------------------------------------------------------------------------------------------

        DO I=1,natom
           DO K=1,3
             accel(I,K) = force(I,K)/massAr                  ! Acceleration based on force
             vel_old(I,K) = vel(I,K)                        ! Store old velocity
             vel(I,K) = vel(I,K) + DT*accel(I,K)              ! Get new velocity (at t+0.tdt)
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
        ! Calculate velocity autocorrelation
        !---------------------------------------------------------------------
        IF (is_NVT == 0) THEN
          CVV_I1 = mod(time_loop-start_NVE-1,CVV_size)+1
          ! Store velocity data
          DO I=1,natom
            DO K=1,3
              vel_store(CVV_I1,I,K)=vel(I,K)
            end DO
          END DO

          ! if we have a full set of data, calculate CVV
          if (time_loop-start_NVE > CVV_size) then
            CVV_count = CVV_count + 1
            CVV_I1 = mod(time_loop-start_NVE,CVV_size)+1
            DO I=1,natom
              DO K=1,3
                DO CVV_loop = 1,CVV_size
                  CVV_I2 = mod(CVV_loop-1+time_loop-start_NVE, 
     :CVV_size)+1
                  CVV_data(CVV_loop)=CVV_data(CVV_loop)
     :+ vel_store(CVV_I1,I,K)*vel_store(CVV_I2,I,K)
                END DO
              END DO
            END DO 
          END IF

          !---------------------------------------------------------------------
          ! Calculate msd
          !---------------------------------------------------------------------
          MSD_I1=MOD(time_loop-start_NVE-1,num_mds_steps)+1
          ! Store positions
          DO I=1,natom
              DO K=1,3
                 p(MSD_I1,I,K)=pos(I,K)
              END DO
           END DO

           ! if we have enough positions stored, calculate values
           IF (time_loop-start_NVE > num_mds_steps) THEN
             MSD_Count= MSD_Count + 1
             MSD_I1=MOD(time_loop-start_NVE,num_mds_steps)+1
             DO I=1,natom
                DO K=1,3
                    DO MSD_loop=1,num_mds_steps
                       MSD_I2=MOD(MSD_loop-1+time_loop
     :-start_NVE,num_mds_steps)+1
                       dist_ij=(p(MSD_I2,I,K)-p(MSD_I1,I,K))**2
                       msd(MSD_loop)=msd(MSD_loop)+ dist_ij
                    END DO
                 END DO
              END DO
           END IF
        END IF ! end IF (is_NVT == 0)


        !---------------------------------------------------------------------
        ! If we're in NVT, scale velocity to keep set temperature
        !---------------------------------------------------------------------
        IF (is_NVT == 1) THEN
            
            KE = kinetic_energy(natom, vel_half,massAr)
            temp = (KE*KE_Temp)
            is_NVT_scale = sqrt(Tref/temp)  ! Scale for velocity

            ! Loop over all atoms and scale velocity
            DO I=1,natom
              DO K=1,3
                vel(I,K)=vel(I,K)*is_NVT_scale
              end do
            end do
        
          !---------------------------------------------------------------------
          ! Every set number of time points, shift com of veloctiy to be 0
          !---------------------------------------------------------------------
          IF(MOD(time_loop,comShiftTime).EQ.0) THEN
            vel_sum(1) = 0.0
            vel_sum(2) = 0.0
            vel_sum(3) = 0.0

            ! Loop over all atoms and find center of mass of velocties
            DO I=1,natom
               DO K=1,3                           ! Loop over x,y,z components
                  vel_sum(K)=vel_sum(K)+vel(I,K)  ! Calculate total velocity of system
               END DO
            END DO

            ! Calculate how much system needs to shift to keep from flying ice cube
            vel_sum(1) = vel_sum(1)/natom
            vel_sum(2) = vel_sum(2)/natom
            vel_sum(3) = vel_sum(3)/natom
            print*,vel_sum(1),vel_sum(2),vel_sum(3)
            ! shift velocity to keep whole system from moving in a direction
            DO I=1,natom                          ! Loop over all atoms
               DO K=1,3                           ! Loop over x,y,z components
                 vel(I,K) = vel(I,K)-vel_sum(K) ! shift velocity
               END DO
            END DO
          END IF
        end if

        !---------------------------------------------------------------------
        ! Store system data if we're in a save time
        !---------------------------------------------------------------------
        IF(MOD(time_loop,nsave).EQ.0) THEN
           IF(is_NVT ==0) THEN
              
              ! Add data to speed distribution
              DO I=1,natom
                this_vel = 0.0
                DO K=1,3
                  this_vel = this_vel+vel_half(I,K)*vel_half(I,K)
                end do
                this_vel = sqrt(this_vel)
                vel_bin=(FLOOR((this_vel-vel_min)/v_bin_width))+1
                IF (vel_bin.LE.vnbin) THEN
                  vel_hist(vel_bin)=vel_hist(vel_bin)+1
                END IF
              end do
          

C              ! Write coordinates to gromac file
C              WRITE(91,*)'After step ',time_loop
C              WRITE(91,*)natom
C              DO I=1,natom
C                WRITE(91,31)I,resname,atomname,I,
C      :(pos(I,K)*1.0E9,K=1,3),
C      :(vel(I,K)*1.0E-3,K=1,3)
C               END DO
C               WRITE(91,*)Length*1.0E9,Length*1.0E9,Length*1.0E9
             END IF
        END IF


        ! Write energies to file
        WRITE(92,*)time,V_tot,KE,V_tot+KE,temp
      END DO ! End of time loop

      !---------------------------------------------------------------------
      ! After time is done, write data to files
      !---------------------------------------------------------------------
      ! Write g_of_r data
      DO I=1,nGrBins
        WRITE(97,*)bin_ends(I),g_of_r(I)
      END DO

      ! Write velocity distribution
      DO I=1,vnbin
        WRITE(98,*)vel_hist(I)
      END DO

      ! Write CVV data
      DO I=1,CVV_size
        WRITE(95,*)CVV_data(I)
      END DO
      WRITE(95,*)CVV_count

      ! Write mds data
      DO I=1,num_mds_steps
        WRITE(96,*)msd(I)
      END DO
      WRITE(96,*)MSD_Count

  31  FORMAT(i5,2a5,i5,3f8.4,3f8.4) ! format of input file

      CLOSE(92) ! Energy
      CLOSE(91) ! positions
      CLOSE(94) ! Temp
      CLOSE(95) ! CVV
      CLOSE(96) ! MSD
      CLOSE(97) ! g(r)
      CLOSE(98) ! velocity distribution
      END
