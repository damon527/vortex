module Prtcl_System
  use MPI
  use m_Timer
  use m_TypeDef
  use Prtcl_Comm
  use Prtcl_Property
  use Prtcl_Geometry
  use Prtcl_decomp_2d
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_IOAndVisu
  use Prtcl_DumpPrtcl
  use Prtcl_Parameters
  use Prtcl_Integration
  use Prtcl_ContactSearch
  use Prtcl_ContactSearchPW
  use m_Decomp2d,only:nrank
  implicit none
  private
    
  !// DEMSystem class 
  type DEMSystem
    integer :: iterNumber   = 0  ! iteration number 
        
    !// timers
    type(timer):: m_total_timer
    type(timer):: m_pre_iter_timer
    type(timer):: m_comm_cs_timer
    type(timer):: m_CSCF_PP_timer
    type(timer):: m_CSCF_PW_timer
    type(timer):: m_Acceleration_timer
    type(timer):: m_integration_timer
    type(timer):: m_write_prtcl_timer
    type(timer):: m_comm_exchange_timer
  contains
    procedure:: Initialize => DEMS_Initialize
    
    ! iterating simulation for n time steps
    procedure:: iterate     => DEMS_iterate
    
    ! performing pre-iterations 
    procedure:: preIteration    => DEMS_preIteration
    
  end type DEMSystem
  type(DEMSystem),public::DEM
  
  integer::iCountACM
contains

!********************************************************************************
!   Initializing DEMSystem object with particles which are inserted from a 
!   predefined plane 
!********************************************************************************
  subroutine  DEMS_Initialize(this,chDEMPrm)
    implicit none
    class(DEMSystem)::this
    character(*),intent(in)::chDEMPrm
    
    ! locals
    integer::ierror
    character(256)::chStr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    !// Initializing main log info
    iCountACM=0
    if(DEM_Opt%RestartFlag) iCountACM=10
    this%IterNumber=DEM_Opt%ifirst-1
    write(chStr,"(A)") 'mkdir -p '//DEM_Opt%ResultsDir//' '//DEM_Opt%RestartDir//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call DEMLogInfo%InitLog(DEM_Opt%ResultsDir,DEM_Opt%RunName,DEM_Opt%LF_file_lvl,DEM_Opt%LF_cmdw_lvl)
    if(nrank==0) call DEMLogInfo%CreateFile(DEM_Opt%RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call DEMLogInfo%OpenFile()
    if(nrank==0) call Write_DEM_Opt_to_Log()

    ! Step1: Physical property
    call DEMProperty%InitPrtclProperty(chDEMPrm)
    call DEMProperty%InitWallProperty(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call DEMLogInfo%OutInfo("Physical properties contains "// trim( num2str(DEM_Opt%numPrtcl_Type ) ) // &
                              " particle types and "//trim( num2str(DEM_Opt%numWall_type ) )// " wall types.",2)
    endif

    ! Step2: set the geometry
    call DEMGeometry%MakeGeometry(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step2: Geometry is set", 1 )
      call DEMLogInfo%OutInfo("Geometry Contains "//trim(num2str(DEMGeometry%num_pWall))//" Plane walls.", 2)
    endif

    ! Step3: initilize all the particle variables
    call GPrtcl_list%AllocateAllVar()
    call DEM_IO%Init_visu(chDEMPrm,1)
    t_restart1=MPI_WTIME()
    if(.not.DEM_Opt%RestartFlag) then
      if(DEM_Opt%numPrtcl>0) call DEM_IO%ReadInitialCoord()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are MAKING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_Opt%numPrtcl)),2)
      endif
      DEM_Opt%np_InDomain = DEM_Opt%numPrtcl   
    else
      if(DEM_Opt%numPrtcl>0) then
        call DEM_IO%Read_Restart()
      else
        DEM_Opt%np_InDomain=0
      endif
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are READING from the Restarting file ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(DEM_Opt%np_InDomain)),2)
      endif
    endif
    if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1

    ! Step4: Initializing visu
    call DEM_IO%Init_visu(chDEMPrm,2)
    call Initialize_DumpPrtcl(chDEMPrm)

    ! Step5: initialize the inter-processors communication
    call DEM_Comm%InitComm()
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    endif

    ! Step6: Initializing contact list and contact force
    t_restart1=MPI_WTIME()
    call GPPW_CntctList%InitContactList()
    if(DEM_Opt%RestartFlag) then
      if(DEM_Opt%np_InDomain>0)call DEM_IO%RestartCL()
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1+t_res_tot
    if(nrank==0 .and. DEM_Opt%RestartFlag) call DEMLogInfo%OutInfo("Restart time [sec] :"//trim(num2str(t_res_tot)),2)
    
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step5: Initializing contact list and contact force models . . . ", 1 )
      if(DEM_Opt%CF_Type == ACM_LSD ) then
        write(chStr,"(A)") "Adaptive linear spring-dashpot model"
      elseif(DEM_Opt%CF_Type == ACM_nLin ) then
        write(chStr,"(A)") "Adaptive non-linear visco-elastic model"
      elseif(DEM_Opt%CF_Type == DEM_LSD ) then
        write(chStr,"(A)") "Typical linear spring-dashpot model"
      elseif(DEM_Opt%CF_Type == DEM_nLin ) then
        write(chStr,"(A)") "Typical non-linear visco-elastic model"
      endif
      call DEMLogInfo%OutInfo("Contact force model is "//trim(chStr), 2)

      if(DEM_Opt%PI_Method==PIM_FE) then
        write(chStr,"(A)") "Forward Euler             "
      elseif(DEM_Opt%PI_Method==PIM_AB2) then
        write(chStr,"(A)") "Adams Bashforth: 2nd Order"
      elseif(DEM_Opt%PI_Method==PIM_AB3) then
        write(chStr,"(A)") "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Linear   movement Integration scheme is : "//trim(chStr),2)

      if(DEM_Opt%PRI_Method==PIM_FE) then
        write(chStr,"(A)") "Forward Euler             "
      elseif(DEM_Opt%PRI_Method==PIM_AB2) then
        write(chStr,"(A)") "Adams Bashforth: 2nd Order"
      elseif(DEM_Opt%PRI_Method==PIM_AB3) then
        write(chStr,"(A)") "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Rotating movement Integration scheme is : "//trim(chStr),2)
    endif

    ! Step7: Initializing contact search method
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step6: Initializing contact search method . . . ", 1)
      call DEMLogInfo%OutInfo("Particle-Particle contact search intialization...",2)
      call DEMLogInfo%OutInfo("Particle-Wall contact search intialization...",2)
    endif
    call DEMContactSearch%InitContactSearch()
    call DEMContactSearchPW%InitContactSearchPW()
    
    ! Step8: timers for recording the execution time of different parts of program
    if(nrank==0) call DEMLogInfo%OutInfo("Step7: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_pre_iter_timer%reset()
    call this%m_comm_cs_timer%reset()
    call this%m_CSCF_PP_timer%reset()
    call this%m_CSCF_PW_timer%reset()
    call this%m_Acceleration_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_comm_exchange_timer%reset()
    call this%m_write_prtcl_timer%reset()
    
    call DEM_IO%dump_visu((DEM_Opt%ifirst-1)/icouple)            
  end subroutine DEMS_Initialize

  !********************************************************************************
  !   iterating over time 
  !   calls all the required methods to do numIter iterations in the DEM system
  !********************************************************************************
  subroutine DEMS_iterate(this,itime)
    implicit none
    class(DEMSystem) this
    integer,intent(in)::itime

    ! locals
    character(256)::chLine
    integer::Consv_Cont(2),Consv_Cont1(2),ierror,npwcs(4)

  IF(UpdateACMflag) THEN
    ! body
    call this%m_total_timer%start()

    ! pre-iteration adjustments 
    call this%m_pre_iter_timer%start()
    call this%preIteration()
    call this%m_pre_iter_timer%finish()

    ! inter-processor commucation for contact search ( ghost particle )
    call this%m_comm_cs_timer%start()
    call DEM_Comm%Comm_For_Cntct()
    call this%m_comm_cs_timer%finish()

    ! finding contacts among particels, and then calculating contact forces
    call this%m_CSCF_PP_timer%start()
    call DEMContactSearch%FindContacts()
    call this%m_CSCF_PP_timer%finish()

    ! finding contacts between particles and walls, and then calculating contact forces
    call this%m_CSCF_PW_timer%start()
    call DEMContactSearchPW%FindContactsPW()
    call this%m_CSCF_PW_timer%finish()

    ! correcting position and velocities 
    call this%m_integration_timer%start()
    iCountACM=iCountACM+1
    call Prtcl_Integrate(iCountACM)
    call this%m_integration_timer%finish()
   
    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call DEM_Comm%Comm_For_Exchange()
    call GPPW_CntctList%RemvReleased()
    call this%m_comm_exchange_timer%finish()
  ENDIF
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, DEM_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    if( mod(this%IterNumber,DEM_Opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime/icouple)
    if( mod(this%IterNumber,DumpPrtclFreq)== 0)      call WriteDumpCache(itime)
    if( mod(this%IterNumber,DEM_Opt%BackupFreq)== 0 .or. itime==DEM_Opt%ilast) then
      call DEM_IO%Write_Restart(itime)
      call DEM_IO%Delete_Prev_Restart(itime)
      call PrtclVarDump(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()
            
  IF(UpdateACMflag) THEN    
    ! output to log file and terminal/command window
    IF((this%IterNumber==DEM_Opt%ifirst .or. mod(this%IterNumber,DEM_Opt%Cmd_LFile_Freq)==0) ) THEN
      Consv_Cont1 =  DEMContactSearch%get_numContact()
      call MPI_REDUCE(Consv_Cont1, Consv_Cont,       2,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(GPPW_CntctList%numCntcts,npwcs,4,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
      if(nrank/=0) return
    
      ! command window and log file output
      call DEMLogInfo%OutInfo("DEM performed "//trim(num2str(this%IterNumber))//" iterations up to here!",1)
      call DEMLogInfo%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(this%m_total_timer%tot_time))//", "// &
      trim(num2str(this%m_total_timer%last_time ))//", "//trim(num2str(this%m_total_timer%average())),2)

      call DEMLogInfo%OutInfo("PreItertion time [tot, ave]        : "//trim(num2str(this%m_pre_iter_timer%tot_time))//", "// &
      trim(num2str(this%m_pre_iter_timer%average())),3)

      call DEMLogInfo%OutInfo("Comm_For_Contact [tot, ave]        : "//trim(num2str(this%m_comm_cs_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_cs_timer%average())), 3)

      call DEMLogInfo%OutInfo("CS and CF P-P time [tot, ave]      : "//trim(num2str(this%m_CSCF_PP_timer%tot_time))//", "// &
      trim(num2str(this%m_CSCF_PP_timer%average())), 3)

      call DEMLogInfo%OutInfo("CS and CF P-W time [tot, ave]      : "//trim(num2str(this%m_CSCF_PW_timer%tot_time))//", "// &
      trim(num2str(this%m_CSCF_PW_timer%average())), 3)
      call DEMLogInfo%OutInfo("Integration time [tot, ave]        : "//trim(num2str(this%m_integration_timer%tot_time))//", "// &
      trim(num2str(this%m_integration_timer%average())), 3)
      call DEMLogInfo%OutInfo("Comm_For_Exchange [tot, ave]       : "//trim(num2str(this%m_comm_exchange_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_exchange_timer%average())), 3)
      call DEMLogInfo%OutInfo("Write to file time [tot, ave]      : "//trim(num2str(this%m_write_prtcl_timer%tot_time))//", "// &
      trim(num2str(this%m_write_prtcl_timer%average())), 3)
      write(chLine,"(A)") "Particle number in  domain:  "//trim(num2str(DEM_Opt%np_InDomain))
      call DEMLogInfo%OutInfo(chLine, 2)
      call DEMLogInfo%OutInfo("Contact information", 2)
      write(chLine,"(A)") "No. consrvtv. contacts, same level | cross level: "// trim(num2str(Consv_Cont(1)))//" | "//trim(num2str(Consv_Cont(2)))
      call DEMLogInfo%OutInfo(chLine, 3)
      write(chLine,"(A)") "No. exact contacts P-P | P-GP | P-FP | P-W : "//trim(num2str(npwcs(1)))//" | "//trim(num2str(npwcs(2)))//" | "//trim(num2str(npwcs(3)))//" | "//trim(num2str(npwcs(4)))
      call DEMLogInfo%OutInfo(chLine, 3)
    ENDIF
  ENDIF
  end subroutine DEMS_iterate
    
  !**********************************************************************
  ! DEMS_preIteration
  !**********************************************************************
  subroutine DEMS_preIteration(this)
    implicit none
    class(DEMSystem)::this
        
    ! update wall neighbor list if necessary
    call DEMContactSearchPW%UpdateNearPrtclsPW(this%iterNumber)
    GPrtcl_cntctForce =zero_r3
    GPrtcl_torque= zero_r3
    GPrtcl_HighSt= "N"
    call GPPW_CntctList%PreIteration()
  end subroutine DEMS_preIteration

  !**********************************************************************
  ! Write_DEM_Opt_to_Log
  !**********************************************************************
  subroutine Write_DEM_Opt_to_Log()
    implicit none

    ! locals
    logical::RestartFlag
    real(RK)::dtDEM,Wall_neighbor_ratio,Prtcl_cs_ratio
    character(64):: RunName,ResultsDir,RestartDir,Geom_Dir
    integer::Wall_max_update_iter,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl
    integer::numPrtcl,numPrtclFix,CS_Method,CF_Type,PI_Method,PRI_Method,GeometrySource
    integer::numPrtcl_Type,numWall_type,CntctList_Size,CS_numlvls,Base_wall_id
    
    logical,dimension(3)::IsPeriodic
    type(real3)::gravity,minpoint,maxpoint
    integer::ifirstDEM,ilastDEM,BackupFreqDEM,SaveVisuDEM
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl,numPrtclFix,dtDEM,gravity,minpoint,maxpoint,CS_Method,CF_Type, &
                          PI_Method,PRI_Method,numPrtcl_Type,numWall_type,CS_numlvls,CntctList_Size,RunName,  &
                          Wall_max_update_iter,Wall_neighbor_ratio,ResultsDir,RestartDir,BackupFreqDEM,       &
                          SaveVisuDEM,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,GeometrySource,Geom_Dir,ifirstDEM, &
                          ilastDEM,Prtcl_cs_ratio,IsPeriodic
#ifdef CFDACM
    NAMELIST/CFDACMCoupling/UpdateACMflag,icouple,nForcingExtra,IBM_Scheme,Klub_pp,Klub_pw,Lub_ratio, &
                            Ndt_coll,IsDryColl,St_Crit,IsAddFluidPressureGradient
#endif

    RestartFlag = DEM_Opt%RestartFlag 
    numPrtcl    = DEM_Opt%numPrtcl   
    numPrtclFix = DEM_Opt%numPrtclFix  
    dtDEM       = DEM_Opt%dt       
    ifirstDEM   = DEM_Opt%ifirst   
    ilastDEM    = DEM_Opt%ilast    
    gravity     = DEM_Opt%gravity  
    minpoint    = DEM_Opt%SimDomain_min 
    maxpoint    = DEM_Opt%SimDomain_max 
    IsPeriodic  = DEM_Opt%IsPeriodic 
           
    Prtcl_cs_ratio=  DEM_Opt%Prtcl_cs_ratio 
    CS_Method  = DEM_Opt%CS_Method 
    CF_Type    = DEM_Opt%CF_Type
    PI_Method  = DEM_Opt%PI_Method  
    PRI_Method = DEM_Opt%PRI_Method
           
    numPrtcl_Type = DEM_Opt%numPrtcl_Type
    numWall_type  = DEM_Opt%numWall_type 
    
    CntctList_Size = DEM_Opt%CntctList_Size 
    CS_numlvls     = DEM_Opt%CS_numlvls    
           
    Base_wall_id  =  DEM_Opt%Base_wall_id 
    Wall_max_update_iter = DEM_Opt%Wall_max_update_iter
    Wall_neighbor_ratio  = DEM_Opt%Wall_neighbor_ratio 
           
    write(RunName,"(A)") DEM_Opt%RunName 
    write(ResultsDir,"(A)") DEM_Opt%ResultsDir 
    write(RestartDir,"(A)")DEM_Opt%RestartDir 
    BackupFreqDEM = DEM_Opt%BackupFreq 
    SaveVisuDEM   = DEM_Opt%SaveVisu 
    Cmd_LFile_Freq= DEM_Opt%Cmd_LFile_Freq 
    LF_file_lvl   = DEM_Opt%LF_file_lvl 
    LF_cmdw_lvl   = DEM_Opt%LF_cmdw_lvl 
           
    GeometrySource=  DEM_Opt%GeometrySource 
    write(Geom_Dir,"(A)") DEM_Opt%Geom_Dir 
    write(DEMLogInfo%nUnit, nml=DEMOptions)
#ifdef CFDACM
    write(DEMLogInfo%nUnit, nml=CFDACMCoupling)
#endif
  end subroutine Write_DEM_Opt_to_Log
end module Prtcl_System
module Prtcl_Integration
  use m_TypeDef
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Parameters
  use m_Parameters,only:gravity,PrGradData,IsUxConst
  implicit none
  private
  real(RK),parameter,dimension(2):: AB2C = [1.5_RK,-0.5_RK]
  real(RK),parameter,dimension(3):: AB3C = [23.0_RK,-16.0_RK,5.0_RK]/12.0_RK
  
  public::Prtcl_Integrate
#ifdef ObliqueWallTest
  logical::IsRotate=.false.
#endif
contains

  !******************************************************************
  ! Prtcl_Integrate
  !******************************************************************
  subroutine Prtcl_Integrate(iCountACM)
    implicit none
    integer,intent(in)::iCountACM
    
    ! locals
    integer::pid,itype
    type(real3)::linVel1,rotVel1,FpLinAcc,FpRotAcc,GravityToTal
    real(RK)::dt,dth,MassEff,InertiaEff,MassInFluid,TimeIntCoe(2)
#ifdef ObliqueWallTest
    real(RK)::rMagnitude,rxDir,ryDir
#endif
    
    if(iCountACM==1) then
      TimeIntCoe(1)=1.0_RK
      TimeIntCoe(2)=0.0_RK
    else
      TimeIntCoe(1)=AB2C(1)
      TimeIntCoe(2)=AB2C(2)
    endif
    GravityToTal=DEM_Opt%Gravity
    if(IsAddFluidPressureGradient .and. IsUxConst) GravityToTal%x= PrGradData(2)
    
    dt=DEM_opt%dt;   dth=dt*0.5_RK
    DO pid =1,GPrtcl_list%nlocal
      itype  = GPrtcl_pType(pid)
      InertiaEff = 1.0_RK/PrtclIBMProp(itype)%InertiaEff
      MassInFluid= PrtclIBMProp(itype)%MassinFluid
    
      ! linear velocity and position
      linVel1= GPrtcl_linVel(1,pid)
      if(GPrtcl_HighSt(pid)=="Y") then   ! Turn off fluid forces for large St collisions
        MassEff = 1.0_RK/DEMProperty%Prtcl_PureProp(itype)%Mass
        FpLinAcc= (MassInFluid*MassEff)*GravityToTal
      else
        MassEff = 1.0_RK/PrtclIBMProp(itype)%MassEff
        FpLinAcc= MassEff*GPrtcl_FpForce(pid)+ (MassInFluid*MassEff)*GravityToTal
      endif
      GPrtcl_linAcc(1,pid)= MassEff*GPrtcl_cntctForce(pid)
      GPrtcl_linVel(1,pid)= linVel1 + dt*(FpLinAcc+ TimeIntCoe(1)*GPrtcl_linAcc(1,pid)+ TimeIntCoe(2)*GPrtcl_linAcc(2,pid))
      GPrtcl_linVel(2,pid)= linVel1
      GPrtcl_linAcc(2,pid)= GPrtcl_linAcc(1,pid)

      ! rotate position
      rotVel1= GPrtcl_rotVel(1,pid)
      GPrtcl_rotAcc(1,pid) = InertiaEff*GPrtcl_torque(pid)
      FpRotAcc = InertiaEff*GPrtcl_FpTorque(pid)
      GPrtcl_rotVel(1,pid)=rotVel1 + dt*(FpRotAcc+ TimeIntCoe(1)*GPrtcl_rotAcc(1,pid)+ TimeIntCoe(2)*GPrtcl_rotAcc(2,pid))
      GPrtcl_rotVel(2,pid)=rotVel1
      GPrtcl_rotAcc(2,pid)=GPrtcl_rotAcc(1,pid)

#ifdef RotateOnly
      GPrtcl_linVel(1,pid)=zero_r3
#endif
#ifdef ObliqueWallTest
      GPrtcl_linVel(1,pid)%z=0.0_RK
      if(.not.IsRotate) then
        GPrtcl_rotVel(1,pid)=zero_r3
        rMagnitude= 1.0_RK/sqrt(DEM_Opt%Gravity%x*DEM_Opt%Gravity%x+DEM_Opt%Gravity%y*DEM_Opt%Gravity%y)
        rxDir= rMagnitude*DEM_Opt%Gravity%x
        ryDir= rMagnitude*DEM_Opt%Gravity%y
        rMagnitude=sqrt(GPrtcl_linVel(1,pid)%x*GPrtcl_linVel(1,pid)%x+GPrtcl_linVel(1,pid)%y*GPrtcl_linVel(1,pid)%y)
        GPrtcl_linVel(1,pid)%x= rMagnitude*rxDir
        GPrtcl_linVel(1,pid)%y= rMagnitude*ryDir
      endif
#endif
      GPrtcl_PosR(pid)=GPrtcl_PosR(pid)+dth*(linVel1+ GPrtcl_linVel(1,pid))
#ifdef ObliqueWallTest    
      if(GPrtcl_PosR(pid)%y<2.0_RK*GPrtcl_PosR(pid)%w) IsRotate=.true.
#endif
    ENDDO

  end subroutine Prtcl_Integrate

end module Prtcl_Integration
#define DEM_ncv_Allowed 1000
module Prtcl_CL_and_CF
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Property
  use Prtcl_Geometry
  use Prtcl_Variables
  use Prtcl_Parameters
#ifdef CFDACM
  use m_Parameters,only:dtMax
  use m_Decomp2d,only: nrank,nproc
#elif CFDDEM
  use m_Decomp2d,only: nrank,nproc
#else
  use Prtcl_decomp_2d,only: nrank,nproc
#endif
  implicit none
  private 
  real(RK),parameter::END_OF_PRTCL = 142857.428571_RK   ! particle IO end flag
    
  integer,dimension(:),allocatable:: Bucket
  integer,dimension(:),allocatable:: id_j
  integer,dimension(:),allocatable:: Next
  integer,dimension(:),allocatable:: CntctStatus
  type(real3),dimension(:),allocatable:: TanDelta
  real(RK),dimension(:),allocatable::VelRel_Init

  integer,dimension(:),allocatable:: id_i
  integer,dimension(:),allocatable:: Head_Cp  ! cp: counterpart
  integer,dimension(:),allocatable:: Next_Cp
    
  type ContactList
    integer:: mBucket
    integer:: numCntcts(4)
    integer:: max_numCntcts
    integer:: NextInsert
  contains
    procedure:: InitContactList => CL_InitContactList
    procedure:: reallocateCL    => CL_reallocateCL
    procedure:: AddContactPP    => CL_AddContactPP
    procedure:: AddContactPPG   => CL_AddContactPPG
    procedure:: AddContactPPFix => CL_AddContactPPFix
    procedure:: AddContactPW    => CL_AddContactPW
    procedure,private:: Find_Insert
    procedure:: PreIteration    => CL_PreIteration
    procedure:: RemvReleased    => CL_RemvReleased
    procedure:: copy            => CL_copy
    procedure:: IsCntct         => CL_IsCntct
    procedure:: getPrtcl_nlink  => CL_getPrtcl_nlink
    procedure:: Gather_Cntctlink=> CL_Gather_Cntctlink
    procedure:: Add_Cntctlink   => CL_Add_Cntctlink
    procedure:: printCL         => CL_printCL
  
    procedure:: Get_numCntcts   => CL_Get_numCntcts
    procedure:: Prepare_Restart => CL_Prepare_Restart
    procedure:: GetNextTanDel_Un=> CL_GetNextTanDel_Un
    procedure:: Add_RestartCntctlink=> CL_Add_RestartCntctlink 
    procedure:: Count_Cntctlink     => CL_Count_Cntctlink
    procedure:: Resemble_Cntctlink  => CL_Resemble_Cntctlink
#ifdef CFDACM
    procedure:: AddLubForcePP    => CL_AddLubForcePP
    procedure:: AddLubForcePPG   => CL_AddLubForcePPG
    procedure:: AddLubForcePPFix => CL_AddLubForcePPFix
    procedure:: AddLubForcePW    => CL_AddLubForcePW
#endif
  end type ContactList
  type(ContactList):: GPPW_CntctList 
    
  public::END_OF_PRTCL,GPPW_CntctList
contains
    
  !**********************************************************************
  ! Initializing the contact list for ContactList class
  !**********************************************************************
  subroutine CL_InitContactList(this)
    implicit none
    class(ContactList) :: this
    integer:: max_numCntcts,i,iErrSum
    integer:: iErr1,iErr2,iErr3,iErr4,iErr5,iErr6,iErr7,iErr8,iErr9
    
    this%numCntcts= 0
    max_numCntcts = DEM_opt%cntctList_Size * GPrtcl_list%mlocal

    this%Max_numCntcts = max_numCntcts
    this%mBucket = GPrtcl_list%mlocal
        
    ! initializing the linked list to stores id pairs
    allocate(Bucket(this%mBucket),    Stat=iErr1)
    allocate(id_j(max_numCntcts),     Stat=iErr2)
    allocate(CntctStatus(max_numCntcts), Stat=iErr3)
    allocate(Next(max_numCntcts),     Stat=iErr4)
    allocate(TanDelta(max_numCntcts), Stat=iErr5)
    allocate(VelRel_Init(max_numCntcts),Stat=iErr6)

    allocate(id_i(max_numCntcts),    Stat=iErr7)
    allocate(Head_Cp(this%mBucket),  Stat=iErr8)
    allocate(Next_cp(max_numCntcts), Stat=ierr9)
    iErrSum=abs(iErr1)+abs(iErr2)+abs(iErr3)+abs(iErr4)+abs(iErr5)+abs(iErr6)+abs(iErr7)+abs(iErr8)+abs(iErr9)
    if(iErrSum/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"CL_InitContactList","Allocation failed ")
        
    Bucket = 0
    do i = 1, max_numCntcts
      Next(i)=-i-1
      CntctStatus(i)=-1
    enddo
    Next(max_numCntcts) = 0
    this%NextInsert = 1
    Head_Cp = -1
    Next_cp = -1
  end subroutine CL_InitContactList

  !**********************************************************************
  ! reallocate contact list
  !**********************************************************************
  subroutine CL_reallocateCL(this,nCL_new)
    implicit none
    class(ContactList):: this
    integer,intent(in):: nCL_new

    ! locals
    integer::sizep,sizen,ierrTmp,ierror=0
    integer,dimension(:),allocatable:: IntVec

    sizep= this%mBucket
    sizen= GPrtcl_list%mlocal    ! 2021-09-08, Zheng Gong
    this%mBucket = sizen

    ! ======= integer vector part =======
    call move_alloc(Bucket, IntVec)
    allocate(Bucket(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    Bucket(1:sizep)=IntVec
    Bucket(sizeP+1:sizen)=0    ! Added at 11:24, 2020-09-11, Gong Zheng

    call move_alloc(Head_Cp,IntVec)
    allocate(Head_Cp(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    Head_Cp(1:sizep)=IntVec
    Head_Cp(sizep+1:sizen)=-1  ! Added at 11:24, 2020-09-11, Gong Zheng
    deallocate(IntVec)
    
    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," CL_reallocateCL"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif   
   ! here nothing is done for id_j, Next, CntctStatus, TanDelta, VelRel_Init
   ! considering that Max_numCntcts is big enough.
   ! If not, a fatal error will occur in CL_AddContactPP/CL_AddContactPW/CL_AddContactPPG/CL_AddContactPPFix
  end subroutine CL_reallocateCL

  !**********************************************************************
  ! Adding a contact pair to the contact list (particle & particle)
  !**********************************************************************
  subroutine CL_AddContactPP(this,id1,id2,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    gid2=GPrtcl_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(1) = this%numCntcts(1) + 1
      call ContactForce_PP(id1,id2,item2,ovrlp)

      id_i(item2)= GPrtcl_id(id1)
      Next_Cp(item2)=Head_Cp(id2)
      Head_Cp(id2)=item2

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPP","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPP","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPP

  !**********************************************************************
  ! Adding a contact pair to the contact list (particle & ghost particle) 
  !**********************************************************************
  subroutine CL_AddContactPPG(this,id1,id2,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    gid2=GhostP_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(2) = this%numCntcts(2) + 1
      call ContactForce_PPG(id1,id2,item2,ovrlp)

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPG","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPG","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPPG

  !**********************************************************************
  ! Adding a contact pair to the contact list (particle & fixed particle) 
  !**********************************************************************
  subroutine CL_AddContactPPFix(this,id1,id2,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    gid2=GPFix_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(3) = this%numCntcts(3) + 1
      call ContactForce_PPFix(id1,id2,item2,ovrlp)

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPFix","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPFix","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPPFix

  !**********************************************************************
  ! Adding a contact pair to the contact list 
  !**********************************************************************
  subroutine CL_AddContactPW(this,id1,id2,ovrlp,nv)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    integer::wid,item1,item2
    real(RK)::ovrlp
    type(real3)::nv

    wid= DEMGeometry%pWall(id2)%wall_id
    call this%Find_Insert(id1,wid,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(4) = this%numCntcts(4) + 1
      call ContactForce_PW(id1,id2,item2,ovrlp,nv)
    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPW","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPW","The container is full and there is no space for new item")      
    endif
  end subroutine CL_AddContactPW
    
  !***********************************************************************************************
  !* searching for an item, if it exists, it would return (/1, index of container that contain the item/),
  !     if it is new, pushing the item into the list and return (/2, index of container that
  !     contains the new item /), otherwise (-2,-2). 
  !********************************************************************************************
  subroutine Find_insert(this,bkt_id,value,item1,item2)
    implicit none
    class(ContactList) this
    integer,intent(in) :: bkt_id,value
    integer,intent(out):: item1, item2

    ! locals
    integer::n,NextI
        
    ! The bucket id is not in the range and must return (-1,-1), this is an error for the linked list
    if(bkt_id>GPrtcl_list%mlocal) then
      item1 =-1; item2 = -1
      return
    end if
        
     n = Bucket(bkt_id) 
     ! the item already exists in the list; returning the proper code and the index of container in which the item exists
    do while(n>0)
      if(id_j(n)==value) then
        item1 = 1; item2 = n
        return
      end if
      n = Next(n)
    enddo
        
    !the container is full and there is no more space for this item
    if(this%NextInsert==0) then
      item1 =-2; item2 = -2
      return
    endif
        
    ! the item is new and it should be pushed into the list, inserting new item in the list 
    NextI = this%NextInsert     
    this%NextInsert = -Next(NextI)
    id_j(NextI)=value
    Next(NextI)=Bucket(bkt_id)
    Bucket(bkt_id) = NextI
    TanDelta(NextI)= zero_r3
    VelRel_Init(NextI)=-2.0_RK
    item1 =2; item2 = NextI
  end subroutine Find_insert

  !**********************************************************************
  ! CL_getPrtcl_nlink
  !**********************************************************************
  function CL_getPrtcl_nlink(this,pid) result(res)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid

    ! locals
    integer::res,n

    res = 0
    n = Bucket(pid)
    do while(n>0)
      res = res + 1
      n = Next(n)
    enddo

    n = Head_cp(pid)
    do while(n .ne. -1)
      res = res + 1
      n = Next_Cp(n)
    enddo
  end function CL_getPrtcl_nlink

  !**********************************************************************
  ! CL_IsCntct
  !**********************************************************************
  function CL_IsCntct(this,pid) result(res)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid

    ! locals
    integer::res

    if(Bucket(pid)>0 .or. Head_cp(pid)>0) then
      res= 1
    else
      res= 0
    endif
  end function CL_IsCntct

  !**********************************************************************
  ! CL_Prepare_Restart
  !**********************************************************************
  subroutine CL_Prepare_Restart(this,nlink_ind)
    implicit none
    class(ContactList)::this
    integer,intent(in)::nlink_ind

    ! locals
    integer::i,CLid

    CLid = nlink_ind
    do i=1,this%Max_numCntcts
      if(CntctStatus(i)>0) then
        CLid = CLid + 1
        CntctStatus(i) = CLid 
      endif
    enddo
  end subroutine CL_Prepare_Restart

  !**********************************************************************
  ! CL_Get_numCntcts
  !**********************************************************************
  subroutine CL_Get_numCntcts(this,nCntct,nTanDel)
    implicit none
    class(ContactList)::this
    integer,intent(inout):: nCntct,nTanDel

    ! locals
    integer:: pid,n

    nCntct=0;nTanDel=0
    DO pid=1,GPrtcl_list%nlocal
      n=Bucket(pid)
      do while(n>0)
        nTanDel=nTanDel+1; n=Next(n)
      enddo
     
      n=Head_cp(pid)
      do while(n>0)
        nCntct=nCntct+1; n=Next_cp(n)
      enddo
    ENDDO
    nCntct=nCntct+nTanDel
  end subroutine CL_Get_numCntcts

  !**********************************************************************
  ! CL_GetNextTanDel_Un
  !**********************************************************************
  subroutine CL_GetNextTanDel_Un(this,TanDel_Un,prev,now)
    implicit none
    class(ContactList)::this 
    type(real4),intent(out)::TanDel_Un
    integer,intent(in)::prev
    integer,intent(out)::now

    ! locals
    integer::i
    
    do i=prev,this%Max_numCntcts
      if(CntctStatus(i)>0) then
        TanDel_Un   = TanDelta(i)
        TanDel_Un%w = VelRel_Init(i)
        now=i;exit
      endif
    enddo
  end subroutine CL_GetNextTanDel_Un

  !**********************************************************************
  ! CL_Gather_Cntctlink
  !**********************************************************************
  subroutine CL_Gather_Cntctlink(this,pid,buf_send,m)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    real(RK),dimension(*),intent(out)::buf_send
    integer,intent(inout)::m

    ! locals
    integer::n

    n = Bucket(pid)
    do while(n>0)
      CntctStatus(n)=3  ! 2021-08-25,added
      buf_send(m)=real(id_j(n));  m=m+1
      buf_send(m)=TanDelta(n)%x;  m=m+1
      buf_send(m)=TanDelta(n)%y;  m=m+1
      buf_send(m)=TanDelta(n)%z;  m=m+1
      buf_send(m)=VelRel_Init(n); m=m+1
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      buf_send(m)=real(id_i(n));  m=m+1
      buf_send(m)=TanDelta(n)%x;  m=m+1
      buf_send(m)=TanDelta(n)%y;  m=m+1
      buf_send(m)=TanDelta(n)%z;  m=m+1  
      buf_send(m)=VelRel_Init(n); m=m+1   
      n = Next_Cp(n)
    enddo
  end subroutine CL_Gather_Cntctlink

  !**********************************************************************
  ! CL_copy
  !**********************************************************************
  subroutine CL_copy(this,id1,id2)
    implicit none
    class(ContactList)::this
    integer,intent(in)::id1,id2

    bucket(id1) = bucket(id2)
    Head_Cp(id1)= Head_Cp(id2)
    bucket(id2) = 0
    Head_Cp(id2)=-1
  end subroutine CL_copy

  !**********************************************************************
  ! CL_Add_Cntctlink
  !**********************************************************************
  subroutine CL_Add_Cntctlink(this,pid,buf_recv,m)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    integer,intent(inout)::m
    real(RK),dimension(*),intent(in)::buf_recv

    ! locals
    integer::NextI
    real(RK)::realt

    do
      realt = buf_recv(m); m=m+1
      if(abs(realt-END_OF_PRTCL)<1.00E-10_RK) return

      NextI = this%NextInsert
      this%NextInsert = -Next(NextI)
      if(this%NextInsert==0) then
        call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Add_Cntctlink","The container is full and there is no space for new item") 
      endif
      id_j(NextI) = int(realt+0.2)
      Next(NextI) = Bucket(pid)
      Bucket(pid) = NextI
      CntctStatus(NextI)=2

      TanDelta(NextI)%x = buf_recv(m); m=m+1
      TanDelta(NextI)%y = buf_recv(m); m=m+1
      TanDelta(NextI)%z = buf_recv(m); m=m+1
      VelRel_Init(NextI)= buf_recv(m); m=m+1
    enddo
  end subroutine CL_Add_Cntctlink

  !**********************************************************************
  ! Removing all released contacts (those which are not in contact in this time step) 
  ! from contact list.
  !**********************************************************************  
  subroutine CL_RemvReleased(this)
    class(ContactList)::this
    
    ! locals
    integer::pid,n,prev,NextI
    
    do pid=1,GPrtcl_list%nlocal ! Modify this part by Zheng Gong, at 2020-09-09
      prev=0
      n=Bucket(pid)
      do while(n>0)
        NextI=Next(n)
        if(CntctStatus(n)==-2) then
          if(prev==0) then
            Bucket(pid)=NextI
          else
            Next(prev)=NextI
          endif
          CntctStatus(n)=-1
        
          Next(n)=-this%nextInsert
          this%nextInsert=n
        else
          prev=n
        endif
        n=NextI
      enddo

      n= Head_cp(pid)
      do while(n>0)
        if(CntctStatus(n)==3) then
          Next(n)=Bucket(pid)
          Bucket(pid)=n
          id_j(n)=id_i(n)
          CntctStatus(n)=2
        endif
        n=Next_cp(n)
      enddo
    enddo
    
    do pid=1,this%Max_numCntcts
      if(CntctStatus(pid)==3) then
        CntctStatus(pid)=-1
        Next(pid)=-this%nextInsert
        this%nextInsert=pid
      endif
    enddo

#ifdef DEBUG_RemvReleased
    ! CntctStatus(pid)==-2 will only appear when some particle escapes the domain
    do pid=1,this%Max_numCntcts
      if(CntctStatus(pid)==-2) then
         print*,nrank,pid,CntctStatus(pid),'############';stop     
      endif
    enddo 
#endif
  end subroutine CL_RemvReleased

  !**********************************************************************
  ! CL_PreIteration
  !**********************************************************************
  subroutine CL_PreIteration(this)
    class(ContactList) this

    ! locals
    integer::i

    this%numCntcts = 0
    do i=1,this%Max_numCntcts
       if(CntctStatus(i)>0) CntctStatus(i) = -2  ! flag previous contact
    enddo
    Head_Cp = -1
    Next_Cp = -1
  end subroutine CL_PreIteration

  !**********************************************************************
  ! CL_Add_RestartCntctlink
  !**********************************************************************
  subroutine CL_Add_RestartCntctlink(this,pid1,ncv,CntctVec,TanDel_Un)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid1,ncv
    integer,dimension(:),intent(in)::CntctVec
    type(real4),dimension(:),intent(in)::TanDel_Un

    ! locals
    real(RK)::dx,dy,dz,dist
    integer::i,j,nlocal,pid2,gid1,gid2,NextI
  
    nlocal=GPrtcl_list%nlocal
    gid1= GPrtcl_id(pid1)
    DO j=1,ncv
      pid2=0;  gid2=CntctVec(j)
      DO i=1,nlocal
        if(GPrtcl_id(i)==gid2) then
          pid2=i; exit
        endif
      ENDDO
      IF(pid2>0) THEN
        dx= GPrtcl_PosR(pid1)%x- GPrtcl_PosR(pid2)%x
        dy= GPrtcl_PosR(pid1)%y- GPrtcl_PosR(pid2)%y
        dz= GPrtcl_PosR(pid1)%z- GPrtcl_PosR(pid2)%z
        dist= sqrt(dx*dx+dy*dy+dz*dz)-1.5_RK*(GPrtcl_PosR(pid1)%w+GPrtcl_PosR(pid2)%w)
        if(gid1>gid2 .and. dist<=0.0_RK) cycle
      ENDIF

      NextI = this%NextInsert     
      this%NextInsert = -Next(NextI)
      if(this%NextInsert==0) then
        call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Add_RestartCntctlink","The container is full and there is no space for new item") 
      endif        
      id_j(NextI) = gid2 ! here gid2 can be a particle_within_this_processor/ghost_particle/wall
      Next(NextI) = Bucket(pid1)
      Bucket(pid1)= NextI
      CntctStatus(NextI)=2
      TanDelta(NextI)= TanDel_Un(j)
      VelRel_Init(NextI)= TanDel_Un(j)%w
      if(pid2>0 .and. gid1<gid2 .and. dist<=0.0_RK) then  ! gid2 is also within this processor
        id_i(NextI)= gid1
        Next_Cp(NextI)=Head_Cp(pid2)
        Head_Cp(pid2)=NextI
      endif
    ENDDO
  end subroutine CL_Add_RestartCntctlink

  !**********************************************************************
  ! CL_Count_Cntctlink
  !**********************************************************************
  subroutine CL_Count_Cntctlink(this,pid,ncv)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    integer,intent(out)::ncv

    ! locals
    integer::n
    
    ncv =0
    n = Bucket(pid)
    do while(n>0)
      ncv = ncv + 1
      ! for monosize particles, a particle can contact with NO MORE THAN 12 neighbor particles.
      if(ncv > DEM_ncv_Allowed) call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Count_Cntctlink","so big ncv")
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      ncv = ncv + 1
      if(ncv > DEM_ncv_Allowed) call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Count_Cntctlink","so big ncv") 
      n = Next_Cp(n)
    enddo
  end subroutine CL_Count_Cntctlink

  !**********************************************************************
  ! CL_Resemble_Cntctlink
  !**********************************************************************
  subroutine CL_Resemble_Cntctlink(this,pid,ncv,CntctVec)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    integer,intent(out)::ncv
    integer,dimension(:),intent(out)::CntctVec

    ! locals
    integer::n
    
    ncv =0
    n = Bucket(pid)
    do while(n>0)
      ncv = ncv + 1
      CntctVec(2*ncv-1) =  id_j(n)
      CntctVec(2*ncv)   =  CntctStatus(n)
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      ncv = ncv + 1
      CntctVec(2*ncv-1) =  id_i(n)
      CntctVec(2*ncv)   =  CntctStatus(n)  
      n = Next_Cp(n)
    enddo
  end subroutine CL_Resemble_Cntctlink

  !**********************************************************************
  ! CL_printCL
  !********************************************************************** 
  subroutine CL_printCL(this,itime)
    implicit none
    class(ContactList)::this
    integer,intent(in)::itime
   
    ! locals
    character(len=128)::filename
    integer::i,n,pid,ierror,nUnit
    
    write(filename,'(A,I10.10,A)')"cntctlist",itime,".txt"
    if(nrank==0) then  
      open(newunit=nUnit,file=filename,status='replace',form='formatted')
      close(nUnit)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    do i=0,nproc-1
      if(nrank==i) then
        open(newunit=nUnit,file=filename,status='old',position='append',form='formatted')
        DO pid=1,GPrtcl_list%nlocal
          n = Bucket(pid)
          do while(n>0)
            write(nUnit,'(3I10,10ES24.15)')nrank,GPrtcl_id(pid),id_j(n),TanDelta(n),VelRel_Init(n)
            n = Next(n)
          enddo
     
          n = Head_cp(pid)
          do while(n>0)
            write(nUnit,'(3I10,10ES24.15)')nrank,GPrtcl_id(pid),id_i(n),TanDelta(n),VelRel_Init(n)
            n = Next_cp(n)
          enddo
        ENDDO
        close(nUnit)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    enddo
  end subroutine CL_printCL
 
  !**********************************************************************
  ! ContactForce_PP
  !**********************************************************************
  subroutine ContactForce_PP(pid,pjd,ind,ovrlp)
    implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    ! locals
    type(real4)::Posi,Posj
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Norm_v,Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif

    Prop_ij=DEMProperty%Prtcl_BnryProp(GPrtcl_pType(pid), GPrtcl_pType(pjd))
    Veli= GPrtcl_linVel(1,pid)
    Velj= GPrtcl_linVel(1,pjd)
    Rvei= GPrtcl_rotVel(1,pid)
    Rvej= GPrtcl_rotVel(1,pjd)
    Posi= GPrtcl_PosR(pid)
    Posj= GPrtcl_PosR(pjd)
    ri= Posi%w; rj= Posj%w
    Norm_v= Posj.nv.Posi  ! Normal vector, Posj-Posi

#define ContactForce_PP
#ifdef CFDACM
#include "ACM_ContactForce_inc.f90"
#else
#include "Prtcl_ContactForce_inc.f90"
#endif
#undef  ContactForce_PP
  end subroutine ContactForce_PP

  !**********************************************************************
  ! ContactForce_PPG
  !**********************************************************************
  subroutine ContactForce_PPG(pid,pjd,ind,ovrlp)
    implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    ! locals
    type(real4)::Posi,Posj
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Norm_v,Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif

    Prop_ij=DEMProperty%Prtcl_BnryProp(GPrtcl_pType(pid), GhostP_pType(pjd))
    if(GPrtcl_id(pid)<GhostP_id(pjd)) then
      Veli= GPrtcl_linVel(1,pid)
      Velj= GhostP_linVel(pjd)
      Rvei= GPrtcl_rotVel(1,pid)
      Rvej= GhostP_rotVel(pjd)
      Posi= GPrtcl_PosR(pid)
      Posj= GhostP_PosR(pjd)
      ri= Posi%w; rj= Posj%w
      Norm_v= Posj.nv.Posi  ! Normal vector, Posj-Posi
#define ContactForce_PPG
#ifdef CFDACM
#include "ACM_ContactForce_inc.f90"
#else
#include "Prtcl_ContactForce_inc.f90"
#endif
#undef  ContactForce_PPG
    else
      Veli= GhostP_linVel(pjd)
      Velj= GPrtcl_linVel(1,pid)
      Rvei= GhostP_rotVel(pjd)
      Rvej= GPrtcl_rotVel(1,pid)
      Posi= GhostP_PosR(pjd)
      Posj= GPrtcl_PosR(pid)
      ri= Posi%w; rj= Posj%w
      Norm_v= Posj.nv.Posi  ! Normal vector, Posj-Posi
#define ContactForce_PGP
#ifdef CFDACM
#include "ACM_ContactForce_inc.f90"
#else
#include "Prtcl_ContactForce_inc.f90"
#endif
#undef  ContactForce_PGP
    endif
  end subroutine ContactForce_PPG

  !**********************************************************************
  ! ContactForce_PPFix
  !**********************************************************************
  subroutine ContactForce_PPFix(pid,pjd,ind,ovrlp)
    implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    ! locals
    type(real4)::Posi,Posj
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Norm_v,Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif

    Prop_ij=DEMProperty%Prtcl_BnryProp(GPrtcl_pType(pid), GPFix_pType(pjd))
    Veli= GPrtcl_linVel(1,pid)
    Velj= zero_r3
    Rvei= GPrtcl_rotVel(1,pid)
    Rvej= zero_r3
    Posi= GPrtcl_PosR(pid)
    Posj= GPFix_PosR(pjd)
    ri= Posi%w; rj= Posj%w
    Norm_v= Posj.nv.Posi  ! Normal vector, Posj-Posi

#define ContactForce_PPFix_W
#ifdef CFDACM
#include "ACM_ContactForce_inc.f90"
#else
#include "Prtcl_ContactForce_inc.f90"
#endif
#undef  ContactForce_PPFix_W
  end subroutine ContactForce_PPFix

  !**********************************************************************
  ! ContactForce_PW
  !**********************************************************************
  subroutine ContactForce_PW(pid,mwi,ind,ovrlp,Norm_v )
    implicit none
    integer,intent(in)::pid,mwi,ind
    real(RK),intent(in):: ovrlp
    type(real3),intent(inout):: Norm_v
    
    ! locals
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif
        
    Prop_ij = DEMProperty%PrtclWall_BnryProp(GPrtcl_pType(pid),DEMGeometry%pWall(mwi)%wall_Type)
    Veli= GPrtcl_linVel(1,pid)
    Velj= DEMGeometry%pWall(mwi)%trans_vel
    Rvei= GPrtcl_rotVel(1,pid)
    Rvej= zero_r3
    ri= GPrtcl_PosR(pid)%w
    rj= 1.00E20_RK*ri
    ! since the normal vector points from particle i to j we must negate the normal vector
    Norm_v = (-1.0_RK)*Norm_v
    
#define ContactForce_PPFix_W
#ifdef CFDACM
#include "ACM_ContactForce_inc.f90"
#else
#include "Prtcl_ContactForce_inc.f90"
#endif
#undef  ContactForce_PPFix_W
  end subroutine ContactForce_PW

#ifdef CFDACM
  !**********************************************************************
  ! Adding lubrication force to the "contact list" (particle & particle) 
  !**********************************************************************
  subroutine CL_AddLubForcePP(this,pid,pjd,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: pid,pjd
    real(RK),intent(in)::ovrlp
    
#ifdef DriftKissTumbleBreugem
    ! locals
    integer::pt_i,pt_j
    type(real3)::Norm_v,LubForce
    type(BinaryProperty)::Prop_ij
    real(RK)::Gravity_Norm,Mass,EpsValue,dlubDist,ratioD
    
    pt_i = GPrtcl_pType(pid)
    pt_j = GPrtcl_pType(pjd)
    dlubDist=dlub_pp(pt_i, pt_j)  
    Prop_ij=DEMProperty%Prtcl_BnryProp(pt_i, pt_j)
    
    EpsValue=1.0E-4
    Gravity_Norm=norm(DEM_Opt%Gravity)
    Mass=2.0_RK*Prop_ij%MassEff
    ratioD=ovrlp/dlubDist-1.0_RK
    ratioD=ratioD*ratioD
    
    ! normal vector, Posj-Posi
    Norm_v = (GPrtcl_PosR(pjd)) .nv. (GPrtcl_PosR(pid))  

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= (Mass*Gravity_Norm/EpsValue)*ratioD*Norm_v
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
    GPrtcl_cntctForce(pjd) = GPrtcl_cntctForce(pjd) + LubForce    
#else
    ! locals
    integer:: pt_i,pt_j
    type(real3)::Norm_v,Vij_n,LubForce

    ! normal vector, Posj-Posi
    pt_i = GPrtcl_pType(pid)
    pt_j = GPrtcl_pType(pjd)
    Norm_v = (GPrtcl_PosR(pjd)) .nv. (GPrtcl_PosR(pid))  

    ! normal relative velocity vectors
    Vij_n = ((GPrtcl_linVel(1,pid)-GPrtcl_linVel(1,pjd)) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pp(pt_i,pt_j)*Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
    GPrtcl_cntctForce(pjd) = GPrtcl_cntctForce(pjd) + LubForce
#endif
  end subroutine CL_AddLubForcePP

  !**********************************************************************
  ! Adding lub force to the "contact list" (particle & ghost particle) 
  !**********************************************************************
  subroutine CL_AddLubForcePPG(this,pid,gjd,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: pid,gjd
    real(RK),intent(in)::ovrlp

#ifdef DriftKissTumbleBreugem
    ! locals
    integer::pt_i,pt_j
    type(real3)::Norm_v,LubForce
    type(BinaryProperty)::Prop_ij
    real(RK)::Gravity_Norm,Mass,EpsValue,dlubDist,ratioD
    
    pt_i = GPrtcl_pType(pid)
    pt_j = GhostP_pType(gjd)
    dlubDist=dlub_pp(pt_i, pt_j)  
    Prop_ij=DEMProperty%Prtcl_BnryProp(pt_i, pt_j)
    
    EpsValue=1.0E-4
    Gravity_Norm=norm(DEM_Opt%Gravity)
    Mass=2.0_RK*Prop_ij%MassEff
    ratioD=ovrlp/dlubDist-1.0_RK
    ratioD=ratioD*ratioD
    
    ! normal vector, Posj-Posi
    Norm_v = (GhostP_PosR(gjd)) .nv. (GPrtcl_PosR(pid))  

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= (Mass*Gravity_Norm/EpsValue)*ratioD*Norm_v
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
#else
    ! locals
    integer:: pt_i,pt_j
    type(real3)::Norm_v,Vij_n,LubForce

    ! normal vector, Posj-Posi
    pt_i = GPrtcl_pType(pid)
    pt_j = GhostP_pType(gjd)
    Norm_v = (GhostP_PosR(gjd)) .nv. (GPrtcl_PosR(pid))  

    ! normal relative velocity vectors
    Vij_n = ((GPrtcl_linVel(1,pid)-GhostP_linVel(gjd)) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pp(pt_i,pt_j)*Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
#endif
  end subroutine CL_AddLubForcePPG

  !**********************************************************************
  ! Adding lub force to the "contact list" (particle & fixed particle) 
  !**********************************************************************
  subroutine CL_AddLubForcePPFix(this,pid,fid,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: pid,fid
    real(RK),intent(in)::ovrlp

    ! locals
    integer:: pt_i,pt_j
    type(real3)::Norm_v,Vij_n,LubForce

    ! normal vector, Posj-Posi
    pt_i = GPrtcl_pType(pid)
    pt_j = GPFix_pType(fid)
    Norm_v = (GPFix_PosR(fid) .nv. GPrtcl_PosR(pid))  

    ! normal relative velocity vectors
    Vij_n = (GPrtcl_linVel(1,pid) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pp(pt_i,pt_j)*Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) -LubForce
  end subroutine CL_AddLubForcePPFix

  !**********************************************************************
  ! Adding lubrication force to the "contact list" (particle & wall) 
  !**********************************************************************
  subroutine CL_AddLubForcePW(this,pid,mwi,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in)::pid,mwi
    real(RK),intent(in):: ovrlp
    
    ! locals
    type(real3):: Vij_n,LubForce,Norm_v
    
    ! normal and tangential velocities vectors 
    Norm_v= DEMGeometry%pWall(mwi)%n
    Vij_n = ((GPrtcl_linVel(1,pid)- DEMGeometry%pWall(mwi)%trans_vel) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pw(GPrtcl_pType(pid)) *Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) -LubForce
  end subroutine CL_AddLubForcePW
#endif
end module Prtcl_CL_and_CF
#ifdef DEM_ncv_Allowed
#undef DEM_ncv_Allowed
#endif
module Prtcl_Comm
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Property
  use Prtcl_Decomp_2d
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
  use Prtcl_ContactSearchPW
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private
#define xm_axis 1
#define xp_axis 2
#define ym_axis 3
#define yp_axis 4
#define zm_axis 5
#define zp_axis 6

  type(real3)::simLen
  logical,dimension(3)::pbc
  real(RK),dimension(6)::dx_pbc
  real(RK),dimension(6)::dy_pbc
  real(RK),dimension(6)::dz_pbc
  real(RK):: xst0_cs,xst1_cs,xst2_cs
  real(RK):: xed0_cs,xed1_cs,xed2_cs
  real(RK):: yst0_cs,yst1_cs,yst2_cs
  real(RK):: yed0_cs,yed1_cs,yed2_cs
  real(RK):: zst0_cs,zst1_cs,zst2_cs
  real(RK):: zed0_cs,zed1_cs,zed2_cs

  type Prtcl_Comm_info
    integer :: msend
    integer :: GhostCS_size
    integer :: Prtcl_Exchange_size
    real(RK):: LenForCS
  contains
    procedure:: InitComm            => PC_InitComm
    procedure:: Comm_For_Cntct      => PC_Comm_For_Cntct
    procedure:: pack_cntct          => PC_pack_cntct
    procedure:: unpack_cntct        => PC_unpack_cntct
    procedure:: Comm_For_Exchange   => PC_Comm_For_Exchange
    procedure:: pack_Exchange       => PC_pack_Exchange
    procedure:: unpack_Exchange     => PC_unpack_Exchange
    procedure:: ISInThisProc        => PC_ISInThisProc
    procedure:: reallocate_sendlist => PC_reallocate_sendlist
    procedure:: reallocate_ghost_for_Cntct => PC_reallocate_ghost_for_Cntct

    procedure,private:: Comm_For_Cntct_fixed  => PC_Comm_For_Cntct_fixed
  end type Prtcl_Comm_info
  type(Prtcl_Comm_info),public::DEM_Comm
#ifdef CFDACM
  integer,allocatable,dimension(:),public::sendlist
#else
  integer,allocatable,dimension(:)::sendlist
#endif
  integer,allocatable,dimension(:):: GhostPFix_id
  integer,allocatable,dimension(:):: GhostPFix_pType
  type(real4),allocatable,dimension(:):: GhostPFix_PosR

contains

  !**********************************************************************
  ! PC_InitComm
  !**********************************************************************
  subroutine PC_InitComm(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    real(RK)::lx,ly,lz,le_cs,vol1,vol2,vol3
    real(RK)::maxD,vol_tot,vol_ghost_cs,vol_sendlist
    integer::numPrtcl,mCS,msend,mFixedCS,ierrTmp,ierror=0
        
    pbc=DEM_Opt%IsPeriodic
    simLen = DEM_Opt%SimDomain_max-DEM_Opt%SimDomain_min

    maxD = 2.0_RK*maxval(DEMProperty%Prtcl_PureProp%Radius)
    lx=DEM_decomp%xEd-DEM_decomp%xSt
    ly=DEM_decomp%yEd-DEM_decomp%ySt
    lz=DEM_decomp%zEd-DEM_decomp%zSt
#ifdef CFDACM
    this%LenForCS= DEM_Opt%Prtcl_cs_ratio*maxD+maxval(dlub_pp)
#else
    this%LenForCS= DEM_Opt%Prtcl_cs_ratio*maxD
#endif
    if(1.02_RK*this%LenForCS>min(min(lx,ly),lz) ) call DEMLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","so big Diameter")
    le_cs = 2.0_RK*this%LenForCS

    ! (id 1) +(ptype 1) +(PosR 3) +(linvel 3) +(rotvel 3)=11
    this%GhostCS_size   = 11

    ! (id 1) +(ptype 1) +(Mark 1) +(PosR   3) +(linvel 3*tsize) +(linAcc 3*tsize) + &
    !                              (theta  3) +(rotVel 3*rsize) +(rotAcc 3*rsize) = 9+6*(tsize+rsize)
    this%Prtcl_Exchange_size = 9 + 6*(GPrtcl_list%tsize + GPrtcl_list%rsize)
#ifdef CFDDEM
    ! (GPrtcl_FpForce 3)  +(GPrtcl_FpForce_old 3) +(GPrtcl_Vfluid 3+3) +(GPrtcl_linVelOld 3)=15
    this%Prtcl_Exchange_size = this%Prtcl_Exchange_size + 15
    if(is_clc_Basset) this%Prtcl_Exchange_size= this%Prtcl_Exchange_size+ 3*GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    ! FpForce(3)+ FpTorque(3)+ FluidIntOld( 2*3)= 6+ 6 =12
    this%Prtcl_Exchange_size= this%Prtcl_Exchange_size +12
    ! PosOld(3)
    if(IBM_Scheme==2) this%Prtcl_Exchange_size= this%Prtcl_Exchange_size +3
#endif

    xst0_cs = DEM_decomp%xSt-this%LenForCS
    xst1_cs = DEM_decomp%xSt
    xst2_cs = DEM_decomp%xSt+this%LenForCS    
    xed0_cs = DEM_decomp%xEd-this%LenForCS
    xed1_cs = DEM_decomp%xEd
    xed2_cs = DEM_decomp%xEd+this%LenForCS

    yst0_cs = DEM_decomp%ySt-this%LenForCS
    yst1_cs = DEM_decomp%ySt
    yst2_cs = DEM_decomp%ySt+this%LenForCS    
    yed0_cs = DEM_decomp%yEd-this%LenForCS
    yed1_cs = DEM_decomp%yEd
    yed2_cs = DEM_decomp%yEd+this%LenForCS 

    zst0_cs = DEM_decomp%zSt-this%LenForCS
    zst1_cs = DEM_decomp%zSt
    zst2_cs = DEM_decomp%zSt+this%LenForCS    
    zed0_cs = DEM_decomp%zEd-this%LenForCS
    zed1_cs = DEM_decomp%zEd
    zed2_cs = DEM_decomp%zEd+this%LenForCS

    dx_pbc=0.0_RK; dy_pbc=0.0_RK; dz_pbc=0.0_RK
    IF(DEM_decomp%Prtcl_Pencil==x_axis)THEN
      if(pbc(1)) then
        dx_pbc(xm_axis)= simLen%x
        dx_pbc(xp_axis)=-simLen%x     
      endif
      if(pbc(2)) then
        if(DEM_decomp%coord1==0)                 dy_pbc(ym_axis)= simLen%y
        if(DEM_decomp%coord1==DEM_decomp%prow-1) dy_pbc(yp_axis)=-simLen%y
      endif
      if(pbc(3)) then
        if(DEM_decomp%coord2==0)                 dz_pbc(zm_axis)= simLen%z
        if(DEM_decomp%coord2==DEM_decomp%pcol-1) dz_pbc(zp_axis)=-simLen%z
      endif
    ELSEIF(DEM_decomp%Prtcl_Pencil==y_axis)THEN
      if(pbc(1)) then
        if(DEM_decomp%coord1==0)                 dx_pbc(xm_axis)= simLen%x
        if(DEM_decomp%coord1==DEM_decomp%prow-1) dx_pbc(xp_axis)=-simLen%x 
      endif
      if(pbc(2)) then
        dy_pbc(ym_axis)= simLen%y
        dy_pbc(yp_axis)=-simLen%y
      endif
      if(pbc(3)) then
        if(DEM_decomp%coord2==0)                 dz_pbc(zm_axis)= simLen%z
        if(DEM_decomp%coord2==DEM_decomp%pcol-1) dz_pbc(zp_axis)=-simLen%z
      endif
    ELSEIF(DEM_decomp%Prtcl_Pencil==z_axis)THEN
      if(pbc(1)) then
        if(DEM_decomp%coord1==0)                 dx_pbc(xm_axis)= simLen%x
        if(DEM_decomp%coord1==DEM_decomp%prow-1) dx_pbc(xp_axis)=-simLen%x 
      endif
      if(pbc(2)) then
        if(DEM_decomp%coord2==0)                 dy_pbc(ym_axis)= simLen%y
        if(DEM_decomp%coord2==DEM_decomp%pcol-1) dy_pbc(yp_axis)=-simLen%y
      endif
      if(pbc(3)) then
        dz_pbc(zm_axis)= simLen%z
        dz_pbc(zp_axis)=-simLen%z
      endif
    ENDIF

    vol_tot= SimLen%x * SimLen%y * SimLen%z
    vol_ghost_cs=(lx+le_cs)*(ly+le_cs)*(lz+le_cs)-lx*ly*lz

    numPrtcl = DEM_opt%numPrtcl
    mCS= int(numPrtcl*real(vol_ghost_cs/vol_tot,RK))
    mCS= 2*mCS
    mCS= min(mCS, numPrtcl)
    mCS= max(mCS, 10)
    GPrtcl_list%mGhost_CS = mCS

    mFixedCS= int(DEM_opt%numPrtclFix*real(vol_ghost_cs/vol_tot,RK))
    mFixedCS= 2*mFixedCS
    mFixedCS= min(mFixedCS, DEM_opt%numPrtclFix)
    mFixedCS= max(mFixedCS, 10)
    GPrtcl_list%nGhostFix_CS= mFixedCS

    vol1=(lx+le_cs)*(ly+le_cs)*this%LenForCS
    vol2=(lx+le_cs)*(lz+le_cs)*this%LenForCS
    vol3=(ly+le_cs)*(lz+le_cs)*this%LenForCS
    vol_sendlist=max(max(vol1,vol2),vol3)
    msend=int(DEM_opt%numPrtclFix*real(vol_sendlist/vol_tot,RK))
    msend=2*msend
    msend=min(msend,DEM_opt%numPrtclFix)
    msend=max(msend,10)
    this%msend=msend
    if(allocated(sendlist)) deallocate(sendlist);
    allocate(GhostPFix_id(mFixedCS),   Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostPFix_pType(mFixedCS),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostPFix_PosR(mFixedCS), Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(sendlist(msend),          Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","Allocation failed-1")
    call this%Comm_For_Cntct_fixed()

    msend=int(numPrtcl*real(vol_sendlist/vol_tot,RK))
    msend=2*msend
    msend=min(msend,numPrtcl)
    msend=max(msend,10)
    this%msend=msend
    if(allocated(sendlist)) deallocate(sendlist);
    allocate(GhostP_id(mCs),    Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_pType(mCs), Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_PosR(mCs),  Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_linVel(mCs),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_rotVel(mCS),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(sendlist(msend),   Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","Allocation failed-2")
  end subroutine PC_InitComm

  !**********************************************************************
  ! Prtcl_Comm_For_Cntct
  !**********************************************************************
  subroutine PC_Comm_For_Cntct(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    real(RK)::px,py,pz
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp
    integer::i,ierror,request(4),SRstatus(MPI_STATUS_SIZE)
    
    ng=0
    SELECT CASE(DEM_decomp%Prtcl_Pencil)
    CASE(x_axis)     ! ccccccccccccccccccccccccc  x-axis  ccccccccccccccccccccccccccc

      ! step1: Handle x-dir
      IF(pbc(1)) THEN
        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo
      ENDIF

      ! step2: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_pType(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(y_axis)     ! ccccccccccccccccccccccccc  y-axis  ccccccccccccccccccccccccccc

      ! step1: Handle y-dir
      IF(pbc(2)) THEN
        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_pType(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(z_axis)     ! ccccccccccccccccccccccccc  z-axis  ccccccccccccccccccccccccccc

      ! step1: Handle z-dir
      IF(pbc(3)) THEN
        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_pType(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py>=yed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>=GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)
    END SELECT
    GPrtcl_list%nGhost_CS = ng
  end subroutine PC_Comm_For_Cntct

  !**********************************************************************
  ! PC_pack_cntct
  !**********************************************************************
  subroutine PC_pack_cntct(this,buf_send,nsendg,nsend,dir)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsendg,nsend,dir

    ! locals
    integer::i,id,m
    real(RK)::dx,dy,dz

    m=1
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)
    do i=1,nsendg
      id=sendlist(i)
      buf_send(m)=real(GhostP_id(id));    m=m+1 ! 01
      buf_send(m)=real(GhostP_pType(id)); m=m+1 ! 02
      buf_send(m)=GhostP_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GhostP_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GhostP_PosR(id)%z+dz;   m=m+1 ! 05
      buf_send(m)=GhostP_linVel(id)%x;    m=m+1 ! 06
      buf_send(m)=GhostP_linVel(id)%y;    m=m+1 ! 07
      buf_send(m)=GhostP_linVel(id)%z;    m=m+1 ! 08
      buf_send(m)=GhostP_rotVel(id)%x;    m=m+1 ! 09
      buf_send(m)=GhostP_rotVel(id)%y;    m=m+1 ! 10
      buf_send(m)=GhostP_rotVel(id)%z;    m=m+1 ! 11 
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(GPrtcl_id(id));    m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id)); m=m+1 ! 02
      buf_send(m)=GPrtcl_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%z+dz;   m=m+1 ! 05
      buf_send(m)=GPrtcl_linVel(1,id)%x;  m=m+1 ! 06
      buf_send(m)=GPrtcl_linVel(1,id)%y;  m=m+1 ! 07
      buf_send(m)=GPrtcl_linVel(1,id)%z;  m=m+1 ! 08
      buf_send(m)=GPrtcl_rotVel(1,id)%x;  m=m+1 ! 09
      buf_send(m)=GPrtcl_rotVel(1,id)%y;  m=m+1 ! 10
      buf_send(m)=GPrtcl_rotVel(1,id)%z;  m=m+1 ! 11 
    enddo
  end subroutine PC_pack_cntct

  !**********************************************************************
  ! PC_unpack_cntct
  !**********************************************************************
  subroutine PC_unpack_cntct(this,buf_recv,n1,n2)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      GhostP_id(i)       =nint(buf_recv(m)); m=m+1 ! 01 
      itype              =nint(buf_recv(m)); m=m+1 ! 02
      GhostP_PosR(i)%x   =buf_recv(m);       m=m+1 ! 03
      GhostP_PosR(i)%y   =buf_recv(m);       m=m+1 ! 04
      GhostP_PosR(i)%z   =buf_recv(m);       m=m+1 ! 05
      GhostP_linVel(i)%x =buf_recv(m);       m=m+1 ! 06 
      GhostP_linVel(i)%y =buf_recv(m);       m=m+1 ! 07 
      GhostP_linVel(i)%z =buf_recv(m);       m=m+1 ! 08 
      GhostP_rotVel(i)%x =buf_recv(m);       m=m+1 ! 09 
      GhostP_rotVel(i)%y =buf_recv(m);       m=m+1 ! 10 
      GhostP_rotVel(i)%z =buf_recv(m);       m=m+1 ! 11
      GhostP_pType(i)    =itype;                
      GhostP_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine PC_unpack_cntct

  !**********************************************************************
  ! PC_Comm_For_Exchange
  !**********************************************************************
  subroutine PC_Comm_For_Exchange(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    integer::i,ierror,request(4),nlocal,nlocalp
    integer::nsend(2),nrecv(2),nlink
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer,dimension(MPI_STATUS_SIZE) :: SRstatus
    real(RK)::px,py,pz

    nlocal=GPrtcl_list%nlocal
    SELECT CASE(DEM_decomp%Prtcl_Pencil)
    CASE(x_axis)     ! ccccccccccccccccccccccccc  x-axis  ccccccccccccccccccccccccccc

      ! step1: Handle x-dir
      IF(pbc(1)) THEN
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            GPrtcl_PosR(i)%x=px - simLen%x
          elseif(px < xst1_cs) then
            GPrtcl_PosR(i)%x=px + simLen%x
          endif
        enddo
      ELSE
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs .or. px < xst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo
      ENDIF

      ! step2: send to yp_axis, and receive from ym_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            GPrtcl_PosR(i)%y=py - simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(3), 9, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,yp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,yp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to ym_axis, and receive from yp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            GPrtcl_PosR(i)%y=py + simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(4),11, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,ym_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,ym_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            GPrtcl_PosR(i)%z=pz - simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(1),13, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            GPrtcl_PosR(i)%z=pz + simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(2),15, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(y_axis)     ! ccccccccccccccccccccccccc  y-axis  ccccccccccccccccccccccccccc

      ! step1: Handle y-dir
      IF(pbc(2)) THEN
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            GPrtcl_PosR(i)%y=py - simLen%y
          elseif(py < yst1_cs) then
            GPrtcl_PosR(i)%y=py + simLen%y
          endif
        enddo
      ELSE
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs .or. py < yst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            GPrtcl_PosR(i)%x=px - simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(3), 9, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            GPrtcl_PosR(i)%x=px + simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(4),11, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            GPrtcl_PosR(i)%z=pz - simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(1),13, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            GPrtcl_PosR(i)%z=pz + simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(2),15, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(z_axis)     ! ccccccccccccccccccccccccc  z-axis  ccccccccccccccccccccccccccc

      ! step1: Handle z-dir
      IF(pbc(3)) THEN
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            GPrtcl_PosR(i)%z=pz - simLen%z
          elseif(pz < zst1_cs) then
            GPrtcl_PosR(i)%z=pz + simLen%z
          endif
        enddo
      ELSE
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs .or. pz < zst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            GPrtcl_PosR(i)%x=px - simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(3), 9, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            GPrtcl_PosR(i)%x=px + simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(4),11, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to yp_axis, and receive from ym_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            GPrtcl_PosR(i)%y=py - simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(1),13, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,yp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,yp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to ym_axis, and receive from yp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            GPrtcl_PosR(i)%y=py + simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(2),15, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,ym_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,ym_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)
    END SELECT
    GPrtcl_list%nlocal = nlocal
  end subroutine PC_Comm_For_Exchange

  !**********************************************************************
  ! PC_pack_Exchange
  !**********************************************************************
  subroutine PC_pack_Exchange(this,buf_send,nsend,nlocalp,dir)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsend,dir,nlocalp

    ! locals
    real(RK)::dx,dy,dz
    integer::i,j,id,m,nlocal

    m=1
    nlocal = nlocalp
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)

#ifdef CFDDEM
    IF(is_clc_Basset)  THEN
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;     m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%x; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%y; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%z; m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%x;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%y;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%z;   m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%z;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%z;    m=m+1
        do j=1, GPrtcl_BassetSeq%nDataLen
          buf_send(m)=GPrtcl_BassetData(j,id)%x; m=m+1
          buf_send(m)=GPrtcl_BassetData(j,id)%y; m=m+1
          buf_send(m)=GPrtcl_BassetData(j,id)%z; m=m+1
        enddo
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    ELSE
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;     m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%x; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%y; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%z; m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%x;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%y;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%z;   m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%z;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%z;    m=m+1
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    ENDIF
#elif CFDACM
    if(IBM_Scheme==2) then
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;       m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%x;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%y;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%z;      m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%z; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%z; m=m+1
        buf_send(m)=GPrtcl_PosOld(id)%x+dx;     m=m+1
        buf_send(m)=GPrtcl_PosOld(id)%y+dy;     m=m+1
        buf_send(m)=GPrtcl_PosOld(id)%z+dz;     m=m+1
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    else
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;       m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%x;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%y;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%z;      m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%z; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%z; m=m+1
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    endif
#else
   DO i=1,nsend
      id = sendlist(i)
      buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
      buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
      buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
      buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
      buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
      buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
      do j=1,GPrtcl_list%tsize
        buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
        buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
      enddo
      do j=1,GPrtcl_list%rsize
        buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
      enddo
      call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
      buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
    ENDDO
#endif

    DO i=nsend,1,-1
      id = sendlist(i)
      call GPrtcl_list%copy(id,nlocal)
      call GPPW_CntctList%copy(id,nlocal)
      call DEMContactSearchPW%copy(id,nlocal)
      nlocal = nlocal -1
    ENDDO
  end subroutine PC_pack_Exchange

  !**********************************************************************
  ! PC_unpack_Exchange
  !**********************************************************************
  subroutine PC_unpack_Exchange(this,buf_recv,nrecv,nlocal,dir)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::nrecv,dir
    integer,intent(inout)::nlocal

    ! locals
    integer::i,j,id,m,itype

    m=1
    id=nlocal+1
#ifdef CFDDEM
    If (is_clc_Basset) then
      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = nint(buf_recv(m)); m=m+1 ! 01
        itype             = nint(buf_recv(m)); m=m+1 ! 02
        GPrtcl_usrMark(id)= nint(buf_recv(m)); m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);       m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);       m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);       m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);       m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);       m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m); m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m); m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m); m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m); m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m); m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m); m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m); m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m); m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x    =buf_recv(m);  m=m+1
        GPrtcl_FpForce(id)%y    =buf_recv(m);  m=m+1
        GPrtcl_FpForce(id)%z    =buf_recv(m);  m=m+1
        GPrtcl_FpForce_old(id)%x=buf_recv(m);  m=m+1
        GPrtcl_FpForce_old(id)%y=buf_recv(m);  m=m+1
        GPrtcl_FpForce_old(id)%z=buf_recv(m);  m=m+1
        GPrtcl_linVelOld(id)%x  =buf_recv(m);  m=m+1
        GPrtcl_linVelOld(id)%y  =buf_recv(m);  m=m+1
        GPrtcl_linVelOld(id)%z  =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(1,id)%x   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(1,id)%y   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(1,id)%z   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(2,id)%x   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(2,id)%y   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(2,id)%z   =buf_recv(m);  m=m+1
        do j=1,GPrtcl_BassetSeq%nDataLen
          GPrtcl_BassetData(j,id)%x= buf_recv(m); m=m+1
          GPrtcl_BassetData(j,id)%y= buf_recv(m); m=m+1
          GPrtcl_BassetData(j,id)%z= buf_recv(m); m=m+1
        enddo
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius
        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    ELSE

      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = nint(buf_recv(m)); m=m+1 ! 01
        itype             = nint(buf_recv(m)); m=m+1 ! 02
        GPrtcl_usrMark(id)= nint(buf_recv(m)); m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);       m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);       m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);       m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);       m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);       m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m); m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m); m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m); m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m); m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m); m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m); m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m); m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m); m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m); m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x    =buf_recv(m);  m=m+1
        GPrtcl_FpForce(id)%y    =buf_recv(m);  m=m+1
        GPrtcl_FpForce(id)%z    =buf_recv(m);  m=m+1
        GPrtcl_FpForce_old(id)%x=buf_recv(m);  m=m+1
        GPrtcl_FpForce_old(id)%y=buf_recv(m);  m=m+1
        GPrtcl_FpForce_old(id)%z=buf_recv(m);  m=m+1
        GPrtcl_linVelOld(id)%x  =buf_recv(m);  m=m+1
        GPrtcl_linVelOld(id)%y  =buf_recv(m);  m=m+1
        GPrtcl_linVelOld(id)%z  =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(1,id)%x   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(1,id)%y   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(1,id)%z   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(2,id)%x   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(2,id)%y   =buf_recv(m);  m=m+1
        GPrtcl_Vfluid(2,id)%z   =buf_recv(m);  m=m+1
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
 
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    ENDIF
#elif CFDACM
    if(IBM_Scheme==2) then
      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = nint(buf_recv(m));  m=m+1 ! 01
        itype             = nint(buf_recv(m));  m=m+1 ! 02
        GPrtcl_usrMark(id)= nint(buf_recv(m));  m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);        m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);        m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);        m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);        m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);        m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);        m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m);  m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m);  m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m);  m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m);  m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m);  m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m);  m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m);  m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m);  m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x = buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%y = buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%z = buf_recv(m);     m=m+1
        GPrtcl_FpTorque(id)%x= buf_recv(m);     m=m+1
        GPrtcl_FpTorque(id)%y= buf_recv(m);     m=m+1
        GPrtcl_FpTorque(id)%z= buf_recv(m);     m=m+1
        GPrtcl_FluidIntold(1,id)%x=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(1,id)%y=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(1,id)%z=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(2,id)%x=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(2,id)%y=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(2,id)%z=buf_recv(m); m=m+1
        GPrtcl_PosOld(id)%x=buf_recv(m);        m=m+1
        GPrtcl_PosOld(id)%y=buf_recv(m);        m=m+1
        GPrtcl_PosOld(id)%z=buf_recv(m);        m=m+1
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
 
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    else
      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = nint(buf_recv(m));  m=m+1 ! 01
        itype             = nint(buf_recv(m));  m=m+1 ! 02
        GPrtcl_usrMark(id)= nint(buf_recv(m));  m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);        m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);        m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);        m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);        m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);        m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);        m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m);  m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m);  m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m);  m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m);  m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m);  m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m);  m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m);  m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m);  m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m);  m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x = buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%y = buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%z = buf_recv(m);     m=m+1
        GPrtcl_FpTorque(id)%x= buf_recv(m);     m=m+1
        GPrtcl_FpTorque(id)%y= buf_recv(m);     m=m+1
        GPrtcl_FpTorque(id)%z= buf_recv(m);     m=m+1
        GPrtcl_FluidIntold(1,id)%x=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(1,id)%y=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(1,id)%z=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(2,id)%x=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(2,id)%y=buf_recv(m); m=m+1
        GPrtcl_FluidIntold(2,id)%z=buf_recv(m); m=m+1
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
 
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    endif
#else
    DO i=1,nrecv
      if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
      GPrtcl_id(id)     = nint(buf_recv(m)); m=m+1 ! 01
      itype             = nint(buf_recv(m)); m=m+1 ! 02
      GPrtcl_usrMark(id)= nint(buf_recv(m)); m=m+1 ! 03
      GPrtcl_PosR(id)%x = buf_recv(m);       m=m+1 ! 04
      GPrtcl_PosR(id)%y = buf_recv(m);       m=m+1 ! 05
      GPrtcl_PosR(id)%z = buf_recv(m);       m=m+1 ! 06
      GPrtcl_theta(id)%x= buf_recv(m);       m=m+1 ! 07
      GPrtcl_theta(id)%y= buf_recv(m);       m=m+1 ! 08
      GPrtcl_theta(id)%z= buf_recv(m);       m=m+1 ! 09
      do j=1,GPrtcl_list%tsize
        GPrtcl_linVel(j,id)%x = buf_recv(m); m=m+1 ! 6* tsize
        GPrtcl_linVel(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_linVel(j,id)%z = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%x = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%z = buf_recv(m); m=m+1
      enddo
      do j=1,GPrtcl_list%rsize
        GPrtcl_rotVel(j,id)%x = buf_recv(m); m=m+1
        GPrtcl_rotVel(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_rotVel(j,id)%z = buf_recv(m); m=m+1
        GPrtcl_rotAcc(j,id)%x = buf_recv(m); m=m+1
        GPrtcl_rotAcc(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_rotAcc(j,id)%z = buf_recv(m); m=m+1 ! 6* rsize
      enddo
      GPrtcl_pType(id)=itype
      GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

      call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
      call DEMContactSearchPW%InsertNearPW(id)
 
      id = id + 1
      nlocal =nlocal + 1
    ENDDO
#endif
  end subroutine PC_unpack_Exchange

  !**********************************************************************
  ! ISInThisProc
  !**********************************************************************
  function PC_ISInThisProc(this,buf_recv,m,dir) result(res)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(inout)::m
    integer,intent(in)::dir
    logical:: res

    !local
    integer:: mp

    res=.true.
    SELECT CASE(dir)
    CASE(xp_axis)
      mp = m+3
      if(buf_recv(mp)>=xed1_cs) res = .false.
    CASE(xm_axis)
      mp = m+3
      if(buf_recv(mp)< xst1_cs) res = .false.
    CASE(yp_axis)
      mp = m+4
      if(buf_recv(mp)>=yed1_cs) res = .false. 
    CASE(ym_axis)
      mp = m+4
      if(buf_recv(mp)< yst1_cs) res = .false.
    CASE(zp_axis)
      mp = m+5
      if(buf_recv(mp)>=zed1_cs) res = .false.
    CASE(zm_axis)
      mp = m+5
      if(buf_recv(mp)< zst1_cs) res = .false.
    END SELECT
    if(res) return
    
    call DEMLogInfo%CheckForError(ErrT_Pass," PC_ISInThisProc"," The following particle is deleted: ")
    call DEMLogInfo%OutInfo(" Exchange direction  is :"//trim(num2str(dir)  ),3)
    call DEMLogInfo%OutInfo("   The particle id   is :"//trim(num2str( nint(buf_recv(m)  ))),3)
    call DEMLogInfo%OutInfo("       particle type is :"//trim(num2str( nint(buf_recv(m+1)))),3)
    call DEMLogInfo%OutInfo("       x-coordinate  is :"//trim(num2str( buf_recv(m+3) )), 3)
    call DEMLogInfo%OutInfo("       y-coordinate  is :"//trim(num2str( buf_recv(m+4) )), 3)
    call DEMLogInfo%OutInfo("       z-coordinate  is :"//trim(num2str( buf_recv(m+5) )), 3)
    call DEMLogInfo%OutInfo(" Present processor   is :"//trim(num2str(nrank)), 3)

    m = m + this%Prtcl_Exchange_size
    do while(abs(buf_recv(m)-END_OF_PRTCL)>=1.00E-10_RK)
      m=m+1
    enddo
    m=m+1
  end function PC_ISInThisProc

  !**********************************************************************
  ! PC_reallocate_ghost_for_Cntct
  !**********************************************************************  
  subroutine PC_reallocate_ghost_for_Cntct(this,ng)
    implicit none
    class(Prtcl_Comm_info)::this
    integer,intent(in)::ng

    ! locals
    integer::sizep,sizen,ierrTmp,ierror=0
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec

    sizep= GPrtcl_list%mGhost_CS
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,DEM_Opt%numPrtcl)
    sizen= max(sizen,ng+1)
    GPrtcl_list%mGhost_CS=sizen  ! NOTE HERE, sometimes GPrtcl_list%mGhost_CS CAN bigger than DEM_Opt%numPrtcl

    ! ======= integer vector part =======
    call move_alloc(GhostP_id,IntVec)
    allocate(GhostP_id(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_id(1:sizep)=IntVec

    call move_alloc(GhostP_pType,IntVec)
    allocate(GhostP_pType(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_pType(1:sizep)=IntVec
    deallocate(IntVec)  

    ! ======= real3 vercor part =======
    call move_alloc(GhostP_linVel,Real3Vec)
    allocate(GhostP_linVel(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_linVel(1:sizep)=Real3Vec
 
    call move_alloc(GhostP_rotVel,Real3Vec)
    allocate(GhostP_rotVel(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_rotVel(1:sizep)=Real3Vec
    deallocate(Real3Vec)

    ! ======= real4 vercor part =======
    call move_alloc(GhostP_PosR,Real4Vec)
    allocate(GhostP_PosR(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_PosR(1:sizep)=Real4Vec
    deallocate(Real4Vec)

    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," PC_reallocate_ghost_for_Cntct"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif   
    !call DEMLogInfo%CheckForError(ErrT_Pass," reallocate_ghost_for_Cntct"," Need to reallocate Ghost variables")
    !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine PC_reallocate_ghost_for_Cntct

  !**********************************************************************
  ! PC_reallocate_sendlist
  !**********************************************************************  
  subroutine PC_reallocate_sendlist(this,ns)
    implicit none
    class(Prtcl_Comm_info)::this
    integer,intent(in)::ns

    ! locals
    integer::sizep,sizen,ierror
    integer,dimension(:),allocatable:: IntVec

    sizep=this%msend
    sizen=int(1.2_RK*real(sizep,kind=RK))
    sizen=min(sizen,DEM_Opt%numPrtcl)
    sizen=max(sizen,ns+1)
    this%msend=sizen
   
    call move_alloc(sendlist, IntVec)
    allocate(sendlist(sizen),stat=ierror)
    sendlist(1:sizep)=IntVec
    deallocate(IntVec)

    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," PC_reallocate_sendlist"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif   
    !call DEMLogInfo%CheckForError(ErrT_Pass," reallocate_sendlist"," Need to reallocate sendlist")
    !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine PC_reallocate_sendlist

  !**********************************************************************
  ! Prtcl_Comm_For_Cntct_fixed
  !**********************************************************************
  subroutine PC_Comm_For_Cntct_fixed(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    real(RK)::px,py,pz
    integer,dimension(:),allocatable:: IntVec
    type(real4),dimension(:),allocatable::Real4Vec
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::i,ierror,request(4),SRstatus(MPI_STATUS_SIZE)
    integer::nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp,GhostFixedCS_size
    
    ! (id 1) +(ptype 1) +(PosR 3) =5
    GhostFixedCS_size   = 5

    ng=0
    SELECT CASE(DEM_decomp%Prtcl_Pencil)
    CASE(x_axis)     ! ccccccccccccccccccccccccc  x-axis  ccccccccccccccccccccccccccc

      ! step1: Handle x-dir
      IF(pbc(1)) THEN
        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo
      ENDIF

      ! step2: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py>=yed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(y_axis)     ! ccccccccccccccccccccccccc  y-axis  ccccccccccccccccccccccccccc

      ! step1: Handle y-dir
      IF(pbc(2)) THEN
        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(z_axis)     ! ccccccccccccccccccccccccc  z-axis  ccccccccccccccccccccccccccc

      ! step1: Handle z-dir
      IF(pbc(3)) THEN
        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py>=yed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>=GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)
    END SELECT

    if(ng>0) then  ! 2021-02-11, Gong Zheng
      if(GPrtcl_list%mlocalFix>0) then
        call move_alloc(GPFix_id,IntVec)
        allocate(GPFix_id(GPrtcl_list%mlocalFix +ng))
        GPFix_id(1:GPrtcl_list%mlocalFix)= IntVec
        GPFix_id(GPrtcl_list%mlocalFix+1:GPrtcl_list%mlocalFix +ng)   = GhostPFix_id(1:ng)

        call move_alloc(GPFix_pType,IntVec)
        allocate(GPFix_pType(GPrtcl_list%mlocalFix +ng))
        GPFix_pType(1:GPrtcl_list%mlocalFix)= IntVec
        GPFix_pType(GPrtcl_list%mlocalFix+1:GPrtcl_list%mlocalFix +ng)= GhostPFix_pType(1:ng)
        deallocate(IntVec)

        call move_alloc(GPFix_PosR,Real4Vec)
        allocate(GPFix_PosR(GPrtcl_list%mlocalFix +ng))
        GPFix_PosR(1:GPrtcl_list%mlocalFix)= Real4Vec
        GPFix_PosR(GPrtcl_list%mlocalFix+1:GPrtcl_list%mlocalFix +ng) = GhostPFix_PosR(1:ng)
        deallocate(Real4Vec)
      else
        if(allocated(GPFix_id))   deallocate(GPFix_id);   allocate(GPFix_id(ng));   GPFix_id= GhostPFix_id(1:ng)
        if(allocated(GPFix_PosR)) deallocate(GPFix_PosR); allocate(GPFix_PosR(ng)); GPFix_PosR=GhostPFix_PosR(1:ng)
        if(allocated(GPFix_pType))deallocate(GPFix_pType);allocate(GPFix_pType(ng));GPFix_pType=GhostPFix_pType(1:ng)      
      endif
    endif
    deallocate(GhostPFix_id, GhostPFix_pType, GhostPFix_PosR)
    GPrtcl_list%nGhostFix_CS = ng
  end subroutine PC_Comm_For_Cntct_fixed

  !**********************************************************************
  ! pack_cntct_fixed
  !**********************************************************************
  subroutine pack_cntct_fixed(buf_send,nsendg,nsend,dir)
    implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsendg,nsend,dir

    ! locals
    integer::i,id,m
    real(RK)::dx,dy,dz

    m=1
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)
    do i=1,nsendg
      id=sendlist(i)
      buf_send(m)=real(GhostPFix_id(id));    m=m+1 ! 01
      buf_send(m)=real(GhostPFix_pType(id)); m=m+1 ! 02
      buf_send(m)=GhostPFix_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GhostPFix_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GhostPFix_PosR(id)%z+dz;   m=m+1 ! 05
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(GPFix_id(id));    m=m+1 ! 01
      buf_send(m)=real(GPFix_pType(id)); m=m+1 ! 02
      buf_send(m)=GPFix_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GPFix_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GPFix_PosR(id)%z+dz;   m=m+1 ! 05
    enddo
  end subroutine pack_cntct_fixed

  !**********************************************************************
  ! unpack_cntct_fixed
  !**********************************************************************
  subroutine unpack_cntct_fixed(buf_recv,n1,n2)
    implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      GhostPFix_id(i)       =nint(buf_recv(m)); m=m+1 ! 01 
      itype                 =nint(buf_recv(m)); m=m+1 ! 02
      GhostPFix_PosR(i)%x   =buf_recv(m);       m=m+1 ! 03
      GhostPFix_PosR(i)%y   =buf_recv(m);       m=m+1 ! 04
      GhostPFix_PosR(i)%z   =buf_recv(m);       m=m+1 ! 05
      GhostPFix_pType(i)    =itype;                
      GhostPFix_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine unpack_cntct_fixed

  !**********************************************************************
  ! reallocate_ghostFix_for_Cntct
  !**********************************************************************  
  subroutine reallocate_ghostFix_for_Cntct(ng)
    implicit none
    integer,intent(in)::ng

    ! locals
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec
    type(real4),dimension(:),allocatable::Real4Vec

    sizep= GPrtcl_list%nGhostFix_CS
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,DEM_Opt%numPrtclFix)
    sizen= max(sizen,ng+1)
    GPrtcl_list%nGhostFix_CS=sizen

    ! ======= integer vector part =======
    call move_alloc(GhostPFix_id,IntVec)
    allocate(GhostPFix_id(sizen))
    GhostPFix_id(1:sizep)=IntVec

    call move_alloc(GhostPFix_pType,IntVec)
    allocate(GhostPFix_pType(sizen))
    GhostPFix_pType(1:sizep)=IntVec
    deallocate(IntVec)  

    ! ======= real4 verctor part =======
    call move_alloc(GhostPFix_PosR,Real4Vec)
    allocate(GhostPFix_PosR(sizen))
    GhostPFix_PosR(1:sizep)=Real4Vec
    deallocate(Real4Vec)
  end subroutine reallocate_ghostFix_for_Cntct

#undef xm_axis
#undef xp_axis
#undef ym_axis
#undef yp_axis
#undef zm_axis
#undef zp_axis
end module Prtcl_Comm
    Vel_w=((ri*Rvei+rj*Rvej).cross.Norm_v)
    Vrij = Veli-Velj + Vel_w ! Relative velocity
    vrn  = Vrij .dot. Norm_v
    if(VelRel_Init(ind)<-1.0_RK)VelRel_Init(ind)=abs(vrn)

    Vij_n= vrn*Norm_v        ! Normal  relative velocity
    Vij_t= Vrij-Vij_n        ! Tangent relative velocity
    
    ! Tangential overlap Vector
    Ovlp_t= TanDelta(ind)
    normTan1=norm(Ovlp_t)
    Ovlp_t= Ovlp_t-(Ovlp_t .dot. Norm_v)*Norm_v
    normTan2=norm(Ovlp_t)
    if(normTan2>1.0E-10_RK) then
      Ovlp_t=(normTan1/normTan2)*Ovlp_t
    else
      Ovlp_t=zero_r3
    endif
    Ovlp_t= (Vij_t*DEM_opt%dt)+Ovlp_t
   
    ! Computing the normal and tangential contact forces
    select case(DEM_opt%CF_Type)
    case(DEM_LSD)
      Vel_in=max(VelRel_Init(ind),1.0E-5_RK); Vel_in=Vel_in**0.2_RK
      k_n= -Prop_ij%StiffnessCoe_n*Vel_in*Vel_in
      d_n= -Prop_ij%DampingCoe_n*Vel_in
      fn =  k_n*ovrlp +d_n* vrn
      k_t= -Prop_ij%StiffnessCoe_t*Vel_in*Vel_in
      d_t= -Prop_ij%DampingCoe_t*Vel_in
      Ftij= k_t*Ovlp_t +d_t*Vij_t
    case(DEM_nLin)
      k_n= -Prop_ij%StiffnessCoe_n
      d_n= -Prop_ij%DampingCoe_n
      fn =  k_n*ovrlp**1.5_RK +d_n*(ovrlp**0.25_RK)*vrn  ! 2.62, P39
      k_t= -Prop_ij%StiffnessCoe_t*sqrt(ovrlp)
      d_t=  0.0_RK  ! No equation is considered for tangential damping yet
      Ftij= k_t*Ovlp_t +d_t*Vij_t                        ! 2.72, P44
    end select
    Fnij= fn*Norm_v

    ! Coulomb's friction law
    ft = norm(Ftij)
    ft_fric = Prop_ij%FrictionCoe_s*abs(fn)
    if(ft>ft_fric) then
      ft_fric = Prop_ij%FrictionCoe_k* abs(fn)
      Ftij =(ft_fric/ft)*Ftij
      Ovlp_t = (1.0_RK/k_t)*Ftij          ! Caution Here !
    endif
    
    ! Computing rolling resistance torque acting on spheres
    !W_hat = Rvei-Rvej
    !w_hat_mag= norm(W_hat) ! 2.123, p56
    !if(w_hat_mag > 1.0E-10_RK) then
    !  W_hat= W_hat/w_hat_mag
    !else
    !  W_hat= zero_r3
    !endif
    !if(DEM_opt%CT_Model == CTM_ConstantTorque) then 
    !  Mrij=(-Prop_ij%FrictionCoe_Roll*abs(fn)*Prop_ij%RadEff)*W_hat  ! 2.122, p56
    !else
    !  Mrij=(-Prop_ij%FrictionCoe_Roll*abs(fn)*Prop_ij%RadEff*norm(Vel_w))*W_hat  ! 2.124, p57
    !endif

    ! Setting the updated contact info pair into the contact list 
    TanDelta(ind) = Ovlp_t

    ! Updating the contact force and torques of particles i and j
    Moment= Norm_v .cross. Ftij
#ifdef ContactForce_PP
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) +(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +ri*Moment!+Mrij
    GPrtcl_cntctForce(pjd)= GPrtcl_cntctForce(pjd) -(Fnij+Ftij)
    GPrtcl_torque(pjd)= GPrtcl_torque(pjd) +rj*Moment!-Mrij
#endif
#ifdef ContactForce_PPG
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) +(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +ri*Moment!+Mrij
#endif
#ifdef ContactForce_PGP
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) -(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +rj*Moment!-Mrij
#endif
#ifdef ContactForce_PPFix_W
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) +(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +ri*Moment!+Mrij
#endif
module Prtcl_ContactSearch
  use Prtcl_Parameters
  use Prtcl_NBS_Munjiza
  use Prtcl_Hrchl_Munjiza
  implicit none
  private
    
  type::ContactSearch
  contains
    procedure:: InitContactSearch => CS_InitContactSearch
    procedure:: FindContacts      => CS_FindContacts
    procedure:: get_numContact    => CS_get_numContact
  end type ContactSearch
  type(ContactSearch),public:: DEMContactSearch
    
contains

  !**********************************************************************
  ! Initializing particle-particle contact search object
  !**********************************************************************
  subroutine CS_InitContactSearch(this )
    implicit none
    class(ContactSearch)::this
    
    SELECT CASE(DEM_opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      allocate( m_NBS_Munjiza )
      call m_NBS_Munjiza%Init_NBSM()
    CASE( CSM_NBS_Munjiza_Hrchl )
      allocate( m_NBS_Munjiza_Hrchl)
      call m_NBS_Munjiza_Hrchl%Init_Munjiza_Hrchl()
    END SELECT
  end subroutine CS_InitContactSearch                            

  !**********************************************************************
  ! finding contact pairs of particles
  !**********************************************************************
  subroutine CS_FindContacts( this )
    implicit none
    class(ContactSearch)::this

    SELECT CASE( DEM_opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      call m_NBS_Munjiza%ContactSearch()
    CASE( CSM_NBS_Munjiza_Hrchl )
      call m_NBS_Munjiza_Hrchl%ContactSearch()
    END SELECT
  end subroutine CS_FindContacts

  !**********************************************************************
  ! CS_get_numContact
  !**********************************************************************
  function CS_get_numContact( this ) result (res)
    implicit none
    class(ContactSearch)::this
    integer,dimension(2):: res
    
    res = 0
    SELECT CASE( DEM_opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      res(1) = m_NBS_Munjiza%num_Cnsv_cntct
    CASE( CSM_NBS_Munjiza_Hrchl )
      res(1) = m_NBS_Munjiza_Hrchl%num_Cnsv_cntct
      res(2) = m_NBS_Munjiza_Hrchl%lvl_num_cnsv_cntct
    END SELECT
  end function CS_get_numContact

end module Prtcl_ContactSearch
module Prtcl_ContactSearchPW
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Geometry
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
  use Prtcl_Decomp_2d,only:int_type,real_type
#else
  use Prtcl_Decomp_2d,only:nrank,int_type,real_type
#endif
  implicit none
  private
    
  integer,dimension(:),allocatable:: Bucket_PWCS ! Head for plane-wall contact search
  integer,dimension(:),allocatable:: id_Wall     ! particle id
  integer,dimension(:),allocatable:: Next_PWCS

  type::ContactSearchPW
    integer :: mHead
    integer :: NextInsert
    real(RK):: MaxWallVel  ! new added
    real(RK):: DeltaXmax1  ! new added
    real(RK):: DeltaXmax2  ! new added
    integer :: Next_iter_update = 0
    integer :: max_nearPrtcl_pWall
    integer :: num_nearPrtcl_pWall = 0
  contains
    procedure:: InitContactSearchPW
    procedure:: FindContactsPW
    procedure:: UpdateNearPrtclsPW
    procedure:: Init_PWCS_List
    procedure:: Reallocate_PWCS_List
    procedure:: reallocate_Bucket
    procedure:: InsertNearPW
    procedure:: copy => CSPW_copy
  end type ContactSearchPW
  type(ContactSearchPW),public::DEMContactSearchPW
    
contains

  !**********************************************************************
  ! Initializing the object
  !**********************************************************************
  subroutine InitContactSearchPW(this )
    implicit none
    class(ContactSearchPW)::this
    type(real3):: wallvel
    integer:: i,nw,iErr1,iErr2,iErr3,iErrSum,ierror
    real(RK)::MaxWallVel,maxvel1
        
    this%max_nearPrtcl_pWall = int(1.5_RK*GPrtcl_list%mlocal)+1
    this%mHead=GPrtcl_list%mlocal
    allocate(Bucket_PWCS(this%mHead), Stat=iErr1)
    allocate(id_Wall(this%max_nearPrtcl_pWall ), Stat=iErr2)
    allocate(Next_PWCS(this%max_nearPrtcl_pWall ), Stat=iErr3)
    iErrSum = abs(iErr1)+abs(iErr2)+abs(iErr3)
    if(iErrSum/=0 ) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize Contact List","Allocation failed ")
    call this%Init_PWCS_List()
        
    MaxWallVel = 0.0_RK
    nw=DEMGeometry%nPW_local
    do i=1,nw
      wallvel=DEMGeometry%pWall(i)%trans_vel
      MaxWallVel = max(MaxWallVel,norm(wallvel)) 
    enddo
    call MPI_ALLREDUCE(MaxWallVel, maxvel1,1, real_type,MPI_MAX,MPI_COMM_WORLD,ierror)

    this%MaxWallVel= maxvel1
    this%DeltaXmax2= maxval(DEMProperty%Prtcl_PureProp%Radius)*(DEM_opt%Wall_neighbor_ratio+1.0_RK)
#ifndef CFDACM
    this%DeltaXmax1= this%DeltaXmax2-maxval(DEMProperty%Prtcl_PureProp%Radius)
#else
    this%DeltaXmax1= this%DeltaXmax2-maxval(DEMProperty%Prtcl_PureProp%Radius)-maxval(dlub_pw)
#endif
  end subroutine InitContactSearchPW
    
  !**********************************************************************
  ! Init_PWCS_List
  !**********************************************************************
  subroutine Init_PWCS_List(this)
    implicit none
    class(ContactSearchPW) :: this 
    integer::i
       
    Bucket_PWCS = 0
    do i=1,this%max_nearPrtcl_pWall 
      Next_PWCS(i)=-i-1
    enddo
    Next_PWCS(this%max_nearPrtcl_pWall) = 0
    this%NextInsert = 1
  end subroutine Init_PWCS_List
  !**********************************************************************
  ! Reallocate_Bucket
  !********************************************************************** 
  subroutine Reallocate_Bucket(this,nB_new)
    implicit none
    class(ContactSearchPW)::this 
    integer,intent(in):: nB_new

    ! locals
    integer::sizep,sizen
    integer,dimension(:),allocatable:: IntVec

    sizep= this%mHead
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, nB_new+1)
    sizen= min(sizen,DEM_Opt%numPrtcl)
    this%mHead = sizen

    call move_alloc(Bucket_PWCS, IntVec)
    allocate(Bucket_PWCS(sizen))
    Bucket_PWCS(1:sizep)=IntVec
    Bucket_PWCS(sizep+1:sizen)=0     ! Added at 22:37, 2020-11-06, Gong Zheng
    deallocate(IntVec)

  end subroutine Reallocate_Bucket 
    
  !**********************************************************************
  ! Reallocate_PWCS_List
  !**********************************************************************    
  subroutine Reallocate_PWCS_List(this)
    implicit none
    class(ContactSearchPW)::this 
      
    ! locals
    integer::i
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec 
      
    sizep=this%max_nearPrtcl_pWall
    sizen=int(1.2_RK*real(sizep,kind=RK))
    this%max_nearPrtcl_pWall=sizen
      
    call move_alloc(id_Wall,IntVec)
    allocate(id_Wall(sizen))
    id_Wall(1:sizep)=IntVec
      
    call move_alloc(Next_PWCS,IntVec)
    allocate(Next_PWCS(sizen))
    Next_PWCS(1:sizep)=IntVec
    deallocate(IntVec)
    
    do i=sizep+1,sizen
      Next_PWCS(i)=-i-1
    enddo
    Next_PWCS(sizen)=0
    this%NextInsert=sizep+1
     
    call DEMLogInfo%CheckForError(ErrT_Pass,"Reallocate_PWCS_List","Need to reallocate particle variables")
    call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)    
  end subroutine Reallocate_PWCS_List
    
  !**********************************************************************
  ! Finding particles which are near all walls
  !**********************************************************************
  subroutine UpdateNearPrtclsPW(this,iterNumber)
    implicit none
    class(ContactSearchPW) :: this
    integer,intent(in):: iterNumber

    ! locals
    real(RK)::max_v,max_a,t,maxreal 
    integer:: nw,i,numIter,wid,nNear,nextI,ierror
    type(real3):: min_point, max_point, pmin_point, pmax_point

    ! checks if in this iteration the neighbour list should be updated
    if(iterNumber<this%Next_iter_update) return
    call this%Init_PWCS_List()
                
    nNear=0
    nw = DEMGeometry%nPW_local
    do wid = 1,nw
      min_point= DEMGeometry%pWall(wid)%min_point
      max_point= DEMGeometry%pWall(wid)%max_point
      pmin_point = min_point - this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK)
      pmax_point = max_point + this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK)

      do i=1,GPrtcl_list%nlocal
        IF(GPrtcl_PosR(i)%x >= pmin_point%x .and. GPrtcl_PosR(i)%x <= pmax_point%x .and. &
           GPrtcl_PosR(i)%y >= pmin_point%y .and. GPrtcl_PosR(i)%y <= pmax_point%y .and. &
           GPrtcl_PosR(i)%z >= pmin_point%z .and. GPrtcl_PosR(i)%z <= pmax_point%z) THEN
          nNear=nNear+1
          if(this%NextInsert==0) call this%Reallocate_PWCS_List()
          nextI=this%NextInsert
          this%NextInsert= -Next_PWCS(nextI)
          id_Wall(nextI)= wid
          Next_PWCS(nextI)=Bucket_PWCS(i)
          Bucket_PWCS(i)=nextI
        ENDIF
      enddo
    enddo
    this%num_nearPrtcl_pWall = nNear
    
    ! calculating number of iterations which should be performed
    ! untill the next update of the list of particles near the wall  
    
    ! calculating the maximum velocity in the system
    maxreal = 0.0_RK
    do i=1,GPrtcl_list%nlocal
      maxreal = max(maxreal, norm(GPrtcl_linVel(1,i)))
    enddo
    call MPI_ALLREDUCE(maxreal,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
    max_v = max_v + this%MaxWallVel
    max_v = max(max_v, 0.3_RK)
    max_v = max_v*1.5_RK
    
    !calculating the maximum acceleration in the system
    maxreal = 0.0_RK
    do i=1,GPrtcl_list%nlocal
      maxreal = max(maxreal, norm(GPrtcl_linAcc(1,i)))
    enddo
    call MPI_ALLREDUCE(maxreal,max_a,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
    max_a = max(max_a, norm(DEM_Opt%gravity))
    max_a = max_a * 1.5_RK
    
    ! solving the equation of motion for t
    ! or  0.5*a*t^2 + v*t - dx = 0 
    ! This equation has two roots, a positive and a negative
    ! we need the positive root
    if(nrank==0) then
      t = (sqrt(max_v**2+ 2.0_RK*max_a*this%DeltaXmax1)-max_v)/max_a
      numIter = max(int(t/DEM_opt%dt), 1)
      numIter = min(numIter, DEM_opt%Wall_max_update_iter)
    endif
    call MPI_BCAST(numIter,1,int_type,0,MPI_COMM_WORLD,ierror)
    this%Next_iter_update =  iterNumber + numIter
    
    if(nrank==0 .and. numIter<DEM_opt%Wall_max_update_iter ) then
      call DEMLogInfo%OutInfo("Neighbor particles of all walls are updated", 4)
      call DEMLogInfo%OutInfo("This process is repeated in Iteration :"//trim(num2str(DEMContactSearchPW%Next_iter_update)),4)
    endif    
  end subroutine UpdateNearPrtclsPW

  !**********************************************************************
  ! CSPW_copy
  !**********************************************************************
  subroutine CSPW_copy(this,id1,id2)
    implicit none
    class(ContactSearchPW):: this
    integer,intent(in)::id1,id2

    ! locals
    integer::n,NextI
   
    n=Bucket_PWCS(id1)
    do while(n>0)
      NextI = Next_PWCS(n)
      Next_PWCS(n) = -this%nextInsert
      this%nextInsert = n      
      n = NextI
    enddo

    Bucket_PWCS(id1) = Bucket_PWCS(id2)
    Bucket_PWCS(id2) = 0

  end subroutine CSPW_copy

  !**********************************************************************
  ! InsertNearPW
  !**********************************************************************
  subroutine InsertNearPW(this,pid)
    implicit none
    class(ContactSearchPW):: this
    integer,intent(in)::pid

    ! locals
    integer:: nw,wid,nextI
    type(real3):: min_point, max_point, pmin_point, pmax_point

    nw = DEMGeometry%nPW_local
    do wid= 1,nw
      min_point= DEMGeometry%pWall(wid)%min_point
      max_point= DEMGeometry%pWall(wid)%max_point
      pmin_point= min_point- this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK)
      pmax_point= max_point+ this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK) 

      IF(GPrtcl_PosR(pid)%x >=pmin_point%x .and. GPrtcl_PosR(pid)%x <=pmax_point%x .and. &
         GPrtcl_PosR(pid)%y >=pmin_point%y .and. GPrtcl_PosR(pid)%y <=pmax_point%y .and. &
         GPrtcl_PosR(pid)%z >=pmin_point%z .and. GPrtcl_PosR(pid)%z <=pmax_point%z) THEN

        if(this%NextInsert==0) call this%Reallocate_PWCS_List()
        nextI=this%NextInsert
        this%NextInsert= -Next_PWCS(nextI)
        id_Wall(nextI)= wid
        Next_PWCS(nextI)=Bucket_PWCS(pid)
        Bucket_PWCS(pid)=nextI
      ENDIF
    enddo
  end subroutine InsertNearPW

#ifdef CFDACM
  !**********************************************************************
  ! Performing contact search to determine particle-wall contacts 
  !**********************************************************************
  subroutine FindContactsPW(this)
    implicit none
    class(ContactSearchPW):: this

    ! locals
    logical:: Iscntct,clcLubFlag
    integer:: i,nid,wid
    real(RK):: ovrlp
    type(real3):: nv

    ! this is a convention, the particle id should be the first item in the contact pair (particle & wall)
    DO i=1,GPrtcl_list%nlocal
       nid=Bucket_PWCS(i)
       do while(nid>0)
         wid=id_Wall(nid)
         Iscntct= DEMGeometry%pWall(wid)%isInContact(GPrtcl_PosR(i),ovrlp,nv,clcLubFlag)
         if(Iscntct)then
           call GPPW_CntctList%AddContactPW(i,wid,ovrlp,nv)
         elseif(clcLubFlag) then
           ovrlp= -ovrlp
           if(ovrlp<=dlub_pw(GPrtcl_pType(i))) call GPPW_CntctList%AddLubForcePW(i,wid,ovrlp) 
         endif 
         nid=Next_PWCS(nid)
       enddo
    ENDDO
  end subroutine FindContactsPW
#else
  !**********************************************************************
  ! Performing contact search to determine particle-wall contacts 
  !**********************************************************************
  subroutine FindContactsPW(this)
    implicit none
    class(ContactSearchPW):: this

    ! locals
    integer:: i,nid,wid
    real(RK):: ovrlp
    type(real3):: nv

    ! this is a convention, the particle id should be the first item in the contact pair (particle & wall)
    DO i=1,GPrtcl_list%nlocal
       nid=Bucket_PWCS(i)
       do while(nid>0)
         wid=id_Wall(nid)
         if(DEMGeometry%pWall(wid)%isInContact(GPrtcl_PosR(i),ovrlp,nv))then
           call GPPW_CntctList%AddContactPW(i,wid,ovrlp,nv) 
         endif 
         nid=Next_PWCS(nid)
       enddo
    ENDDO
  end subroutine FindContactsPW
#endif
    
end module Prtcl_ContactSearchPW
module Prtcl_Decomp_2d
  use MPI
  use m_TypeDef
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_MeshAndMetries,only:dx,dz
  use m_Parameters,only:p_row,p_col
  use m_Decomp2d,only:nrank,nproc,y1start,y1end
#endif
  implicit none
  private

#if !defined(CFDDEM) && !defined(CFDACM)
  integer,public:: nrank  ! local MPI rank 
  integer,public:: nproc  ! total number of processors
#endif

  integer,public:: int_type,real_type,real3_type,real4_type
  integer,public:: int_byte,real_byte,real3_byte,real4_byte

  TYPE Prtcl_DECOMP_INFO

    ! define neighboring blocks
    ! second dimension 8 neighbour processors:
    !        1:4, 4 edge neighbours; 5:8, 4 cornor neighbours; 0, current processor(if any)
    integer::Prtcl_Pencil
    integer::prow
    integer::pcol
    integer::coord1
    integer::coord2
    integer,dimension(4):: ProcNgh
    real(RK)::xSt,ySt,zSt !min domain of the current and neighbor processors 
    real(RK)::xEd,yEd,zEd !max domain of the current and neighbor processors
  contains
    procedure:: Init_DECOMP => PDI_Init_DECOMP
  end type Prtcl_DECOMP_INFO
  type(Prtcl_DECOMP_INFO),public :: DEM_decomp
contains

  !**********************************************************************
  ! PDI_Init_DECOMP
  !**********************************************************************
#if !defined(CFDDEM) && !defined(CFDACM)
  subroutine PDI_Init_DECOMP(this,chFile,RowCol)
    implicit none
    class(Prtcl_DECOMP_INFO)::this
    character(*),intent(in)::chFile
    integer,intent(in),optional::RowCol(2)
#else
  subroutine PDI_Init_DECOMP(this)
    implicit none
    class(Prtcl_DECOMP_INFO)::this
#endif
   
    ! locals
    type(real3)::pMin,Pmax,SimLen
    integer::row,col,coord1,coord2,idTop,idBottom,idLeft,idRight
#if !defined(CFDDEM) && !defined(CFDACM)
    character::pencil
    real(RK)::xpart,ypart,zpart
    integer::nUnitFile,ierror
    NAMELIST/PrtclDomainDecomp/ row,col,pencil

    open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror )
    if(ierror /= 0 .and. nrank==0) then
      print*, "  PDI_Init_DECOMP: Cannot open file"//trim(adjustl(chFile))
      STOP
    endif
    read(nUnitFile, nml=PrtclDomainDecomp)
    close(nUnitFile,IOSTAT=ierror)
    if(present(RowCol)) then
      row=RowCol(1); col=RowCol(2)
    endif
    
    if(row*col /= nproc .and. nrank==0) then
      STOP "  PDI_Init_DECOMP: Invalid 2D processor grid- nproc/= row*col"
    endif
    if(pencil == "x" .or. pencil == "X") then
      this%Prtcl_Pencil = x_axis
    elseif(pencil == "y" .or. pencil == "Y") then
      this%Prtcl_Pencil = y_axis
    elseif(pencil == "z" .or. pencil == "Z") then
      this%Prtcl_Pencil = z_axis
    else
      if(nrank==0) STOP "  PDI_Init_DECOMP: Invalid pencil sign"
    endif
#else
    row=p_row; col=p_col
    this%Prtcl_Pencil = y_axis    
#endif
    this%prow= row
    this%pcol= col
    call Init_Prctl_MPI_TYPE()
  
    ! ------------------------- neigbor information begins------------------------
    ! 
    ! 2D domain decomposition method in LPT_MPI(from left to right:  x1-pencil, y1-pencil, z1-pencil)
    ! 
    !   If we have 6 processors = 3 row * 2 col
    ! 
    !     the arrangement of the subdomains(nrank) is as follow:
    !       y               x               x
    !       |        4 5    |        4 5    |        4 5   
    !       |        2 3    |        2 3    |        2 3
    !       |_ _ _z  0 1    |_ _ _z  0 1    |_ _ _y  0 1
    !
    !     the arrangement of the coord1 is as follow:
    !       y               x               x
    !       |        2 2    |        2 2    |        2 2
    !       |        1 1    |        1 1    |        1 1
    !       |_ _ _z  0 0    |_ _ _z  0 0    |_ _ _y  0 0
    ! 
    !     the arrangement of the coord2 is as follow:
    !       y               x               x
    !       |        0 1    |        0 1    |        0 1
    !       |        0 1    |        0 1    |        0 1
    !       |_ _ _ z 0 1    |_ _ _z  0 1    |_ _ _y  0 1
    ! 
    !     neighbor index:
    ! 
    !       y               x               x
    !       |        6 3 5  |        6 3 5  |        6 3 5
    !       |        2 0 1  |        2 0 1  |        2 0 1
    !       |_ _ _ z 7 4 8  |_ _ _z  7 4 8  |_ _ _y  7 4 8
    ! 
    !        Here 0 means the center subdomain, and 1-8 stands for the relative location of the eight neighbors

    coord1 = int ( nrank / col)
    coord2 = mod ( nrank,  col)
    this%coord1 = coord1
    this%coord2 = coord2

    ! Firstly, all the boundaries are assumed to be periodic
    idTop     = mod(coord1+1,    row)  ! top
    idBottom  = mod(coord1+row-1,row)  ! bottom 
    idLeft    = mod(coord2+col-1,col)  ! left
    idRight   = mod(coord2+1,    col)  ! right
    this%ProcNgh(1) = coord1   * col + idRight 
    this%ProcNgh(2) = coord1   * col + idLeft  
    this%ProcNgh(3) = idTop    * col + coord2  
    this%ProcNgh(4) = idBottom * col + coord2  
   !this%ProcNgh(5) = idTop    * col + idRight
   !this%ProcNgh(6) = idTop    * col + idLeft
   !this%ProcNgh(7) = idBottom * col + idLeft
   !this%ProcNgh(8) = idBottom * col + idRight

   ! Secondly, modify the edge neighbour ids
    IF(coord1==0) THEN
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(4)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(4)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(4)=MPI_PROC_NULL
    ENDIF
    IF(coord1==row-1) THEN
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(3)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(3)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(3)=MPI_PROC_NULL
    ENDIF
    IF(coord2==0) THEN
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(2)=MPI_PROC_NULL  
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(2)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(2)=MPI_PROC_NULL
    ENDIF
    IF(coord2==col-1) THEN
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(1)=MPI_PROC_NULL 
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(1)=MPI_PROC_NULL  
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(1)=MPI_PROC_NULL
    ENDIF

    pMin = DEM_Opt%SimDomain_min
    pMax = DEM_Opt%SimDomain_max
    SimLen = pMax - pMin
#if !defined(CFDDEM) && !defined(CFDACM)
    if(this%Prtcl_Pencil ==x_axis) then
      ypart = SimLen%y/real(row,kind=RK)
      zpart = SimLen%z/real(col,kind=RK)
      this%xSt= pMin%x
      this%xEd= pMax%x
      this%ySt= pMin%y + ypart*real(coord1,  kind=RK)
      this%yEd= pMin%y + ypart*real(coord1+1,kind=RK)
      this%zSt= pMin%z + zpart*real(coord2,  kind=RK)
      this%zEd= pMin%z + zpart*real(coord2+1,kind=RK)
    elseif(this%Prtcl_Pencil ==y_axis) then
      xpart = SimLen%x/real(row,kind=RK)
      zpart = SimLen%z/real(col,kind=RK)
      this%xSt= pMin%x + xpart*real(coord1,  kind=RK)
      this%xEd= pMin%x + xpart*real(coord1+1,kind=RK)
      this%ySt= pMin%y
      this%yEd= pMax%y
      this%zSt= pMin%z + zpart*real(coord2,  kind=RK)
      this%zEd= pMin%z + zpart*real(coord2+1,kind=RK)
    elseif(this%Prtcl_Pencil ==z_axis) then
      xpart = SimLen%x/real(row,kind=RK)
      ypart = SimLen%y/real(col,kind=RK)
      this%xSt= pMin%x + xpart*real(coord1,  kind=RK)
      this%xEd= pMin%x + xpart*real(coord1+1,kind=RK)
      this%ySt= pMin%y + ypart*real(coord2,  kind=RK)
      this%yEd= pMin%y + ypart*real(coord2+1,kind=RK)
      this%zSt= pMin%z
      this%zEd= pMax%z
    endif
#else
    this%xSt= real(y1start(1)-1, kind=RK)*dx
    this%xEd= real(y1end(1),     kind=RK)*dx
    this%ySt= pMin%y
    this%yEd= pMax%y
    this%zSt= real(y1start(3)-1, kind=RK)*dz
    this%zEd= real(y1end(3),     kind=RK)*dz
#endif
  end subroutine PDI_Init_DECOMP
