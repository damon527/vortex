
  !**********************************************************************
  ! Init_Prctl_MPI_TYPE
  !**********************************************************************
  subroutine Init_Prctl_MPI_TYPE()
    implicit none
    integer::ierror
    integer,dimension(4)::disp,blocklen,blocktype
  
    ! integer
    int_type = MPI_INTEGER
    call MPI_TYPE_SIZE(int_type,   int_byte,   ierror)

    ! real
    if(RK==4) then
      real_type = MPI_REAL
    else
      real_type = MPI_DOUBLE_PRECISION
    endif
    call MPI_TYPE_SIZE(real_type,  real_byte,  ierror)

    ! real3 type
    blocklen(1:3)=1
    blocktype(1:3)=real_type
    disp(1)=0
    disp(2)=disp(1)+real_byte
    disp(3)=disp(2)+real_byte
    call MPI_TYPE_STRUCT(3,blocklen(1:3),disp(1:3),blocktype(1:3),real3_type,ierror)
    call MPI_TYPE_COMMIT(real3_type,ierror)
    call MPI_TYPE_SIZE(real3_type, real3_byte, ierror)

    ! real4 type
    blocklen(1:4)=1
    blocktype(1:4)=real_type
    disp(1)=0
    disp(2)=disp(1)+real_byte
    disp(3)=disp(2)+real_byte
    disp(4)=disp(3)+real_byte
    call MPI_TYPE_STRUCT(4,blocklen(1:4),disp(1:4),blocktype(1:4),real4_type,ierror)
    call MPI_TYPE_COMMIT(real4_type,ierror)
    call MPI_TYPE_SIZE(real4_type, real4_byte, ierror)
  end subroutine Init_Prctl_MPI_TYPE
end module Prtcl_Decomp_2d
module Prtcl_DEMSystem
  use MPI
  use m_Timer
  use m_TypeDef
  use Prtcl_Comm
  use Prtcl_Property
  use Prtcl_Geometry
  use Prtcl_IOAndVisu
  use Prtcl_Variables
  use Prtcl_decomp_2d
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
  use Prtcl_Integration
  use Prtcl_ContactSearch
  use Prtcl_ContactSearchPW
#ifdef CFDDEM
  use Prtcl_DumpPrtcl
  use m_Decomp2d,only: nrank
#endif  
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
   
  integer::iCountDEM
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
    integer:: ierror
    character(256):: chStr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    !// Initializing main log info and visu
    iCountDEM=0
    if(DEM_Opt%RestartFlag) iCountDEM=10
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

    ! Step3: initilize all the particle variables and IO
    call GPrtcl_list%AllocateAllVar()
    call DEM_IO%Init_visu(chDEMPrm,1)
    t_restart1=MPI_WTIME()
#ifdef CFDDEM
    if(.not.DEM_Opt%RestartFlag) then
      if(DEM_Opt%numPrtcl>0) call DEM_IO%ReadInitialCoord()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Initial Particle coordinates are READING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_Opt%numPrtcl)),2)
      endif
      DEM_Opt%np_InDomain = DEM_Opt%numPrtcl
      if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
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
      if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedRestart()
    endif
#else
    if(.not.DEM_Opt%RestartFlag) then
      call GPrtcl_list%MakingAllPrtcl(chDEMPrm)
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are MAKING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_Opt%numPrtcl)),2)
      endif
      DEM_Opt%np_InDomain = DEM_Opt%numPrtcl
    else
      call DEM_IO%Read_Restart()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(DEM_Opt%np_InDomain)),2)
      endif
    endif
    if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
#endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1

    call DEM_IO%Init_visu(chDEMPrm,2)
#ifdef CFDDEM
    call Initialize_DumpPrtcl(chDEMPrm)
#endif

    ! Step4: initialize the inter-processors communication
    call DEM_Comm%InitComm()
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    endif

    ! Step5: Initializing contact list and contact force 
    call  GPPW_CntctList%InitContactList()
    t_restart1=MPI_WTIME()
    if(DEM_Opt%RestartFlag) then
      if(DEM_Opt%np_InDomain>0)call DEM_IO%RestartCL()
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1+t_res_tot
    if(nrank==0 .and. DEM_Opt%RestartFlag) call DEMLogInfo%OutInfo("Restart time [sec] :"//trim(num2str(t_res_tot)),2)
    
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step5: Initializing contact list and contact force models . . . ", 1 )
      if( DEM_Opt%CF_Type == DEM_LSD ) then
        write(chStr,"(A)") "linear spring-dashpot with limited tangential displacement"
      elseif( DEM_Opt%CF_Type == DEM_nLin ) then
        write(chStr,"(A)") "non-linear visco-elastic model with limited tangential displacement"
      endif
      call DEMLogInfo%OutInfo("Contact force model is "//trim(chStr), 2 )

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

    ! Step6: Initializing contact search method
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step6: Initializing contact search method . . . ", 1)
      call DEMLogInfo%OutInfo("Particle-Particle contact search intialization...",2)
      call DEMLogInfo%OutInfo("Particle-Wall contact search intialization...",2)
    endif
    call DEMContactSearch%InitContactSearch()
    call DEMContactSearchPW%InitContactSearchPW()
    
    ! Step7: timers for recording the execution time of different parts of program
    if(nrank==0) call DEMLogInfo%OutInfo("Step7: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_pre_iter_timer%reset()
    call this%m_comm_cs_timer%reset()
    call this%m_CSCF_PP_timer%reset()
    call this%m_CSCF_PW_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_write_prtcl_timer%reset()
    call this%m_comm_exchange_timer%reset()
    
#ifdef CFDDEM
    call DEM_IO%dump_visu((DEM_Opt%ifirst-1)/icouple)            
#else
    call DEM_IO%dump_visu(DEM_Opt%ifirst-1)
#endif
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
    character(256):: chLine
    integer:: Consv_Cont(2),Consv_Cont1(2),ierror,npwcs(4)

#ifdef CFDDEM
  IF(UpdateDEMflag) THEN
#endif
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

    ! calculate linear and angular accelerations, position and velocities 
    call this%m_integration_timer%start()
    iCountDEM=iCountDEM+1
    call Prtcl_Integrate(iCountDEM)
    call this%m_integration_timer%finish()
   
    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call DEM_Comm%Comm_For_Exchange()
    call GPPW_CntctList%RemvReleased()
    call this%m_comm_exchange_timer%finish()
#ifdef CFDDEM
  ENDIF
#endif
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, DEM_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
#ifdef CFDDEM
    if( mod(this%IterNumber,DEM_Opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime/icouple)
    if( mod(this%IterNumber,DumpPrtclFreq)== 0)      call WriteDumpCache(itime)
#else
    if( mod(this%IterNumber,DEM_Opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime)
#endif
    if( mod(this%IterNumber,DEM_Opt%BackupFreq)== 0 .or. itime==DEM_Opt%ilast) then
      call DEM_IO%Write_Restart(itime)
#ifdef CFDDEM
      call DEM_IO%WriteFixedRestart(itime)
      call PrtclVarDump(itime)
#endif
      call DEM_IO%Delete_Prev_Restart(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()

#ifdef CFDDEM
  IF(UpdateDEMflag) THEN
#endif    
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
#ifdef CFDDEM
  ENDIF
#endif
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
#ifdef CFDDEM
    NAMELIST/CFDDEMCoupling/icouple,UpdateDEMflag,is_clc_Lift,is_clc_Basset,is_clc_Basset_fixed,is_clc_ViscousForce,&
                            is_clc_PressureGradient,is_clc_ViscousForce,is_clc_FluidAcc,FluidAccCoe,SaffmanConst,   &
                            RatioSR,IsAddFluidPressureGradient
    NAMELIST/BassetOptions/ mWinBasset, mTailBasset, BassetAccuracy, BassetTailType
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
           
    Base_wall_id   = DEM_Opt%Base_wall_id 
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
#ifdef CFDDEM
    write(DEMLogInfo%nUnit, nml=CFDDEMCoupling)
    write(DEMLogInfo%nUnit, nml=BassetOptions)
#endif
  end subroutine Write_DEM_Opt_to_Log

end module Prtcl_DEMSystem
module Prtcl_DumpPrtcl
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
  use m_Decomp2d,only: nrank
  implicit none
  private
#define RKP_Dump 4
  
  character(128)::DumpPrtclDir
  logical::DumpPrtclFlag,ResetDumpFlag
  integer,  dimension(:,:),allocatable::DumpInteMat
  real(RKP_Dump),dimension(:,:),allocatable::DumpRealMat
  integer::nDumpPrtclSize,mDumpPrtclSize,DumpPrtclFreq,iDump
  
  public::DumpPrtclFreq,Initialize_DumpPrtcl,WriteDumpCache,PrtclVarDump
contains
!#define OnlyDumpFpForce
#define nDumpPrtclInte 4

#ifdef OnlyDumpFpForce
#define nDumpPrtclReal 9
#else
#ifdef CFDACM
#define nDumpPrtclReal 21
#else
#define nDumpPrtclReal 9
#endif
#endif

#define Prtcl_Dump_Flag 93

  !**********************************************************************
  ! Initialize_DumpPrtcl
  !**********************************************************************
  subroutine Initialize_DumpPrtcl(chFile)
    implicit none
    character(*),intent(in)::chFile

    ! locals
    real(RK)::yDump
    character(256)::chStr
    integer:: pid,nUnitFile,ierror,ierrTmp
    namelist/DumpPrtclOptions/DumpPrtclFlag,ResetDumpFlag,yDump,mDumpPrtclSize,DumpPrtclDir,DumpPrtclFreq
  
    yDump=0.0_RK
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","Cannot open file:"//trim(adjustl(chFile)))
    read(nUnitFile, nml=DumpPrtclOptions)
    close(nUnitFile,IOSTAT=ierror)

    if(.not. DumpPrtclFlag) return
    if(ResetDumpFlag) then
      do pid=1,GPrtcl_list%nlocal
        if(GPrtcl_PosR(pid)%y >= yDump) then
          GPrtcl_usrMark(pid)=Prtcl_Dump_Flag
        else
          GPrtcl_usrMark(pid)=1
        endif
      enddo
    endif

    nDumpPrtclSize=0; iDump=0; ierror=0
    if(mDumpPrtclSize<10000 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","So small mDumpPrtclSize:"//trim(num2str(mDumpPrtclSize)))
    endif
    allocate(DumpInteMat(nDumpPrtclInte,mDumpPrtclSize),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(DumpRealMat(nDumpPrtclReal,mDumpPrtclSize),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","Allocation failed")

    write(chStr,"(A)") 'mkdir -p '//trim(adjustl(DumpPrtclDir))//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))

#ifdef OnlyDumpFpForce
    if(nrank==0) call DEMLogInfo%OutInfo("Choose to only dump Fluid-particle force",2)
#else
    if(nrank==0) call DEMLogInfo%OutInfo("Choose to dump full particle information",2)
#endif
  end subroutine Initialize_DumpPrtcl

  !**********************************************************************
  ! WriteDumpCache
  !**********************************************************************
  subroutine WriteDumpCache(itime)
    implicit none
    integer,intent(in)::itime

    ! locals
    integer::k,pid,nlocal

    if(.not. DumpPrtclFlag) return
    nlocal=GPrtcl_list%nlocal
    do pid=1,nlocal
      if(GPrtcl_usrMark(pid)/=Prtcl_Dump_Flag) cycle
      k=1
      nDumpPrtclSize=nDumpPrtclSize+1
      DumpInteMat(1,nDumpPrtclSize)=itime
      DumpInteMat(2,nDumpPrtclSize)=GPrtcl_id(pid)
      DumpInteMat(3,nDumpPrtclSize)=GPrtcl_pType(pid)
      DumpInteMat(4,nDumpPrtclSize)=GPPW_CntctList%IsCntct(pid)
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%x, RKP_Dump);       k=k+1 ! 01
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%y, RKP_Dump);       k=k+1 ! 02         
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%z, RKP_Dump);       k=k+1 ! 03         
#ifdef OnlyDumpFpForce
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%x, RKP_Dump);    k=k+1 ! 04 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%y, RKP_Dump);    k=k+1 ! 05 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%z, RKP_Dump);    k=k+1 ! 06       
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%x, RKP_Dump);   k=k+1 ! 07
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%y, RKP_Dump);   k=k+1 ! 08
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%z, RKP_Dump);   k=k+1 ! 09
#else
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%x, RKP_Dump);   k=k+1 ! 04 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%y, RKP_Dump);   k=k+1 ! 05      
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%z, RKP_Dump);   k=k+1 ! 06      
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%x, RKP_Dump);   k=k+1 ! 07     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%y, RKP_Dump);   k=k+1 ! 08    
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%z, RKP_Dump);   k=k+1 ! 09    
#ifdef CFDACM     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%x, RKP_Dump);    k=k+1 ! 10     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%y, RKP_Dump);    k=k+1 ! 11 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%z, RKP_Dump);    k=k+1 ! 12       
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%x, RKP_Dump);   k=k+1 ! 13
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%y, RKP_Dump);   k=k+1 ! 14
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%z, RKP_Dump);   k=k+1 ! 15
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%x, RKP_Dump); k=k+1 ! 16
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%y, RKP_Dump); k=k+1 ! 17
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%z, RKP_Dump); k=k+1 ! 18
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%x, RKP_Dump);     k=k+1 ! 19
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%y, RKP_Dump);     k=k+1 ! 20
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%z, RKP_Dump);     k=k+1 ! 21
#endif
#endif
      if(nDumpPrtclSize==mDumpPrtclSize) then
        call PrtclVarDump(itime)
        iDump=iDump+1
      endif
    enddo
    iDump=0
  end subroutine WriteDumpCache

  !**********************************************************************
  ! PrtclVarDump
  !**********************************************************************
  subroutine PrtclVarDump(itime)
    implicit none
    integer,intent(in)::itime

    ! locals
    integer::ierror,nUnit
    character(128)::chFile

    if(.not.DumpPrtclFlag)return
    if(nDumpPrtclSize==0) return
    write(chFile,'(A,I5.5,A,I10.10,A,I2.2)')trim(DumpPrtclDir)//'rank',nrank,'_',itime/icouple,'_',iDump

    open(newunit=nUnit,file=trim(chFile),status='replace',form='unformatted',access='stream',IOSTAT=ierror)
    IF(ierror/=0) THEN
      call DEMLogInfo%CheckForError(ErrT_Pass,"PrtclVarDump","Cannot open file: "//trim(chFile))
    ELSE
      write(nUnit)nDumpPrtclSize ! Added by Zheng Gong, 2023-05-04
      write(nUnit)DumpInteMat(:,1:nDumpPrtclSize)
      write(nUnit)DumpRealMat(:,1:nDumpPrtclSize)
    ENDIF
    close(nUnit,IOSTAT=ierror)

    nDumpPrtclSize=0
  end subroutine PrtclVarDump
end module Prtcl_DumpPrtcl

#undef nDumpPrtclInte
#undef nDumpPrtclReal
#undef Prtcl_Dump_Flag

#ifdef OnlyDumpFpForce
#undef OnlyDumpFpForce
#endif

#undef RKP_Dump
!********************************************************************
! 
!  Purpose::
!    1) Provide a much cut down module to partitinon the   
!       surface of a sphere into regions of equal area.
! 
!  References:
!    1) Leopardi P. A partition of the unit sphere into regions of 
!       equal area and small diameter[J]. Electronic Transactions on 
!       Numerical Analysis Etna, 2006, 25(1):2006.
!    2) http://eqsp.sourceforge.net/
!    3) https://sourceforge.net/projects/eqsp/
! 
!  Author: Zheng Gong
!  Date:   2021-08-13
! 
!********************************************************************
module Prtcl_EqualSphere
  use m_TypeDef,only:RK
  implicit none
  private
  public::eq_sphere
contains
  !******************************************************************
  ! eq_sphere
  !******************************************************************
  subroutine eq_sphere(Point)
    implicit none
    real(RK),dimension(:,:),intent(out)::Point

    ! locals
    integer,dimension(:),allocatable::n_regions
    integer::k,m,nMarker,nColumn,n_collars,nCount,n_top,n_bot
    real(RK),parameter::PI=3.141592653589793238462643383279502884_RK
    real(RK)::rx,ry,rz,rnorm,r_regions,a_top,a_bot,area_tot,Psi,Phi,aTemp
    real(RK)::area_cap,c_polar,discrepancy,a_fitting,area_top,area_bot,offset

    nMarker=size(Point,1)
    nColumn=size(Point,2)
    if(nMarker<1 .or. nColumn/=3) then
      print*,"Error in eq_sphere, nMarker or nColumn wrong:",nMarker,nColumn; stop
    endif
    Point(1,:)=(/0.0_RK, 0.0_RK, 1.0_RK/);       if(nMarker==1)return
    Point(2,:)=(/1.0_RK, 0.0_RK, 0.0_RK/); 
    Point(nMarker,:)=(/0.0_RK, 0.0_RK,-1.0_RK/); if(nMarker< 4)return

    area_cap=4.0_RK*PI/real(nMarker,RK)
    c_polar =2.0_RK*asin(sqrt(area_cap/(4.0_RK*PI)))
    n_collars= max(1,nint((PI-2.0_RK*c_polar)/sqrt(area_cap)))
    allocate(n_regions(n_collars+2))
    n_regions(1)=1;n_regions(n_collars+2)=1
    discrepancy = 0.0_RK; area_top=area_cap
    a_fitting = (PI-2.0_RK*c_polar)/real(n_collars,RK)
    do k=1,n_collars
      area_bot=c_polar+real(k,RK)*a_fitting
      area_bot=sin(0.5_RK*area_bot)
      area_bot=4*PI*area_bot*area_bot
      r_regions=(area_bot-area_top)/area_cap
      area_top=area_bot
      n_regions(k+1)=nint(r_regions+discrepancy)
      discrepancy=discrepancy+r_regions-n_regions(k+1)
    enddo

    nCount= 2; offset= 0.0_RK
    a_top = c_polar; area_tot=area_cap
    do k=1,n_collars
      n_top=n_regions(k+1)
      n_bot=n_regions(k+2)
      area_tot=area_tot+n_top*area_cap
      a_bot= 2.0_RK*asin(sqrt(area_tot/(4.0_RK*PI)))
      Psi= 0.5_RK*(a_top+a_bot); a_top=a_bot
      do m=1,n_top
        aTemp=real(2*m-1,RK)*PI/real(n_top,RK) +2.0_RK*PI*offset
        Phi = aTemp -2.0*PI*floor(aTemp/(2.0*PI))
        rx=sin(Psi)*cos(Phi)
        ry=sin(Psi)*sin(Phi)
        rz=cos(Psi)
        rnorm=sqrt(rx*rx+ry*ry+rz*rz)
        Point(nCount,:)=(/rx,ry,rz/)/rnorm
        nCount=nCount+1
      enddo
      offset=offset+real(n_top-n_bot+gcd(n_top,n_bot),RK)/real(2*n_top*n_bot,RK)
      offset=offset-floor(offset)
    enddo
    if(nMarker /= nCount) then
      print*,'Error in eq_sphere: nMarker /= nCount',nMarker,nCount; stop
    endif
  end subroutine eq_sphere

  !******************************************************************
  ! gcd: Greatest Common Divisor
  !******************************************************************
  function gcd(m1,n1) result(r)
    implicit none
    integer,intent(in)::m1,n1
    integer::m,n,r
    m=m1; n=n1
    do 
      r=mod(m,n); m=n; n=r
      if(r==0)exit
    enddo
    r=m
  end function gcd
end module Prtcl_EqualSphere

#ifdef test_EqualSphere
!******************************************************************
! Main program
!******************************************************************
program main
  use Prtcl_EqualSphere
  implicit none
  integer::k,nMarker
  real(8),dimension(:,:),allocatable::Point

  print*,'nMarker='
  read*,nMarker
  allocate(Point(nMarker,3))
  call eq_sphere(Point)
  print*,'    nid,      x,      y,      z'
  do k=1,nMarker
    print*,k,Point(k,:)
  enddo
end program main
#endif
module prtcl_Geometry
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Property
  use Prtcl_decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private
    
  integer:: MaxWallSize=6
  type(real3):: pmin, pmax             ! global domain
  type(real3):: pmin_local, pmax_local ! local domain

  type PlaneWall
    type(real3):: P1       ! first point   
    type(real3):: P2       ! second point  
    type(real3):: P3       ! third point 
    type(real3):: P4       ! fourth point
    type(real3):: trans_vel = zero_r3 ! translational velocity
    integer:: user_id      ! user supplied wall id
    integer:: wall_Type    ! property type of wall material
        
    logical:: bothSide     ! checking if both sides are active 
    logical:: isInfinite
    
    integer:: wall_id      ! program generated wall id
    real (RK):: d          ! d in the implicit equation: ax+by+cz+d = 0   
    type(real3):: n        ! normal vector
    type(real3):: min_point
    type(real3):: max_point
  contains
    procedure:: isInContact => PW_isInContact
    procedure:: IsInDomain  => PW_IsInDomain
  end type PlaneWall 
    
  type Geometry
    integer:: num_pWall                              ! number of plane walls (total)
    integer:: nPW_local                              ! number of plane walls (local) 
    type(PlaneWall),allocatable,dimension(:):: pWall ! a vector that stores all plane walls 
  contains
    procedure:: InitAllocate  => G_InitAllocate
    procedure:: add_PlaneWall => G_add_PlaneWall
    procedure:: MakeGeometry  => G_MakeGeometry
  end type Geometry
  type(Geometry),public :: DEMGeometry
contains

#ifndef CFDACM
  !**********************************************************************
  ! determining if an sphere (box) has a contact with this wall
  !**********************************************************************
  logical function PW_isInContact(this, box,ovrlp,nv )
    implicit none
    class(PlaneWall)::this
    type(real4),intent(in):: box
    real(RK),intent(out):: ovrlp
    type(real3),intent(out)::nv
        
    !locals
    real(RK)::Radius, dist,t
    type(real3):: p,cp
    
    PW_isInContact = .false.
    Radius = box%w
    p = box

    dist = (this%n .dot. p) + this%d
    if(this%bothside) then
      ovrlp = Radius - abs(dist)
    else
      if( dist < 0.0_RK ) return
      ovrlp = Radius - dist
    end if
    if(ovrlp<0.0_RK) return
        
    IF(this%isInfinite) THEN
      PW_isInContact = .true.
      nv = sign(1.0_RK, dist) *this%n
      return
    ELSE
      t = -dist
      cp = t * this%n + p
      if(PW_IsInPlane(this%p1,this%p2,this%p3, this%p4,cp)) then
        PW_isInContact = .true.
        nv = sign(1.0_RK, dist) *this%n
        return
      endif
               
      if(Line_point_check(this%p1,this%p2, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p2,this%p3, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
      if(Line_point_check(this%p3,this%p4, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p4,this%p1, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
    ENDIF
  end function PW_isInContact

#else
  !**********************************************************************
  ! determining if an sphere (box) has a contact with this wall
  !**********************************************************************
  logical function PW_isInContact(this, box,ovrlp,nv,clcLubFlag )
    implicit none
    class(PlaneWall)::this
    type(real4),intent(in):: box
    real(RK),intent(out):: ovrlp
    type(real3),intent(out)::nv
    logical,intent(out)::clcLubFlag
        
    !locals
    type(real3):: p,cp
    real(RK)::Radius, dist,t
    
    clcLubFlag= .false.
    PW_isInContact = .false.
    Radius = box%w
    p = box

    dist = (this%n .dot. p) + this%d
    if(this%bothside) then
      ovrlp = Radius - abs(dist)
    else
      if(dist < 0.0_RK) return
      ovrlp = Radius - dist
    endif
    if(ovrlp<0.0_RK) then
      IF(this%isInfinite) THEN
        clcLubFlag= .true.
      ELSE
        t = -dist
        cp= t * this%n + p
        if(PW_IsInPlane(this%p1,this%p2,this%p3, this%p4,cp)) then
          clcLubFlag= .true.
        endif
      ENDIF
      return
    endif
        
    IF(this%isInfinite) THEN
      PW_isInContact = .true.
      nv = sign(1.0_RK, dist) *this%n
      return
    ELSE
      t  = -dist
      cp = t * this%n + p
      if(PW_IsInPlane(this%p1,this%p2,this%p3, this%p4,cp)) then
        PW_isInContact = .true.
        nv = sign(1.0_RK, dist) *this%n
        return
      endif
               
      if(Line_point_check(this%p1,this%p2, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p2,this%p3, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
      if(Line_point_check(this%p3,this%p4, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p4,this%p1, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
    ENDIF
  end function PW_isInContact
#endif

  !**********************************************************************
  ! checking if the plane have some overlpa region with the domain
  !**********************************************************************
  function PW_IsInDomain(this,dpmin, dpmax) result(res)
    implicit none
    class(PlaneWall)::this
    type(real3),intent(in):: dpmin, dpmax
    logical:: res
      
    ! locals
    real(RK)::t
    type(real3),dimension(8):: DomPoint
    integer,dimension(8)::intarr
    integer::i,sumintarr
    type(real3)::fp1,fp2,fp3,fp4,sp1,sp2,sp3,sp4
      
    res= .false.
    DomPoint(1)=real3(dpmin%x, dpmin%y, dpmin%z)
    DomPoint(2)=real3(dpmax%x, dpmin%y, dpmin%z)
    DomPoint(3)=real3(dpmax%x, dpmax%y, dpmin%z)
    DomPoint(4)=real3(dpmin%x, dpmax%y, dpmin%z)
    DomPoint(5)=real3(dpmin%x, dpmin%y, dpmax%z)
    DomPoint(6)=real3(dpmax%x, dpmin%y, dpmax%z)
    DomPoint(7)=real3(dpmax%x, dpmax%y, dpmax%z)
    DomPoint(8)=real3(dpmin%x, dpmax%y, dpmax%z)
  
    ! Firstly, check whether there are some points within the domain or not
    if( PointIsInDomain(this%p1, dpmin, dpmax) .or. PointIsInDomain(this%p2, dpmin, dpmax)  .or. &
      PointIsInDomain(this%p3, dpmin, dpmax) .or. PointIsInDomain(this%p4, dpmin, dpmax)) then
      res= .true.
      return
    endif
         
    ! Secondly, check whether the reflections overlap or not
    ! reflection in x-y plane
    fp1= real3(this%p1%x,this%p1%y,0.0_RK); fp2= real3(this%p2%x,this%p2%y,0.0_RK)
    fp3= real3(this%p3%x,this%p3%y,0.0_RK); fp4= real3(this%p4%x,this%p4%y,0.0_RK)
    sp1= real3(dpmin%x,  dpmin%y,  0.0_RK); sp2= real3(dpmax%x,  dpmin%y,  0.0_RK)
    sp3= real3(dpmax%x,  dpmax%y,  0.0_RK); sp4= real3(dpmin%x,  dpmax%y,  0.0_RK) 
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return
      
    ! reflection in x-z plane
    fp1= real3(this%p1%x,0.0_RK,this%p1%z); fp2= real3(this%p2%x,0.0_RK,this%p2%z)
    fp3= real3(this%p3%x,0.0_RK,this%p3%z); fp4= real3(this%p4%x,0.0_RK,this%p4%z) 
    sp1= real3(dpmin%x,  0.0_RK,dpmin%z);   sp2= real3(dpmin%x,  0.0_RK,dpmax%z)
    sp3= real3(dpmax%x,  0.0_RK,dpmax%z);   sp4= real3(dpmax%x,  0.0_RK,dpmin%z)
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return
      
    ! reflection in y-z plane
    fp1= real3(0.0_RK,this%p1%y,this%p1%z); fp2= real3(0.0_RK,this%p2%y,this%p2%z)
    fp3= real3(0.0_RK,this%p3%y,this%p3%z); fp4= real3(0.0_RK,this%p4%y,this%p4%z) 
    sp1= real3(0.0_RK,dpmin%y,  dpmin%z);   sp2= real3(0.0_RK,dpmax%y,  dpmin%z)
    sp3= real3(0.0_RK,dpmax%y,  dpmax%z);   sp4= real3(0.0_RK,dpmin%y,  dpmax%z)
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return 
      
    ! Thirdly, check whether the domain is totally located in one side of the plane 
    intarr=-99
    do i=1,8
      t = ((this%n .dot. DomPoint(i))+ this%d )
      if(t>1.0E-8_RK) then
        intarr(i)=1
      elseif ( t< -1.0E-8_RK ) then
        intarr(i)=0
      endif
    enddo
    sumintarr=sum(intarr)
    if(sumintarr==0 .or. sumintarr==8) return

    res= .true.
  end function PW_IsInDomain
  
  !**********************************************************************
  ! checking if the point lays within the boundaries of the plane
  !**********************************************************************
  logical function PW_IsInPlane( p1,p2,p3,p4, cp )
    implicit none
  type(real3),intent(in)::p1,p2,p3,p4, cp
        
    !// locals
    real(RK):: p1p3, p2p4,p1p4,p2p3,p1p2,p3p4,p2p2
    type(real3) p1p, p2p, p3p, p4p
        
    !// body    
    ! Here (p1,p2,p3,p4) are the four point of the plane.
    PW_IsInPlane = .false.    
    p1p = P1-cp
  p2p = P2-cp
  p3p = P3-cp
    p4p = P4-cp
        
    ! first condition u.w<0
    ! u.w = [(p1-p)x(p2-p)].[(p3-p)x(p4-p)] = (p1p.p3p)(p2p.p4p) - (p1p.p4p)(p2p.p3p)
  p1p3 = p1p .dot. p3p
    p2p4 = p2p .dot. p4p
    p1p4 = p1p .dot. p4p
    p2p3 = p2p .dot. p3p
  if(p1p3*p2p4-p1p4*p2p3<0.0_RK) return
  
  ! second condition v.x < 0
    ! v.x = [(p2-p)x(p3-p)].[(p4-p)x(p1-p)] = (p1p.p3p)(p2p.p4p) - (p1p.p2p)(p3p.p4p)
  p1p2 = p1p .dot. p2p
  p3p4 = p3p .dot. p4p
  if(p1p3*p2p4-p1p2*p3p4<0.0_RK) return

  ! third condition u.v < 0
  ! u.v = [(p1-p)x(p2-p)].[(p2-p)x(p3-p)] = (p1p.p2p)(p2p.p3p) - (p1p.p3p)(p2p.p2p)
  p2p2 = p2p .dot. p2p
  if(p1p2*p2p3-p1p3*p2p2<0.0_RK) return 
  
    PW_IsInPlane = .true.
  end function PW_IsInPlane  
    
  !**********************************************************************
  ! checking if the point lays within the line
  !**********************************************************************
  logical function Line_point_check( lp1,lp2, dpos, nv,ovrlp )
    implicit none
    type(real3),intent(in) :: lp1,lp2
    type(real4),intent(in) :: dpos
    type(real3),intent(out):: nv
    real(RK),intent(out)::ovrlp
    
    real(RK):: t, r,length
    type(real3):: w,v,pos,cp
    
    pos = dpos
    w = pos-lp1 
    v = lp2-lp1
    length = (lp1.dist.lp2)
    r = dpos%w
    t = (w.dot.v )/(v.dot.v)
    
    Line_point_check = .false.
    if( t>= 0.0_RK .and. t<= 1.0_RK ) then
      cp = (v * t) + lp1
    elseif( t >= (-r/length)  .and. t <0.0_RK )then
      cp = lp1
    elseif( t> 1.0_RK  .and. t>= (1.0_RK+r/length) ) then
      cp = lp2
    else
      Line_point_check = .false.
      return
    endif
    
    ovrlp = r - (cp .dist. pos)
    if( ovrlp >= 0.0_RK )then
      nv = (pos .nv. cp)   
      Line_point_check = .true.
      return 
    endif
  end function Line_point_check

  !**********************************************************************
  ! checking whether the point lays within the domain
  !**********************************************************************    
  function PointIsInDomain(point, dpmin, dpmax) result(res)
    implicit none
    type(real3),intent(in)::point, dpmin, dpmax
    logical::res
      
    !locals
    real(RK)::realeps=1.000E-10_RK
      
    ! here I slightly expand the domain, to make sure the point  on the domain surface can be regarded as "PointIsInDomain"
    res = .false.
    if(point%x<dpmin%x-realeps .or.  point%x>dpmax%x+realeps) return
    if(point%y<dpmin%y-realeps .or.  point%y>dpmax%y+realeps) return
    if(point%z<dpmin%z-realeps .or.  point%z>dpmax%z+realeps) return
    res = .true.
  end function PointIsInDomain
    
  !**********************************************************************
  ! checking whether the two reflection plane have some overlap region
  !**********************************************************************
  function IsReflectionOvlp(fp1,fp2,fp3,fp4,sp1,sp2,sp3,sp4) result(res)
    implicit none
    type(real3),intent(in)::fp1,fp2,fp3,fp4,sp1,sp2,sp3,sp4
    logical::res
       
    ! Here (fp1,fp2,fp3,fp4) represents first plane
    ! And (sp1,sp2,sp3,sp4) represent  second plane
    res = .true.
       
    ! Firstly, check whether there are some points from one plane located in the inner region of the other plane  
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp1)) return
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp2)) return
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp3)) return
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp4)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp1)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp2)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp3)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp4)) return
       
    ! secondly, check whether there are some inner insection points among the lines from different plane
    if(IsLineInsertInner(fp1,fp2,sp1,sp2)) return
    if(IsLineInsertInner(fp1,fp2,sp2,sp3)) return
    if(IsLineInsertInner(fp1,fp2,sp3,sp4)) return
    if(IsLineInsertInner(fp1,fp2,sp4,sp1)) return
       
    if(IsLineInsertInner(fp2,fp3,sp1,sp2)) return
    if(IsLineInsertInner(fp2,fp3,sp2,sp3)) return
    if(IsLineInsertInner(fp2,fp3,sp3,sp4)) return
    if(IsLineInsertInner(fp2,fp3,sp4,sp1)) return
     
    if(IsLineInsertInner(fp3,fp4,sp1,sp2)) return
    if(IsLineInsertInner(fp3,fp4,sp2,sp3)) return
    if(IsLineInsertInner(fp3,fp4,sp3,sp4)) return
    if(IsLineInsertInner(fp3,fp4,sp4,sp1)) return

    if(IsLineInsertInner(fp4,fp1,sp1,sp2)) return
    if(IsLineInsertInner(fp4,fp1,sp2,sp3)) return
    if(IsLineInsertInner(fp4,fp1,sp3,sp4)) return
    if(IsLineInsertInner(fp4,fp1,sp4,sp1)) return
    res = .false.
        
  end function IsReflectionOvlp
    
  !**********************************************************************
  ! checking whether the two line have inner insertion
  !**********************************************************************    
  function IsLineInsertInner(L1_p1, L1_p2, L2_p1, L2_p2) result(res)
    implicit none
    type(real3),intent(in)::L1_p1, L1_p2, L2_p1, L2_p2
    logical :: res
      
    ! locals
    real(RK)::a,b,c,d,e,dist,t1,t2
    type(real3):: v1,v2,L21_p1,v1_crs_v2, L21_crs_v1,L21_crs_v2,cp
      
    ! Here I use L1_p1 and L1_p2 to express the starting and ending points of the line L1
    !            L2_p1 and L2_p2 to express the starting and ending points of the line L2
    ! And
    !       direction vector of line L1 : v1 = L1_p2 - L1_p1
    !       direction vector of line L2 : v2 = L2_p2 - L2_p1
    ! So any point on line L1 and L2 can be determined by the following parameter equation:
    !       line L1:   P1 = L1_p1 + a * v1, where 'a' is the parameter
    !       line L2:   P2 = L2_p1 + b * v2, where 'b' is the parameter
    ! The inner insert point of line L1 and line L2 should satisfy the following expression:
    !       P_insertion = P1 = P2
    ! i.e.
    !       L1_p1 + a * v1  = L2_p1 + b * v2                                 (1)
    ! We use the cross production of v1 to both sides of Eq.(1):
    !       (L1_p1 + a * v1) x v1 = (L2_p1 + b * v2) x v1
    ! i.e.
    !       (L1_p1 - L2_p1) x v1 =  b* v2 x v1                               (2)    
    ! S
    !       b = sign_b*norm((L2_p1 - L1_p1) x v1 ) /norm( v1 x v2 )          (3)
    ! Where
    !       sign_b= 1,  if ((L2_p1 - L1_p1) x v1 ) .dot. ( v1 x v2 ) )>0     (4)
    !             =-1,  if ((L2_p1 - L1_p1) x v1 ) .dot. ( v1 x v2 ) )<0
    ! Similarly, we can get  that:
    !       a = sign_a*norm((L2_p1 - L1_p1) x v2 ) / norm( v1 x v2 )         (5)
    !       sign_a= 1,  if ((L2_p1 - L1_p1) x v2 ) .dot. ( v1 x v2 ) )>0
    !             =-1,  if ((L2_p1 - L1_p1) x v2 ) .dot. ( v1 x v2 ) )<0
    ! If    0 =<a <= 1, and 0 =<b <= 1, L1 and L2 have inner insection point.
    ! In this function, I use the following notes:
    !     L21_p1 = L2_p1 - L1_p1
    !     v1_crs_v2 = v1 x v2
    !     L21_crs_v1 = (L2_p1 - L1_p1) x v1
    !     L21_crs_v2 = (L2_p1 - L1_p1) x v2
    !     c = norm( (L2_p1 - L1_p1) x v1 ) =norm(L21_crs_v1)
    !     d = norm( (L2_p1 - L1_p1) x v2 ) =norm(L21_crs_v2)
    !     e = norm( v1 x v2 ) =norm(v1_crs_v2)
      
    v1 = L1_p2 - L1_p1
    v2 = L2_p2 - L2_p1
    L21_p1 = L2_p1 - L1_p1
      
    v1_crs_v2 = v1 .cross. v2
    L21_crs_v1= L21_p1 .cross. v1
    L21_crs_v2= L21_p1 .cross. v2
      
    c= norm(L21_crs_v1)
    d= norm(L21_crs_v2)
    e= norm(v1_crs_v2)
      
    ! L1 and L2 are parallel or they two are on the same line.
    IF(e<1.0E-10_RK) THEN
      t1 = (L21_p1 .dot. v1 )/ (v1 .dot. v1)
      cp = (v1 * t1) + L1_p1
      dist = cp .dist. L2_p1  ! the distance between two lines
      if(dist >1.0E-6_RK) then
        res = .false.
      else
        t2 = ((L2_p2-L1_p1).dot.v1)/(v1.dot.v1)
        if((t1<0.0_RK.and.t2<0.0_RK) .or. (t1>1.0_RK .and. t2>1.0_RK) ) then
          res = .false.
        else
          res = .true.
        endif
      endif
      return
    ENDIF
      
    if((L21_crs_v2 .dot. v1_crs_v2)>0.0_RK) then
      a = d/e
    else
      a= - d/e
    endif
    if((L21_crs_v1 .dot. v1_crs_v2)>0.0_RK) then
      b = c/e
    else
      b= - c/e
    endif      
    if(a>=-1.00E-10_RK .and. a<=1.0_RK+1.00E-10_RK   .and. &
       b>=-1.00E-10_RK .and. b<=1.0_RK+1.00E-10_RK)  then
      res = .true.
    else   
      res = .false.
    endif
  end function IsLineInsertInner
  
  !**********************************************************************
  ! Initializing the geometry object 
  !**********************************************************************
  subroutine G_InitAllocate(this)
    implicit none
    class(Geometry):: this
    integer:: nUnitFile,ierror
    character(128) :: chFile
        
    ! locals
    real(RK):: LenExp
         
    LenExp=  1.2_RK*maxval( DEMProperty%Prtcl_PureProp%Radius )
#ifdef CFDACM
    LenExp= LenExp+maxval(dlub_pw)
#endif
    pmin = DEm_opt%SimDomain_min - LenExp *real3(1.0_RK,1.0_RK,1.0_RK)
    pmax = DEm_opt%SimDomain_max + LenExp *real3(1.0_RK,1.0_RK,1.0_RK)
    pmin_local =real3(DEM_decomp%xSt,DEM_decomp%ySt,DEM_decomp%zSt)-LenExp*real3(1.0_RK,1.0_RK,1.0_RK)
    pmax_local =real3(DEM_decomp%xEd,DEM_decomp%yEd,DEM_decomp%zEd)+LenExp*real3(1.0_RK,1.0_RK,1.0_RK)
        
    this%num_pWall = 0
    this%nPW_local = 0
    allocate(this%pWall( MaxWallSize ))
        
    if(nrank/=0) return
    write(chFile,"(A)") trim(DEM_opt%ResultsDir)//"WallsFor"//trim(DEM_opt%RunName)//".backup"
    open(newunit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0.and.nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_InitAllocate","Cannot open file: "//trim(chFile))
    close(nUnitFile,IOSTAT=ierror)       
  end subroutine G_InitAllocate

  !**********************************************************************
  ! Add a plane to the geometry object 
  !**********************************************************************
  subroutine G_add_PlaneWall(this, p1, p2, p3, p4 , user_id ,prop_type,  both, infinite, t_vel_t)
    implicit none
    class(Geometry) this
    type(real3),intent(in) :: p1, p2, p3, p4   ! corner points 
    integer,intent(in) :: user_id, prop_type
    logical,optional,intent(in) :: both        !  both side active status
    logical,optional,intent(in) :: infinite
    type(real3),optional,intent(in) :: t_vel_t ! translational velocity
        
    ! locals
    type(PlaneWall)::wall
    type(real3)::t_vel,ln
    character(128)::chFile
    integer::nUnitFile,ierror
    logical::lboth,linfinite
    type(PlaneWall),dimension(:),allocatable:: wall_temp
    
    IF(this%nPW_local == MaxWallSize) THEN
      MaxWallSize  = int( real(MaxWallSize, RK) *1.2_RK) +1
      call move_alloc(this%pWall, wall_temp)
      allocate(this%pWall(MaxWallSize),Stat=ierror) 
      if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Reallocation failed, 2 ")
      this%pWall(1: this%nPW_local) = wall_temp
      deallocate(wall_temp)
    ENDIF        
        
    ! checking for both side spec. 
    lboth = .false.   !.false.
    linfinite=.true.
    t_vel=zero_r3
    if(present(both))    lboth = both
    if(present(infinite))linfinite=infinite
    if(present(t_vel_t)) t_vel = t_vel_t
    
    ! assignments 
    ln = (p2-p1).cross.(p3-p1)
    wall%p1 = p1
    wall%p2 = p2
    wall%p3 = p3
    wall%p4 = p4
    wall%min_point%x= min(p4%x,min(p3%x,min(p2%x,p1%x)))
    wall%min_point%y= min(p4%y,min(p3%y,min(p2%y,p1%y)))
    wall%min_point%z= min(p4%z,min(p3%z,min(p2%z,p1%z)))
    wall%max_point%x= max(p4%x,max(p3%x,max(p2%x,p1%x)))
    wall%max_point%y= max(p4%y,max(p3%y,max(p2%y,p1%y)))
    wall%max_point%z= max(p4%z,max(p3%z,max(p2%z,p1%z)))
    wall%user_id = user_id
    wall%wall_Type = prop_type
    wall%bothSide = lboth
    wall%IsInfinite=linfinite
    wall%n = (1.0_RK/norm(ln))*ln
    wall%d = -(wall%n .dot. p1)
    wall%trans_vel = t_vel
    if(abs((wall%n .dot. p4) +wall%d) >= 0.00001_RK .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot create a plane wall, wall No. "//num2str(wall%wall_id))
    endif
        
    IF(wall%IsInDomain(pmin, pmax)) THEN
      this%num_pWall= this%num_pWall +1
      wall%wall_id  = this%num_pWall + DEM_opt%Base_wall_id
      if(wall%IsInDomain(pmin_local,pmax_local)) then
        this%nPW_local = this%nPW_local +1
        this%pWall(this%nPW_local) = wall 
      endif
        
      if(nrank/=0) return
      write(chFile,"(A)") trim(DEM_opt%ResultsDir)//"WallsFor"//trim(DEM_opt%RunName)//".backup"
      open(newunit=nUnitFile, file=chFile, status='old',position='append',form='formatted',IOSTAT=ierror )
      if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot open file: "//trim(chFile))
      write(nUnitFile,* ) p1
      write(nUnitFile,* ) p2
      write(nUnitFile,* ) p3
      write(nUnitFile,* ) p4
      close(nUnitFile,IOSTAT=ierror)
    ELSE
      if(nrank/=0) return
      call DEMLogInfo%CheckForError(ErrT_Pass,"G_add_PlaneWall","  The following plane ISNOT within the simulation domain: ")
      call DEMLogInfo%OutInfo("   It will be skipped :",3, .true.)
      call DEMLogInfo%OutInfo("   Point 1: "//trim(num2str(p1%x))//'  '//trim(num2str(p1%y))//'  '//trim(num2str(p1%z)), 3, .true.)
      call DEMLogInfo%OutInfo("   Point 2: "//trim(num2str(p2%x))//'  '//trim(num2str(p2%y))//'  '//trim(num2str(p2%z)), 3, .true.)
      call DEMLogInfo%OutInfo("   Point 3: "//trim(num2str(p3%x))//'  '//trim(num2str(p3%y))//'  '//trim(num2str(p3%z)), 3, .true.)
      call DEMLogInfo%OutInfo("   Point 4: "//trim(num2str(p4%x))//'  '//trim(num2str(p4%y))//'  '//trim(num2str(p4%z)), 3, .true.)
    ENDIF
  end subroutine G_add_PlaneWall

  !**********************************************************************
  ! MakeGeometry
  !**********************************************************************     
  subroutine G_MakeGeometry(this,chFile)
    implicit none
    class(Geometry)::this
    character(*),intent(in)::chFile
        
    !locals
    integer:: i,nplane,nUnitFile,ierror
    integer,allocatable,dimension(:):: user_id, wall_Type
    logical,allocatable,dimension(:):: BothSide, IsInfinite
    type(real3):: p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12
    type(real3):: p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24
    type(real3),allocatable,dimension(:)::Point1,Point2,Point3,Point4,TraVel
    NAMELIST/GeometryMakingNumPlane/nplane
    NAMELIST/GeometryMakingParam/Point1,Point2,Point3,Point4,TraVel,user_id, wall_Type,bothSide, IsInfinite
        
    call this%InitAllocate()
    if(DEM_Opt%GeometrySource ==0) then  ! add the geometry directly

      !sandbox
      p01= real3(-0.075_RK,  0.30_RK, -0.075_RK)
      p02= real3( 0.075_RK,  0.30_RK, -0.075_RK)
      p03= real3( 0.075_RK,  0.30_RK,  0.075_RK)
      p04= real3(-0.075_RK,  0.30_RK,  0.075_RK)

      p05= real3(-0.075_RK,  0.14_RK, -0.075_RK)
      p06= real3( 0.075_RK,  0.14_RK, -0.075_RK)
      p07= real3( 0.075_RK,  0.14_RK,  0.075_RK)
      p08= real3(-0.075_RK,  0.14_RK,  0.075_RK)
  
      p09= real3(-0.015_RK,  0.05_RK, -0.015_RK)
      p10= real3( 0.015_RK,  0.05_RK, -0.015_RK)
      p11= real3( 0.015_RK,  0.05_RK,  0.015_RK)
      p12= real3(-0.015_RK,  0.05_RK,  0.015_RK)

      p13= real3(-0.015_RK, -0.05_RK, -0.015_RK)
      p14= real3( 0.015_RK, -0.05_RK, -0.015_RK)
      p15= real3( 0.015_RK, -0.05_RK,  0.015_RK)
      p16= real3(-0.015_RK, -0.05_RK,  0.015_RK)

      p17= real3(-0.075_RK, -0.14_RK, -0.075_RK)
      p18= real3( 0.075_RK, -0.14_RK, -0.075_RK)
      p19= real3( 0.075_RK, -0.14_RK,  0.075_RK)
      p20= real3(-0.075_RK, -0.14_RK,  0.075_RK)

      p21= real3(-0.075_RK, -0.30_RK, -0.075_RK)
      p22= real3( 0.075_RK, -0.30_RK, -0.075_RK)
      p23= real3( 0.075_RK, -0.30_RK,  0.075_RK)
      p24= real3(-0.075_RK, -0.30_RK,  0.075_RK)            
      call this%add_PlaneWall( p01, p05, p06, p02, 1, 1, infinite=.false. ) !
      call this%add_PlaneWall( p05, p09, p10, p06, 1, 1, infinite=.false. ) !   
      call this%add_PlaneWall( p09, p13, p14, p10, 1, 1, infinite=.false. ) 
      call this%add_PlaneWall( p13, p17, p18, p14, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p17, p21, p22, p18, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p01, p04, p08, p05, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p05, p08, p12, p09, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p09, p12, p16, p13, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p13, p16, p20, p17, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p17, p20, p24, p21, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p03, p07, p08, p04, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p07, p11, p12, p08, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p11, p15, p16, p12, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p15, p19, p20, p16, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p19, p23, p24, p20, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p02, p06, p07, p03, 1, 1, infinite=.false. )    
      call this%add_PlaneWall( p06, p10, p11, p07, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p10, p14, p15, p11, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p14, p18, p19, p15, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p18, p22, p23, p19, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p01, p02, p03, p04, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p21, p24, p23, p22, 1, 1, infinite=.false. )
          
    elseif(DEM_Opt%GeometrySource ==1) then  ! add the geometry from the NAMELIST "&GeometryMakingParam"
                
      open(newunit=nUnitFile, file=chFile,status='old', form='formatted', IOSTAT=ierror)
      if(ierror /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_MakeGeometry","Cannot open file:"//trim(chFile))
      read(nUnitFile, nml=GeometryMakingNumPlane)
      if(nplane<1) then
        close(nUnitFile, IOSTAT=ierror);return
      endif
      allocate(user_id(nplane), wall_Type(nplane),bothSide(nplane),IsInfinite(nplane))
      allocate(Point1(nplane),Point2(nplane),Point3(nplane),Point4(nplane),TraVel(nplane))
      rewind(nUnitFile)
      read(nUnitFile, nml=GeometryMakingParam)
      close(nUnitFile, IOSTAT=ierror)
      do i=1,nplane
        call this%add_PlaneWall(Point1(i),Point2(i),Point3(i),Point4(i),user_id(i),wall_Type(i),bothSide(i),IsInfinite(i),TraVel(i))
      enddo
    elseif(DEM_Opt%GeometrySource ==2) then  ! add the geometry from the external STL file
            
    endif
  end subroutine G_MakeGeometry
 
end module prtcl_Geometry
module Prtcl_Hrchl_Munjiza
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Property
  use Prtcl_CL_and_CF
  use Prtcl_Variables
  use Prtcl_decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private

  integer::nlocal,nlocalp,nghost
  real(RK)::xst_cs,yst_cs,zst_cs
    
  type(integer3),dimension(:),allocatable:: box_index       ! integer coordinate of box
  integer,dimension(:), allocatable :: NextX
  integer,dimension(:), allocatable :: NextY
  integer,dimension(:), allocatable :: NextZ  
    
  type::lvl_Munjiza
    integer:: lvl              ! level number
    integer:: lvl_multiple     ! for lvl=1, lvl_multiple=1
    integer:: numPrtcl_lvl = 0 ! number of  particles in this level        
    real(RK):: minD_lvl        ! the minimum diameter of bonnding boxes in this level
    real(RK):: maxD_lvl        ! the maximum diameter of bounding boxes in this level

    real(RK):: cell_len_lvl      ! length of cell
    integer::  nx_lvl      ! number of divisions in x direction
    integer::  ny_lvl      ! number of divisions in y direction
    integer::  nz_lvl      ! number of divisions in z direction
        
    integer,dimension(:), allocatable :: HeadY
    integer,dimension(:), allocatable :: HeadX0
    integer,dimension(:), allocatable :: HeadX
    integer,dimension(:), allocatable :: HeadX2        
    integer:: curr_xList_ind
    integer:: curr_xList2_ind        
        
    integer,dimension(:,:), allocatable :: HeadZ0 ! head list for (iy-1), lower row
    integer,dimension(:,:), allocatable :: HeadZ  ! head list for iy, current row 
    integer,dimension(:,:), allocatable :: HeadZ2
    integer,dimension(0:2):: curr_zList0_ind        
    integer,dimension(0:2):: curr_zList_ind
    integer,dimension(0:2):: curr_zList2_ind
  end type lvl_Munjiza
  type(lvl_Munjiza),allocatable,dimension(:):: lvls_Munjiza ! level

  type::NBS_Munjiza_Hrchl
    integer:: mbox
    integer:: num_lvls = 1            ! number of levels
    integer:: num_Cnsv_cntct = 0      ! number of conservative contacts in the broad search phase
    integer:: lvl_num_cnsv_cntct = 0  ! number of conservative contact in this level        
  contains
    procedure:: Init_Munjiza_Hrchl
        
    ! performing a full contact search 
    procedure:: ContactSearch => NBSMH_ContactSearch
    ! calculating integer coordinates of all boxes
    procedure:: clcBoxIndexAndBuildYList   => NBSMH_clcBoxIndex_and_BuildYList 
            
    procedure:: BuildXList    => NBSMH_BuildXList
    procedure:: BuildXList2   => NBSMH_BuildXList2
    procedure:: BuildZList0   => NBSMH_BuildZList0
    procedure:: BuildZList    => NBSMH_BuildZList
    procedure:: BuildZList2   => NBSMH_BuildZList2
        
    procedure:: LoopNBSMask   => NBSMH_LoopNBSMask
    procedure:: Loop_CrossMask=> NBSMH_Loop_CrossMask
    procedure:: FineSearch    => NBSMH_FineSearch
        
    procedure,private:: OneLevelBroadSearch
    procedure:: Grow_Box_And_Next => NBSMH_Grow_Box_And_Next        
  end type NBS_Munjiza_Hrchl
  type(NBS_Munjiza_Hrchl),public,allocatable :: m_NBS_Munjiza_Hrchl
    
contains
 
  !*********************************************************************
  ! NBS_Munjiza_Hrchl
  !*********************************************************************
  subroutine Init_Munjiza_Hrchl(this )
    implicit none
    class(NBS_Munjiza_Hrchl):: this
        
    integer:: numLevels, idh, nx,ny,nz, dmax_min,numPrtcl_lvl,i, id_level,nLevel_final
    integer:: iErr1, iErr2, iErr3, iErr4, iErr5, iErr6, iErr7,iErrSum
    real(RK):: minD, maxD, minD_lvl, maxD_lvl, cell_len,Diam
    real(RK):: xed_cs,yed_cs,zed_cs
    type(integer3)::numCell
        
    maxD   = 2.0_RK*maxval( DEMProperty%Prtcl_PureProp%Radius )
    minD   = 2.0_RK*minval( DEMProperty%Prtcl_PureProp%Radius )
    dmax_min = int(maxD/minD)
    if( DEM_opt%CS_numlvls <=0) then
      select case( dmax_min )
      case(0:1)
        numLevels = 1
      case(2:3)
        numLevels = 2
      case(4:7)
        numLevels = 3
      case(8:15)
        numLevels = 4
      case(16:31)
        numLevels = 5
      case default
        numLevels = 6
      end select
    else
      numLevels = DEM_opt%CS_numlvls            
    endif

    maxD_lvl = maxD
    nLevel_final = 0
    DO idh = 1, numLevels
      numPrtcl_lvl = 0
      minD_lvl = maxD_lvl/2.0_RK; if(idh==numLevels) minD_lvl=0.0_RK

      do i = 1, DEM_opt%numPrtcl_Type
        Diam = 2.0_RK * DEMProperty%Prtcl_PureProp(i)%Radius
        if(Diam >minD_lvl .and. Diam <= maxD_lvl) then
          numPrtcl_lvl = numPrtcl_lvl + DEMProperty%nPrtcl_in_Bin(i)
        endif
      enddo
      if(numPrtcl_lvl>0) nLevel_final=nLevel_final+1
      maxD_lvl = minD_lvl
    ENDDO
    this%num_lvls = nLevel_final
    this%mbox = GPrtcl_list%mlocal + GPrtcl_list%mGhost_CS
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Contact search method is NBS Munjiza Hierarchy", 3 )
      call DEMLogInfo%OutInfo(" Number of levels is :" // trim( num2str(nLevel_final)),3)
    endif

    allocate( lvls_Munjiza(nLevel_final),  Stat = iErr1 )
    allocate( box_index(this%mbox), Stat = iErr2 ) 
    allocate( NextX( this%mbox ),   STAT = iErr3 )
    allocate( NextY( this%mbox ),   STAT = iErr4 )
    allocate( NextZ( this%mbox ),   STAT = iErr5 ) 
    iErrSum =  abs(iErr1) + abs(iErr2) + abs(iErr3) + abs(iErr4) + abs(iErr5) 
    if(iErrSum /= 0) then
      call DEMLogInfo%CheckForError( ErrT_Abort, "Init_Munjiza_Hrchl", "Allocations failed 1" )
    endif
    NextX = -1
    NextY = -1
    NextZ = -1

    cell_len =  DEM_Opt%Prtcl_cs_ratio * maxD
    xst_cs = DEM_decomp%xSt - cell_len*1.05_RK
    yst_cs = DEM_decomp%ySt - cell_len*1.05_RK
    zst_cs = DEM_decomp%zSt - cell_len*1.05_RK
    xed_cs = DEM_decomp%xEd + cell_len*1.05_RK
    yed_cs = DEM_decomp%yEd + cell_len*1.05_RK 
    zed_cs = DEM_decomp%zEd + cell_len*1.05_RK
    nx = int((xed_cs-xst_cs)/cell_len)+1
    ny = int((yed_cs-yst_cs)/cell_len)+1
    nz = int((zed_cs-zst_cs)/cell_len)+1

    id_level = 0    
    maxD_lvl = maxD ! setting maximum diameter of the first level equal to maximum diameter of bounding boxes
    DO idh = 1, numLevels
      ! the minimum diameter of the level is half of the maximum diameter
      ! modifying the minimum diameter of the last level and sets it to a very small value
      minD_lvl = maxD_lvl/2.0_RK 
      if( idh == numLevels) minD_lvl = 0.0_RK

      numPrtcl_lvl = 0
      do i = 1, DEM_opt%numPrtcl_Type
        Diam = 2.0_RK * DEMProperty%Prtcl_PureProp(i)%Radius
        if(Diam >minD_lvl .and. ((idh==1) .or. (idh>1 .and. Diam <= maxD_lvl))) then
          numPrtcl_lvl = numPrtcl_lvl + DEMProperty%nPrtcl_in_Bin(i)
        endif
      enddo
            
      if(numPrtcl_lvl>0) then
        id_level = id_level + 1
        do i = 1, DEM_opt%numPrtcl_Type
          Diam = 2.0_RK * DEMProperty%Prtcl_PureProp(i)%Radius
          if(Diam >minD_lvl .and. ((idh==1) .or. (idh>1 .and. Diam <= maxD_lvl))) then
            DEMProperty%CS_Hrchl_level(i) = id_level
          endif
        enddo
        lvls_Munjiza(id_level)%lvl = id_level
        lvls_Munjiza(id_level)%lvl_multiple = idh
        lvls_Munjiza(id_level)%nx_lvl = nx
        lvls_Munjiza(id_level)%ny_lvl = ny
        lvls_Munjiza(id_level)%nz_lvl = nz
        lvls_Munjiza(id_level)%cell_len_lvl = cell_len
        lvls_Munjiza(id_level)%numPrtcl_lvl = numPrtcl_lvl

        allocate(lvls_Munjiza(id_level)%HeadY( 0:ny+1), STAT= iErr1)
        allocate(lvls_Munjiza(id_level)%HeadX0(0:nx+1), STAT= iErr2)
        allocate(lvls_Munjiza(id_level)%HeadX( 0:nx+1), STAT= iErr3)
        allocate(lvls_Munjiza(id_level)%HeadX2(0:nx+1), STAT= iErr4)
        allocate(lvls_Munjiza(id_level)%HeadZ0(0:2,0:nz+1), STAT= iErr5)
        allocate(lvls_Munjiza(id_level)%HeadZ( 0:2,0:nz+1), STAT= iErr6)
        allocate(lvls_Munjiza(id_level)%HeadZ2(0:2,0:nz+1), STAT= iErr7)
        iErrSum=abs(iErr1)+abs(iErr2)+abs(iErr3)+abs(iErr4)+abs(iErr5)+ abs(iErr6)+abs(iErr7)
        if(iErrSum/= 0) then
          call DEMLogInfo%CheckForError(ErrT_Abort,"Init_Munjiza_Hrchl","Allocation failed 2")
        endif
                
        lvls_Munjiza(id_level)%HeadY  = -1
        lvls_Munjiza(id_level)%HeadX0 = -1
        lvls_Munjiza(id_level)%HeadX  = -1
        lvls_Munjiza(id_level)%HeadX2 = -1
                
        lvls_Munjiza(id_level)%curr_xList_ind  = -1
        lvls_Munjiza(id_level)%curr_XList2_ind = -1
                
        lvls_Munjiza(id_level)%HeadZ0 = -1
        lvls_Munjiza(id_level)%HeadZ  = -1
        lvls_Munjiza(id_level)%HeadZ2 = -1
        lvls_Munjiza(id_level)%curr_zList0_ind = -1
        lvls_Munjiza(id_level)%curr_zList_ind  = -1
        lvls_Munjiza(id_level)%curr_ZList2_ind = -1

        !>>> log file
        numCell = integer3(nx,ny,nz)
        call DEMLogInfo%OutInfo("Level "//trim( num2str(id_level)), 4)
        call DEMLogInfo%OutInfo("   Cell size is [m]: "// trim(num2str(cell_len)), 4)
        call DEMLogInfo%OutInfo("   number of cells (x,y,z) :"//trim(num2str(numCell)),4,.true.)
        call DEMLogInfo%OutInfo("   number of particles : "//trim(num2str(numPrtcl_lvl)),4,.true.)
      endif
              
      ! halving the maximum diameter of the level to be the maximum diameter  for the next level
      nx= nx*2
      ny= ny*2
      nz= nz*2            
      maxD_lvl = minD_lvl
      cell_len = cell_len/2.0_RK
    ENDDO

  end subroutine Init_Munjiza_Hrchl
    
  !******************************************************************
  ! NBSMH_Grow_Box_And_Next
  !******************************************************************
  subroutine NBSMH_Grow_Box_And_Next(this,nbox)
    implicit none
    class(NBS_Munjiza_Hrchl):: this 
    integer,intent(in)::nbox
    
    ! lcoals
    integer:: sizen
    
    sizen= int(1.2_RK*real(this%mbox,kind=RK))
    sizen= max(sizen, nbox+1)

    deallocate(box_index); allocate(box_index(sizen))
    deallocate(NextX); allocate(NextX(sizen))
    deallocate(NextY); allocate(NextY(sizen))
    deallocate(NextZ); allocate(NextZ(sizen))

    call DEMLogInfo%CheckForError(ErrT_Pass," NBSMH_Grow_Box_And_Next"," Need to reallocate Box_And_Next")
    call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),    3)
    call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(this%mbox)),3)
    call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),    3)

    this%mbox=sizen
  end subroutine NBSMH_Grow_Box_And_Next

  !*********************************************************************
  ! a full contact search for all levels
  !*********************************************************************
  subroutine NBSMH_ContactSearch(this)
    implicit none
    class(NBS_Munjiza_Hrchl) :: this
    integer :: idh
    
    nlocal=GPrtcl_list%nlocal
    if(nlocal == 0) return
    nlocalp = nlocal + 1
    nghost=GPrtcl_list%nGhost_CS
    this%num_Cnsv_cntct = 0
    this%lvl_num_cnsv_cntct = 0

    ! grid index and Ylists of all levels
    call this%clcBoxIndexAndBuildYList()
    
    ! contact search of all levels
    do idh = this%num_lvls, 1, -1
      call this%OneLevelBroadSearch( lvls_Munjiza(idh) )            
    end do
  end subroutine NBSMH_ContactSearch

  !******************************************************************
  ! calculating integer coordinates of all boxes
  !******************************************************************    
  subroutine NBSMH_clcBoxIndex_and_BuildYList(this)
    implicit none
    class(NBS_Munjiza_Hrchl)::this
    integer::i,idh,m,n12, iy
    real(RK):: rpdx,rpdy, rpdz,cell_len

    n12=nlocal + nghost
    if(n12>this%mbox)call this%Grow_Box_And_Next(n12)

    ! nullifying list Y
    do idh=1,this%num_lvls
      lvls_Munjiza(idh)%HeadY = -1
    enddo 
        
    do i= 1,nlocal
      rpdx = GPrtcl_PosR(i)%x - xst_cs
      rpdy = GPrtcl_PosR(i)%y - yst_cs
      rpdz = GPrtcl_PosR(i)%z - zst_cs
      idh= DEMProperty%CS_Hrchl_level(GPrtcl_pType(i))
      cell_len = lvls_Munjiza(idh)%cell_len_lvl
      box_index(i)%x = floor(rpdx/cell_len)+1
      iy = floor(rpdy/cell_len)+1; box_index(i)%y = iy
      box_index(i)%z = floor(rpdz/cell_len)+1

      NextY(i) = lvls_Munjiza(idh)%HeadY(iy)
      lvls_Munjiza(idh)%HeadY(iy) = i
    enddo   
 
    m=1
    do i=nlocal+1,n12
      rpdx = GhostP_PosR(m)%x - xst_cs
      rpdy = GhostP_PosR(m)%y - yst_cs
      rpdz = GhostP_PosR(m)%z - zst_cs
      idh= DEMProperty%CS_Hrchl_level(GhostP_pType(m))
      cell_len = lvls_Munjiza(idh)%cell_len_lvl
      box_index(i)%x = floor(rpdx/cell_len)+1
      iy = floor(rpdy/cell_len)+1; box_index(i)%y = iy
      box_index(i)%z = floor(rpdz/cell_len)+1

      NextY(i) = lvls_Munjiza(idh)%HeadY(iy)
      lvls_Munjiza(idh)%HeadY(iy) = i
      m=m+1    
    enddo

  end subroutine NBSMH_clcBoxIndex_and_BuildYList
 
  !******************************************************************
  ! constructing Xlist of row iy 
  !******************************************************************
  subroutine NBSMH_BuildXList(this,base_lvl,iy)
  implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza):: base_lvl
    integer,intent(in)::iy ! row index
  integer:: n,ix
        
    ! first checking if the xlist of row iy has been constructed previously 
    if(iy==base_lvl%curr_xList_ind) return

    ! nullifying the xlist of current row but keeps the previous raw
    base_lvl%HeadX= -1
    base_lvl%curr_xList_ind = iy
        
    n = base_lvl%HeadY(iy)
    do while (n.ne.-1)
    ix = box_index(n)%x
    NextX( n ) = base_lvl%HeadX(ix)
    base_lvl%HeadX(ix) = n  
    n = NextY(n)
  enddo
  end subroutine NBSMH_BuildXList
    
  !******************************************************************
  ! constructing Xlist of row iy+1
  !******************************************************************
  subroutine NBSMH_BuildXList2(this,base_lvl,iy)
  implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in)::iy ! row index
  integer:: n,ix
        
  ! first checking if the xlist of row iy has been constructed previously 
    if(iy == base_lvl%curr_xList2_ind) return

    ! nullifying the xlist of current row but keeps the previous raw
  base_lvl%HeadX2= -1
    base_lvl%curr_xList_ind = iy
        
  n = base_lvl%HeadY(iy)
  do while(n.ne.-1)
    ix = box_index(n)%x
    NextX(n) = base_lvl%HeadX2(ix)
    base_lvl%HeadX2(ix) = n  
    n=NextY(n)
  enddo
  end subroutine NBSMH_BuildXList2    

  !*********************************************************************
  !   Constructing the ZList0 of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSMH_BuildZList0(this,base_lvl, ix,m,lcheck)
    implicit none
    class(NBS_Munjiza_Hrchl ):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in):: ix ! column index
    integer,intent(in):: m  ! the column location with respect to ix  
    logical,optional,intent(in):: lcheck
    integer:: n,iz
  
  if(present(lcheck)) then
      if(lcheck.and.ix == base_lvl%curr_zList0_ind(m)) return
    endif
        
    ! 0: previous (left) column, 1: current column, 2: next (right) column
  base_lvl%HeadZ0(m,:) = -1
  base_lvl%curr_ZList0_ind(m) = ix
  n = base_lvl%HeadX0(ix+m-1) ! reading from the row below iy (or iy-1)
  do while (n.ne.-1)
      iz = box_index(n)%z
    NextZ(n) = base_lvl%HeadZ0(m,iz)
      base_lvl%HeadZ0(m,iz) = n
      n = NextX(n)
  enddo
  end subroutine NBSMH_BuildZList0    

  !*********************************************************************
  !   Constructing the ZList of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSMH_BuildZList(this, base_lvl,ix,m,lcheck)
    implicit none
    class(NBS_Munjiza_Hrchl ):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in):: ix ! col index
    integer,intent(in):: m  ! the column location with resect to ix 
    logical,optional,intent(in):: lcheck
    integer::n,iz
      
    if(present(lcheck)) then
      if(lcheck.and.ix == base_lvl%curr_zList_ind(m)) return
    endif
        
  ! nullifying the zlist of current col ix, but keeps the previous and next cols
    ! 0: previous (left) column, 1: current column, 2: next (right) column
  base_lvl%HeadZ(m,:) = -1
  base_lvl%curr_ZList_ind(m) = ix
  n = base_lvl%HeadX(ix+m-1) ! reading from current row iy
        
    do while (n.ne.-1)
      iz = box_index(n)%z
    NextZ(n) = base_lvl%HeadZ(m,iz)
      base_lvl%HeadZ(m,iz) = n
      n = NextX(n)
  enddo
  end subroutine NBSMH_BuildZList

  !*********************************************************************
  !   Constructing the ZList2 of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSMH_BuildZList2(this,base_lvl,ix, m,lcheck )
    implicit none
    class(NBS_Munjiza_Hrchl ):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in) :: ix ! col index
    integer,intent(in) :: m  ! the column location with resect to ix 
    logical,optional,intent(in):: lcheck
    integer::n,iz
      
    if(present(lcheck)) then
      if(lcheck.and.ix == base_lvl%curr_zList2_ind(m)) return
    endif
        
  ! nullifying the zlist of current col ix, but keeps the previous and next cols
    ! 0: previous (left) column, 1: current column, 2: next (right) column
  base_lvl%HeadZ2(m,:) = -1
  base_lvl%curr_ZList2_ind(m) = ix
  n = base_lvl%HeadX2(ix+m-1) ! reading from current row iy
        
    do while (n.ne.-1)
      iz = box_index(n)%z
    NextZ(n) = base_lvl%HeadZ2(m,iz)
      base_lvl%HeadZ2(m,iz) = n
      n = NextX(n)
  enddo
  end subroutine NBSMH_BuildZList2
    
  !******************************************************************************
  !One level contact search (intra-level and cross level with next levels)
  !******************************************************************************
  subroutine OneLevelBroadSearch(this,base_lvl)
    implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza) base_lvl
    integer:: idh,lvlCoe
    integer:: ix,iy,iz,crs_indx,crs_indy,crs_indz
    
    ! same level
    base_lvl%HeadX0(:) = -1
    ! starting the loop over all rows in base level
    DO iy=1,base_lvl%ny_lvl
      call this%BuildXList(base_lvl, iy)
      ! constructing the xLists of lvl at above rows
      DO idh=base_lvl%lvl-1, 1, -1
        lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
        crs_indy=(iy -1)/lvlCoe + 1
        if(crs_indy == 1) then
          ! for the first row, the xLists of below and top rows should be constructed
          lvls_Munjiza(idh)%HeadX0(:) = -1
          call this%BuildXlist(lvls_Munjiza(idh), 1)
        endif
        call this%BuildXList2(lvls_Munjiza(idh),  crs_indy+1 )
      ENDDO
                        
      ! if row is non-empty
      IF(base_lvl%HeadY(iy) .ne. -1 ) THEN
        base_lvl%HeadZ(0,:) = -1
        base_lvl%HeadZ0(0,:)= -1
               
        ! Creating the zlist of the lower row and current ix (column)
        call this%BuildZList0(base_lvl, 1, 1)
        do ix = 1,base_lvl%nx_lvl
          call this%BuildZList(base_lvl, ix, 1)
          call this%BuildZList0(base_lvl,ix, 2)
                
          ! constructing the zLists of the lvl at the right column
          do idh =base_lvl%lvl-1, 1, -1
                
            lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
            crs_indx = (ix-1)/lvlCoe + 1
            if(crs_indx == 1 ) then
              ! for the first column, the zLists of current and left columns should be constructed
              lvls_Munjiza(idh)%HeadZ0(0,:) = -1
              lvls_Munjiza(idh)%HeadZ(0,:)  = -1
              lvls_Munjiza(idh)%HeadZ2(0,:) = -1
              call this%BuildZList0(lvls_Munjiza(idh),crs_indx,1, .true. )    
              call this%BuildZList(lvls_Munjiza(idh), crs_indx,1, .true. )    
              call this%BuildZList2(lvls_Munjiza(idh),crs_indx,1, .true. )
            endif
                        
            call this%BuildZList0(lvls_Munjiza(idh),crs_indx,2, .true. )    
            call this%BuildZList(lvls_Munjiza(idh), crs_indx,2, .true. )    
            call this%BuildZList2(lvls_Munjiza(idh),crs_indx,2, .true. )  
          enddo                        
                              
          if(base_lvl%Headx(ix).ne.-1) then
                    
            do iz= 1, base_lvl%nz_lvl
              ! same level NBS mask check
              call this%LoopNBSMask(base_lvl, iz)
              do idh = base_lvl%lvl-1, 1, -1
                lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
                crs_indz = (iz-1)/lvlCoe + 1
                call this%Loop_CrossMask(base_lvl,iz,crs_indz,lvls_Munjiza(idh))
              enddo
            enddo
          endif
                
          ! same row, subs
          base_lvl%HeadZ(0,:) = base_lvl%HeadZ(1,:)
            
          ! lower row, subs
          base_lvl%HeadZ0(0,:) = base_lvl%HeadZ0(1,:)
          base_lvl%HeadZ0(1,:) = base_lvl%HeadZ0(2,:)
                
          do idh = base_lvl%lvl-1, 1, -1
            lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
            if(mod(ix, lvlCoe) == 0)then
              ! swap zlists
                        
              lvls_Munjiza(idh)%HeadZ0(0,:) = lvls_Munjiza(idh)%HeadZ0(1,:)
              lvls_Munjiza(idh)%HeadZ0(1,:) = lvls_Munjiza(idh)%HeadZ0(2,:)
              lvls_Munjiza(idh)%curr_ZList0_ind(0) = lvls_Munjiza(idh)%curr_ZList0_ind(1)
              lvls_Munjiza(idh)%curr_ZList0_ind(1) = lvls_Munjiza(idh)%curr_ZList0_ind(2)
    
              lvls_Munjiza(idh)%HeadZ(0,:) = lvls_Munjiza(idh)%HeadZ(1,:)
              lvls_Munjiza(idh)%HeadZ(1,:) = lvls_Munjiza(idh)%HeadZ(2,:)
              lvls_Munjiza(idh)%curr_ZList_ind(0) = lvls_Munjiza(idh)%curr_ZList_ind(1)
              lvls_Munjiza(idh)%curr_ZList_ind(1) = lvls_Munjiza(idh)%curr_ZList_ind(2)
    
              lvls_Munjiza(idh)%HeadZ2(0,:) = lvls_Munjiza(idh)%HeadZ2(1,:)
              lvls_Munjiza(idh)%HeadZ2(1,:) = lvls_Munjiza(idh)%HeadZ2(2,:)
              lvls_Munjiza(idh)%curr_ZList2_ind(0) = lvls_Munjiza(idh)%curr_ZList2_ind(1)
              lvls_Munjiza(idh)%curr_ZList2_ind(1) = lvls_Munjiza(idh)%curr_ZList2_ind(2)
                        
            endif 
          enddo
        enddo
      ENDIF
        
      base_lvl%Headx0(:) = base_lvl%Headx(:)
      DO idh= base_lvl%lvl-1, 1, -1
        lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
        if(mod(iy, lvlCoe) == 0)then
          ! swap xlists
                
          lvls_Munjiza(idh)%HeadX0 = lvls_Munjiza(idh)%HeadX
          lvls_Munjiza(idh)%HeadX = lvls_Munjiza(idh)%HeadX2
          lvls_Munjiza(idh)%curr_xList_ind = lvls_Munjiza(idh)%curr_xList2_ind
        endif
      ENDDO
    ENDDO
    
  end subroutine OneLevelBroadSearch

  !**********************************************************************
  ! finding contacts between particles in the target cell and particles in
  ! cells determined by NBS mask.  
  !**********************************************************************
  subroutine NBSMH_LoopNBSMask( this, base_lvl, iz)
    implicit none
    class(NBS_Munjiza_Hrchl) :: this
    type(lvl_Munjiza):: base_lvl
    integer,intent(in)  :: iz
    integer m, n, i, lx
  
    m = base_lvl%HeadZ(1,iz)
    DO WHILE( m .ne. -1 )

      IF(m<nlocalp) THEN

        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          call this%FineSearch(n,m)
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = base_lvl%HeadZ(1,iz-1)
        DO WHILE (n .ne. -1)
          call this%FineSearch(n,m)
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = base_lvl%HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            call this%FineSearch(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = base_lvl%HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              call this%FineSearch(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              n = NextZ(n)
            ENDDO
          enddo
        enddo
        m= NextZ(m)
      ELSE
        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nlocalp) then
            call this%FineSearch(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = base_lvl%HeadZ(1,iz-1)
        DO WHILE (n .ne. -1)
          if(n<nlocalp) then
            call this%FineSearch(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = base_lvl%HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nlocalp) then
              call this%FineSearch(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            endif
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = base_lvl%HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nlocalp) then
                call this%FineSearch(n,m)
                this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              endif
              n = NextZ(n)
            ENDDO
          enddo
        enddo
        m= NextZ(m)
      ENDIF
    ENDDO
  end subroutine  NBSMH_LoopNBSMask

  !********************************************************************** 
  ! particle  fine search
  !**********************************************************************    
  subroutine NBSMH_FineSearch(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza_Hrchl):: this
    integer,intent(in)  :: pid1, pid2
    integer::gid
    real(RK):: dx,dy,dz,dr,d2sum,dr2,ovrlp

    IF(pid1>nlocal) THEN
      gid=pid1-nlocal
      dr= GhostP_PosR(gid)%w + GPrtcl_PosR(pid2)%w
      dr2= dr*dr
      dx= GhostP_PosR(gid)%x - GPrtcl_PosR(pid2)%x
      d2sum=dx*dx;             if(d2sum>dr2) return
      dy= GhostP_PosR(gid)%y - GPrtcl_PosR(pid2)%y
      d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
      dz= GhostP_PosR(gid)%z - GPrtcl_PosR(pid2)%z
      d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
      ovrlp = dr-sqrt(d2sum)  
      call GPPW_CntctList%AddContactPPG(pid2,gid,ovrlp)
    ELSEIF(pid2>nlocal) THEN
      gid=pid2-nlocal
      dr= GPrtcl_PosR(pid1)%w + GhostP_PosR(gid)%w
      dr2= dr*dr
      dx= GPrtcl_PosR(pid1)%x - GhostP_PosR(gid)%x
      d2sum=dx*dx;             if(d2sum>dr2) return
      dy= GPrtcl_PosR(pid1)%y - GhostP_PosR(gid)%y
      d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
      dz= GPrtcl_PosR(pid1)%z - GhostP_PosR(gid)%z
      d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
      ovrlp = dr-sqrt(d2sum)  
      call GPPW_CntctList%AddContactPPG(pid1,gid,ovrlp)
    ELSE
      dr= GPrtcl_PosR(pid1)%w + GPrtcl_PosR(pid2)%w
      dr2= dr*dr
      dx= GPrtcl_PosR(pid1)%x - GPrtcl_PosR(pid2)%x
      d2sum=dx*dx;             if(d2sum>dr2) return
      dy= GPrtcl_PosR(pid1)%y - GPrtcl_PosR(pid2)%y
      d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
      dz= GPrtcl_PosR(pid1)%z - GPrtcl_PosR(pid2)%z
      d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
      ovrlp = dr-sqrt(d2sum)        

      ! this is a convention, the lower id should be the first item in the contact pair (particle & particle)
      if(GPrtcl_id(pid1) < GPrtcl_id(pid2) ) then
         call GPPW_CntctList%AddContactPP(pid1,pid2,ovrlp)
      else
         call GPPW_CntctList%AddContactPP(pid2,pid1,ovrlp)
      endif
    ENDIF
        
  end subroutine NBSMH_FineSearch 
    
  !******************************************************************************
  !   Performing a cross-level contact search between the current level (this)
  ! and CrossMunjiza level. A cross level contact search performed between 
  ! particles from target cell in current level and 27 neighbor cells from
  ! CrossMunjiza level
  !******************************************************************************
  subroutine NBSMH_Loop_CrossMask(this, base_lvl, iz, crs_iz, CrossMunjiza)
    implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza)::base_lvl
    integer, intent(in):: iz,crs_iz
    type(lvl_Munjiza),intent(in):: CrossMunjiza
    integer::m,n,lx,l
    
    m = base_lvl%HeadZ(1,iz)
    DO WHILE(m.ne.-1)
      IF(m<nlocalp) THEN
        ! first, looping all 9 cells of cross level located in the row below of the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ0(lx,crs_iz+l)
            do while( n .ne. -1)
              ! performing the narrow phase search
              call this%FineSearch(n,m)
              this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              n = NextZ(n)
            enddo
          enddo
        enddo
         
        ! second, looping all 9 cells of cross level located in the same row as the target cell
        do lx=0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              call this%FineSearch(n,m)
              this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              n = NextZ(n)
            enddo
          enddo
        enddo
          
        ! third, looping all 9 cells of cross level located in the row above the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ2(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              call this%FineSearch(n,m)
              this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              n = NextZ(n)
            enddo
          enddo
        enddo
        m = NextZ(m)

      ELSE
        ! first, looping all 9 cells of cross level located in the row below of the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ0(lx,crs_iz+l)
            do while( n .ne. -1)
              ! performing the narrow phase search
              if(n<nlocalp .or. m<nlocalp) then
                call this%FineSearch(n,m)
                this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              endif
              n = NextZ(n)
            enddo
          enddo
        enddo
         
        ! second, looping all 9 cells of cross level located in the same row as the target cell
        do lx=0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              if(n<nlocalp .or. m<nlocalp) then
                call this%FineSearch(n,m)
                this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              endif
              n = NextZ(n)
            enddo
          enddo
        enddo
          
        ! third, looping all 9 cells of cross level located in the row above the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ2(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              if(n<nlocalp .or. m<nlocalp) then
                call this%FineSearch(n,m)
                this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              endif
              n = NextZ(n)
            enddo
          enddo
        enddo
        m = NextZ(m)

      ENDIF
    ENDDO
  end subroutine NBSMH_Loop_CrossMask

end module Prtcl_Hrchl_Munjiza
module Prtcl_Integration
  use m_TypeDef
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Parameters
#ifdef CFDDEM
  use m_Parameters,only:gravity,PrGradData,IsUxConst
#endif
  implicit none
  private 
  real(RK),parameter,dimension(2):: AB2C = [1.5_RK,-0.5_RK]
  real(RK),parameter,dimension(3):: AB3C = [23.0_RK,-16.0_RK,5.0_RK]/12.0_RK
    
  public::Prtcl_Integrate
contains
  
  !**********************************************************************
  ! calculating acceleration of particles (linear and angular) 
  !**********************************************************************
#ifdef CFDDEM
  subroutine clc_Acceleration()
    implicit none

    ! locals
    integer:: i,itype
    type(real3):: Fpforce,GravityToTal
    real(RK)::Mass,Inertia,MassInFluid,MassTot,PrtclDensity
    
    GravityToTal=DEM_Opt%Gravity
    if(IsAddFluidPressureGradient .and. IsUxConst) GravityToTal%x= PrGradData(2)
    
    do i=1,GPrtcl_list%nlocal
      itype   = GPrtcl_pType(i)
      Mass    = DEMProperty%Prtcl_PureProp(itype)%Mass
      
      MassTot    = DEMProperty%Prtcl_PureProp(itype)%MassOfFluid*0.5_RK +Mass
      MassInFluid= DEMProperty%Prtcl_PureProp(itype)%MassInFluid
      Inertia    = DEMProperty%Prtcl_PureProp(itype)%Inertia
      Fpforce = 1.50_RK*GPrtcl_FpForce(i)-0.50_RK*GPrtcl_FpForce_old(i)
      
      GPrtcl_linAcc(1,i) = (1.0_RK/MassTot)*(GPrtcl_cntctForce(i)+Fpforce)+ MassInFluid/MassTot*GravityToTal
      GPrtcl_rotAcc(1,i) = (1.0_RK/Inertia)*GPrtcl_torque(i)
    enddo
  end subroutine clc_Acceleration

#else
  subroutine clc_Acceleration()
    implicit none

    ! locals
    integer:: i,itype
    real(RK)::Mass, Inertia
    
    do i=1,GPrtcl_list%nlocal
      itype   = GPrtcl_pType(i)
      Mass    = DEMProperty%Prtcl_PureProp(itype)%Mass
      Inertia = DEMProperty%Prtcl_PureProp(itype)%Inertia
      GPrtcl_linAcc(1,i) = (1.0_RK/Mass)*(GPrtcl_cntctForce(i))+DEM_opt%gravity
      GPrtcl_rotAcc(1,i) = (1.0_RK/Inertia)*(GPrtcl_torque(i))
    enddo
  end subroutine clc_Acceleration
#endif

  !**********************************************************************
  ! calculating acceleration of particles (linear and angular) 
  !**********************************************************************
  subroutine Prtcl_Integrate(iCountDEM)
    implicit none
    integer,intent(in)::iCountDEM
    
    ! locals
    integer::i,nlocal
    real(RK)::dt,TimeIntCoe(3)
    type(real3)::linVel1,linVel2,rotVel1,rotVel2
    
    call clc_Acceleration()
    
    dt=DEM_opt%dt
    nlocal = GPrtcl_list%nlocal


    
    ! linear position
    if(DEM_Opt%PI_Method==PIM_FE) then
      DO i = 1,nlocal 
        GPrtcl_PosR(i) = GPrtcl_PosR(i) + GPrtcl_linVel(1,i) *dt
        GPrtcl_linVel(1,i) = GPrtcl_linVel(1,i) + GPrtcl_linAcc(1,i) *dt
      ENDDO

    elseif(DEM_Opt%PI_Method==PIM_AB2) then
      if(iCountDEM==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
      else
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
      endif
      DO i = 1,nlocal 
        linVel1=GPrtcl_linVel(1,i)
        GPrtcl_PosR(i)=GPrtcl_PosR(i)+(TimeIntCoe(1)*linVel1 + TimeIntCoe(2)*GPrtcl_linVel(2,i))*dt
        GPrtcl_linVel(1,i)=linVel1+(TimeIntCoe(1)*GPrtcl_linAcc(1,i)+TimeIntCoe(2)*GPrtcl_linAcc(2,i))*dt
        GPrtcl_linVel(2,i)=linVel1
        GPrtcl_linAcc(2,i)=GPrtcl_linAcc(1,i)
      ENDDO

    elseif(DEM_Opt%PI_Method==PIM_AB3 ) then
      if(iCountDEM==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
        TimeIntCoe(3)=0.0_RK
      elseif(iCountDEM==2) then
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
        TimeIntCoe(3)=0.0_RK
      else
        TimeIntCoe(1)=AB3C(1)
        TimeIntCoe(2)=AB3C(2)
        TimeIntCoe(3)=AB3C(3)  
      endif
      DO i=1,nlocal
        linVel1=GPrtcl_linVel(1,i)
        linVel2=GPrtcl_linVel(2,i)
                
        GPrtcl_PosR(i)=GPrtcl_PosR(i)+(TimeIntCoe(1)*linVel1+TimeIntCoe(2)*linVel2+TimeIntCoe(3)*GPrtcl_linVel(3,i))*dt
        GPrtcl_linVel(1,i) =linVel1+(TimeIntCoe(1)*GPrtcl_linAcc(1,i)+TimeIntCoe(2)*GPrtcl_linAcc(2,i)+ TimeIntCoe(3)*GPrtcl_linAcc(3,i))*dt

        GPrtcl_linVel(3,i) = linVel2
        GPrtcl_linVel(2,i) = linVel1
        GPrtcl_linAcc(3,i) = GPrtcl_linAcc(2,i)
        GPrtcl_linAcc(2,i) = GPrtcl_linAcc(1,i)
      ENDDO
    endif
    
    ! rotate position
    if(DEM_Opt%PRI_Method==PIM_FE) then
      DO i=1,nlocal
        GPrtcl_theta(i)= GPrtcl_theta(i)+ GPrtcl_rotVel(1,i)*dt
        GPrtcl_rotVel(1,i) = GPrtcl_rotVel(1,i) + GPrtcl_rotAcc(1,i) *dt
      ENDDO        

    elseif(DEM_Opt%PRI_Method==PIM_AB2) then
      DO i=1,nlocal 
        rotVel1=GPrtcl_linVel(1,i)
                
        GPrtcl_theta(i)=GPrtcl_theta(i)+(TimeIntCoe(1)*rotVel1+TimeIntCoe(2)*GPrtcl_rotVel(2,i))*dt
        GPrtcl_rotVel(1,i)=rotVel1+(TimeIntCoe(1)*GPrtcl_rotAcc(1,i)+TimeIntCoe(2)*GPrtcl_rotAcc(2,i))*dt
                
        GPrtcl_rotVel(2,i)=rotVel1
        GPrtcl_rotAcc(2,i)=GPrtcl_rotAcc(1,i)
      ENDDO        

    elseif(DEM_Opt%PRI_Method==PIM_AB3 ) then
      DO i=1,nlocal 
        rotVel1=GPrtcl_rotVel(1,i)
        rotVel2=GPrtcl_rotVel(2,i)
        GPrtcl_theta(i)=GPrtcl_theta(i)+(TimeIntCoe(1)*rotVel1+TimeIntCoe(2)*rotVel2+TimeIntCoe(3)*GPrtcl_rotVel(3,i))*dt
        GPrtcl_rotVel(1,i)=rotVel1+(TimeIntCoe(1)*GPrtcl_rotAcc(1,i)+TimeIntCoe(2)*GPrtcl_rotAcc(2,i)+ &
                                    TimeIntCoe(3)*GPrtcl_rotAcc(3,i))*dt
                
        GPrtcl_rotVel(3,i) = rotVel2
        GPrtcl_rotVel(2,i) = rotVel1
        GPrtcl_rotAcc(3,i) = GPrtcl_rotAcc(2,i)
        GPrtcl_rotAcc(2,i) = GPrtcl_rotAcc(1,i)
      ENDDO
    endif
    
  end subroutine Prtcl_Integrate

end module Prtcl_Integration
module Prtcl_IOAndVisu
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Comm
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank,nproc
#endif
  implicit none
  private

  integer,parameter:: IK = 4
  integer::Prev_BackUp_itime= 53456791
  logical::saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel
  logical::save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque
#ifdef CFDACM
  logical::save_HighSt
#endif

  type::part_io_size_vec
    integer,dimension(1)::sizes
    integer,dimension(1)::subsizes
    integer,dimension(1)::starts
  end type part_io_size_vec
  type::part_io_size_mat
    integer,dimension(2)::sizes
    integer,dimension(2)::subsizes
    integer,dimension(2)::starts
  end type part_io_size_mat

  type:: Prtcl_IO_Visu
  contains
    procedure:: Init_visu     =>  PIO_Init_visu
    procedure:: Final_visu    =>  PIO_Final_visu
    procedure:: Dump_visu     =>  PIO_Dump_visu
    procedure:: Read_Restart  =>  PIO_Read_Restart
    procedure:: ReadFixEdCoord=>  PIO_ReadFixEdCoord
    procedure:: RestartCL     =>  PIO_RestartCL
    procedure:: Write_Restart =>  PIO_Write_Restart
    procedure:: Delete_Prev_Restart =>  PIO_Delete_Prev_Restart
    procedure,private:: Write_XDMF  =>  PIO_Write_XDMF
#ifdef CFDDEM
    procedure:: ReadInitialCoord  =>  PIO_ReadInitialCoord
    procedure:: ReadFixEdRestart  =>  PIO_ReadFixEdRestart
    procedure:: WriteFixEdRestart =>  PIO_WriteFixEdRestart
#endif
#ifdef CFDACM
    procedure:: ReadInitialCoord  =>  PIO_ReadInitialCoord
#endif
  end type Prtcl_IO_Visu
  type(Prtcl_IO_Visu),public:: DEM_IO

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_vector,  Prtcl_dump_int_matrix
    module procedure Prtcl_dump_real_vector, Prtcl_dump_real3_vector, Prtcl_dump_real3_matrix
  end interface Prtcl_dump
contains

  !**********************************************************************
  ! PIO_Init_visu
  !**********************************************************************
  subroutine PIO_Init_visu(this,chFile,iStage)
    implicit none
    class(Prtcl_IO_Visu)::this
    character(*),intent(in)::chFile
    integer,intent(in)::iStage
    
    ! locals
    character(128)::XdmfFile
    integer::nUnitFile,ierror,indent,nflds,ifld
#ifdef CFDACM
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel,   &
                               save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque,save_HighSt
#else
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel,   &
                               save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque
#endif
  
    if(iStage==1) then
      open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
      if(ierror/=0)call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu", "Cannot open file: "//trim(chFile))
      read(nUnitFile, nml=PrtclVisuOption)
      if(nrank==0)write(DEMLogInfo%nUnit, nml=PrtclVisuOption)
      close(nUnitFile,IOSTAT=ierror)
      return
    endif

    ! initialize the XDMF/XDF file
    if(nrank/=0) return
    write(XdmfFile,"(A)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu","Cannot open file:  "//trim(XdmfFile))
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! Time series
    indent =  4
    nflds = (DEM_Opt%ilast - DEM_Opt%ifirst +1)/DEM_Opt%SaveVisu  + 1
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    indent = indent + 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Time TimeType="List">'
    indent = indent + 4
    write(nUnitFile,'(A,I6,A)')repeat(' ',indent)//'<DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no') repeat(' ',indent)
    do ifld = 1,nflds
#if defined(CFDDEM) || defined(CFDACM)
      write(nUnitFile,'(I9)',advance='no') ((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)/icouple
#else
      write(nUnitFile,'(I9)',advance='no')  (ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1
#endif
    enddo
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,fmt='(A)')repeat(' ',indent)//'</Time>'
    close(nUnitFile,IOSTAT=ierror)
    if( .not. saveXDMFOnce) return

    do ifld = 1,nflds
#if defined(CFDDEM) || defined(CFDACM)
      call this%Write_XDMF(((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)/icouple)
#else
      call this%Write_XDMF((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)
#endif
    enddo

    ! XDMF/XMF Tail
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)
  end subroutine PIO_Init_visu

  !**********************************************************************
  ! PIO_Delete_Prev_Restart
  !**********************************************************************
  subroutine PIO_Delete_Prev_Restart(this,itime)
    implicit none
    class(Prtcl_IO_Visu)::this
    integer:: itime

    ! locals
    integer::nUnit,ierror
    character(128)::chFile

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if(nrank/=0) return
    write(chFile,"(A,I10.10)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName),Prev_BackUp_itime 
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)

#ifdef CFDDEM
    if(Is_clc_FluidAcc .or. (is_clc_Basset .and. is_clc_Basset_fixEd)) then
      write(chFile,"(A,I10.10)") trim(DEM_opt%RestartDir)//"FixEdSpheresRestart",Prev_BackUp_itime
      open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
      close(unit=nUnit,status='delete',IOSTAT=ierror)
    endif
    Prev_BackUp_itime = itime/icouple
#elif CFDACM
    Prev_BackUp_itime = itime/icouple
#else
    Prev_BackUp_itime = itime
#endif
  end subroutine PIO_Delete_Prev_Restart

  !**********************************************************************
  ! PIO_ReadFixEdCoord
  !**********************************************************************
  subroutine PIO_ReadFixEdCoord(this)
    implicit none
    class(Prtcl_IO_visu)::this
 
    ! locals
    type(real3)::real3t
    character(128)::chFile
    integer,parameter::NumRead=100
    integer(kind=8)::byte_total,disp
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd,radius,diam
    type(real4),allocatable,dimension(:)::real4Vec
    integer,allocatable,dimension(:)::nP_in_bin,nP_in_bin_reduce,IntVec
    integer::i,pid,ierror,nLeft,nRead,numPrtclFix,nfix,nfix_sum,pType,nfixNew,pbyte,nUnit
    real(RK),dimension(5,NumRead)::FixEdPrtcl

    ! FixEdPrtcl_type
    pbyte=real_byte*5

    numPrtclFix = DEM_opt%numPrtclFix
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
#ifdef ChanBraunJFM2011
    ySt=-yEd
#endif
#ifdef JiaYan
    if(DEM_decomp%ProcNgh(1) == MPI_PROC_NULL)zEd=zEd+2.0E-3
    if(DEM_decomp%ProcNgh(2) == MPI_PROC_NULL)zSt=zSt-2.0E-3   
    if(DEM_decomp%ProcNgh(3) == MPI_PROC_NULL)xEd=xEd+2.0E-3
    if(DEM_decomp%ProcNgh(4) == MPI_PROC_NULL)xSt=xSt-2.0E-3
#endif

    ! The data storage sequence in file "FixedSpheresCoord.dat" is as follow:
    ! Particle1: Position(real3 type), Diameter(real type), Prtcl_Type(real type)
    ! Particle2: Position(real3 type), Diameter(real type), Prtcl_Type(real type) ..
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"FixedSpheresCoord.dat"
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord","Cannot open file: "//trim(chFile))
    inquire(unit=nUnit,Pos=disp); disp=disp-1_8
    rewind(unit=nUnit,IOSTAT=ierror)
    byte_total= int(pbyte,8)*int(numPrtclFix,8)
    if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord","file byte wrong")

    nfix=0; nLeft=numPrtclFix; pid=0; disp= 1_8
    allocate(nP_in_bin(DEM_opt%numPrtcl_Type));        nP_in_bin=0
    allocate(nP_in_bin_reduce(DEM_opt%numPrtcl_Type)); nP_in_bin_reduce=0
    DO 
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp,IOSTAT=ierror)FixEdPrtcl(:,1:nRead)
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x=FixEdPrtcl(1,i)
        real3t%y=FixEdPrtcl(2,i)
        real3t%z=FixEdPrtcl(3,i)
        diam    =FixEdPrtcl(4,i)
        pType   =int(FixEdPrtcl(5,i)+0.2)
        radius= DEMProperty%Prtcl_PureProp(pType)%Radius
        if( abs(2.0_RK*radius/diam -1.0_RK)>1.0E-6 .and. nrank==0 )then
          call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord","Diameter not coordinate")
        endif
        if(real3t%x< xSt .or. real3t%y< ySt .or. real3t%z< zSt .or. &
           real3t%x>=xEd .or. real3t%y>=yEd .or. real3t%z>=zEd) cycle
        if(nfix>=GPrtcl_list%mlocalFix)  then
          nfixNew= int(1.2_RK*real(nfix,kind=RK))
          nfixNew= max(nfixNew, nfix+1)
          nfixNew= min(nfixNew,DEM_Opt%numPrtclFix)
          GPrtcl_list%mlocalFix= nfixNew

          call move_alloc(GPFix_id,IntVec)
          allocate(GPFix_id(nfixNew))
          GPFix_id(1:nfix)=IntVec
          call move_alloc(GPFix_pType,IntVec)
          allocate(GPFix_pType(nfixNew))
          GPFix_pType(1:nfix)=IntVec
          deallocate(IntVec)
          call move_alloc(GPFix_PosR,real4Vec)
          allocate(GPFix_PosR(nfixNew))
          GPFix_PosR(1:nfix)=real4Vec
          deallocate(real4Vec)
        endif
        nfix=nfix+1
         
        GPFix_id(nfix)     = pid  + DEM_opt%numPrtcl         ! NOTE HERE
        GPFix_pType(nfix)  = pType
        GPFix_PosR(nfix)   = real3t
        GPFix_PosR(nfix)%w = radius
        nP_in_bin(pType)   = nP_in_bin(pType)+1
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    close(nUnit,IOSTAT=ierror)

    call MPI_ALLREDUCE(nP_in_bin, nP_in_bin_reduce,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    DEMProperty%nPrtcl_in_Bin= DEMProperty%nPrtcl_in_Bin+ nP_in_bin_reduce
    deallocate(nP_in_bin,nP_in_bin_reduce)

    if(nfix>0) then
      call move_alloc(GPFix_id,IntVec)
      allocate(GPFix_id(nfix))
      GPFix_id=IntVec(1:nfix)
      call move_alloc(GPFix_pType,IntVec)
      allocate(GPFix_pType(nfix))
      GPFix_pType=IntVec(1:nfix)
      deallocate(IntVec)
      call move_alloc(GPFix_PosR,real4Vec)
      allocate(GPFix_PosR(nfix))
      GPFix_PosR=real4Vec(1:nfix)
      deallocate(real4Vec)
#ifdef CFDDEM
      allocate(GPFix_VFluid(2,nfix));GPFix_VFluid=zero_r3
#endif
    else
      deallocate(GPFix_id,GPFix_pType,GPFix_PosR)
    endif

    GPrtcl_list%mlocalFix = nfix
    call MPI_REDUCE(nfix,nfix_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nfix_sum/= numPrtclFix .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord"," nfix_sum/= numPrtclFix " )
    endif
  end subroutine PIO_ReadFixEdCoord

  !**********************************************************************
  ! PIO_RestartCL
  !**********************************************************************
  subroutine PIO_RestartCL(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(24)::ch
    type(real4)::real4t
    character(128)::chFile
    integer,parameter::NumRead=2000    
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    type(real3),allocatable,dimension(:)::PosVec
    integer,allocatable,dimension(:)::ncvVec,iCLvec
    integer::itime,nUnit,ierror,nlocal,np,nreal3,pbyte,ncvMax
    integer::intvec(2),i,j,ncv,nCntct,nCntctTotal,nLeft,nRead
    integer(kind=8)::disp,disp_pos,disp_ncv,disp_CL,disp_TanStart,disp_Tan
    type(real4),dimension(:),allocatable::TanDel_Un
    integer,dimension(:),allocatable::CntctVec

    itime = DEM_Opt%ifirst - 1
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
#if defined(CFDDEM) || defined(CFDACM)
    write(ch,'(I10.10)')itime/icouple
#else
    write(ch,'(I10.10)')itime
#endif
    
    ! Begin to read Restart_Contact_List
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_RestartCL","Cannot open file: "//trim(chFile))
    np= DEM_Opt%np_InDomain; disp=1_8 + int_byte  ! firstly skip the np_InDomain
    read(nUnit,pos=disp,IOSTAT=ierror)nCntctTotal; disp=disp+int_byte;
    
    nreal3 = 2*(1+GPrtcl_list%tsize+GPrtcl_list%rsize)
#ifdef CFDDEM
    nreal3=nreal3 +3      !(GPrtcl_FpForce, GPrtcl_linVelOld, GPrtcl_VFluid(1,:))
    if(Is_clc_Basset) nreal3= nreal3+ GPrtcl_BassetSeq%nDataLen
    disp=disp+2*int_byte; !skip Is_clc_Basset and HistoryStage
#endif
#ifdef CFDACM
    nreal3=nreal3+2
#endif
    pbyte=int_byte*3 + nreal3*real3_byte  ! corresponding to the subroutine 'PIO_Write_Restart'    

    allocate(PosVec(NumRead),ncvVec(NumRead))
    
    ! Determine the maxinum of ncvVec
    nLeft=np; ncvMax=0
    disp_ncv = disp + pbyte*np
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp_ncv,IOSTAT=ierror)ncvVec(1:nRead)
      disp_ncv=disp_ncv+int(int_byte,8)*int(nRead,8)
      do i=1,nRead
        if(ncvMax<ncvVec(i)) ncvMax=ncvVec(i)
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    if(ncvMax==0) ncvMax=1
    allocate(TanDel_Un(ncvMax))
    allocate(CntctVec(ncvMax))
        
    nlocal=0; nLeft=np
    disp_pos = disp
    disp_ncv = disp + pbyte*np
    disp_CL  = disp_ncv+ int_byte*np
    disp_TanStart= disp_CL + 2*nCntctTotal*int_byte
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp_pos,IOSTAT=ierror)PosVec(1:nRead)
      read(nUnit,pos=disp_ncv,IOSTAT=ierror)ncvVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      disp_ncv=disp_ncv+int(int_byte,8)*int(nRead,8)
      do i=1,nRead
        ! ncv: number of particles/walls which have overlap with this particle
        ncv=ncvVec(i)
        if(PosVec(i)%x< xSt .or. PosVec(i)%y< ySt .or. PosVec(i)%z< zSt  .or. &
           PosVec(i)%x>=xEd .or. PosVec(i)%y>=yEd .or. PosVec(i)%z>=zEd) then
          disp_CL=disp_CL+int_byte*2*ncv
        else
          nlocal= nlocal+ 1
          do j=1,ncv
            read(nUnit,pos=disp_CL,IOSTAT=ierror)intvec(1:2); disp_CL=disp_CL+int_byte*2
            disp_Tan = disp_TanStart + real4_byte*(intvec(2)-1)
            read(nUnit,pos=disp_Tan,IOSTAT=ierror)real4t
            CntctVec(j) = intvec(1)
            TanDel_Un(j)= real4t
          enddo
          if(ncv>0) call GPPW_CntctList%Add_RestartCntctlink(nlocal,ncv,CntctVec,TanDel_Un)
        endif
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    deallocate(PosVec, ncvVec)
    deallocate(TanDel_Un, CntctVec)
    close(nUnit,IOSTAT=ierror)
  end subroutine PIO_RestartCL

#if defined(CFDDEM) || defined(CFDACM)
  !**********************************************************************
  ! PIO_ReadInitialCoord
  !**********************************************************************
  subroutine PIO_ReadInitialCoord(this)
    implicit none
    class(Prtcl_IO_visu)::this
 
    ! locals
    type(real3)::real3t
    character(128)::chFile
    integer,parameter::NumRead=100
    integer(kind=8)::byte_total,disp
    integer,allocatable,dimension(:)::nP_in_bin
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd,radius,diam
    integer::i,pid,ierror,nLeft,nRead,numPrtcl,nlocal,nlocal_sum,pType,nUnit,pbyte
    real(RK),dimension(5,NumRead)::InitPrtclIn

    ! InitPrtcl_type
    pbyte=real_byte*5

    numPrtcl = DEM_opt%numPrtcl
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd

    ! The data storage sequence in file "SpheresCoord.dat" is as follow:
    ! Particle1: Position(real3 type), Diameter(real type), Prtcl_Type(real type)
    ! Particle2: Position(real3 type), Diameter(real type), Prtcl_Type(real type) ...
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"SpheresCoord.dat"
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","Cannot open file: "//trim(chFile))
    inquire(unit=nUnit,Pos=disp); disp=disp-1_8; 
    rewind(unit=nUnit,IOSTAT=ierror)
    byte_total=int(pbyte,8)*int(numPrtcl,8)
    if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","file byte wrong")

    allocate(nP_in_bin(DEM_opt%numPrtcl_Type))
    nlocal=0; nLeft=numPrtcl; pid=0; disp= 1_8; nP_in_bin=0
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp,IOSTAT=ierror)InitPrtclIn(:,1:nRead)
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x =InitPrtclIn(1,i)
        real3t%y =InitPrtclIn(2,i)
        real3t%z =InitPrtclIn(3,i)
        diam     =InitPrtclIn(4,i)
        pType=int(InitPrtclIn(5,i)+0.2)
        radius= DEMProperty%Prtcl_PureProp(pType)%Radius
        if(abs(2.0_RK*radius/diam -1.0_RK)>1.0E-6 .and. nrank==0 )then
          call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","Diameter not coordinate")
        endif
        if(real3t%x< xSt .or. real3t%y< ySt .or. real3t%z< zSt .or. &
           real3t%x>=xEd .or. real3t%y>=yEd .or. real3t%z>=zEd) cycle
        if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
        nlocal=nlocal+1
       
        GPrtcl_id(nlocal)     = pid
        GPrtcl_pType(nlocal)  = pType
        GPrtcl_PosR(nlocal)   = real3t
        GPrtcl_PosR(nlocal)%w = radius
        nP_in_bin(pType)= nP_in_bin(pType)+1      
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    close(nUnit,IOSTAT=ierror)
    GPrtcl_list%nlocal = nlocal

    call MPI_ALLREDUCE(nP_in_bin, DEMProperty%nPrtcl_in_Bin,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    deallocate(nP_in_bin)

    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","nlocal_sum/= numPrtcl " )
    endif
  end subroutine PIO_ReadInitialCoord
#endif
#ifdef CFDDEM
  !**********************************************************************
  ! PIO_WriteFixEdRestart
  !**********************************************************************
  subroutine PIO_WriteFixEdRestart(this,itime)
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    ! locals
    character(128)::chFile
    logical::IsWriteBasset
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    integer,parameter::NumWrite=100
    real(RK),dimension(5,NumWrite)::FixEdPrtclOut
    integer(kind=MPI_OFFSET_KIND)::disp,disp_pos,FileSize
    integer::bgn_ind,color,key,ierror,PrtclFix_WORLD,fh,nLeft,nWrite,i,pbyte,pid,intVec(2)

    if(DEM_opt%numPrtclFix<1) return
    IsWriteBasset=.false.
    if(is_clc_Basset .and. is_clc_Basset_fixEd) IsWriteBasset=.true.
    if((.not.Is_clc_FluidAcc) .and. (.not.IsWriteBasset)) return

    ! Create and initialize file
    if(IsWriteBasset) then
      intVec=[1,GPrtcl_BassetSeq%HistStageFix]
    else
      intVec=[0,0]
    endif
    write(chFile,"(A,I10.10)") trim(DEM_opt%RestartDir)//"FixEdSpheresRestart",itime/icouple
    if(nrank==0) then
      open(newunit=fh,file=trim(chFile),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
      write(fh,pos=1_8,IOSTAT=ierror) intVec(1:2)
      close(fh,IOSTAT=ierror)
    endif
    disp=int_byte*2_8
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    
    ! Create the Prtcl_GROUP
    bgn_ind= clc_bgn_ind(GPrtcl_list%mlocalFix)
    color = 1; key=nrank
    if(GPrtcl_list%mlocalFix<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,PrtclFix_WORLD,ierror)
    if(color==2) return

    ! PosD, id, pType
    pbyte= real_byte*5_8
    FileSize=disp+int(pbyte,8)*int(DEM_opt%numPrtclFix,8)
    call MPI_FILE_OPEN(PrtclFix_WORLD, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_BARRIER(PrtclFix_WORLD,ierror)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(PrtclFix_WORLD,ierror)
      
    nLeft=GPrtcl_list%mlocalFix;  pid=0
    disp_pos=disp+int(pbyte,8)*int(bgn_ind,8)
    DO 
      nWrite=min(nLeft,NumWrite)
      do i=1,nWrite
        pid=pid+1
        FixEdPrtclOut(1,i)=GPFix_PosR(pid)%x
        FixEdPrtclOut(2,i)=GPFix_PosR(pid)%y
        FixEdPrtclOut(3,i)=GPFix_PosR(pid)%z
        FixEdPrtclOut(4,i)=real(GPFix_id(pid),RK)
        FixEdPrtclOut(5,i)=real(GPFix_pType(pid),RK)
      enddo
      call MPI_FILE_WRITE_AT(fh,disp_pos,FixEdPrtclOut,5*nWrite,real_type,MPI_STATUS_IGNORE, ierror)
      disp_pos=disp_pos+int(pbyte,8)*int(nWrite,8)
      nLeft=nLeft-nWrite
      if(nLeft==0)exit
    ENDDO
    disp=disp+int(pbyte,8)*int(DEM_opt%numPrtclFix,8)
    call MPI_BARRIER(PrtclFix_WORLD,ierror)

    ! Write GPFix_Vfluid and GPFix_BassetData
    if(Is_clc_FluidAcc) then
      pvsize%sizes(1)   = DEM_opt%numPrtclFix
      pvsize%subsizes(1)= GPrtcl_list%mlocalFix
      pvsize%starts(1)  = bgn_ind
      call Prtcl_dump(fh, disp, GPFix_Vfluid(1,:), pvsize)
    endif
    if(IsWriteBasset) then
      pmsize%sizes(1)   = GPrtcl_BassetSeq%nDataLen;  pmsize%sizes(2)   = DEM_opt%numPrtclFix
      pmsize%subsizes(1)= GPrtcl_BassetSeq%nDataLen;  pmsize%subsizes(2)= GPrtcl_list%mlocalFix
      pmsize%starts(1)  = 0 ;                         pmsize%starts(2)  = bgn_ind
      call Prtcl_dump(fh, disp, GPFix_BassetData, pmsize)
    endif
    call MPI_FILE_CLOSE(fh, ierror) 
    call MPI_COMM_FREE(PrtclFix_WORLD,ierror)
  end subroutine PIO_WriteFixEdRestart

  !**********************************************************************
  ! PIO_ReadFixEdRestart
  !**********************************************************************
  subroutine PIO_ReadFixEdRestart(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    type(real3)::real3t
    logical::IsReadBasset
    character(128)::chFile
    integer,parameter::NumRead=500
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    integer(kind=8)::byte_total,disp,disp_Fluid
    type(real3),allocatable,dimension(:)::Real3Vec
    integer,dimension(:),allocatable::nP_in_bin,nP_in_bin_reduce,idVec,IntVec
    integer::i,pbyte,pid,iPos,pType,nUnit,ierror,numPrtclFix,nfix,nfixNew,nRead,nLeft,nfix_sum,int_t(2)
    real(RK)::FixEdPrtclIn(5,NumRead),FixEdPrtclOne(5)

    IsReadBasset=.false.
    if(is_clc_Basset .and. is_clc_Basset_fixEd) IsReadBasset=.true.
    if((.not.Is_clc_FluidAcc) .and. (.not.IsReadBasset)) then
      call this%ReadFixEdCoord()
      return
    endif

    numPrtclFix = DEM_opt%numPrtclFix
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd

    ! PosD, id, pType
    pbyte=real_byte*5

    write(chFile,'(I10.10)')(DEM_Opt%ifirst - 1)/icouple
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"FixEdSpheresRestart"//trim(adjustl(chFile))
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart","Cannot open file:"//trim(chFile))

    byte_total= int_byte*2+pbyte*numPrtclFix ! 'int_byte*2' corresponds to HistStageFix in  PIO_WriteFixEdRestart
    if(is_clc_FluidAcc) byte_total= byte_total+ real3_byte*numPrtclFix
    if(IsReadBasset) byte_total= byte_total+ real3_byte*numPrtclFix*GPrtcl_BassetSeq%nDataLen
    inquire(unit=nUnit,Pos=disp); disp=disp-1_8
    rewind(unit=nUnit,IOSTAT=ierror)
    if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoordRestart","file byte wrong")

    disp=1_8
    read(nUnit,pos=disp,IOSTAT=ierror)int_t(1:2); disp=disp+int_byte*2
    if(IsReadBasset) then
      if(int_t(1)/= 1) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart"," Is_clc_Basset_fixEd Wrong 1" )
      GPrtcl_BassetSeq%HistStageFix= int_t(2)
    else
      if(int_t(1)/= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart"," Is_clc_Basset_fixEd Wrong 2" )
    endif   

    nfix=0; nLeft=numPrtclFix; pid=0
    allocate(idVec(GPrtcl_list%mlocalFix))           
    allocate(nP_in_bin(DEM_opt%numPrtcl_Type));        nP_in_bin=0
    allocate(nP_in_bin_reduce(DEM_opt%numPrtcl_Type)); nP_in_bin_reduce=0    
    DO 
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp,IOSTAT=ierror)FixEdPrtclIn(:,1:nRead)
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x=FixEdPrtclIn(1,i)
        real3t%y=FixEdPrtclIn(2,i)
        real3t%z=FixEdPrtclIn(3,i)
        if(real3t%x< xSt .or. real3t%y< ySt .or. real3t%z< zSt .or. &
           real3t%x>=xEd .or. real3t%y>=yEd .or. real3t%z>=zEd) cycle
        if(nfix>=GPrtcl_list%mlocalFix)  then
          nfixNew= int(1.2_RK*real(nfix,kind=RK))
          nfixNew= max(nfixNew, nfix+1)
          nfixNew= min(nfixNew,numPrtclFix)
          GPrtcl_list%mlocalFix= nfixNew
          call move_alloc(idVec,IntVec)
          allocate(idVec(nfixNew))
          idVec(1:nfix)=IntVec
          deallocate(IntVec)
        endif
        nfix=nfix+1
        idVec(nfix)=pid
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    call MPI_ALLREDUCE(nP_in_bin, nP_in_bin_reduce,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    DEMProperty%nPrtcl_in_Bin= DEMProperty%nPrtcl_in_Bin+ nP_in_bin_reduce
    deallocate(nP_in_bin,nP_in_bin_reduce)

    GPrtcl_list%mlocalFix = nfix
    call MPI_REDUCE(nfix,nfix_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nfix_sum/= numPrtclFix .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart"," nfix_sum/= numPrtclFix " )
    endif

    if(nfix>0) then
      call move_alloc(idVec,IntVec)
      allocate(idVec(nfix))
      idVec=IntVec(1:nfix)
      deallocate(IntVec)
      deallocate(GPFix_id);   allocate(GPFix_id(nfix))
      deallocate(GPFix_pType);allocate(GPFix_pType(nfix))
      deallocate(GPFix_PosR); allocate(GPFix_PosR(nfix))
      allocate(GPFix_VFluid(2,nfix));GPFix_VFluid=zero_r3
      if(IsReadBasset) then
        allocate(GPFix_BassetData(GPrtcl_BassetSeq%nDataLen, nfix), Stat= ierror)
        if(ierror/=0) call DEMLogInfo%checkForError(ErrT_Abort,"PIO_ReadFixEdRestart","allocate wrong2")
        GPFix_BassetData=zero_r3
      endif
    else
      if(allocated(IntVec))deallocate(IntVec)
      deallocate(GPFix_id,GPFix_pType,GPFix_PosR)
    endif

    allocate(Real3Vec(GPrtcl_BassetSeq%nDataLen))
    do pid=1,nfix
      iPos=idVec(pid)-1
      disp=1_8+int_byte*2+int(iPos,8)*int(pbyte,8)
      read(nUnit,pos=disp,IOSTAT=ierror)FixEdPrtclOne
      GPFix_id(pid)   = int(FixEdPrtclOne(4)+0.2)
      pType = int(FixEdPrtclOne(5)+0.2)
      GPFix_pType(pid)= pType
      GPFix_PosR(pid)%x= FixEdPrtclOne(1)
      GPFix_PosR(pid)%y= FixEdPrtclOne(2)
      GPFix_PosR(pid)%z= FixEdPrtclOne(3)
      GPFix_PosR(pid)%w= DEMProperty%Prtcl_PureProp(pType)%Radius
      if(Is_clc_FluidAcc) then
        disp_Fluid= 1_8+int_byte*2+int(numPrtclFix,8)*int(pbyte,8) +int(iPos,8)*int(real3_byte,8)
        read(nUnit,pos=disp_Fluid,IOSTAT=ierror)real3t
        GPFix_Vfluid(1,pid)=real3t
      endif
      if(IsReadBasset) then
        disp_Fluid= 1_8+int_byte*2+int(numPrtclFix,8)*int(pbyte+real3_byte,8) +int(iPos,8)*int(real3_byte*GPrtcl_BassetSeq%nDataLen,8)
        read(nUnit,pos=disp_Fluid,IOSTAT=ierror)Real3Vec(1:GPrtcl_BassetSeq%nDataLen)
        GPFix_BassetData(:,pid)= Real3Vec
      endif
    enddo
    deallocate(idVec,Real3Vec)
    close(nUnit,IOSTAT=ierror)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end subroutine PIO_ReadFixEdRestart
#endif

  !**********************************************************************
  ! PIO_Read_Restart
  !**********************************************************************
  subroutine PIO_Read_Restart(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(24)::ch
    character(128)::chFile
    integer,parameter::NumRead=2000
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    integer,allocatable,dimension(:):: nP_in_bin
    type(real3),allocatable,dimension(:)::real3Vec,PosVec
    integer::nUnit,tsize,rsize,nlocal_sum,nreal3,itype
    integer::itime,ierror,nlocal,i,k,np,nLeft,nRead,int_t(3)
    integer(kind=MPI_OFFSET_KIND)::disp,disp_pos,disp_int,disp_real3
    
    itime = DEM_Opt%ifirst - 1
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
#if defined(CFDDEM) || defined(CFDACM)
    write(ch,'(I10.10)')itime/icouple
#else
    write(ch,'(I10.10)')itime
#endif

    ! Begin to read Restart_file
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart","Cannot open file: "//trim(chFile))
    disp =1_8; read(nUnit,pos=disp,IOSTAT=ierror)np; disp=disp+int_byte
    if(np>DEM_Opt%numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart: "," np_InDomain > numPrtcl " )
    endif
    DEM_Opt%np_InDomain = np
    disp=disp+int_byte   ! skip nCntctTotal
#ifdef CFDDEM
    read(nUnit,pos=disp,IOSTAT=ierror)int_t(1:2); disp=disp+2*int_byte;
    if(Is_clc_Basset) then
      if(int_t(1)/= 1) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," Is_clc_Basset Wrong 1" )
      GPrtcl_BassetSeq%HistoryStage= int_t(2)
    else
      if(int_t(1)/= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," Is_clc_Basset Wrong 2" )
    endif
#endif

    tsize=GPrtcl_list%tsize; rsize=GPrtcl_list%rsize;
    nreal3 = 2*(1+tsize+rsize)
#ifdef CFDDEM
    nreal3=nreal3+3
    if(Is_clc_Basset) nreal3=nreal3+GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    nreal3=nreal3+2
#endif
    nreal3=nreal3-1   ! "-1" corresponds to GPrtcl_Pos
    allocate(real3Vec(nreal3),PosVec(NumRead))
    allocate(nP_in_bin(DEM_opt%numPrtcl_Type)); nP_in_bin=0

    nlocal=0; nLeft=np
    disp_pos  = disp
    disp_int  = disp_pos+ real3_byte*np
    disp_real3= disp_int+ int_byte*np*3
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp_pos,IOSTAT=ierror)PosVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      do i=1,nRead
        if(PosVec(i)%x>=xSt .and. PosVec(i)%x< xEd .and. PosVec(i)%y>=ySt .and. &
           PosVec(i)%y< yEd .and. PosVec(i)%z>=zSt .and. PosVec(i)%z<zEd) then
          if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1

          read(nUnit,pos=disp_int,IOSTAT=ierror)int_t(1:3)
          GPrtcl_id(nlocal)=int_t(1)      ! id
          itype=int_t(2)
          GPrtcl_pType(nlocal)=itype      ! pType
          nP_in_bin(itype)= nP_in_bin(itype)+1
          GPrtcl_UsrMark(nlocal)=int_t(3) ! Usr_Mark

          GPrtcl_PosR(nlocal)= PosVec(i)  ! PosR
          GPrtcl_PosR(nlocal)%w=DEMProperty%Prtcl_PureProp(itype)%Radius
          k=0;
          read(nUnit,pos=disp_real3,IOSTAT=ierror)real3Vec(1:nreal3)
          GPrtcl_LinVel(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinVec
          GPrtcl_LinAcc(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinAcc
          GPrtcl_theta(nlocal)           =real3Vec(k+1);         k=k+1     ! Theta
          GPrtcl_RotVel(1:rsize,nlocal)  =real3Vec(k+1:k+rsize); k=k+rsize ! RotVel
          GPrtcl_RotAcc(1:rsize,nlocal)  =real3Vec(k+1:k+rsize); k=k+rsize ! RotAcc
#ifdef CFDACM        
          GPrtcl_FluidIntOld(1:2,nlocal)= real3Vec(k+1:k+2);     k=k+2     ! FluidIntOld
#endif
#ifdef CFDDEM     
          GPrtcl_FpForce(nlocal)         =real3Vec(k+1);         k=k+1    
          GPrtcl_linVelOld(nlocal)       =real3Vec(k+1);         k=k+1
          GPrtcl_VFluid(1,nlocal)        =real3Vec(k+1);         k=k+1 
          if(Is_clc_Basset) then
            itype=GPrtcl_BassetSeq%nDataLen
            GPrtcl_BassetData(1:itype,nlocal)=real3Vec(k+1:k+itype);     k=k+itype
          endif
#endif
        endif
        disp_int  = disp_int  + int_byte*3
        disp_real3= disp_real3+ real3_byte*nreal3       
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    deallocate(PosVec,Real3Vec)
    call MPI_ALLREDUCE(nP_in_bin, DEMProperty%nPrtcl_in_Bin,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    deallocate(nP_in_bin)

    GPrtcl_list%nlocal = nlocal
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= np .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart: "," nlocal_sum/= np_InDomain " )
    endif
    close(nUnit,IOSTAT=ierror)
  end subroutine PIO_Read_Restart

  !**********************************************************************
  ! PIO_Write_Restart
  !**********************************************************************
  subroutine PIO_Write_Restart(this,itime)
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime
   
    ! locals 
    character(24)::ch
    character(128)::chFile
    type(real4)::TanDel_Un
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    integer,parameter::NumRestart=100
    integer,allocatable,dimension(:)::IntVec
    integer,allocatable,dimension(:,:)::IntMat
    type(real3),allocatable,dimension(:)::real3Vec
    type(real4),allocatable,dimension(:)::real4Vec
    integer(kind=MPI_OFFSET_KIND)::disp,bgn_byte,FileSize
    integer::pid,i,j,k,nlocal,bgn_ind,nreal3,rsize,nRestart
    integer::nCntct,nCntctTotal,nCnt_ind,nTanDel,nTan_ind,nLeft
    integer::ierror,fh,color,key,Prtcl_WORLD1,Prtcl_WORLD2,ncv,prev,now,iCntct,tsize
    integer,dimension(:),allocatable::CntctVec

    if(DEM_Opt%np_InDomain<1) return
    nlocal = GPrtcl_list%nlocal
    call GPPW_CntctList%Get_numCntcts(nCntct,nTanDel)    
#if defined(CFDDEM) || defined(CFDACM)
    write(ch,'(I10.10)')itime/icouple
#else
    write(ch,'(I10.10)')itime
#endif
  
    ! Calculate the bgn_ind and nlink_ind
    bgn_ind = clc_bgn_ind(nlocal)
    nCnt_ind= clc_bgn_ind(nCntct)
    nTan_ind= clc_bgn_ind(nTanDel)
    call MPI_ALLREDUCE(nCntct,nCntctTotal,1,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    call GPPW_CntctList%Prepare_Restart(nTan_ind)

    ! Create and empty file, Write DEM_Opt%np_InDomain,nCntctTotal, in the begining of the Restart file
    disp = 1_8
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    if(nrank==0) then
      open(newunit=fh,file=trim(chFile),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
      write(unit=fh,pos=disp,IOSTAT=ierror) DEM_Opt%np_InDomain,nCntctTotal; disp=disp+int_byte*2
#ifdef CFDDEM
      if(Is_clc_Basset) then
        i=1; j=GPrtcl_BassetSeq%HistoryStage
      else
        i=0; j=0
      endif
      write(unit=fh,pos=disp,IOSTAT=ierror) i,j
#endif
      close(fh,IOSTAT=ierror)
    endif
#ifdef CFDDEM
    disp = int_byte*4
#else
    disp = int_byte*2
#endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    
    ! Create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD1,ierror)
    if(color==2) return

    ! Begin to write Restart file    
    call MPI_FILE_OPEN(Prtcl_WORLD1, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    pvsize%sizes(1)   = DEM_Opt%np_InDomain
    pvsize%subsizes(1)= nlocal
    pvsize%starts(1)  = bgn_ind
    allocate(real3Vec(nlocal))
    do pid=1,nlocal
      real3Vec(pid)=GPrtcl_PosR(pid)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)

    pmsize%sizes(1)     = 3;   pmsize%sizes(2)     = DEM_Opt%np_InDomain
    pmsize%subsizes(1)  = 3;   pmsize%subsizes(2)  = nlocal
    pmsize%starts(1)    = 0;   pmsize%starts(2)    = bgn_ind
    allocate(IntMat(3,nlocal))
    do pid=1,nlocal
      IntMat(1,pid)= GPrtcl_id(pid)
      IntMat(2,pid)= GPrtcl_pType(pid)
      IntMat(3,pid)= GPrtcl_UsrMark(pid)
    enddo
    call Prtcl_dump(fh,disp,IntMat(1:3,1:nlocal),  pmsize)
    deallocate(IntMat)
    call MPI_FILE_CLOSE(fh,ierror)

    tsize=GPrtcl_list%tsize
    rsize=GPrtcl_list%rsize
    nreal3 = 2*(1+GPrtcl_list%tsize+GPrtcl_list%rsize)
#ifdef CFDDEM
    nreal3=nreal3+3
    if(Is_clc_Basset) nreal3=nreal3+GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    nreal3=nreal3+2
#endif
    nreal3=nreal3-1   ! "-1" corresponds to GPrtcl_Pos
    FileSize=disp+int(nreal3*real3_byte,8)*int(DEM_Opt%np_InDomain,8)
    call MPI_FILE_OPEN(Prtcl_WORLD1, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_BARRIER(Prtcl_WORLD1,ierror)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(Prtcl_WORLD1,ierror)      
      
    allocate(real3Vec(NumRestart*nreal3))
    nLeft=nlocal; pid=0
    bgn_byte=disp+int(nreal3*real3_byte,8)*int(bgn_ind,8)
    DO
      nRestart=min(nLeft,NumRestart)
      k=0
      do i=1,nRestart
        pid=pid+1
        real3Vec(k+1:k+tsize)=GPrtcl_LinVel(1:tsize,pid); k=k+tsize
        real3Vec(k+1:k+tsize)=GPrtcl_LinAcc(1:tsize,pid); k=k+tsize
        real3Vec(k+1)        =GPrtcl_theta(pid);          k=k+1
        real3Vec(k+1:k+rsize)=GPrtcl_RotVel(1:rsize,pid); k=k+rsize
        real3Vec(k+1:k+rsize)=GPrtcl_RotAcc(1:rsize,pid); k=k+rsize
#ifdef CFDDEM
        real3Vec(k+1)        =GPrtcl_FpForce(pid);        k=k+1
        real3Vec(k+1)        =GPrtcl_linVelOld(pid);      k=k+1
        real3Vec(k+1)        =GPrtcl_VFluid(1,pid);       k=k+1
        if(Is_clc_Basset) then
          j=GPrtcl_BassetSeq%nDataLen
          real3Vec(k+1:k+j) =GPrtcl_BassetData(1:j,pid);  k=k+j
        endif
#endif
#ifdef CFDACM
        real3Vec(k+1:k+2)=GPrtcl_FluidIntOld(1:2,pid);    k=k+2
#endif
      enddo
      call MPI_FILE_WRITE_AT(fh,bgn_byte,real3Vec,k,real3_type,MPI_STATUS_IGNORE,ierror)
      bgn_byte=bgn_byte+int(nreal3*real3_byte,8)*int(nRestart,8)
      nLeft=nLeft-nRestart
      if(nLeft==0)exit
    ENDDO
    deallocate(real3Vec)
    disp=disp+int(nreal3*real3_byte,8)*int(DEM_Opt%np_InDomain,8)
    call MPI_BARRIER(Prtcl_WORLD1,ierror)
      
    ! ncv: number of particles/walls which have overlap with this particle
    iCntct=0
    allocate(IntVec(nlocal))
    DO pid=1,nlocal
      call GPPW_CntctList%Count_Cntctlink(pid,ncv)
      IntVec(pid)=ncv
      if(iCntct<ncv) iCntct=ncv
    ENDDO
    if(iCntct==0) iCntct=1
    call Prtcl_dump(fh,disp,IntVec(1:nlocal),pvsize)
    allocate( CntctVec(2*iCntct) )
    deallocate(IntVec)
    call MPI_FILE_CLOSE(fh,ierror)

    ! Begin to write Contact List file
    color = 1; key=nrank
    if(nTanDel<=0) color=2
    call MPI_COMM_SPLIT(Prtcl_WORLD1,color,key,Prtcl_WORLD2,ierror)
    call MPI_COMM_FREE( Prtcl_WORLD1, ierror)
    if(color==2) return
    call MPI_FILE_OPEN(Prtcl_WORLD2, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    pvsize%sizes(1)   = 2*nCntctTotal
    pvsize%subsizes(1)= 2*nCntct     
    pvsize%starts(1)  = 2*nCnt_ind   
    allocate(IntVec(max(2*nCntct,1)))
    iCntct=0
    DO pid=1,nlocal
      call GPPW_CntctList%Resemble_Cntctlink(pid,ncv,CntctVec)
      if(ncv>0) then
        IntVec(iCntct+1:iCntct+2*ncv)=CntctVec(1:2*ncv)
        iCntct=iCntct+2*ncv
      endif
    ENDDO
    call Prtcl_dump(fh,disp,IntVec,pvsize)
    deallocate(IntVec)
    deallocate(CntctVec)
    call MPI_FILE_CLOSE(fh,ierror)

    call MPI_ALLREDUCE(nTanDel,nCntctTotal,1,int_type,MPI_SUM,Prtcl_WORLD2,ierror)
    FileSize=disp+int(real4_byte,8)*int(nCntctTotal,8)    
    call MPI_FILE_OPEN(Prtcl_WORLD2, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_BARRIER(Prtcl_WORLD2,ierror)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(Prtcl_WORLD2,ierror)
      
    prev=1; j=0; k=0
    disp = disp + nTan_ind * real4_byte
    allocate(real4Vec(nlocal))
    DO
      call GPPW_CntctList%GetNextTanDel_Un(TanDel_Un,prev,now);  prev=now+1
      j=j+1; k=k+1; real4Vec(j)=TanDel_Un
      if(k==nTanDel) then
        call MPI_FILE_WRITE_AT(fh,disp, real4Vec, j, real4_type, MPI_STATUS_IGNORE, ierror)
        disp = disp +j*real4_byte; exit
      endif
      if(j==nlocal) then
        call MPI_FILE_WRITE_AT(fh,disp, real4Vec, j, real4_type, MPI_STATUS_IGNORE, ierror)
        disp = disp +j*real4_byte; j=0
      endif
    ENDDO
    deallocate(real4Vec)
    call MPI_FILE_CLOSE(fh, ierror)
    call MPI_COMM_FREE(Prtcl_WORLD2, ierror)
  end subroutine PIO_Write_Restart

  !**********************************************************************
  ! Purpose:
  !   Create a xdmf/xmf file in order to view the simulation results
  !     by Paraview directly
  ! 
  ! Original Author: 
  !   Pedro Costa
  ! 
  ! Modified by:
  !   Zheng Gong
  ! 
  ! Original Source file is downloaded from ( April 2020 ):
  !   https://github.com/p-costa/gen_xdmf_particles
  !
  !**********************************************************************
  subroutine PIO_Write_XDMF(this,itime) 
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    ! locals
    character(128)::chFile
    integer(kind=MPI_OFFSET_KIND)::disp
    integer:: indent,nUnitFile,ierror,np,dims,iprec

    if(nrank/=0) return 
    np=DEM_Opt%np_InDomain
    write(chFile,"(A)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=chFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Write_XDMF","Cannot open file: "//trim(chFile))
    endif

    indent = 8; disp = 0_MPI_OFFSET_KIND
    write(chFile,"(A)") "PartVisuFor"//trim(DEM_opt%RunName)
    dims=3; iprec=RK
    write(nUnitFile,'(A,I10.10,A)')repeat(' ',indent)//'<Grid Name="T',itime,'" GridType="Uniform">'
    indent = indent + 4
    write(nUnitFile,'(A,I9,A)')repeat(' ',indent)//'<Topology TopologyType="Polyvertex" NodesPerElement="',np,'"/>'
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Geometry GeometryType="'//"XYZ"//'">'
    indent = indent + 4
    write(nUnitFile,'(A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
          ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    indent = indent + 4
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(chFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Geometry>'

    IF(save_ID) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"ID","Scalar","Int",disp)
    ENDIF
    IF(save_Diameter) THEN
      dims=1; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Diameter","Scalar","Float",disp)
    ENDIF
    IF(save_Type) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Type","Scalar","Int",disp)
    ENDIF
    IF(save_UsrMark) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"UsrMark","Scalar","Int",disp)
    ENDIF
    IF(save_LinVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"LinVel","Vector","Float",disp)
    ENDIF
    IF(save_LinAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"LinAcc","Vector","Float",disp)
    ENDIF
    IF(save_Theta) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Theta","Vector","Float",disp)
    ENDIF
    IF(save_RotVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"RotVel","Vector","Float",disp)
    ENDIF
    IF(save_RotAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"RotAcc","Vector","Float",disp)
    ENDIF
    IF(save_CntctForce) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"CntctForce","Vector","Float",disp)
    ENDIF
    IF(save_Torque) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Torque","Vector","Float",disp)
    ENDIF
#ifdef CFDACM
    IF(save_HighSt) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"IsHighSt","Scalar","Int",disp)
    ENDIF
#endif    
    write(nUnitFile,'(A)')'        </Grid>'
    close(nUnitFile,IOSTAT=ierror)
  end subroutine PIO_Write_XDMF

  !**********************************************************************
  ! Write_XDMF_One
  !**********************************************************************
  subroutine Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,chName,chAttribute,chDataType,disp)
    implicit none
    integer,intent(in)::nUnitFile,dims,iprec,np,itime
    character(*),intent(in)::chFile,chName,chAttribute,chDataType
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    
    ! locals
    integer:: indent
    indent = 12

    write(nUnitFile,'(A)')repeat(' ',indent)//'<Attribute Type="'//trim(chAttribute)//'" Center="Node" Name="'//trim(chName)//'">'
    indent = indent + 4
    write(nUnitFile,'(3A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
          ' DataType="',trim(chDataType),'" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    indent = indent + 4
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(chFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Attribute>'

  end subroutine Write_XDMF_One

  !**********************************************************************
  ! PIO_Final_visu
  !**********************************************************************
  subroutine PIO_Final_visu(this)
    implicit none 
    class(Prtcl_IO_Visu)::this

    ! locals
    integer::nUnitFile,ierror
    character(128)::XdmfFile

    if(nrank/=0 .or. saveXDMFOnce) return
    write(xdmfFile,"(A)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Final_visu","Cannot open file: "//trim(XdmfFile))
    endif
    ! XDMF/XMF Tail
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile)

  end subroutine PIO_Final_visu

  !**********************************************************************
  ! PIO_Dump_visu
  !**********************************************************************
  subroutine PIO_Dump_visu(this, itime)
    implicit none
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    ! locals
    character(128)::chFile
    type(part_io_size_vec)::pvsize
    integer(kind=MPI_OFFSET_KIND)::disp
    integer,allocatable,dimension(:)::intVec
    real(RK),allocatable,dimension(:)::realVec
    type(real3),allocatable,dimension(:)::real3Vec
    integer::ierror,fh,i,color,key,Prtcl_WORLD,nlocal,bgn_ind

    ! write xdmf file first
    if(.not.saveXDMFOnce) call this%Write_XDMF(itime)

    ! update the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind=clc_bgn_ind(nlocal)

    ! Create and empty file
    write(chFile,"(A,I10.10)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName),itime
    if(nrank==0) then
      open(newunit=fh,file=trim(chFile),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
      close(fh,IOSTAT=ierror)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    
    ! create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierror)
    if(color==2) return
    
    ! begin to dump
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND
    pvsize%sizes(1)     = DEM_Opt%np_InDomain
    pvsize%subsizes(1)  = nlocal
    pvsize%starts(1)    = bgn_ind
    if(nlocal<=0) return

    allocate(real3Vec(nlocal))
    do i=1,nlocal
      real3Vec(i)=GPrtcl_PosR(i)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)
    if(save_ID) call Prtcl_dump(fh,disp, GPrtcl_id(1:nlocal),  pvsize)
    if(save_Diameter) then
      allocate(realVec(nlocal))
      do i=1,nlocal
        realVec(i)= 2.0_RK*GPrtcl_PosR(i)%w
      enddo   
      call Prtcl_dump(fh,disp, realVec(1:nlocal),  pvsize)
      deallocate(realVec)
    endif
    if(save_Type)       call Prtcl_dump(fh,disp, GPrtcl_pType(1:nlocal),     pvsize)
    if(save_UsrMark)    call Prtcl_dump(fh,disp, GPrtcl_UsrMark(1:nlocal),   pvsize)
    if(save_LinVel)     call Prtcl_dump(fh,disp, GPrtcl_LinVel(1,1:nlocal),  pvsize)
    if(save_LinAcc)     call Prtcl_dump(fh,disp, GPrtcl_LinAcc(1,1:nlocal),  pvsize)
    if(save_Theta)      call Prtcl_dump(fh,disp, GPrtcl_Theta(1:nlocal),     pvsize)
    if(save_RotVel)     call Prtcl_dump(fh,disp, GPrtcl_RotVel(1,1:nlocal),  pvsize)
    if(save_RotAcc)     call Prtcl_dump(fh,disp, GPrtcl_RotAcc(1,1:nlocal),  pvsize)
    if(save_CntctForce) call Prtcl_dump(fh,disp, GPrtcl_CntctForce(1:nlocal),pvsize)
    if(save_Torque)     call Prtcl_dump(fh,disp, GPrtcl_Torque(1:nlocal),    pvsize)
#ifdef CFDACM
    if(save_HighSt)  then
      allocate(IntVec(nlocal))
      do i=1,nlocal
        if(GPrtcl_HighSt(i)=="N") then
          IntVec(i)=0
        else
          IntVec(i)=1
        endif
      enddo   
      call Prtcl_dump(fh,disp, IntVec(1:nlocal),  pvsize)
      deallocate(IntVec)
    endif
#endif     
    call MPI_FILE_CLOSE(fh, ierror) 
    call MPI_COMM_FREE( Prtcl_WORLD, ierror)
  end subroutine PIO_Dump_visu

  !**********************************************************************
  ! clc_bgn_ind
  !**********************************************************************
  function clc_bgn_ind(nlocal) result(bgn_ind)
    implicit none
    integer,intent(in)::nlocal
    integer::bgn_ind

    ! locals
    integer::end_ind,ierror,SRstatus(MPI_STATUS_SIZE)
  
    bgn_ind=0
    if(nproc<=1) return
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    IF(nrank==0)THEN
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,int_type, nrank+1,0,MPI_COMM_WORLD,ierror)
    ELSEIF(nrank /= nproc-1) THEN
      call MPI_RECV(bgn_ind, 1,int_type, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierror)
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,int_type, nrank+1,0,MPI_COMM_WORLD,ierror)
    ELSE
      call MPI_RECV(bgn_ind, 1,int_type, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierror)
    ENDIF
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end function clc_bgn_ind

  !**********************************************************************
  ! Prtcl_dump_int_vector
  !**********************************************************************
  subroutine Prtcl_dump_int_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    integer,dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer:: ierror,newtype
    integer,dimension(1) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant integer vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, int_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,int_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),int_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * int_byte
  end subroutine Prtcl_dump_int_vector

  !**********************************************************************
  ! Prtcl_dump_int_matrix
  !**********************************************************************
  subroutine Prtcl_dump_int_matrix(fh,disp,var,pmsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_mat),intent(in)::pmsize    
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    integer,dimension(1:pmsize%subsizes(1),1:pmsize%subsizes(2)),intent(in)::var

    ! locals
    integer :: ierror,newtype
    integer, dimension(2) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes     = pmsize%sizes
    subsizes  = pmsize%subsizes
    starts    = pmsize%starts

    ! write the particle relevant real matrix
    call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, int_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,int_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1)*subsizes(2),int_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * sizes(2) * int_byte
  end subroutine Prtcl_dump_int_matrix

  !**********************************************************************
  ! Prtcl_dump_real_vector
  !**********************************************************************
  subroutine Prtcl_dump_real_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    real(RK),dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer::ierror,newtype
    integer,dimension(1):: sizes,subsizes,starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant real vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * real_byte
  end subroutine Prtcl_dump_real_vector

  !**********************************************************************
  ! Prtcl_dump_real3_vector
  !**********************************************************************
  subroutine Prtcl_dump_real3_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    type(real3),dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer::ierror,newtype
    integer,dimension(1)::sizes,subsizes,starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant real3 vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real3_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp+sizes(1)*real3_byte
  end subroutine Prtcl_dump_real3_vector

  !**********************************************************************
  ! Prtcl_dump_real3_matrix
  !**********************************************************************
  subroutine Prtcl_dump_real3_matrix(fh,disp,var,pmsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_mat),intent(in)::pmsize    
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    type(real3),dimension(1:pmsize%subsizes(1),1:pmsize%subsizes(2)),intent(in)::var

    ! locals
    integer::ierror,newtype
    integer,dimension(2)::sizes,subsizes,starts

    ! calculate sizes, subsizes and starts
    sizes     = pmsize%sizes
    subsizes  = pmsize%subsizes
    starts    = pmsize%starts

    ! write the particle relevant real matrix
    call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1)*subsizes(2),real3_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp+sizes(1)*sizes(2)*real3_byte
  end subroutine Prtcl_dump_real3_matrix
end module Prtcl_IOAndVisu
module Prtcl_NBS_Munjiza
  use m_TypeDef
  use m_LogInfo
  use Prtcl_CL_and_CF
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Parameters
  use Prtcl_decomp_2d
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private
  
  real(RK)::xst_cs,yst_cs,zst_cs
  integer::nFixed,nFixedP,nlocal,nFixedAndlocal,nFixedAndlocalP,nNeedCSTotoal
  
  type(integer3),dimension(:),allocatable:: box_index ! integer coordinate of box
  integer,dimension(:),allocatable :: NextX
  integer,dimension(:),allocatable :: NextY
  integer,dimension(:),allocatable :: NextZ
  integer,dimension(:),allocatable :: HeadY
  integer,dimension(:),allocatable :: HeadX
  integer,dimension(:),allocatable :: HeadX0
  integer,dimension(:,:),allocatable :: HeadZ  ! head list for iy, current row 
  integer,dimension(:,:),allocatable :: HeadZ0 ! head list for (iy-1), lower row
    
  type::NBS_Munjiza
    integer:: mbox
    real(RK)::maxDiam
    real(RK)::cell_len   ! length of cell
    integer:: nx         ! number of divisions in x direction
    integer:: ny         ! number of divisions in y direction
    integer:: nz         ! number of divisions in z direction       
    integer:: num_Cnsv_cntct = 0 !number of conservative contacts in the broad search phase
  contains
    procedure:: Init_NBSM
        
    ! performing contact search (includes all steps)
    procedure:: ContactSearch=> NBSM_ContactSearch
    procedure:: clcBoxIndex  => NBSM_clcBoxIndex ! calculating integer coordinates of all boxes
    procedure:: BuildYList   => NBSM_BuildYList  ! constructing YList
    procedure:: BuildXList   => NBSM_BuildXList  ! constructing XList
    procedure:: BuildZList   => NBSM_BuildZList  ! constructing ZList
    procedure:: BuildZList0  => NBSM_BuildZList0
    procedure:: FineSearch1  => NBSM_FineSearch1
    procedure:: FineSearch2  => NBSM_FineSearch2
    procedure:: FineSearch3  => NBSM_FineSearch3
    procedure:: LoopNBSMask  => NBSM_LoopNBSMask
  end type NBS_Munjiza
  type(NBS_Munjiza),public,allocatable:: m_NBS_Munjiza
    
contains
  !******************************************************************
  ! Initializing NBS_Munjiza object
  !******************************************************************
  subroutine Init_NBSM(this)
    implicit none
    class(NBS_Munjiza)::this

    ! locals
    type(integer3)::numCell
    real(RK)::xed_cs,yed_cs,zed_cs
    integer::i,iErr1,iErr2,iErr3,iErr4,iErr5,iErr6,iErr7,iErr8,iErr9,iErrSum

    this%maxDiam = 2.0_RK*maxval( DEMProperty%Prtcl_PureProp%Radius )
#ifndef CFDACM
    this%cell_len= DEM_Opt%Prtcl_cs_ratio*this%maxDiam
#else
    this%cell_len= DEM_Opt%Prtcl_cs_ratio*this%maxDiam +maxval(dlub_pp)
#endif

    xst_cs = DEM_decomp%xSt - this%cell_len*1.05_RK
    yst_cs = DEM_decomp%ySt - this%cell_len*1.05_RK
    zst_cs = DEM_decomp%zSt - this%cell_len*1.05_RK
    xed_cs = DEM_decomp%xEd + this%cell_len*1.05_RK
    yed_cs = DEM_decomp%yEd + this%cell_len*1.05_RK 
    zed_cs = DEM_decomp%zEd + this%cell_len*1.05_RK
    this%nx = int((xed_cs-xst_cs)/this%cell_len)+1
    this%ny = int((yed_cs-yst_cs)/this%cell_len)+1
    this%nz = int((zed_cs-zst_cs)/this%cell_len)+1
        
    numcell = integer3(this%nx, this%ny, this%nz)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Contact search method is NBS Munjiza", 3 )
      call DEMLogInfo%OutInfo("Cell size is [m]: "// trim(num2str(this%cell_len)), 4)
      call DEMLogInfo%OutInfo("Number of cells considered is (x,y,z) :"//trim(num2str(numCell)),4) 
    endif       
        
    nFixed= GPrtcl_list%mlocalFix+ GPrtcl_list%nGhostFix_CS; nFixedP=nFixed+1
    this%mbox = nFixed+ GPrtcl_list%mlocal + GPrtcl_list%mGhost_CS

    allocate(box_index(this%mbox),   STAT=iErr1) 
    allocate(HeadY(this%ny),         STAT=iErr2)
    allocate(HeadX(this%nx),         STAT=iErr3)
    allocate(HeadX0(this%nx),        STAT=iErr4)
    allocate(HeadZ(0:1,0:this%nz+1), STAT=iErr5)
    allocate(HeadZ0(0:2,0:this%nz+1),STAT=iErr6)
    allocate(NextY(this%mbox),       STAT=iErr7)
    allocate(NextX(this%mbox),       STAT=iErr8)
    allocate(NextZ(this%mbox),       STAT=iErr9)

    iErrSum=abs(iErr1)+abs(iErr2)+abs(iErr3)+abs(iErr4)+abs(iErr5)+abs(iErr6)+abs(iErr7)+abs(iErr8)+abs(iErr9)
    if(iErrSum/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"Init_NBSM: ","Allocation failed " )
        
    HeadY = -1; NextY = -1
    HeadX = -1; NextX = -1; HeadX0 = -1
    HeadZ = -1; NextZ = -1; HeadZ0 = -1

    ! Firstly, clculate the Box_index for fixed particles(including the relevant ghost fixed particles)
    do i=1,nFixed
      box_index(i)%x = floor(( GPFix_PosR(i)%x - xst_cs )/this%cell_len)+1
      box_index(i)%y = floor(( GPFix_PosR(i)%y - yst_cs )/this%cell_len)+1
      box_index(i)%z = floor(( GPFix_PosR(i)%z - zst_cs )/this%cell_len)+1
    enddo
  end subroutine Init_NBSM

  !******************************************************************
  ! performing a contact search on all particles (includes all steps)
  !******************************************************************
  subroutine NBSM_ContactSearch(this)    
    implicit none
    class(NBS_Munjiza):: this      

    ! locals
    integer::ix,iy,iz

    nlocal=GPrtcl_list%nlocal
    this%num_Cnsv_cntct = 0
    if(nlocal == 0) return
    nFixedAndlocal  = nFixed + nlocal
    nFixedAndlocalP = nFixedAndlocal + 1
    nNeedCSTotoal = nFixedAndlocal + GPrtcl_list%nGhost_CS

    call this%clcBoxIndex()
    call this%BuildYList()

    HeadX0 = -1
    DO iy = 1, this%ny
      call this%BuildXList(iy)
      if(HeadY(iy).ne. -1 ) then
        HeadZ(0,:)  = -1
        HeadZ0(0,:) = -1
        call this%BuildZlist0(1,1)
        do ix = 1, this%nx
          call this%BuildZlist(ix,1)
          call this%BuildZlist0(ix,2)

          if(HeadX(ix).ne.-1) then
            do iz=1,this%nz
              call this%LoopNBSMask(iz)
            enddo
          endif
          HeadZ(0,:) = HeadZ(1,:)  ! same row, subs
          HeadZ0(0,:)= HeadZ0(1,:) ! lower row, subs
          HeadZ0(1,:)= HeadZ0(2,:)
        enddo
      endif
      HeadX0 = HeadX
    ENDDO
  end subroutine NBSM_ContactSearch

  !******************************************************************
  ! calculating integer coordinates of all boxes
  !******************************************************************
  subroutine NBSM_clcBoxIndex(this)
    implicit none
    class(NBS_Munjiza):: this

    ! locals
    real(RK)::rpdx,rpdy,rpdz
    integer::i,m,sizen,ierrTmp,ierror=0
    type(integer3),allocatable,dimension(:)::Int3Vec
  
    if(nNeedCSTotoal>this%mbox) then
      sizen= int(1.2_RK*real(this%mbox,kind=RK))
      sizen= max(sizen, nNeedCSTotoal+1)

      call move_alloc(box_index,Int3Vec)
      allocate(box_index(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      if(nFixed>0)box_index(1:nFixed)= Int3Vec(1:nFixed)
      deallocate(Int3Vec)

      deallocate(NextX); allocate(NextX(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      deallocate(NextY); allocate(NextY(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      deallocate(NextZ); allocate(NextZ(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      if(ierror/=0) then
        call DEMLogInfo%CheckForError(ErrT_Abort," NBSM_clcBoxIndex"," Reallocate wrong!")
        call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
      endif      
      !call DEMLogInfo%CheckForError(ErrT_Pass," NBSM_clcBoxIndex"," Need to reallocate Box_And_Next")
      !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),    3)
      !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(this%mbox)),3)
      !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),    3)
      this%mbox=sizen
    endif

    m=1
    do i=nFixed+1, nFixedAndlocal
      rpdx = GPrtcl_PosR(m)%x - xst_cs
      rpdy = GPrtcl_PosR(m)%y - yst_cs
      rpdz = GPrtcl_PosR(m)%z - zst_cs
      box_index(i)%x = floor(rpdx/this%cell_len)+1
      box_index(i)%y = floor(rpdy/this%cell_len)+1
      box_index(i)%z = floor(rpdz/this%cell_len)+1
      m=m+1
    enddo

    m=1
    do i=nFixedAndlocalP, nNeedCSTotoal
      rpdx = GhostP_PosR(m)%x - xst_cs
      rpdy = GhostP_PosR(m)%y - yst_cs
      rpdz = GhostP_PosR(m)%z - zst_cs
      box_index(i)%x = floor(rpdx/this%cell_len)+1
      box_index(i)%y = floor(rpdy/this%cell_len)+1
      box_index(i)%z = floor(rpdz/this%cell_len)+1
      m=m+1    
    enddo
  end subroutine NBSM_clcBoxIndex
    
  !******************************************************************
  ! constructing Ylist of particles
  !******************************************************************
  subroutine NBSM_BuildYList(this)
    implicit none
    class(NBS_Munjiza) this
    integer:: i,iy

    HeadY = -1  ! nullifying list Y
    do i=1,nNeedCSTotoal  
      iy=box_index(i)%y
      NextY(i) = HeadY(iy) 
      HeadY(iy)= i
    enddo
  end subroutine NBSM_BuildYList

  !******************************************************************
  ! constructing Xlist of row iy 
  !******************************************************************
  subroutine NBSM_BuildXList(this,iy)
    implicit none
    class(NBS_Munjiza)::this
    integer,intent(in)::iy ! row index
    integer::n,ix

    ! nullifying the xlist of current row but keeps the previous raw
    HeadX= -1
        
    n = HeadY(iy)
    do while (n .ne. -1)
      ix = box_index(n)%x
      NextX(n) = HeadX(ix)
      HeadX(ix)= n  
      n = NextY(n)
    enddo
  end subroutine NBSM_BuildXList

  !*********************************************************************
  !   Constructing the ZList of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSM_BuildZList(this,ix,m)
    implicit none
    class(NBS_Munjiza ) this
    integer,intent(in) :: ix ! col index
    integer,intent(in) :: m  ! the column location with resect to ix
    integer:: n,iz
 
    ! nullifying the zlist of current col ix, but keeps the previous and next cols
    ! 0: previous (left) column, 1: current column, 2: next (right) column
    HeadZ(m,:) = -1
    n = HeadX(ix+m-1) ! reading from current row iy
    do while ( n .ne. -1 )
      iz = box_index(n)%z
      NextZ(n) = HeadZ(m,iz)
      HeadZ(m,iz) = n
      n = NextX(n)
    enddo
  end subroutine NBSM_BuildZList
    
  !*********************************************************************
  !   Constructing the ZList0 of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSM_BuildZList0(this,ix,m)
    implicit none
    class(NBS_Munjiza ):: this
    integer,intent(in):: ix ! column index
    integer,intent(in):: m  ! the column location with respect to ix
    integer:: n,iz,ixm
        
    ! 0: previous (left) column, 1: current column, 2: next (right) column
    HeadZ0(m,:) = -1
    ixm=ix+m-1
    if(ixm>this%nx) return
   
    n = HeadX0(ixm) ! reading from the row below iy (or iy-1)
    do while ( n .ne. -1 )
      iz = box_index(n)%z
      NextZ(n) = HeadZ0(m,iz)
      HeadZ0(m,iz) = n
      n = NextX(n)
    enddo
  end subroutine NBSM_BuildZList0

  !**********************************************************************
  ! finding contacts between particles in the target cell and particles in
  ! cells determined by NBS mask.  
  !**********************************************************************
  subroutine NBSM_LoopNBSMask(this, iz)
    implicit none
    class(NBS_Munjiza) this
    integer,intent(in)  :: iz
    integer m, n, i, lx
  
    m = HeadZ(1,iz)
    DO WHILE(m.ne.-1)
      IF(m<nFixedP) THEN               !==================================

        ! over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch1(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch1(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nFixedAndlocalP .and. n>nFixed ) then
              call this%FineSearch1(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            endif
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nFixedAndlocalP .and. n>nFixed ) then
                call this%FineSearch1(n,m)
                this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              endif
              n = NextZ(n)
            ENDDO
          enddo
        enddo

      ELSEIF(m>nFixedAndlocal) THEN    !==================================

        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch3(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch3(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nFixedAndlocalP .and. n>nFixed ) then
              call this%FineSearch3(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            endif
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nFixedAndlocalP .and. n>nFixed ) then
                call this%FineSearch3(n,m)
                this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              endif
              n = NextZ(n)
            ENDDO
          enddo
        enddo

      ELSE                             !==================================

       !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nFixedP) then
            call this%FineSearch1(m,n)
          elseif(n>nFixedAndlocal) then
            call this%FineSearch3(m,n)
          else
            call this%FineSearch2(m,n)
          endif
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n<nFixedP) then
            call this%FineSearch1(m,n)
          elseif(n>nFixedAndlocal) then
            call this%FineSearch3(m,n)
          else
            call this%FineSearch2(m,n)
          endif
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nFixedP) then
              call this%FineSearch1(m,n)
            elseif(n>nFixedAndlocal) then
              call this%FineSearch3(m,n)
            else
              call this%FineSearch2(m,n)
            endif
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nFixedP) then
                call this%FineSearch1(m,n)
              elseif(n>nFixedAndlocal) then
                call this%FineSearch3(m,n)
              else
                call this%FineSearch2(m,n)
              endif
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              n = NextZ(n)
            ENDDO
          enddo
        enddo
      ENDIF                            !==================================

      m = NextZ(m)
    ENDDO
  end subroutine NBSM_LoopNBSMask

#ifdef CFDACM
  !********************************************************************** 
  ! particle fine search (Moving-particles with Fixed-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch1(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp,drlub

    pid1Temp= pid1- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPFix_PosR(pid2)%w
    drlub= dr+ dlub_pp(GPrtcl_pType(pid1Temp),GPFix_pType(pid2))
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GPFix_PosR(pid2)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPFix_PosR(pid2)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPFix_PosR(pid2)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return

    ovrlp = dr-sqrt(d2sum)  
    if(ovrlp>=0.0_RK) then
      call GPPW_CntctList%AddContactPPFix(pid1Temp,pid2,ovrlp)
    else
      call GPPW_CntctList%AddLubForcePPFix(pid1Temp,pid2,-ovrlp) 
    endif
  end subroutine NBSM_FineSearch1 

  !********************************************************************** 
  ! particle fine search(Moving-particles with Moving-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch2(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp,pid2Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,drlub,ovrlp

    pid1Temp= pid1- nFixed
    pid2Temp= pid2- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPrtcl_PosR(pid2Temp)%w
    drlub= dr+ dlub_pp(GPrtcl_pType(pid1Temp),GPrtcl_pType(pid2Temp))
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GPrtcl_PosR(pid2Temp)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPrtcl_PosR(pid2Temp)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPrtcl_PosR(pid2Temp)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)        

    if(ovrlp>=0.0_RK) then
      ! this is a convention, the lower id should be the first item in the contact pair (particle & particle)
      if(GPrtcl_id(pid1Temp) < GPrtcl_id(pid2Temp) ) then
        call GPPW_CntctList%AddContactPP(pid1Temp,pid2Temp,ovrlp)
      else
        call GPPW_CntctList%AddContactPP(pid2Temp,pid1Temp,ovrlp)
      endif
    else
      call GPPW_CntctList%AddLubForcePP(pid1Temp,pid2Temp,-ovrlp)
    endif
  end subroutine NBSM_FineSearch2

  !********************************************************************** 
  ! particle fine search (Moving-particles with Ghost-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch3(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::gid,pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,drlub,ovrlp

    pid1Temp= pid1- nFixed
    gid     = pid2- nFixedAndlocal
    dr= GPrtcl_PosR(pid1Temp)%w + GhostP_PosR(gid)%w
    drlub= dr+ dlub_pp(GPrtcl_pType(pid1Temp),GhostP_pType(gid))
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GhostP_PosR(gid)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GhostP_PosR(gid)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GhostP_PosR(gid)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)
    
    if(ovrlp>=0.0_RK) then
      call GPPW_CntctList%AddContactPPG(pid1Temp,gid,ovrlp)
    else
      call GPPW_CntctList%AddLubForcePPG(pid1Temp,gid,-ovrlp) 
    endif
  end subroutine NBSM_FineSearch3
#else
  !********************************************************************** 
  ! particle fine search (Moving-particles with Fixed-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch1(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp

    pid1Temp= pid1- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPFix_PosR(pid2)%w
    dr2= dr*dr
    
    dx= GPrtcl_PosR(pid1Temp)%x - GPFix_PosR(pid2)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPFix_PosR(pid2)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPFix_PosR(pid2)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    
    ovrlp = dr-sqrt(d2sum)  
    call GPPW_CntctList%AddContactPPFix(pid1Temp,pid2,ovrlp)
  end subroutine NBSM_FineSearch1 

  !********************************************************************** 
  ! particle fine search(Moving-particles with Moving-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch2(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp,pid2Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp

    pid1Temp= pid1- nFixed
    pid2Temp= pid2- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPrtcl_PosR(pid2Temp)%w
    dr2= dr*dr
    dx= GPrtcl_PosR(pid1Temp)%x - GPrtcl_PosR(pid2Temp)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPrtcl_PosR(pid2Temp)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPrtcl_PosR(pid2Temp)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)        

    ! this is a convention, the lower id should be the first item in the contact pair (particle & particle)
    if(GPrtcl_id(pid1Temp) < GPrtcl_id(pid2Temp) ) then
      call GPPW_CntctList%AddContactPP(pid1Temp,pid2Temp,ovrlp)
    else
      call GPPW_CntctList%AddContactPP(pid2Temp,pid1Temp,ovrlp)
    endif
  end subroutine NBSM_FineSearch2 

  !********************************************************************** 
  ! particle fine search (Moving-particles with Ghost-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch3(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::gid,pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp

    pid1Temp= pid1- nFixed
    gid     = pid2- nFixedAndlocal
    dr= GPrtcl_PosR(pid1Temp)%w + GhostP_PosR(gid)%w
    dr2= dr*dr
    dx= GPrtcl_PosR(pid1Temp)%x - GhostP_PosR(gid)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GhostP_PosR(gid)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GhostP_PosR(gid)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)  
    call GPPW_CntctList%AddContactPPG(pid1Temp,gid,ovrlp)
  end subroutine NBSM_FineSearch3 
#endif
end module Prtcl_NBS_Munjiza
module Prtcl_Parameters
  use m_TypeDef
  use m_LogInfo
#ifdef CFDDEM
  use m_Parameters,only: dtMax,ifirst,ilast,BackupFreq,SaveVisu,xlx,yly,zlz,BcOption
#endif
#ifdef CFDACM
  use m_Decomp2d,only: nrank
  use m_Parameters,only: dtMax,ifirst,ilast,BackupFreq,SaveVisu,xlx,yly,zlz,BcOption,ischeme,FI_AB2,FI_RK2,FI_RK3
#endif
  implicit none
  private

  ! Log 
  type(LogType),public::DEMLogInfo
  
  integer,parameter,public:: x_axis = 1
  integer,parameter,public:: y_axis = 2
  integer,parameter,public:: z_axis = 3    
    
  integer,parameter,public:: CSM_NBS_Munjiza       = 1
  integer,parameter,public:: CSM_NBS_Munjiza_Hrchl = 2
  integer,parameter,public:: DEM_LSD = 1
  integer,parameter,public:: DEM_nLin= 2
#ifdef CFDACM
  integer,parameter,public:: ACM_LSD = 3
  integer,parameter,public:: ACM_nLin= 4
#endif
  integer,parameter,public:: PIM_FE  = 1
  integer,parameter,public:: PIM_AB2 = 2
  integer,parameter,public:: PIM_AB3 = 3

#ifdef CFDACM
  real(RK),public::Klub_pp, Klub_pw,Lub_ratio,Ndt_coll,St_Crit
  logical,public::UpdateACMflag,IsDryColl,IsAddFluidPressureGradient
  integer,public::icouple,nForcingExtra,IBM_Scheme,idem_advance_start(3),idem_advance_end(3)
#endif
#ifdef CFDDEM
  integer,public::  icouple
  real(RK),public:: RatioSR                            ! search radius ratio
  real(RK),public:: FluidAccCoe,SaffmanConst
  logical,public::  IsAddFluidPressureGradient,UpdateDEMflag
  logical,public::  is_clc_Lift, is_clc_Basset,is_clc_Basset_fixed, is_clc_ViscousForce, is_clc_PressureGradient,is_clc_FluidAcc
  integer,public::  mWinBasset, mTailBasset, BassetAccuracy, BassetTailType 
  type BassetDataSeq
    integer:: HistoryStage=0
    integer:: HistStageFix=0
    integer:: nDataLen
    integer:: iWindowStart
    integer:: iWindowEnd
    integer:: iTailStart
    integer:: iTailEnd
  end type BassetDataSeq
  type(BassetDataSeq),public:: GPrtcl_BassetSeq
#endif
   
  ! default values  
  type DEMS_Options
    character(64)::RunName  = "DEMRun" ! run name
    character(64)::ResultsDir= "."     ! result directory 
    character(64)::RestartDir="."      ! restart directory
    ! If GeometrySource =2, please give a STL file routine, if not, just ignore this variable      
    character(64):: Geom_Dir ="DEMGeom.stl" 
    
    logical:: RestartFlag=.false.
    logical,dimension(3):: IsPeriodic = .false.
    
    integer:: numPrtcl    = 8000     ! total particle number
    integer:: numPrtclFix = 1000     ! total fixed particle number
    integer:: np_InDomain            ! particle in domain
    integer:: ifirst                 ! first time step
    integer:: ilast                  ! last time step
    integer:: CS_Method = CSM_NBS_Munjiza      ! contact search method        
    integer:: CF_Type   = DEM_LSD              ! contact force type
    integer:: PI_Method = PIM_FE   ! integration scheme for translational motion
    integer:: PRI_Method = PIM_FE  ! integration scheme for rotational motion  
    integer:: numPrtcl_Type=1      ! number of particle type
    integer:: numWall_type =1      ! number of wall type
    ! contact list size, 6 means every particle can contact with 12 neighbour particles/walls in average 
    integer:: CntctList_Size = 6
    ! means default behavior, number of levels in multi-level contact search 
    integer:: CS_numlvls = 0  
    ! Global wall id will start at (Base_wall_id+1). Please make sure that  Base_wall_id≥numPrtcl    
    integer:: Base_wall_id = 1000000000  
    ! Near wall list will be updated no more than every 100 iterations
    integer:: Wall_max_update_iter = 100
    integer:: SaveVisu      = 1000     ! save frequency for visulizing file
    integer:: BackupFreq    = 100000   ! save frequency for restarting file
    integer:: Cmd_LFile_Freq= 500      ! report frequency in the terminal 
    integer:: LF_file_lvl   = 5        ! logfile report level      
    integer:: LF_cmdw_lvl   = 3        ! terminal report level
    ! Where is the geometry from? 0: Added directly in the program, 1: From DEM.prm, 2: From external STL file
    integer:: GeometrySource =0
                            
    real(RK):: dt   =  1.0E-5_RK     ! time step       
    real(RK):: Prtcl_cs_ratio =1.0_RK
    ! The particle withthin 2*MaxRadius, will be considered into the NEAR WALL LIST 
    real(RK):: Wall_neighbor_ratio = 2.0_RK 
        
    type(real3):: gravity = real3(0.0_RK,-9.81_RK,0.0_RK) ! gravity or other constant body forces if any
    type(real3):: SimDomain_min
    type(real3):: SimDomain_max
  contains 
    procedure :: ReadDEMOption => DO_ReadDEMOption
  end type DEMS_Options
  type(DEMS_Options),public::  DEM_opt
contains

  !**********************************************************************
  ! DO_ReadDEMOption
  !**********************************************************************
  subroutine DO_ReadDEMOption(this, chFile)
    implicit none
    class(DEMS_Options):: this
    character(*),intent(in)::chFile
           
    ! locals
    logical::RestartFlag
    character(64)::RunName,ResultsDir,RestartDir,Geom_Dir
    real(RK)::Wall_neighbor_ratio,Prtcl_cs_ratio,gravity(3)
    integer::Wall_max_update_iter,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl
    integer::numPrtcl,numPrtclFix,CS_Method,CF_Type,PI_Method,PRI_Method,GeometrySource
    integer::numPrtcl_Type,numWall_type,CntctList_Size,CS_numlvls,nUnitFile,ierror
        
    logical::IsPeriodic(3)=.false.
    real(RK)::dtDEM,minpoint(3),maxpoint(3)
    integer:: ifirstDEM,ilastDEM,BackupFreqDEM,SaveVisuDEM
#if defined(CFDDEM) || defined(CFDACM)
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl,numPrtclFix,gravity,CS_Method,CF_Type,PI_Method,PRI_Method,   &
                          numPrtcl_Type,numWall_type,CS_numlvls,CntctList_Size,Wall_max_update_iter,RunName, &
                          Wall_neighbor_ratio,ResultsDir,RestartDir,Geom_Dir,Cmd_LFile_Freq,LF_file_lvl,     &
                          LF_cmdw_lvl,GeometrySource,Prtcl_cs_ratio
#else
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl,numPrtclFix,dtDEM,gravity,minpoint,maxpoint,CS_Method,CF_Type, &
                          PI_Method,PRI_Method,numPrtcl_Type,numWall_type,CS_numlvls,CntctList_Size,RunName,  &
                          Wall_max_update_iter,Wall_neighbor_ratio,ResultsDir,RestartDir,BackupFreqDEM,       &
                          SaveVisuDEM,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,GeometrySource,Geom_Dir,ifirstDEM, &
                          ilastDEM,Prtcl_cs_ratio,IsPeriodic
#endif

#ifdef CFDDEM
    NAMELIST/CFDDEMCoupling/icouple,UpdateDEMflag,is_clc_Lift,is_clc_Basset,is_clc_Basset_fixed,is_clc_ViscousForce, &
                            is_clc_PressureGradient,is_clc_ViscousForce,is_clc_FluidAcc,FluidAccCoe,SaffmanConst,    &
                            RatioSR,IsAddFluidPressureGradient
    NAMELIST/BassetOptions/ mWinBasset, mTailBasset, BassetAccuracy, BassetTailType 
#endif
#ifdef CFDACM
    integer::icouple_sub
    NAMELIST/CFDACMCoupling/UpdateACMflag,icouple,nForcingExtra,IBM_Scheme,Klub_pp,Klub_pw,Lub_ratio, &
                            Ndt_coll,IsDryColl,St_Crit,IsAddFluidPressureGradient
#endif
              
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) then
       print*, "Cannot open file: "//trim(adjustl(chFile)); STOP
    endif
    read(nUnitFile, nml=DEMOptions)

#ifdef CFDDEM 
    rewind(nUnitFile)
    read(nUnitFile, nml=CFDDEMCoupling)
    rewind(nUnitFile)
    read(nUnitFile, nml=BassetOptions )
    dtDEM     =  dtMax/real(icouple,kind=RK)
    ifirstDEM =  icouple*(ifirst-1)+1
    ilastDEM  =  icouple* ilast
    BackupFreqDEM = icouple* BackupFreq
    SaveVisuDEM   = icouple* SaveVisu
    minpoint= (/0.0_RK, 0.0_RK,0.0_RK /)
    maxpoint= (/xlx, yly, zlz/)
    if(BcOption(1)==0)IsPeriodic(1)=.true.
    if(BcOption(3)==0)IsPeriodic(2)=.true.
    if(BcOption(5)==0)IsPeriodic(3)=.true.

    ! Basset history part
    GPrtcl_BassetSeq%HistoryStage= 0
    GPrtcl_BassetSeq%nDataLen    = 1+ (mWinBasset+2) + mTailBasset
    GPrtcl_BassetSeq%iWindowStart= 2
    GPrtcl_BassetSeq%iWindowEnd  = GPrtcl_BassetSeq%iWindowStart+ (mWinBasset+2)-1
    GPrtcl_BassetSeq%iTailStart  = GPrtcl_BassetSeq%iWindowEnd  + 1
    GPrtcl_BassetSeq%iTailEnd    = GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    rewind(nUnitFile)
    read(nUnitFile, nml=CFDACMCoupling)
    if(ischeme==FI_RK2 .and. mod(icouple, 2)/=0 .and. nrank==0) then
      print*,"icouple WRONG!!! icouple must be multiple of  2 for FI_RK2 scheme"; STOP
    endif
    if(ischeme==FI_RK3 .and. mod(icouple,15)/=0 .and. nrank==0) then
      print*,"icouple WRONG!!! icouple must be multiple of 15 for FI_RK3 scheme"; STOP
    endif    
    if(ischeme==FI_AB2) then
      idem_advance_start(1)=1
      idem_advance_end(1)  =icouple
    elseif(ischeme==FI_RK2) then
      icouple_sub= icouple/2
      idem_advance_start(1)=1
      idem_advance_start(2)=icouple_sub+1 
      idem_advance_end(1)  =icouple_sub
      idem_advance_end(2)  =icouple_sub*2
    elseif(ischeme==FI_RK3) then
      icouple_sub= icouple/15
      idem_advance_start(1)=1
      idem_advance_start(2)=icouple_sub*8  + 1 
      idem_advance_start(3)=icouple_sub*10 + 1
      idem_advance_end(1)  =icouple_sub*8
      idem_advance_end(2)  =icouple_sub*10
      idem_advance_end(3)  =icouple_sub*15
    endif
    dtDEM     =  dtMax/real(icouple,kind=RK)
    ifirstDEM =  icouple*(ifirst-1)+1
    ilastDEM  =  icouple* ilast
    BackupFreqDEM = icouple* BackupFreq
    SaveVisuDEM   = icouple* SaveVisu
    minpoint= (/0.0_RK, 0.0_RK,0.0_RK /)
    maxpoint= (/xlx, yly, zlz/)
    if(BcOption(1)==0)IsPeriodic(1)=.true.
    if(BcOption(3)==0)IsPeriodic(2)=.true.
    if(BcOption(5)==0)IsPeriodic(3)=.true.
#endif
    close(nUnitFile,IOSTAT=ierror)
           
    this%RestartFlag = RestartFlag
    this%numPrtcl    = numPrtcl
    this%numPrtclFix = numPrtclFix
    this%dt       = dtDEM
    this%ifirst   = ifirstDEM
    this%ilast    = ilastDEM
    this%gravity  = gravity
    this%SimDomain_min = minpoint
    this%SimDomain_max = maxpoint
    this%IsPeriodic = IsPeriodic
           
    this%Prtcl_cs_ratio = Prtcl_cs_ratio
    this%CS_Method = CS_Method
    this%CF_Type   = CF_Type
    this%PI_Method = PI_Method  
    this%PRI_Method= PRI_Method
           
    this%numPrtcl_Type= numPrtcl_Type
    this%numWall_type = numWall_type
           
    this%CntctList_Size= CntctList_Size
    this%CS_numlvls    = CS_numlvls
           
    this%Base_wall_id =  max(numPrtcl+numPrtclFix, 1000000000) 
    this%Wall_max_update_iter= Wall_max_update_iter
    this%Wall_neighbor_ratio = Wall_neighbor_ratio 
           
    write(this%RunName,"(A)") RunName
    write(this%ResultsDir,"(A)") ResultsDir
    write(this%RestartDir,"(A)") RestartDir
    this%BackupFreq = BackupFreqDEM
    this%SaveVisu = SaveVisuDEM
    this%Cmd_LFile_Freq = Cmd_LFile_Freq
    this%LF_file_lvl = LF_file_lvl
    this%LF_cmdw_lvl = LF_cmdw_lvl
    this%GeometrySource = GeometrySource
    if( this%GeometrySource==2 ) write(this%Geom_Dir,"(A)")Geom_Dir

  end subroutine DO_ReadDEMOption
end module Prtcl_Parameters
module Prtcl_Property
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Parameters
#ifdef CFDDEM
  use m_Decomp2d,only: nrank
  use m_Parameters,only: FluidDensity
#elif  CFDACM
  use Prtcl_EqualSphere
  use m_Decomp2d,only: nrank
  use m_Parameters,only: FluidDensity,xnu
  use m_MeshAndMetries,only: dx,dyUniform,dz
#else
  use Prtcl_decomp_2d,only: nrank
#endif
  implicit none
  private
    
  type PureProperty
    real(RK):: Radius = 0.005_RK
    real(RK):: Density= 2500.0_RK
    real(RK):: PoissonRatio= 0.25_RK
    real(RK):: YoungsModulus=5.0E6_RK
    real(RK):: Mass
    real(RK):: Volume
    real(RK):: Inertia
#ifdef CFDDEM
    real(RK):: MassInFluid
    real(RK):: MassOfFluid
#endif
  end type PureProperty

  type,public:: BinaryProperty
    real(RK):: RadEff       ! effective Radiusius
    real(RK):: MassEff      ! effective Mass

    real(RK):: StiffnessCoe_n
    real(RK):: StiffnessCoe_t        
    real(RK):: DampingCoe_n = 0.0_RK
    real(RK):: DampingCoe_t = 0.0_RK
    real(RK):: RestitutionCoe_n = 0.95_RK  ! Normal Resitution Coefficient
    real(RK):: FrictionCoe_s = 0.80_RK     ! Coefficient of static  friction
    real(RK):: FrictionCoe_k = 0.15_RK     ! Coefficient of kinetic friction
#ifdef CFDACM
    real(RK):: Vel_Crit                    ! Turn off fluid forces for large St collsions
    real(RK):: Kn_Grav                
#endif
  end type BinaryProperty
    
  type PhysicalProperty
    integer,allocatable,dimension(:) :: nPrtcl_in_Bin
    integer,allocatable,dimension(:) :: CS_Hrchl_level ! level for NBS-Munjiza-Hierarchy Contact Search Method  
    type(pureProperty),allocatable,dimension(:) :: Prtcl_PureProp
    type(pureProperty),allocatable,dimension(:) :: Wall_PureProp
    type(BinaryProperty),allocatable,dimension(:,:):: Prtcl_BnryProp    
    type(BinaryProperty),allocatable,dimension(:,:):: PrtclWall_BnryProp         
  contains
    procedure:: InitPrtclProperty
    procedure:: InitWallProperty
  end type PhysicalProperty
  type(PhysicalProperty),public::DEMProperty
  
#ifdef CFDACM
  type IBMProperty
    integer::  nPartition
    real(RK):: IBPVolume
    real(RK):: MassIBL
    real(RK):: InertiaIBL
    real(RK):: MassInFluid
    real(RK):: MassEff
    real(RK):: InertiaEff
  end type IBMProperty
  type(IBMProperty),dimension(:),allocatable,public:: PrtclIBMProp
  type(real3),allocatable,dimension(:,:),public::SpherePartitionCoord
  real(RK),allocatable,dimension(:,:),public:: dlub_pp,LubCoe_pp
  real(RK),allocatable,dimension(:),  public:: dlub_pw,LubCoe_pw
#endif

contains
  !*******************************************************
  ! initializing the size distribution with property 
  !*******************************************************
  subroutine InitPrtclProperty(this,chFile)
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK)::FrictionCoe_s_PP,FrictionCoe_k_PP,RestitutionCoe_n_PP
    real(RK),dimension(:),allocatable:: Bin_Divided,Density,Diameter,YoungsModulus_P,PoissonRatio_P
    namelist/ParticlePhysicalProperty/Bin_Divided, Density, Diameter,YoungsModulus_P,PoissonRatio_P, &
                                      FrictionCoe_s_PP,FrictionCoe_k_PP,RestitutionCoe_n_PP
    integer:: i,j,iTV(8),nPType,nUnitFile,ierror,sum_prtcl,bin_pnum,prdiff,bin_id
    real(RK):: sum_divided,rtemp,Radius
    type(PureProperty)::pari,parj
    type(BinaryProperty)::Bnry
#ifdef CFDACM
    integer::nSpSize
    real(RK)::dxyz
    integer,dimension(:),allocatable::  nPartition
    real(RK),dimension(:),allocatable:: RetractionRatio
    real(RK),dimension(:,:),allocatable:: PointTemp
    namelist/ParticleIBMProperty/ nPartition, RetractionRatio
#endif
        
    nPType  = DEM_Opt%numPrtcl_Type
    allocate( Bin_Divided(nPType))
    allocate( Density(nPType))
    allocate( Diameter(nPType))
    allocate( YoungsModulus_P(nPType))
    allocate( PoissonRatio_P(nPType))
        
    allocate( this%CS_Hrchl_level(nPType))
    allocate( this%nPrtcl_in_Bin(nPType))
    allocate( this%Prtcl_PureProp(nPType))
    allocate( this%Prtcl_BnryProp(nPType,nPType)) 
           
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 ) call DEMLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file:"//trim(chFile))
    read(nUnitFile, nml=ParticlePhysicalProperty)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=ParticlePhysicalProperty)
#ifdef CFDACM
    dxyz= (dx*dyUniform*dz)**(0.333333333333333333_RK)
    allocate( nPartition(nPType))
    allocate( RetractionRatio(nPType))
    allocate( PrtclIBMProp(nPType))
    allocate(dlub_pp(nPType,nPType))
    allocate(LubCoe_pp(nPType,nPType))
    rewind(nUnitFile)
    read(nUnitFile, nml=ParticleIBMProperty)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=ParticleIBMProperty)
#endif
    close(nUnitFile,IOSTAT=ierror)
        
    ! calculate this%nPrtcl_in_Bin
    call date_and_time(values=iTV); !iTV=0
    call random_seed(size= i)
    call random_seed(put = iTV(7)*iTV(8)+[(j,j=1,i)])
    sum_divided=0.0_RK
    do i=1,nPType
      sum_divided = sum_divided + Bin_Divided(i)
    enddo
    sum_prtcl = 0
    do i=1,nPType
      bin_pnum = int(DEM_Opt%numPrtcl*Bin_Divided(i)/sum_divided)
      this%nPrtcl_in_Bin(i)= bin_pnum
      sum_prtcl = sum_prtcl + bin_pnum
    enddo
    prdiff = DEM_Opt%numPrtcl - sum_prtcl
    if( prdiff > 0 ) then
      do i=1, prdiff 
        call random_number(rtemp)
        bin_id = int(rtemp*nPType) + 1
        this%nPrtcl_in_Bin(bin_id) = this%nPrtcl_in_Bin(bin_id)  + 1
      enddo
    endif
        
    ! calculate particle properties
    do i = 1, nPType
      Radius= 0.5_RK*Diameter(i)
      this%Prtcl_PureProp(i)%Density= Density(i)
      this%Prtcl_PureProp(i)%Radius = Radius
      this%Prtcl_PureProp(i)%YoungsModulus = YoungsModulus_P(i)
      this%Prtcl_PureProp(i)%PoissonRatio = PoissonRatio_P(i)

      this%Prtcl_PureProp(i)%Volume = 1.333333333333333333_RK*Pi*Radius**3
      this%Prtcl_PureProp(i)%Mass = Density(i)*this%Prtcl_PureProp(i)%Volume
      this%Prtcl_PureProp(i)%Inertia = 0.4_RK*this%Prtcl_PureProp(i)%Mass*Radius**2
#ifdef CFDDEM
      this%Prtcl_PureProp(i)%MassInFluid= (Density(i)-FluidDensity)*this%Prtcl_PureProp(i)%Volume
      this%Prtcl_PureProp(i)%MassOfFluid= FluidDensity*this%Prtcl_PureProp(i)%Volume
#endif
#ifdef CFDACM
      rtemp= Diameter(i)*0.5_RK -RetractionRatio(i)*dxyz
      PrtclIBMProp(i)%nPartition = nPartition(i)
      PrtclIBMProp(i)%IBPVolume= PI*dxyz*(4.0_RK*rtemp*rtemp+ dxyz*dxyz/3.0_RK)/real(nPartition(i),RK)
      PrtclIBMProp(i)%MassIBL  = FluidDensity* PI*dxyz*(4.0_RK*rtemp*rtemp+ dxyz*dxyz/3.0_RK)
      PrtclIBMProp(i)%InertiaIBL=0.6666666666666666667_RK*(rtemp**2)*PrtclIBMProp(i)%MassIBL
      PrtclIBMProp(i)%MassInFluid= (Density(i)-FluidDensity)*this%Prtcl_PureProp(i)%Volume

      select case(IBM_Scheme)
      case(0)   ! 0: Explicit,Uhlmann(2005,JCP)
        PrtclIBMProp(i)%MassEff   = (1.0_RK-FluidDensity/Density(i))*this%Prtcl_PureProp(i)%Mass
        PrtclIBMProp(i)%InertiaEff= (1.0_RK-FluidDensity/Density(i))*this%Prtcl_PureProp(i)%Inertia

      case(1)   ! 1: Explicit,Kempe(2012,JCP)
        PrtclIBMProp(i)%MassEff   = this%Prtcl_PureProp(i)%Mass
        PrtclIBMProp(i)%InertiaEff= this%Prtcl_PureProp(i)%Inertia

      case(2)   ! 2: Semi-implicit,Tschisgale(2017,JCP)
        PrtclIBMProp(i)%MassEff   = this%Prtcl_PureProp(i)%Mass    + PrtclIBMProp(i)%MassIBL
        PrtclIBMProp(i)%InertiaEff= this%Prtcl_PureProp(i)%Inertia + PrtclIBMProp(i)%InertiaIBL

      case default
        call DEMLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","IBM scheme wrong!!!")
      end select
#endif
    enddo
        
    do j=1,nPType
      parj = this%Prtcl_PureProp(j)
      do i=1,nPType
        pari = this%Prtcl_PureProp(i)        
        Bnry= clc_BnryPrtcl_Prop(pari,parj,FrictionCoe_s_PP,FrictionCoe_k_PP,RestitutionCoe_n_PP,.false.)
        this%Prtcl_BnryProp(i,j)= Bnry  
#ifdef CFDACM
        dlub_pp(i,j)= Lub_ratio* dxyz
        LubCoe_pp(i,j)= Klub_pp*xnu*FluidDensity*Bnry%RadEff*Bnry%RadEff/dlub_pp(i,j)
#endif
      enddo
    enddo

#ifdef CFDACM
    ! Sphere partition informations
    nSpSize=maxval(nPartition)
    allocate(SpherePartitionCoord(nSpSize,nPType))
    SpherePartitionCoord=zero_r3
    do i=1,nPType
      allocate(PointTemp(nPartition(i),3))
      call eq_Sphere(PointTemp)
      do j=1,nPartition(i)
        rtemp= Diameter(i)*0.5_RK -RetractionRatio(i)*dxyz
        SpherePartitionCoord(j,i)=real3(PointTemp(j,1),PointTemp(j,2),PointTemp(j,3))*rtemp
      enddo
      deallocate(PointTemp)
    enddo
#endif
  end subroutine InitPrtclProperty

  !*****************************************************************************    
  ! setting the physical property of wall 
  !*****************************************************************************
  subroutine InitWallProperty(this,chFile)
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    type(PureProperty)pari,wall
    integer::i,j,nPType,nWType,nUnitFile,ierror
    real(RK)::FrictionCoe_s_PW,FrictionCoe_k_PW,RestitutionCoe_n_PW
    real(RK),dimension(:),allocatable:: YoungsModulus_W,PoissonRatio_W
    namelist /WallPhysicalProperty/YoungsModulus_W,PoissonRatio_W,FrictionCoe_s_PW,FrictionCoe_k_PW,RestitutionCoe_n_PW
            
    nPType= DEM_Opt%numPrtcl_Type
    nWType= DEM_Opt%numWall_type
    allocate(YoungsModulus_W(nWType))
    allocate(PoissonRatio_W(nWType))
    allocate(this%Wall_PureProp(nWType))
    allocate(this%PrtclWall_BnryProp(nPType,nWType))
#ifdef CFDACM 
    allocate(dlub_pw(nPType))
    allocate(LubCoe_pw(nPType))
#endif
   
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"InitWallProperty: " ,"Cannot open file:"//trim(chFile))
    endif
    read(nUnitFile, nml=WallPhysicalProperty)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=WallPhysicalProperty)
    close(nUnitFile)        
        
    do i=1,nWType
       this%Wall_PureProp(i)%Density = 1.0E50_RK
       this%Wall_PureProp(i)%Radius  = 1.0E50_RK
       this%Wall_PureProp(i)%YoungsModulus= YoungsModulus_W(i)
       this%Wall_PureProp(i)%PoissonRatio = PoissonRatio_W(i)
           
       this%Wall_PureProp(i)%Mass    = 1.0E50_RK
       this%Wall_PureProp(i)%Volume  = 1.0E50_RK
       this%Wall_PureProp(i)%Inertia = 1.0E50_RK
    enddo

    do i = 1, nPType      
      pari = this%Prtcl_PureProp(i)  
      do j = 1, nWType
        wall = this%Wall_PureProp(j)     
        this%PrtclWall_BnryProp(i,j)=clc_BnryPrtcl_Prop(pari,wall,FrictionCoe_s_PW,FrictionCoe_k_PW,RestitutionCoe_n_PW,.true.)
      enddo
#ifdef CFDACM
      dlub_pw(i)  = Lub_ratio*(dx*dyUniform*dz)**(0.333333333333333333_RK)
      LubCoe_pw(i)= Klub_pw*xnu*FluidDensity*(this%Prtcl_PureProp(i)%Radius)**2/dlub_pw(i)
#endif
    enddo       
  end subroutine InitWallProperty

  !*****************************************************************************    
  ! Calculating the binary contact properties
  !*****************************************************************************
  function clc_BnryPrtcl_Prop(pari,parj,FrictionCoe_s,FrictionCoe_k,RestitutionCoe_n,iswall) result(Bnry)
    implicit none
    class(PureProperty),intent(in):: pari, parj
    real(RK),intent(in)::FrictionCoe_s,FrictionCoe_k,RestitutionCoe_n
    logical,intent(in) ::iswall

    ! locals
    type(BinaryProperty):: Bnry
    real(RK)::kappa,pri,prj,ShearModi,ShearModj,YoungsModEff,ShearModEff,Eta,K_hertz
#ifdef CFDACM
    real(RK)::DensityEff,Vel_Crit,Kn_Grav,GravityMag,GRAV_OVERLAP,Lamda,Beta,rA,rB,rC,AlphaTau2,Tau_c0
#endif
    
    if(.not.iswall) then
      Bnry%RadEff = (pari%Radius*parj%Radius)/(pari%Radius +parj%Radius)
      Bnry%MassEff= (pari%Mass  *parj%Mass  )/(pari%Mass   +parj%Mass )
    else
      Bnry%RadEff = pari%Radius
      Bnry%MassEff= pari%Mass         
    endif
    Bnry%RestitutionCoe_n = RestitutionCoe_n
    Bnry%FrictionCoe_s = FrictionCoe_s
    Bnry%FrictionCoe_k = FrictionCoe_k

    pri=pari%PoissonRatio
    prj=parj%PoissonRatio
    ShearModi=pari%YoungsModulus/(2.0_RK*(1.0_RK+pri))
    ShearModj=parj%YoungsModulus/(2.0_RK*(1.0_RK+prj))
    YoungsModEff=1.0_RK/((1.0_RK-pri*pri)/pari%YoungsModulus + (1.0_RK-prj*prj)/parj%YoungsModulus) ! 2.47, p31
    ShearModEff =1.0_RK/((2.0_RK- pri)/ShearModi + (2.0_RK-prj)/ShearModj)                          ! 2.71, p43
    kappa=((1.0_RK-pri)/ShearModi+(1.0_RK-prj)/ShearModj)/((1.0_RK-0.5_RK*pri)/ShearModi+(1.0_RK-0.5_RK*prj)/ShearModj)
    kappa=abs(kappa)
    Eta=log(RestitutionCoe_n); Eta=Eta*Eta
    if(DEM_Opt%CF_Type == DEM_LSD) then                                             
      Bnry%StiffnessCoe_n = 1.2024_RK*(sqrt(Bnry%MassEff)*(YoungsModEff**2)*Bnry%RadEff)**(0.4_RK)        ! 2.46, p31
      Bnry%DampingCoe_n   =-2.0_RK*log(RestitutionCoe_n)*sqrt(Bnry%MassEff*Bnry%StiffnessCoe_n)/sqrt(PI*PI+Eta) ! 2.44, p31
      Bnry%StiffnessCoe_t = Bnry%StiffnessCoe_n*kappa                                                     ! 2.52, p34
      Bnry%DampingCoe_t   = Bnry%DampingCoe_n*sqrt(kappa)
    elseif(DEM_Opt%CF_Type == DEM_nLin) then
      K_hertz =1.333333333333333333_RK*YoungsModEff*sqrt(Bnry%RadEff)               ! 2.62, P39
      Bnry%StiffnessCoe_n= K_hertz
      Bnry%DampingCoe_n  = -2.2664_RK*log(RestitutionCoe_n)*sqrt(Bnry%MassEff*K_hertz)/sqrt(Eta+10.1354_RK) ! 2.66, p40
      Bnry%StiffnessCoe_t= 5.3333333333333_RK*ShearModEff*sqrt(Bnry%RadEff) ! 2.72, P44
      Bnry%DampingCoe_t  = 0.0_RK ! No equation is considered for tangential damping yet
    endif
#ifdef CFDACM
    if(DEM_Opt%CF_Type == ACM_LSD) then      ! Costa et al./Physics Review E 92,053012 (2015)
      Bnry%StiffnessCoe_n= Bnry%MassEff*(PI*PI+Eta)
      Bnry%DampingCoe_n  =-2.0_RK*Bnry%MassEff*log(RestitutionCoe_n)
      Bnry%StiffnessCoe_t= Bnry%StiffnessCoe_n*kappa
      Bnry%DampingCoe_t  = Bnry%DampingCoe_n*sqrt(kappa)
    elseif(DEM_Opt%CF_Type == ACM_nLin) then ! E. Biegert et al./Journal of Computational Physics 340(2017): 105-127
      rA=0.716_RK; rB=0.830_RK; 
      rC=0.744_RK; Tau_c0=3.218_RK
      AlphaTau2=1.111_RK*1.111_RK +3.218_RK*3.218_RK
      Lamda= (sqrt(0.25_RK*rC*rC*Eta*Eta+AlphaTau2*Eta)-0.5_RK*rC*Eta)/AlphaTau2
      Beta= Tau_c0/sqrt(abs(1.0_RK-rA*Lamda-rB*Lamda*Lamda))
      Bnry%StiffnessCoe_n= Bnry%MassEff*(Beta**(2.5_RK))
      Bnry%DampingCoe_n  = 2.0_RK*Bnry%MassEff*Lamda*Beta
      Bnry%StiffnessCoe_t= Bnry%MassEff*(PI*PI+Eta)*kappa
      Bnry%DampingCoe_t  =-2.0_RK*Bnry%MassEff*log(RestitutionCoe_n)*sqrt(kappa)
    endif
    DensityEff= (pari%Density * parj%Density)/(pari%Density + parj%Density)
    Vel_Crit= 4.5_RK*St_Crit*FluidDensity*xnu/(Bnry%RadEff*DensityEff)
    if(IsDryColl) then
      Bnry%Vel_Crit= Vel_Crit
    else
      Bnry%Vel_Crit= 1.0E100_RK ! If IsDryColl=F, we set a very huge Vel_Crit
    endif

    GRAV_OVERLAP=0.001_RK ! E. Biegert et al./Journal of Computational Physics 340(2017): 105-127
    GravityMag=norm(DEM_Opt%gravity)
    if(DEM_Opt%CF_Type == ACM_LSD) then
      kn_Grav= pari%Mass*GravityMag/(GRAV_OVERLAP*pari%Radius)
      if(.not.iswall) kn_Grav= max(parj%Mass*GravityMag/(GRAV_OVERLAP*parj%Radius),kn_Grav)
    elseif(DEM_Opt%CF_Type == ACM_nLin) then
      kn_Grav= pari%Mass*GravityMag/(GRAV_OVERLAP*pari%Radius)**(1.5_RK)
      if(.not.iswall) kn_Grav= max(parj%Mass*GravityMag/(GRAV_OVERLAP*parj%Radius)**(1.5_RK),kn_Grav)
    endif
    Bnry%kn_Grav=kn_Grav
#endif
  end function clc_BnryPrtcl_Prop
end module Prtcl_Property
!********************************************************************!
!*    file name  : Prtcl_Variables.f90                              *!
!*    module name: Prtcl_Variables                                  *!  
!*                                                                  *!
!*    purpose:                                                      *! 
!*      1) All datas required to represent spherical particles      *!
!*      2) Initialize all the particle variables                    *!
!*                                                                  *!
!*  Author: Zheng Gong           Date: 23:Feb:2020                  *!
!*                                                                  *!
!********************************************************************!

module Prtcl_Variables
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Property
  use Prtcl_decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private

  integer,dimension(:),allocatable,public:: GPrtcl_id
  integer,dimension(:),allocatable,public:: GPrtcl_pType
  integer,dimension(:),allocatable,public:: GPrtcl_usrMark
  type(real4),dimension(:),allocatable,public::   GPrtcl_PosR
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linAcc
  type(real3),dimension(:),allocatable,public::   GPrtcl_theta
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_rotVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_rotAcc
  type(real3),dimension(:),allocatable,public::   GPrtcl_cntctForce
  type(real3),dimension(:),allocatable,public::   GPrtcl_torque

  integer,dimension(:),allocatable,public::       GPFix_id
  integer,dimension(:),allocatable,public::       GPFix_pType
  type(real4),dimension(:),allocatable,public::   GPFix_PosR
#ifdef CFDDEM
  type(real3),dimension(:,:),allocatable,public:: GPFix_VFluid 
  type(real3),dimension(:),  allocatable,public:: GPrtcl_FpForce,GPrtcl_FpForce_old,GPrtcl_linVelOld
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_Vfluid
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_BassetData
  type(real3),dimension(:,:),allocatable,public:: GPFix_BassetData
#endif
#ifdef CFDACM
  type(real3),dimension(:),allocatable,public::   GPrtcl_FpForce
  type(real3),dimension(:),allocatable,public::   GPrtcl_FpTorque
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_FluidIntegrate
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_FluidIntOld
  type(real4),dimension(:),allocatable,public::   GPrtcl_PosOld
  character(len=1),dimension(:),allocatable,public:: GPrtcl_HighSt
#endif

  type VarList
    integer:: nlocal     ! number of particles in local processor
    integer:: mlocal     ! the possible maxium # of particles in local processor
    integer:: mlocalFix  ! number of fixed particles in local processor

    integer:: nGhost_CS
    integer:: mGhost_CS
    integer:: nGhostFix_CS

    integer:: tsize  ! size for translational veloctity and acceleration
    integer:: rsize  ! size for rotational veloctity and acceleration
  contains
    procedure:: AllocateAllVar     => GL_AllocateAllVar
#if !defined(CFDDEM) && !defined(CFDACM)
    procedure:: MakingAllPrtcl     => GL_MakingAllPrtcl
#endif
    procedure:: ReallocatePrtclVar => GL_ReallocatePrtclVar
    procedure:: copy               => GL_copy
  end type VarList
  type(VarList),public::GPrtcl_list 

  ! Ghost particle variables
  integer,dimension(:),allocatable,public:: GhostP_id
  integer,dimension(:),allocatable,public:: GhostP_pType
  type(real4),dimension(:),allocatable,public:: GhostP_PosR
  type(real3),dimension(:),allocatable,public:: GhostP_linVel
  type(real3),dimension(:),allocatable,public:: GhostP_rotVel
    
contains

  !**********************************************************************
  ! GL_AllocateAllVar
  !**********************************************************************
  subroutine GL_AllocateAllVar(this)
    implicit none
    class(VarList)::this

    ! locals
    type(real3)::SimLen
    real(RK)::xst,xed,yst,yed,zst,zed,vol_tot,vol_local
    integer:: mlocal,numPrtcl,mlocalFix,numPrtclFix,ierrTmp,ierror=0

    numPrtcl    = DEM_opt%numPrtcl
    numPrtclFix = DEM_opt%numPrtclFix

    ! step0: determine initial misze
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    ! step1: allocating memory for particles
    if(DEM_opt%PI_Method==PIM_FE) then
      this%tsize=1
    elseif (DEM_opt%PI_Method==PIM_AB2) then
      this%tsize=2
    elseif(DEM_opt%PI_Method==PIM_AB3) then
      this%tsize=3
    endif
    if(DEM_opt%PRI_Method==PIM_FE) then
      this%rsize=1
    elseif (DEM_opt%PRI_Method==PIM_AB2) then
      this%rsize=2
    elseif(DEM_opt%PRI_Method==PIM_AB3) then
      this%rsize=3
    endif 

    SimLen = DEM_Opt%SimDomain_max - DEM_Opt%SimDomain_min
    vol_tot= SimLen%x * SimLen%y * SimLen%z
    
    vol_local =(xed-xst)*(yed-yst)*(zed-zst)
    mlocal = int(vol_local/vol_tot*real(numPrtcl,kind=RK))
    mlocal = int(1.5_RK*mlocal)
    mlocal = min(mlocal, numPrtcl)
    mlocal = max(mlocal, 10)
        
    allocate(GPrtcl_id(mlocal),                Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_pType(mlocal),             Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_usrMark(mlocal),           Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_PosR(mlocal),              Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_linVel(this%tsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_linAcc(this%tsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_theta(mlocal),             Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_rotVel(this%rsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_rotAcc(this%rsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_cntctForce(mlocal),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_torque(mlocal),            Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 1")

    mlocalFix = int(vol_local/vol_tot*real(numPrtclFix,kind=RK))
    mlocalFix = int(1.5_RK*mlocalFix)
    mlocalFix = max(mlocalFix, 10)
    mlocalFix = min(mlocalFix, numPrtclFix)
    allocate(GPFix_id(mlocalFix),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPFix_pType(mlocalFix),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPFix_PosR(mlocalFix), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 2")
#ifdef CFDDEM
    allocate(GPrtcl_FpForce(mlocal),    Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FpForce_old(mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_linVelOld(mlocal),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_Vfluid(2,mlocal),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation Fpforce failed")
    GPrtcl_FpForce     = zero_r3
    GPrtcl_FpForce_old = zero_r3
    GPrtcl_linVelOld   = zero_r3
    GPrtcl_Vfluid      = zero_r3
    if(Is_clc_Basset) then
      allocate(GPrtcl_BassetData(GPrtcl_BassetSeq%nDataLen, mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation GPrtcl_BassetData failed")
      GPrtcl_BassetData= zero_r3
    endif
#endif
#ifdef CFDACM
    allocate(GPrtcl_FpForce(mlocal),         Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FpTorque(mlocal),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FluidIntegrate(2,mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FluidIntOld(2,mlocal),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_HighSt(mlocal),          Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(IBM_Scheme==2) allocate(GPrtcl_PosOld(mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"AllocateAllVar","Allocation Fpforce failed")
    GPrtcl_FpForce  = zero_r3
    GPrtcl_FpTorque = zero_r3
    GPrtcl_FluidIntegrate= zero_r3
    GPrtcl_FluidIntOld   = zero_r3
    GPrtcl_HighSt = "N"
    if(IBM_Scheme==2) GPrtcl_PosOld=zero_r4
#endif

    GPrtcl_id = 0
    GPrtcl_pType = 1
    GPrtcl_usrMark = 1
    GPrtcl_PosR = zero_r4
    GPrtcl_linVel  = zero_r3
    GPrtcl_linAcc  = zero_r3
    GPrtcl_theta   = zero_r3
    GPrtcl_rotVel  = zero_r3
    GPrtcl_rotAcc  = zero_r3
    GPrtcl_cntctForce = zero_r3
    GPrtcl_torque     = zero_r3
    GPFix_id = 0
    GPFix_pType = 1
    GPFix_PosR = zero_r4

    ! step2: Initialize this%mlocal,this%nlocal
    this%nlocal    = 0
    this%mlocal    = mlocal
    this%mlocalFix = mlocalFix
  end subroutine GL_AllocateAllVar

#if !defined(CFDDEM) && !defined(CFDACM)
  !**********************************************************************
  ! GL_MakingAllPrtcl
  !**********************************************************************
  subroutine GL_MakingAllPrtcl(this, chFile)
    implicit none
    class(VarList)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK)::Distance_Ratio,VelMag
    character,dimension(3)::Fill_order
    real(RK),dimension(3)::MkPrtclMinpoint,MkPrtclMaxpoint
    type(real3):: MkPrtclMinpoint_real3,MkPrtclMaxpoint_real3
    NAMELIST /ParticleMakingOption/MkPrtclMinpoint,MkPrtclMaxpoint, Fill_order, Distance_Ratio, VelMag
    
    real(RK),dimension(3):: vel_dir
    integer,dimension(:),allocatable:: sum_bin
    real(RK)::xst,xed,yst,yed,zst,zed,maxRad,rvelt,randt
    type(real3)::cntr,l1_vec,l2_vec,l3_vec,lmin_p,lmax_p,dx,pos1,pos2
    integer::i,j,k,iTV(8),bin_id,numPrtcl,nUnitFile,ierror,nlocal,nlocal_sum
   
    open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl", "Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=ParticleMakingOption)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=ParticleMakingOption)
    close(nUnitFile,IOSTAT=ierror)
    MkPrtclMinpoint_real3%x=max(MkPrtclMinpoint(1),DEM_Opt%SimDomain_min%x)
    MkPrtclMinpoint_real3%y=max(MkPrtclMinpoint(2),DEM_Opt%SimDomain_min%y)
    MkPrtclMinpoint_real3%z=max(MkPrtclMinpoint(3),DEM_Opt%SimDomain_min%z)
    MkPrtclMaxpoint_real3%x=min(MkPrtclMaxpoint(1),DEM_Opt%SimDomain_max%x)
    MkPrtclMaxpoint_real3%y=min(MkPrtclMaxpoint(2),DEM_Opt%SimDomain_max%y)
    MkPrtclMaxpoint_real3%z=min(MkPrtclMaxpoint(3),DEM_Opt%SimDomain_max%z)
    numPrtcl = DEM_opt%numPrtcl
        
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    ! step3: Assign Prtcl_Type, Prtcl_id and Prtcl_PosR
    call Fill_Vectors(Fill_order, l1_vec, l2_vec, l3_vec)
    call date_and_time(values=iTV); !iTV=0
    call random_seed(size= i)
    call random_seed(put = iTV(7)*iTV(8)+[(j,j=1,i)])
    maxRad = maxval( DEMProperty%Prtcl_PureProp%Radius )
    lmin_p = MkPrtclMinpoint_real3 + Distance_Ratio* maxRad*real3(1.0_RK,1.0_RK,1.0_RK)
    lmax_p = MkPrtclMaxpoint_real3 - Distance_Ratio* maxRad*real3(1.0_RK,1.0_RK,1.0_RK)
    dx = Distance_Ratio *2.0_RK*maxRad*real3(1.0_RK,1.0_RK,1.0_RK)
        
    allocate(sum_bin(DEM_opt%numPrtcl_Type))
    sum_bin(1)=DEMProperty%nPrtcl_in_Bin(1)
    do j=2, DEM_opt%numPrtcl_Type
      sum_bin(j)=sum_bin(j-1)+DEMProperty%nPrtcl_in_Bin(j)
    enddo

    cntr = lmin_p       ! start point
    nlocal=0
    do i=1,numPrtcl
      do j=1,DEM_opt%numPrtcl_Type
        if(i<=sum_bin(j)) then
          bin_id = j; exit
        endif
      enddo

      if(cntr%x>=xst.and.cntr%x<xed.and.cntr%y>=yst.and.cntr%y<yed.and.cntr%z>=zst.and.cntr%z<zed) then
        if(nlocal>=this%mlocal) call this%ReallocatePrtclVar(nlocal)
        nlocal=nlocal+1
        GPrtcl_id(nlocal) = i
        GPrtcl_pType(nlocal)= bin_id
        GPrtcl_PosR(nlocal) = cntr
        GPrtcl_PosR(nlocal)%w = DEMProperty%Prtcl_PureProp(bin_id)%Radius
      endif

      cntr = cntr + l1_vec * dx
      if((l1_vec.dot.cntr)>=(l1_vec.dot.lmax_p))then
        cntr = (lmin_p*l1_vec) + ((cntr+dx)*l2_vec)+(cntr*l3_vec)
        if((l2_vec.dot.cntr)>=(l2_vec.dot.lmax_p))then
          cntr = (cntr*l1_vec)+(lmin_p*l2_vec)+((cntr+dx)*l3_vec)
          if((l3_vec.dot.cntr) >= (l3_vec.dot.lmax_p) .and. nrank==0) then
            call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Not enough space for positioning" ) 
          endif
        endif
      endif
    enddo

    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: "," nlocal_sum/= numPrtcl " )
    endif
    this%nlocal = nlocal
    deallocate(sum_bin)

    ! randomize the GPrtcl_PosR
    do k=1,3
      do i=1,this%nlocal
        call random_number(randt)
        j=int(randt*this%nlocal)+1
        pos1=GPrtcl_PosR(i)
        pos2=GPrtcl_PosR(j)
        GPrtcl_PosR(i)=pos2
        GPrtcl_PosR(j)=pos1
      enddo
    enddo
        
    ! step3: assign the remaining variables
    do i=1,this%nlocal
      call random_number(vel_dir)
      vel_dir=2.0_RK*vel_dir-1.0_RK
      rvelt=sqrt(vel_dir(1)*vel_dir(1)+vel_dir(2)*vel_dir(2)+vel_dir(3)*vel_dir(3))
      if(rvelt>1.0E-10_RK) then
        vel_dir(1)=vel_dir(1)/rvelt
        vel_dir(2)=vel_dir(2)/rvelt
        vel_dir(3)=vel_dir(3)/rvelt
      else
        vel_dir(1)= -1.0_RK*l3_vec%x
        vel_dir(2)= -1.0_RK*l3_vec%y
        vel_dir(3)= -1.0_RK*l3_vec%z
      endif
      GPrtcl_linVel(1,i)%x=VelMag*vel_dir(1)
      GPrtcl_linVel(1,i)%y=VelMag*vel_dir(2)
      GPrtcl_linVel(1,i)%z=VelMag*vel_dir(3)
    enddo
#ifdef DEMObliqueCollideDry
    do i=1,this%nlocal
      GPrtcl_PosR(i)%x=(xSt+xEd)*0.5_RK
      GPrtcl_PosR(i)%y=ySt+3.0_RK*DEMProperty%Prtcl_PureProp(GPrtcl_pType(i))%Radius
      GPrtcl_PosR(i)%z=(zSt+zEd)*0.5_RK
      GPrtcl_linVel(1,i)= DEM_opt%Gravity
    enddo
    DEM_opt%Gravity=zero_r3
#endif
  end subroutine GL_MakingAllPrtcl

  !**********************************************************************
  ! Fill order Vector 
  !**********************************************************************
  subroutine Fill_Vectors( fill_order , l1 , l2, l3 )
    implicit none
    character,dimension(3) :: fill_order
    type(real3),intent(out):: l1, l2, l3
    
    if(fill_order(1)=="x" .or. fill_order(1)=="X") then
      l1 = real3(1.0_RK,0.0_RK,0.0_RK)
    elseif(fill_order(1)=="y" .or. fill_order(1)=="Y") then
      l1 = real3(0.0_RK,1.0_RK,0.0_RK)
    elseif(fill_order(1)=="z" .or. fill_order(1)=="Z") then
      l1 = real3(0.0_RK,0.0_RK,1.0_RK)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 1 ")
    endif
            
    if(fill_order(2)=="x".or.fill_order(2)=="X") then
      l2 = real3(1.0_RK,0.0_RK,0.0_RK)
    elseif(fill_order(2)=="y".or.fill_order(2)=="Y") then
      l2 = real3(0.0_RK,1.0_RK,0.0_RK)
    elseif(fill_order(2)=="z".or.fill_order(2)=="Z") then
      l2 = real3(0.0_RK,0.0_RK,1.0_RK)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 2 ")
    endif            
            
    if(fill_order(3)=="x".or.fill_order(3)=="X") then
      l3 = real3(1.0_RK,0.0_RK,0.0_RK)
    elseif(fill_order(3)=="y".or.fill_order(3)=="Y") then
      l3 = real3(0.0_RK,1.0_RK,0.0_RK) 
    elseif(fill_order(3)=="z".or.fill_order(3)=="Z") then
      l3 = real3(0.0_RK,0.0_RK,1.0_RK)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 3 ")
    endif
  end subroutine Fill_Vectors
#endif

  !**********************************************************************
  ! copy i2 to i1
  !**********************************************************************
  subroutine GL_copy(this,i1,i2)
    implicit none
    class(VarList)::this
    integer,intent(in)::i1,i2

    if(i1==i2) return
    GPrtcl_id(i1)       = GPrtcl_id(i2)
    GPrtcl_pType(i1)    = GPrtcl_pType(i2)
    GPrtcl_usrMark(i1)  = GPrtcl_usrMark(i2)
    GPrtcl_PosR(i1)     = GPrtcl_PosR(i2)
    GPrtcl_linVel(:,i1) = GPrtcl_linVel(:,i2)
    GPrtcl_linAcc(:,i1) = GPrtcl_linAcc(:,i2)
    GPrtcl_theta(i1)    = GPrtcl_theta(i2)
    GPrtcl_rotVel(:,i1) = GPrtcl_rotVel(:,i2)
    GPrtcl_rotAcc(:,i1) = GPrtcl_rotAcc(:,i2)
#ifdef CFDDEM
    GPrtcl_FpForce(i1)     = GPrtcl_FpForce(i2)
    GPrtcl_FpForce_old(i1) = GPrtcl_FpForce_old(i2)
    GPrtcl_linVelOld(i1)   = GPrtcl_linVelOld(i2)
    GPrtcl_Vfluid(:,i1)    = GPrtcl_Vfluid(:,i2)
    if(is_clc_Basset)  GPrtcl_BassetData(:,i1)= GPrtcl_BassetData(:,i2)
#endif
#ifdef CFDACM
    GPrtcl_FpForce(i1) = GPrtcl_FpForce(i2)
    GPrtcl_FpTorque(i1)= GPrtcl_FpTorque(i2)
    GPrtcl_FluidIntOld(:,i1)= GPrtcl_FluidIntOld(:,i2)
    if(IBM_Scheme==2) GPrtcl_PosOld(i1) = GPrtcl_PosOld(i2)
#endif
  end subroutine GL_copy

  !**********************************************************************
  ! Reallocate particle varaibles 
  !**********************************************************************
  subroutine GL_ReallocatePrtclVar(this,np_new)
    implicit none
    class(VarList)::this
    integer,intent(in)::np_new

    ! locals
    integer:: sizep,sizen,ierrTmp,ierror=0
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec
    type(real3),dimension(:,:),allocatable::Real3Arr
#ifdef CFDACM
    character(len=1),dimension(:),allocatable::ChaVec
#endif

    sizep= this%mlocal
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, np_new+1)
    sizen= min(sizen,DEM_Opt%numPrtcl)
    this%mlocal=sizen

    ! ======= integer vector part =======
    call move_alloc(GPrtcl_id, IntVec)
    allocate(GPrtcl_id(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_id(1:sizep)=IntVec
    GPrtcl_id(sizep+1:sizen)=0

    call move_alloc(GPrtcl_pType, IntVec)
    allocate(GPrtcl_pType(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_pType(1:sizep)=IntVec
    GPrtcl_pType(sizep+1:sizen)=1

    call move_alloc(GPrtcl_usrMark, IntVec)
    allocate(GPrtcl_usrMark(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_usrMark(1:sizep)=IntVec
    GPrtcl_usrMark(sizep+1:sizen)=1
    deallocate(IntVec)

    ! ======= real3 vercor part =======
    call move_alloc(GPrtcl_theta,Real3Vec)
    allocate(GPrtcl_theta(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_theta(1:sizep)=Real3Vec
    GPrtcl_theta(sizep+1:sizen)=zero_r3
#ifdef CFDDEM
    call move_alloc(GPrtcl_FpForce,Real3Vec)
    allocate(GPrtcl_FpForce(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpForce(1:sizep)=Real3Vec
    GPrtcl_FpForce(sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_FpForce_old,Real3Vec)
    allocate(GPrtcl_FpForce_old(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpForce_old(1:sizep)=Real3Vec
    GPrtcl_FpForce_old(sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_linVelOld,Real3Vec)
    allocate(GPrtcl_linVelOld(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_linVelOld(1:sizep)=Real3Vec
    GPrtcl_linVelOld(sizep+1:sizen)=zero_r3
#endif
#ifdef CFDACM
    call move_alloc(GPrtcl_FpForce,Real3Vec)
    allocate(GPrtcl_FpForce(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpForce(1:sizep)=Real3Vec
    GPrtcl_FpForce(sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_FpTorque,Real3Vec)
    allocate(GPrtcl_FpTorque(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpTorque(1:sizep)=Real3Vec
    GPrtcl_FpTorque(sizep+1:sizen)=zero_r3
#endif
    deallocate(Real3Vec)
 
    deallocate(GPrtcl_cntctForce)
    allocate(GPrtcl_cntctForce(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)

    deallocate(GPrtcl_torque)
    allocate(GPrtcl_torque(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    ! ======= real3 matrix part =======
    call move_alloc(GPrtcl_linVel,Real3Arr)
    allocate(GPrtcl_linVel(this%tsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_linVel(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linVel(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_linAcc,Real3Arr)
    allocate(GPrtcl_linAcc(this%tsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_linAcc(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linAcc(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_rotVel,Real3Arr)
    allocate(GPrtcl_rotVel(this%rsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_rotVel(1:this%rsize,1:sizep)=Real3Arr
    GPrtcl_rotVel(1:this%rsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_rotAcc,Real3Arr)
    allocate(GPrtcl_rotAcc(this%rsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_rotAcc(1:this%rsize,1:sizep)=Real3Arr
    GPrtcl_rotAcc(1:this%rsize,sizep+1:sizen)=zero_r3
#ifdef CFDDEM
    call move_alloc(GPrtcl_Vfluid,Real3Arr)
    allocate(GPrtcl_Vfluid(2,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_Vfluid(1:2,1:sizep)=Real3Arr
    GPrtcl_Vfluid(1:2,sizep+1:sizen)=zero_r3

    if(is_clc_Basset) then
      call move_alloc(GPrtcl_BassetData,Real3Arr)
      allocate(GPrtcl_BassetData(GPrtcl_BassetSeq%nDataLen ,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      GPrtcl_BassetData(1:GPrtcl_BassetSeq%nDataLen, 1:sizep)=Real3Arr
      GPrtcl_BassetData(1:GPrtcl_BassetSeq%nDataLen, sizep+1:sizen)=zero_r3
    endif
#endif
#ifdef CFDACM
    call move_alloc(GPrtcl_FluidIntOld,Real3Arr)
    allocate(GPrtcl_FluidIntOld(2,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FluidIntOld(1:2,1:sizep)=Real3Arr    
    GPrtcl_FluidIntOld(1:2,sizep+1:sizen)=zero_r3    

    deallocate(GPrtcl_FluidIntegrate)
    allocate(GPrtcl_FluidIntegrate(2,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
#endif   
    deallocate(Real3Arr) 

    ! ======= real4 vercor part =======
    call move_alloc(GPrtcl_PosR,Real4Vec)
    allocate(GPrtcl_PosR(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_PosR(1:sizep)=Real4Vec
    GPrtcl_PosR(sizep+1:sizen)=zero_r4
#ifdef CFDACM
    if(IBM_Scheme==2) then
      call move_alloc(GPrtcl_PosOld,Real4Vec)
      allocate(GPrtcl_PosOld(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      GPrtcl_PosOld(1:sizep)=Real4Vec
      GPrtcl_PosOld(sizep+1:sizen)=zero_r4
    endif
#endif
    deallocate(Real4Vec)

#ifdef CFDACM
    ! ======= character vercor part =======
    call move_alloc(GPrtcl_HighSt,ChaVec)
    allocate(GPrtcl_HighSt(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_HighSt(1:sizep)=ChaVec
    GPrtcl_HighSt(sizep+1:sizen)="N"
    deallocate(ChaVec)
#endif

    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," GL_ReallocatePrtclVar"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif
    !call DEMLogInfo%CheckForError(ErrT_Pass," ReallocatePrtclVar"," Need to reallocate particle variables")
    !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine GL_ReallocatePrtclVar

end module Prtcl_Variables
