program main_channelLPT
  use MPI
  use m_Timer
  use m_Tools
  use LPT_System
  use m_Decomp2d
  use LPT_Fpforce
  use m_Variables
  use m_IOAndVisu  
  use m_Parameters
  use LPT_decomp_2d
  use LPT_IOAndVisu
  use LPT_Variables
  use LPT_Parameters
  use m_ChannelSystem
  use m_MeshAndMetries
  use m_TypeDef,only:num2str
  use m_Poisson,only:Destory_Poisson_FFT_Plan
  implicit none
  integer::intT,pid,ierror
  type(timer)::CoupleTimer
  character(len=128)::chPrm
  character(len=10)::RowColStr
#ifdef CFDFourthOrder
  integer::BcOption(6)
#endif

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! ================== initialize channel options ==================  
  intT=command_argument_count()
  if((intT/=2 .and. intT/=4) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif
  call get_command_argument(1,chPrm)
  call ReadAndInitParameters(chPrm)
  if(intT==4) then
    call get_command_argument(3,RowColStr)
    read(RowColStr,*) p_row
    call get_command_argument(4,RowColStr)
    read(RowColStr,*) p_col
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  
#ifdef CFDFourthOrder
  if(FlowType==FT_CH) then
    BcOption=(/0,0,-1,-1,0,0/)
  elseif(FlowType==FT_HC) then
    BcOption=(/0,0,-1,-2,0,0/)
  endif
#endif
  call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
  call ChannelInitialize(chPrm)    ! Topest level initialing for Channel body

  ! ================== initialize LPT options ==================
  call get_command_argument(2,chPrm)
  call LPT_opt%ReadLPTOption( chPrm)
  call LPT_decomp%Init_DECOMP()
  call LPT%Initialize(chPrm)       ! Topest level initialing for LPT body

  ! ================== initialize CFD-LPT coupling part ==================
  call InitFpForce(chPrm)

  ! =============== dump initial visulizing files ===============
  asso_Q: associate(Q_vor =>RealArr1)
  call  dump_visu(ifirst-1,ux,uy,uz,pressure,Q_vor)  ! channel3d
  end associate asso_Q
  if(LPT_Opt%RestartFlag)call LPT_IO%dump_visu(ifirst-1)
  print*, nrank,GPrtcl_list%nlocal

  call CoupleTimer%reset()
  do itime=ifirst,ilast

    ! CFD-LPT coupling part
    call CoupleTimer%start()
    call PrepareInterpolation()
    call clc_VelInterpolation(ux,uy,uz)
    if(itime==ifirst .and. (.not. LPT_Opt%RestartFlag)) then
      do pid=1,GPrtcl_list%nlocal
        GPrtcl_linVel(:,pid)=GPrtcl_VFluid(pid)
      enddo
      call LPT_IO%dump_visu(ifirst-1)                    ! LPT
    endif
    call clc_FpForce() !clc_FpForce(ux,uy,uz,pressure)
#ifdef CFDLPT_TwoWay
    call distribute_FpForce()
#endif

    call FinalFpForce()
    call CoupleTimer%finish()

    ! Largrangian Particle Trackiing part
    call LPT%iterate(itime)
    
    ! CFD Iterate
    call ChannelIterate()
    
    if(nrank==0 .and. mod(itime, Cmd_LFile_Freq)==0) then
      call MainLog%OutInfo("Coupling time  [tot, last, ave] [sec]: "//trim(num2str(CoupleTimer%tot_time))//", "// &
          trim(num2str(CoupleTimer%last_time ))//", "//trim(num2str(CoupleTimer%average())),2)
    endif  
  enddo
  call LPT_IO%Final_visu()

  if(nrank==0)call MainLog%OutInfo("Good job! ChannelLPM finished successfully at "//time2str(),1)
  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
end program main_channelLPT

module LPT_Variables
  use MPI
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_decomp_2d
  use LPT_Parameters
  use m_Decomp2d,only: nrank
  implicit none
  private

  integer,dimension(:),allocatable,public:: GPrtcl_id
  integer,dimension(:),allocatable,public:: GPrtcl_pType
  integer,dimension(:),allocatable,public:: GPrtcl_usrMark
  type(real4),dimension(:),allocatable,public::   GPrtcl_PosR
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linAcc
  type(real3),dimension(:),allocatable,public:: GPrtcl_PosOld
  type(real3),dimension(:),allocatable,public:: GPrtcl_VFluid
  
  type(real3),dimension(:),  allocatable,public:: GPrtcl_FpForce

  type VarList
    integer:: nlocal     ! number of particles in local processor
    integer:: mlocal     ! the possible maxium # of particles in local processor
    integer:: tsize      ! size for translational veloctity and acceleration
  contains
    procedure:: AllocateAllVar     => GL_AllocateAllVar
    procedure:: MakingAllPrtcl     => GL_MakingAllPrtcl
    procedure:: ReallocatePrtclVar => GL_ReallocatePrtclVar
    procedure:: copy               => GL_copy
  end type VarList
  type(VarList),public::GPrtcl_list
    
contains

  !**********************************************************************
  ! GL_AllocateAllVar
  !**********************************************************************
  subroutine GL_AllocateAllVar(this)
    implicit none
    class(VarList)::this

    ! locals
    type(real3)::SimLen
    integer::mlocal,numPrtcl,ierrTmp,ierror=0
    real(RK)::xst,xed,yst,yed,zst,zed,vol_tot,vol_local

    numPrtcl    = LPT_opt%numPrtcl

    ! step1: allocating memory for particles
    if (LPT_opt%PI_Method==PIM_AB2) then
      this%tsize=2
    elseif(LPT_opt%PI_Method==PIM_AB3) then
      this%tsize=3
    endif

    ! step0: determine initial msize
    xst=LPT_decomp%xSt; xed=LPT_decomp%xEd
    yst=LPT_decomp%ySt; yed=LPT_decomp%yEd
    zst=LPT_decomp%zSt; zed=LPT_decomp%zEd

    SimLen = LPT_Opt%SimDomain_max - LPT_Opt%SimDomain_min
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
    allocate(GPrtcl_PosOld(mlocal),            Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_VFluid(mlocal),            Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar","Allocation failed 1")

    allocate(GPrtcl_FpForce(mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 2")    
    
    GPrtcl_id = 0
    GPrtcl_pType = 1
    GPrtcl_usrMark = 1
    GPrtcl_PosR = zero_r4
    GPrtcl_linVel  = zero_r3
    GPrtcl_linAcc  = zero_r3
    GPrtcl_PosOld  = zero_r3
    GPrtcl_VFluid  = zero_r3
    GPrtcl_FpForce = zero_r3
    
    ! step2: Initialize this%mlocal,this%nlocal
    this%nlocal  = 0
    this%mlocal  = mlocal
  end subroutine GL_AllocateAllVar

  !**********************************************************************
  ! GL_MakingAllPrtcl
  !**********************************************************************
  subroutine GL_MakingAllPrtcl(this, chFile)
    implicit none
    class(VarList)::this
    character(*),intent(in)::chFile
        
    ! locals     
    logical::IsRandomDist
    character,dimension(3):: Fill_order
    integer,dimension(:),allocatable:: sum_bin    
    integer::i,j,k,iTV(8),bin_id,numPrtcl,nUnitFile,ierror,nlocal,nlocal_sum
    real(RK)::xst,xed,yst,yed,zst,zed,maxRad,randt,Distance_Ratio,MkPrtclMinpoint(3), MkPrtclMaxpoint(3),realt(3)
    type(real3):: MkPrtclMinpoint_real3, MkPrtclMaxpoint_real3,cntr,l1_vec,l2_vec,l3_vec,lmin_p,lmax_p,dx,pos1,pos2
    NAMELIST/ParticleMakingOption/ MkPrtclMinpoint,MkPrtclMaxpoint,Fill_order,Distance_Ratio,IsRandomDist
        
    open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl", "Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=ParticleMakingOption)
    if(nrank==0)write(LPTLogInfo%nUnit, nml=ParticleMakingOption)
    close(nUnitFile,IOSTAT=ierror)
    MkPrtclMinpoint_real3 = MkPrtclMinpoint
    MkPrtclMaxpoint_real3 = MkPrtclMaxpoint
    numPrtcl = LPT_opt%numPrtcl
        
    xst=LPT_decomp%xSt; xed=LPT_decomp%xEd
    yst=LPT_decomp%ySt; yed=LPT_decomp%yEd
    zst=LPT_decomp%zSt; zed=LPT_decomp%zEd

    ! step3: Assign Prtcl_Type, Prtcl_id and Prtcl_PosR
    call Fill_Vectors(Fill_order, l1_vec, l2_vec, l3_vec)
    call date_and_time(values=iTV); !iTV=0
    if(IsRandomDist) iTV=0
    call random_seed(size= i)
    call random_seed(put = iTV(7)*iTV(8)+[(j,j=1,i)])
    if(IsRandomDist) Distance_Ratio=1.0_RK
    maxRad = maxval( LPTProperty%Prtcl_PureProp%Radius )
    lmin_p = MkPrtclMinpoint_real3 + Distance_Ratio* maxRad*real3(1.0_RK,1.0_RK,1.0_RK)
    lmax_p = MkPrtclMaxpoint_real3 - Distance_Ratio* maxRad*real3(1.0_RK,1.0_RK,1.0_RK)
    dx = Distance_Ratio*2.0_RK*maxRad*real3(1.0_RK,1.0_RK,1.0_RK)
        
    allocate(sum_bin(LPT_opt%numPrtcl_Type))
    sum_bin(1)=LPTProperty%nPrtcl_in_Bin(1)
    do j=2, LPT_opt%numPrtcl_Type
      sum_bin(j)=sum_bin(j-1)+LPTProperty%nPrtcl_in_Bin(j)
    enddo

    nlocal=0; bin_id=-1
    IF(IsRandomDist) THEN
      do i=1,numPrtcl
        do j=1,LPT_opt%numPrtcl_Type
          if(i<=sum_bin(j)) then
            bin_id = j; exit
          endif
        enddo
                
        call random_number(realt)
        cntr%x=  lmin_p%x +realt(1)*(lmax_p%x-lmin_p%x)
        cntr%y=  lmin_p%y +realt(2)*(lmax_p%y-lmin_p%y)
        cntr%z=  lmin_p%z +realt(3)*(lmax_p%z-lmin_p%z)         
        if(cntr%x>=xst.and.cntr%x<xed.and.cntr%y>=yst.and.cntr%y<yed.and.cntr%z>=zst.and.cntr%z<zed) then
          if(nlocal>=this%mlocal) call this%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1
          GPrtcl_id(nlocal) = i
          GPrtcl_pType(nlocal)= bin_id
          GPrtcl_PosR(nlocal) = cntr
          GPrtcl_PosR(nlocal)%w = LPTProperty%Prtcl_PureProp(bin_id)%Radius
        endif
      enddo    
    ELSE
      cntr = lmin_p       ! start point
      do i=1,numPrtcl
        do j=1,LPT_opt%numPrtcl_Type
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
          GPrtcl_PosR(nlocal)%w = LPTProperty%Prtcl_PureProp(bin_id)%Radius
        endif

        cntr = cntr + l1_vec * dx
        if((l1_vec.dot.cntr)>=(l1_vec.dot.lmax_p))then
          cntr = (lmin_p*l1_vec) + ((cntr+dx)*l2_vec)+(cntr*l3_vec)
          if((l2_vec.dot.cntr)>=(l2_vec.dot.lmax_p))then
            cntr = (cntr*l1_vec)+(lmin_p*l2_vec)+((cntr+dx)*l3_vec)
            if((l3_vec.dot.cntr) >= (l3_vec.dot.lmax_p) .and. nrank==0) then
              call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Not enough space for positioning" ) 
            endif
          endif
        endif
      enddo
    ENDIF
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: "," nlocal_sum/= numPrtcl " )
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
  end subroutine GL_MakingAllPrtcl
  
  !**********************************************************************
  ! Fill order Vector 
  !**********************************************************************
  subroutine Fill_Vectors(fill_order,l1,l2,l3)
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
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 1 ")
    endif
            
    if(fill_order(2)=="x".or.fill_order(2)=="X") then
      l2 = real3(1.0_RK,0.0_RK,0.0_RK)
    elseif(fill_order(2)=="y".or.fill_order(2)=="Y") then
      l2 = real3(0.0_RK,1.0_RK,0.0_RK)
    elseif(fill_order(2)=="z".or.fill_order(2)=="Z") then
      l2 = real3(0.0_RK,0.0_RK,1.0_RK)
    else
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 2 ")
    endif            
            
    if(fill_order(3)=="x".or.fill_order(3)=="X") then
      l3 = real3(1.0_RK,0.0_RK,0.0_RK)
    elseif(fill_order(3)=="y".or.fill_order(3)=="Y") then
      l3 = real3(0.0_RK,1.0_RK,0.0_RK) 
    elseif(fill_order(3)=="z".or.fill_order(3)=="Z") then
      l3 = real3(0.0_RK,0.0_RK,1.0_RK)
    else
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 3 ")
    endif
  end subroutine Fill_Vectors

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
    GPrtcl_FpForce(i1)  = GPrtcl_FpForce(i2)
  end subroutine GL_copy

  !**********************************************************************
  ! Reallocate particle varaibles 
  !**********************************************************************
  subroutine GL_ReallocatePrtclVar(this,np_new)
    implicit none
    class(VarList)::this
    integer,intent(in)::np_new

    ! locals
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec
    type(real3),dimension(:,:),allocatable::Real3Arr

    sizep= this%mlocal
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, np_new+1)
    sizen= min(sizen,LPT_Opt%numPrtcl)
    this%mlocal=sizen

    ! ======= integer vector part =======
    call move_alloc(GPrtcl_id, IntVec)
    allocate(GPrtcl_id(sizen))
    GPrtcl_id(1:sizep)=IntVec
    GPrtcl_id(sizep+1:sizen)=0

    call move_alloc(GPrtcl_pType, IntVec)
    allocate(GPrtcl_pType(sizen))
    GPrtcl_pType(1:sizep)=IntVec
    GPrtcl_pType(sizep+1:sizen)=1

    call move_alloc(GPrtcl_usrMark, IntVec)
    allocate(GPrtcl_usrMark(sizen))
    GPrtcl_usrMark(1:sizep)=IntVec
    GPrtcl_usrMark(sizep+1:sizen)=1
    deallocate(IntVec)

    ! ======= real3 vector part =======
    deallocate(GPrtcl_PosOld)
    allocate(GPrtcl_PosOld(sizen))
    deallocate(GPrtcl_VFluid)
    allocate(GPrtcl_VFluid(sizen))
    call move_alloc(GPrtcl_FpForce,Real3Vec)
    allocate(GPrtcl_FpForce(sizen))
    GPrtcl_FpForce(1:sizep)=Real3Vec
    GPrtcl_FpForce(sizep+1:sizen)=zero_r3
    
    ! ======= real3 matrix part =======
    call move_alloc(GPrtcl_linVel,Real3Arr)
    allocate(GPrtcl_linVel(this%tsize,sizen))
    GPrtcl_linVel(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linVel(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_linAcc,Real3Arr)
    allocate(GPrtcl_linAcc(this%tsize,sizen))
    GPrtcl_linAcc(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linAcc(1:this%tsize,sizep+1:sizen)=zero_r3
    deallocate(Real3Arr) 

    ! ======= real4 vercor part =======
    call move_alloc(GPrtcl_PosR,Real4Vec)
    allocate(GPrtcl_PosR(sizen))
    GPrtcl_PosR(1:sizep)=Real4Vec
    GPrtcl_PosR(sizep+1:sizen)=zero_r4
    deallocate(Real4Vec)

    call LPTLogInfo%CheckForError(ErrT_Pass," ReallocatePrtclVar"," Need to reallocate particle variables")
    call LPTLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call LPTLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call LPTLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine GL_ReallocatePrtclVar

end module LPT_Variables

module LPT_Comm
  use MPI
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_Decomp_2d
  use LPT_Variables
  use LPT_Parameters
  use LPT_ContactSearchPW
  use m_Decomp2d,only: nrank
  implicit none
  private

  logical,dimension(3)::pbc
  type(real3)::simLen
  integer,parameter::xm_axis=1
  integer,parameter::xp_axis=2
  integer,parameter::ym_axis=3
  integer,parameter::yp_axis=4
  integer,parameter::zm_axis=5
  integer,parameter::zp_axis=6
  real(RK),dimension(6)::dx_pbc
  real(RK),dimension(6)::dy_pbc
  real(RK),dimension(6)::dz_pbc  

  real(RK):: xst1_cs
  real(RK):: xed1_cs
  real(RK):: yst1_cs
  real(RK):: yed1_cs
  real(RK):: zst1_cs
  real(RK):: zed1_cs
  integer,allocatable,dimension(:)::sendlist

  type Prtcl_Comm_info
    integer :: msend
    integer :: Prtcl_Exchange_size
  contains
    procedure:: InitComm            => PC_InitComm
    procedure:: Comm_For_Exchange   => PC_Comm_For_Exchange
    procedure:: pack_Exchange       => PC_pack_Exchange
    procedure:: unpack_Exchange     => PC_unpack_Exchange
    procedure:: ISInThisProc        => PC_ISInThisProc
    procedure:: reallocate_sendlist => PC_reallocate_sendlist
  end type Prtcl_Comm_info
  type(Prtcl_Comm_info),public::LPTComm

contains

  !**********************************************************************
  ! PC_InitComm
  !**********************************************************************
  subroutine PC_InitComm(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    integer::iErr01
        
    pbc=LPT_Opt%IsPeriodic
    simLen = LPT_Opt%SimDomain_max-LPT_Opt%SimDomain_min

    xst1_cs = LPT_decomp%xSt   
    xed1_cs = LPT_decomp%xEd
    yst1_cs = LPT_decomp%ySt
    yed1_cs = LPT_decomp%yEd 
    zst1_cs = LPT_decomp%zSt
    zed1_cs = LPT_decomp%zEd

    dx_pbc=0.0_RK; dy_pbc=0.0_RK; dz_pbc=0.0_RK
    if(pbc(1)) then
      if(LPT_decomp%coord1==0)                 dx_pbc(xm_axis)= simLen%x
      if(LPT_decomp%coord1==LPT_decomp%prow-1) dx_pbc(xp_axis)=-simLen%x 
    endif
    if(pbc(2)) then
      dy_pbc(ym_axis)= simLen%y
      dy_pbc(yp_axis)=-simLen%y
    endif
    if(pbc(3)) then
      if(LPT_decomp%coord2==0)                 dz_pbc(zm_axis)= simLen%z
      if(LPT_decomp%coord2==LPT_decomp%pcol-1) dz_pbc(zp_axis)=-simLen%z
    endif

    ! (id 1) +(ptype 1) +(Mark 1) +(PosR   3) +(linvel 3*tsize) +(linAcc 3*tsize) +(GPrtcl_FpForce 3) 
    ! = 9+6*tsize
    this%Prtcl_Exchange_size = 9 + 6*GPrtcl_list%tsize
    this%msend=GPrtcl_list%mlocal
    allocate(sendlist(this%msend),   Stat=iErr01)
    if(iErr01/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","Allocation failed2")

  end subroutine PC_InitComm

  !**********************************************************************
  ! PC_Comm_For_Exchange
  !**********************************************************************
  subroutine PC_Comm_For_Exchange(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    integer::i,ierror,request(4),nlocal,nlocalp
    integer::nsend(2),nrecv(2)
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer,dimension(MPI_STATUS_SIZE) :: SRstatus
    real(RK)::px,py,pz

    nlocal=GPrtcl_list%nlocal

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
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo
    ENDIF

    ! step2: send to xp_axis, and receive from xm_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        px=GPrtcl_PosR(i)%x
        if(px >= xed1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2) + this%Prtcl_Exchange_size
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(3), 9, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,xp_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif
      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xp_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! step3: send to xm_axis, and receive from xp_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        px=GPrtcl_PosR(i)%x
        if(px < xst1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2)+ this%Prtcl_Exchange_size
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(4),11, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,xm_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif
      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xm_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! step4: send to zp_axis, and receive from zm_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        pz=GPrtcl_PosR(i)%z
        if(pz >= zed1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2)+ this%Prtcl_Exchange_size
        endif
      enddo

    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(1),13, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,zp_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif
      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zp_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! step5: send to zm_axis, and receive from zp_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        pz=GPrtcl_PosR(i)%z
        if(pz < zst1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2)+ this%Prtcl_Exchange_size
        endif
      enddo

    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(2),15, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,zm_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif

      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zm_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)
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
    integer::i,j,id,m,nlocal
    real(RK)::dx,dy,dz

    m=1
    nlocal = nlocalp
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)
   
    DO i=1,nsend
      id = sendlist(i)
      buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
      buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
      buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
      do j=1,GPrtcl_list%tsize
        buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
        buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
      enddo
      buf_send(m)=GPrtcl_FpForce(id)%x;     m=m+1
      buf_send(m)=GPrtcl_FpForce(id)%y;     m=m+1
      buf_send(m)=GPrtcl_FpForce(id)%z;     m=m+1
    ENDDO
    DO i=nsend,1,-1
      id = sendlist(i)
      call GPrtcl_list%copy(id,nlocal)
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
    DO i=1,nrecv
      if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
      GPrtcl_id(id)     = nint(buf_recv(m)); m=m+1 ! 01
      itype             = nint(buf_recv(m)); m=m+1 ! 02
      GPrtcl_usrMark(id)= nint(buf_recv(m)); m=m+1 ! 03
      GPrtcl_PosR(id)%x = buf_recv(m);       m=m+1 ! 04
      GPrtcl_PosR(id)%y = buf_recv(m);       m=m+1 ! 05
      GPrtcl_PosR(id)%z = buf_recv(m);       m=m+1 ! 06
      do j=1,GPrtcl_list%tsize
        GPrtcl_linVel(j,id)%x = buf_recv(m); m=m+1 ! 6* tsize
        GPrtcl_linVel(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_linVel(j,id)%z = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%x = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%z = buf_recv(m); m=m+1
      enddo
      GPrtcl_FpForce(id)%x    = buf_recv(m); m=m+1
      GPrtcl_FpForce(id)%y    = buf_recv(m); m=m+1
      GPrtcl_FpForce(id)%z    = buf_recv(m); m=m+1
      GPrtcl_pType(id)=itype
      GPrtcl_PosR(id)%w =LPTProperty%Prtcl_PureProp(itype)%Radius
 
      id = id + 1
      nlocal =nlocal + 1
    ENDDO
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
    
    call LPTLogInfo%CheckForError(ErrT_Pass," PC_ISInThisProc"," The following particle is deleted: ")
    call LPTLogInfo%OutInfo(" Exchange direction  is :"//trim(num2str(dir)  ),3)
    call LPTLogInfo%OutInfo("   The particle id   is :"//trim(num2str( nint(buf_recv(m)  ))),3)
    call LPTLogInfo%OutInfo("       particle type is :"//trim(num2str( nint(buf_recv(m+1)))),3)
    call LPTLogInfo%OutInfo("       x-coordinate  is :"//trim(num2str( buf_recv(m+3) )), 3)
    call LPTLogInfo%OutInfo("       y-coordinate  is :"//trim(num2str( buf_recv(m+4) )), 3)
    call LPTLogInfo%OutInfo("       z-coordinate  is :"//trim(num2str( buf_recv(m+5) )), 3)
    call LPTLogInfo%OutInfo(" Present processor   is :"//trim(num2str(nrank)), 3)

    m = m + this%Prtcl_Exchange_size
    m=m+1
  end function PC_ISInThisProc

  !**********************************************************************
  ! PC_reallocate_sendlist
  !**********************************************************************  
  subroutine PC_reallocate_sendlist(this,ns)
    implicit none
    class(Prtcl_Comm_info)::this
    integer,intent(in)::ns

    ! locals
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec

    sizep= this%msend
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,LPT_Opt%numPrtcl)
    sizen=max(sizen,ns+1)
    this%msend=sizen
   
    call move_alloc(sendlist, IntVec)
    allocate(sendlist(sizen))
    sendlist(1:sizep)=IntVec
    deallocate(IntVec)

    call LPTLogInfo%CheckForError(ErrT_Pass," reallocate_sendlist"," Need to reallocate sendlist")
    call LPTLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call LPTLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call LPTLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine PC_reallocate_sendlist

end module LPT_Comm

module LPT_ContactSearchPW
  use MPI
  use m_TypeDef
  use LPT_Property
  use LPT_Geometry
  use LPT_Variables
  use LPT_Parameters
  use m_Decomp2d,only: nrank
  use m_Parameters,only:yly,xlx,zlz  
  use LPT_Decomp_2d,only:int_type,real_type
  implicit none
  private

  type::ContactSearchPW
  contains
    procedure:: FindContactsPW
  end type ContactSearchPW
  type(ContactSearchPW),public::LPTContactSearchPW
    
contains

  !**********************************************************************
  ! Performing contact search to determine particle-wall contacts 
  !**********************************************************************
  subroutine FindContactsPW(this)
    implicit none
    class(ContactSearchPW):: this

    ! locals
    integer:: pid
    real(RK)::Posy,radius

#ifdef UseDEMWallContact
    integer::wid
    DO pid=1,GPrtcl_list%nlocal
      do wid=1,LPTGeometry%nPW_local
         if(LPTGeometry%pWall(wid)%isInContact(GPrtcl_PosR(pid),ovrlp,nv))then
           GPrtcl_PosR(pid)%y=  GPrtcl_PosOld(pid)%y
           GPrtcl_linVel(1,pid)%y= - GPrtcl_linVel(2,pid)%y               
         endif 
      enddo
    ENDDO
#else
    DO pid=1,GPrtcl_list%nlocal
      Posy=GPrtcl_PosR(pid)%y
      radius= GPrtcl_PosR(pid)%w
      if(Posy<=radius .or. Posy+radius>=yly) then
        GPrtcl_PosR(pid)%y    =   GPrtcl_PosOld(pid)%y
        GPrtcl_linVel(1,pid)%y= - GPrtcl_linVel(2,pid)%y
      endif
    ENDDO
#endif
  end subroutine FindContactsPW
    
end module LPT_ContactSearchPW

module LPT_Geometry
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_decomp_2d
  use LPT_Parameters
  use m_Decomp2d,only:nrank
  use m_Parameters,only:xlx,yly,zlz
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
    real(RK):: d          ! d in the implicit equation: ax+by+cz+d = 0   
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
    procedure:: MakeGeometry =>G_MakeGeometry
  end type Geometry
    
  type(Geometry),public :: LPTGeometry
    
contains
    
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
      t = -((this%n .dot. p)+ this%d )
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
    if(a>=0.0_RK-1.00E-10_RK .and. a<=1.0_RK+1.00E-10_RK   .and. &
       b>=0.0_RK-1.00E-10_RK .and. b<=1.0_RK+1.00E-10_RK)  then
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
         
    LenExp=  1.2_RK*maxval( LPTProperty%Prtcl_PureProp%Radius )
    pmin = LPT_opt%SimDomain_min - LenExp *real3(1.0_RK,1.0_RK,1.0_RK)
    pmax = LPT_opt%SimDomain_max + LenExp *real3(1.0_RK,1.0_RK,1.0_RK)
    pmin_local =real3(LPT_decomp%xSt,LPT_decomp%ySt,LPT_decomp%zSt)-LenExp*real3(1.0_RK,1.0_RK,1.0_RK)
    pmax_local =real3(LPT_decomp%xEd,LPT_decomp%yEd,LPT_decomp%zEd)+LenExp*real3(1.0_RK,1.0_RK,1.0_RK)
        
    this%num_pWall = 0
    this%nPW_local = 0
    allocate(this%pWall( MaxWallSize ))
        
    if(nrank/=0) return
    write(chFile,"(A)") trim(LPT_opt%ResultsDir)//"WallsFor"//trim(LPT_opt%RunName)//".backup"
    open(newunit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"G_InitAllocate","Cannot open file: "//trim(chFile))
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
    logical ::lboth,linfinite
    type(PlaneWall)::wall
    type(PlaneWall),dimension(:),allocatable:: wall_temp
    type(real3)::t_vel, ln
    integer:: nUnitFile,ierror
    character(128) :: chFile        
 
    IF(this%nPW_local == MaxWallSize) THEN
      MaxWallSize  = int( real(MaxWallSize, RK) *1.2_RK) +1
      call move_alloc(this%pWall, wall_temp)
      allocate(this%pWall(MaxWallSize),Stat=ierror) 
      if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Reallocation failed, 2 ")
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
    wall%d = - (wall%n .dot. p1)
    wall%trans_vel = t_vel
    if(abs((wall%n .dot. p4) +wall%d) >= 0.00001_RK .and. nrank==0) then
      call LPTLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot create a plane wall, wall No. "//num2str(wall%wall_id))
    endif
        
    IF(wall%IsInDomain(pmin, pmax)) THEN
      this%num_pWall= this%num_pWall +1
      wall%wall_id  = this%num_pWall
      if(wall%IsInDomain(pmin_local,pmax_local)) then
        this%nPW_local = this%nPW_local +1
        this%pWall(this%nPW_local) = wall 
      endif
        
      if(nrank/=0) return
      write(chFile,"(A)") trim(LPT_opt%ResultsDir)//"WallsFor"//trim(LPT_opt%RunName)//".backup"
      open(newunit=nUnitFile, file=chFile, status='old',position='append',form='formatted',IOSTAT=ierror )
      if(ierror/=0 .and. nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot open file: "//trim(chFile))
      write(nUnitFile,* ) p1
      write(nUnitFile,* ) p2
      write(nUnitFile,* ) p3
      write(nUnitFile,* ) p4
      close(nUnitFile,IOSTAT=ierror)
    ELSE
      if(nrank/=0) return
      call LPTLogInfo%CheckForError(ErrT_Pass,"G_add_PlaneWall","  The following plane ISNOT within the simulation domain: ")
      call LPTLogInfo%OutInfo("   It will be skipped :",3, .true.)
      call LPTLogInfo%OutInfo("   Point 1: "//trim(num2str(p1%x))//'  '//trim(num2str(p1%y))//'  '//trim(num2str(p1%z)), 3, .true.)
      call LPTLogInfo%OutInfo("   Point 2: "//trim(num2str(p2%x))//'  '//trim(num2str(p2%y))//'  '//trim(num2str(p2%z)), 3, .true.)
      call LPTLogInfo%OutInfo("   Point 3: "//trim(num2str(p3%x))//'  '//trim(num2str(p3%y))//'  '//trim(num2str(p3%z)), 3, .true.)
      call LPTLogInfo%OutInfo("   Point 4: "//trim(num2str(p4%x))//'  '//trim(num2str(p4%y))//'  '//trim(num2str(p4%z)), 3, .true.)
    ENDIF
  end subroutine G_add_PlaneWall

  !**********************************************************************
  ! MakeGeometry
  !**********************************************************************     
  subroutine G_MakeGeometry(this)
    implicit none
    class(Geometry)::this
        
    !locals
    type(real3):: p01,p02,p03,p04
        
    call this%InitAllocate()
    p01= real3( 0.0_RK,  0.0_RK,  0.0_RK)
    p02= real3( 0.0_RK,  0.0_RK,   zlz)
    p03= real3(  xlx,  0.0_RK,   zlz)
    p04= real3(  xlx,  0.0_RK,  0.0_RK)
    call this%add_PlaneWall( p01, p02, p03, p04, 1, 1, infinite=.true. )
    p01= real3( 0.0_RK,  yly,  0.0_RK)
    p02= real3(  xlx,  yly,  0.0_RK)
    p03= real3(  xlx,  yly,   zlz)
    p04= real3( 0.0_RK,  yly,   zlz)
    call this%add_PlaneWall( p01, p02, p03, p04, 1, 1, infinite=.true. )
  end subroutine G_MakeGeometry
 
end module LPT_Geometry

module LPT_decomp_2d
  use MPI
  use m_TypeDef
  use LPT_Parameters
  use m_MeshAndMetries,only:dx,dz
  use m_Parameters,only:p_row,p_col
  use m_Decomp2d,only:nrank,y1start,y1end,myProcNghBC
  implicit none
  private

  integer,public:: int_type,real_type,real3_type
  integer,public:: int_byte,real_byte,real3_byte
  TYPE Prtcl_DECOMP_INFO

    ! define neighboring blocks
    ! second dimension 8 neighbour processors:
    !        1:4, 4 edge neighbours; 5:8, 4 cornor neighbours; 0, current processor(if any)
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
  TYPE(Prtcl_DECOMP_INFO), public :: LPT_decomp

contains

  !**********************************************************************
  ! PDI_Init_DECOMP
  !**********************************************************************
  subroutine PDI_Init_DECOMP(this)
    implicit none
    class(Prtcl_DECOMP_INFO)::this

    ! locals
    integer::i

    this%prow= p_row
    this%pcol= p_col
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

    this%coord1 = int ( nrank / p_col)
    this%coord2 = mod ( nrank,  p_col)
    DO i=1,4
      if(myProcNghBC(2,i)<0) then
        this%ProcNgh(i)= MPI_PROC_NULL
      else
        this%ProcNgh(i)= myProcNghBC(2,i)
      endif
    ENDDO
    this%xSt= real(y1start(1)-1,kind=RK)*dx
    this%xEd= real(y1end(1),    kind=RK)*dx
    this%zSt= real(y1start(3)-1,kind=RK)*dz
    this%zEd= real(y1end(3),    kind=RK)*dz
    this%ySt= LPT_Opt%SimDomain_min%y
    this%yEd= LPT_Opt%SimDomain_max%y
  end subroutine PDI_Init_DECOMP

  !**********************************************************************
  ! Init_Prctl_MPI_TYPE
  !**********************************************************************
  subroutine Init_Prctl_MPI_TYPE()
    implicit none
    integer::ierror
    integer,dimension(3)::disp,blocklen,blocktype
  
    ! integer
    int_type = MPI_INTEGER
    call MPI_TYPE_SIZE(int_type,int_byte,ierror)

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
  end subroutine Init_Prctl_MPI_TYPE
  
end module LPT_decomp_2d

module LPT_Fpforce
  use MPI
  use m_LogInfo
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use LPT_Property
  use LPT_Variables
  use LPT_parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
#ifdef CFDLPT_TwoWay
  use m_Variables,only: FpForce_x,FpForce_y,FpForce_z
#endif
  implicit none
  private

#ifdef CFDFourthOrder
  real(RK),public,allocatable,dimension(:):: xc   ! center coordinate in x-dir
  real(RK),public,allocatable,dimension(:):: zc   ! center coordinate in z-dir
#endif
  real(RK),allocatable,dimension(:)::iDistRatioYp,iDistRatioYc

  integer(kind=2),allocatable,dimension(:,:)::indxyz
  type(real3),allocatable,dimension(:)::RatioYp_interp,RatioYc_interp

  type(HaloInfo):: hi_ux_interp,   hi_uz_interp       ! halo info type for interpolation(velocity)
  integer:: xmp_interp, xep_interp, zmp_interp, zep_interp ! index constraints for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
  integer:: xmc_interp, xec_interp, zmc_interp, zec_interp ! index constraints for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
  
  procedure(),pointer::PrepareInterpolation,clc_VelInterpolation,distribute_FpForce
  public:: InitFpForce,PrepareInterpolation,clc_VelInterpolation
  public:: clc_FpForce,distribute_FpForce,FinalFpForce
contains
 
  !******************************************************************
  ! InitFpForce
  !******************************************************************
  subroutine InitFpForce(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    integer::j,nUnitFile,ierror,InterpAccuracy
    namelist/CFDLPT_interpolation/InterpAccuracy
 
    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(0_2)-20 .or. nyc>huge(0_2)-20 .or. nzc>huge(0_2)-20)) then
      call MainLog%CheckForError(ErrT_Abort,"InitFpForce","kind=2 is not enough for indxyz")
    endif
    
#ifdef CFDFourthOrder
    ! xc,zc, center coordinate interval in x-dir and z-dir
    allocate(xc(0:nxp))
    allocate(zc(0:nzp))
    do j=0,nxp
      xc(j)=dx*(real(j,RK)-0.5_RK)
    enddo
    do j=0,nzp
      zc(j)=dz*(real(j,RK)-0.5_RK)
    enddo    
#endif
  
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitFpForce","Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=CFDLPT_interpolation)
    close(nUnitFile,IOSTAT=ierror)
    
    if(InterpAccuracy==1) then      ! tri-linear interpolation
      PrepareInterpolation => PrepareInterpolation_1
      clc_VelInterpolation => clc_VelInterpolation_1 
#if defined(CFDLPT_TwoWay)
      distribute_FpForce   => distribute_FpForce_1
#endif      
      xmp_interp=y1start(1)
      xep_interp=y1end(1)     
      zmp_interp=y1start(3)
      zep_interp=y1end(3)
      xmc_interp=y1start(1)-1;   if( myProcNghBC(y_pencil,4)<0 ) xmc_interp=1
      xec_interp=y1end(1);       if( myProcNghBC(y_pencil,3)<0 ) xec_interp=nxc -1
      zmc_interp=y1start(3)-1;   if( myProcNghBC(y_pencil,2)<0 ) zmc_interp=1
      zec_interp=y1end(3);       if( myProcNghBC(y_pencil,1)<0 ) zec_interp=nzc -1
             
    elseif(InterpAccuracy==2) then  ! Quadratic interpolation
      PrepareInterpolation => PrepareInterpolation_2
      clc_VelInterpolation => clc_VelInterpolation_2
#if defined(CFDLPT_TwoWay)
      distribute_FpForce   => distribute_FpForce_2
#endif            
      ! ux
      hi_ux_interp%pencil = y_pencil
      hi_ux_interp%xmh=0;  hi_ux_interp%xph=2
      hi_ux_interp%ymh=0;  hi_ux_interp%yph=0
      hi_ux_interp%zmh=0;  hi_ux_interp%zph=0

      ! uz
      hi_uz_interp%pencil = y_pencil
      hi_uz_interp%xmh=0;  hi_uz_interp%xph=0
      hi_uz_interp%ymh=0;  hi_uz_interp%yph=0
      hi_uz_interp%zmh=0;  hi_uz_interp%zph=2

      ! index for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
      xmp_interp=y1start(1)-1;   if( myProcNghBC(y_pencil,4)<0 ) xmp_interp=1
      xep_interp=y1end(1);       if( myProcNghBC(y_pencil,3)<0 ) xep_interp=nxc -1 ! Modified by Zheng Gong,2021-09-23
      zmp_interp=y1start(3)-1;   if( myProcNghBC(y_pencil,2)<0 ) zmp_interp=1
      zep_interp=y1end(3);       if( myProcNghBC(y_pencil,1)<0 ) zep_interp=nzc -1 ! Modified by Zheng Gong,2021-09-23
      xmc_interp=y1start(1)-1;   if( myProcNghBC(y_pencil,4)<0 ) xmc_interp=1
      xec_interp=y1end(1)  -1;   if( myProcNghBC(y_pencil,3)<0 ) xec_interp=nxc -2
      zmc_interp=y1start(3)-1;   if( myProcNghBC(y_pencil,2)<0 ) zmc_interp=1
      zec_interp=y1end(3)  -1;   if( myProcNghBC(y_pencil,1)<0 ) zec_interp=nzc -2
    else
      call MainLog%CheckForError(ErrT_Abort,"InitFpForce","wrong LPT_opt%InterpAccuracy")
    endif
    
    ! calculate inverse distribution retio
    allocate(iDistRatioYp(0:nyp),iDistRatioYc(0:nyp))
    do j=0,nyp
      iDistRatioYp(j)=rdx*rdyc(j)*rdz/FluidDensity ! Note Here
      iDistRatioYc(j)=rdx*rdyp(j)*rdz/FluidDensity 
    enddo
  end subroutine InitFpForce

  !******************************************************************
  ! PrepareInterpolation_1
  !******************************************************************
  subroutine PrepareInterpolation_1()
    implicit none
    
    ! locals
    type(real3)::pos,RatioYp,RatioYc
    integer::pid,ic,jc,kc,nlocal,js,je
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    
    nlocal= GPrtcl_list%nlocal
    if(nlocal > 0) then
      allocate(indxyz(6,nlocal),RatioYp_interp(nlocal),RatioYc_interp(nlocal))
    endif

    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)

      ! if pos%y is within [0,yly), jc will be within [1,nyc]
      js=0
      je=nyp+1
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(pos%y< yp(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      ic= floor(pos%x*rdx)+1; ic=min(ic,y1end(1)); ic=max(ic,y1start(1));
      kc= floor(pos%z*rdz)+1; kc=min(kc,y1end(3)); kc=max(kc,y1start(3));
      if(jc>nyc) call MainLog%CheckForError(ErrT_Abort,"InitDistribute","wrong jc")

      ! index for interpolation as follow, first-order largrange-interpolation:
      idxp_interp=ic
      if(pos%x>xc(ic)) then
        idxc_interp= min(ic,  xec_interp)
      else
        idxc_interp= max(ic-1,xmc_interp)
      endif

      idyp_interp=jc
      if(pos%y>yc(jc)) then
        idyc_interp=min(jc,nyc-1)
      else
        idyc_interp=max(jc-1,1)
      endif

      idzp_interp=kc
      if(pos%z>zc(kc)) then
        idzc_interp= min(kc,  zec_interp)     
      else
        idzc_interp= max(kc-1,zmc_interp)
      endif
      
      indxyz(1,pid)=int(idxc_interp,2)
      indxyz(2,pid)=int(idxp_interp,2)
      indxyz(3,pid)=int(idyc_interp,2)
      indxyz(4,pid)=int(idyp_interp,2)
      indxyz(5,pid)=int(idzc_interp,2)
      indxyz(6,pid)=int(idzp_interp,2)

      RatioYp%x=(yp(idyp_interp+1)-pos%y)*rdyp(idyp_interp) !(y2cord-pos%y)/(y2cord-y1cord)
      RatioYp%y=1.0_RK-RatioYp%x
      RatioYc%x=(yc(idyc_interp+1)-pos%y)*rdyc(idyc_interp+1)
      RatioYc%y=1.0_RK-RatioYc%x
      RatioYp_interp(pid)=RatioYp
      RatioYc_interp(pid)=RatioYc
    ENDDO
  end subroutine PrepareInterpolation_1
    
  !******************************************************************
  ! PrepareInterpolation_2
  !******************************************************************
  subroutine PrepareInterpolation_2()
    implicit none
    
    ! locals
    real(RK)::y1cord,y2cord,y3cord
    type(real3)::pos,RatioYp,RatioYc
    integer::pid,ic,jc,kc,nlocal,js,je
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    
    nlocal= GPrtcl_list%nlocal
    if(nlocal > 0) then
      allocate(indxyz(6,nlocal),RatioYp_interp(nlocal),RatioYc_interp(nlocal))
    endif

    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)

      ! if pos%y is within [0,yly), jc will be within [1,nyc]
      js=0
      je=nyp+1
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(pos%y< yp(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      ic= floor(pos%x*rdx)+1; ic=min(ic,y1end(1)); ic=max(ic,y1start(1));
      kc= floor(pos%z*rdz)+1; kc=min(kc,y1end(3)); kc=max(kc,y1start(3));
      if(jc>nyc) call MainLog%CheckForError(ErrT_Abort,"InitDistribute","wrong jc")

      ! index for interpolation as follow, second-order largrange-interpolation:
      idxc_interp=max(xmc_interp,ic-1)          ! Modified by Zheng Gong,2021-09-23
      idxc_interp=min(xec_interp,idxc_interp)   ! Modified by Zheng Gong,2021-09-23
      if(pos%x>xc(ic)) then
        idxp_interp= min(ic,  xep_interp)       ! Modified by Zheng Gong,2021-09-23, "xep_interp-2" => "xep_interp"
      else
        idxp_interp= max(ic-1,xmp_interp)
      endif

      idyc_interp=max(1,    jc-1)               ! Modified by Zheng Gong,2021-09-23
      idyc_interp=min(nyc-2,idyc_interp)        ! Modified by Zheng Gong,2021-09-23
      if(pos%y>yc(jc)) then
        idyp_interp=min(jc,nyc-1)
      else
        idyp_interp=max(jc-1,1)
      endif

      idzc_interp=max(zmc_interp,kc-1)          ! Modified by Zheng Gong,2021-09-23
      idzc_interp=min(zec_interp,idzc_interp)   ! Modified by Zheng Gong,2021-09-23
      if(pos%z>zc(kc)) then
        idzp_interp= min(kc,  zep_interp)       ! Modified by Zheng Gong,2021-09-23, "zep_interp-2" => "zep_interp"
      else
        idzp_interp= max(kc-1,zmp_interp)
      endif

      indxyz(1,pid)=int(idxc_interp,2)
      indxyz(2,pid)=int(idxp_interp,2)
      indxyz(3,pid)=int(idyc_interp,2)
      indxyz(4,pid)=int(idyp_interp,2)
      indxyz(5,pid)=int(idzc_interp,2)
      indxyz(6,pid)=int(idzp_interp,2)

      y1cord= yp(idyp_interp  )
      y2cord= yp(idyp_interp+1)
      y3cord= yp(idyp_interp+2)
      RatioYp%x=((pos%y-y2cord)*(pos%y-y3cord))/((y1cord-y2cord)*(y1cord-y3cord))
      RatioYp%y=((pos%y-y1cord)*(pos%y-y3cord))/((y2cord-y1cord)*(y2cord-y3cord))
      RatioYp%z=1.0_RK-RatioYp%x-RatioYp%y

      y1cord= yc(idyc_interp  )
      y2cord= yc(idyc_interp+1)
      y3cord= yc(idyc_interp+2)
      RatioYc%x=((pos%y-y2cord)*(pos%y-y3cord))/((y1cord-y2cord)*(y1cord-y3cord))
      RatioYc%y=((pos%y-y1cord)*(pos%y-y3cord))/((y2cord-y1cord)*(y2cord-y3cord))
      RatioYc%z=1.0_RK-RatioYc%x-RatioYc%y

      RatioYp_interp(pid)=RatioYp
      RatioYc_interp(pid)=RatioYc
    ENDDO
  end subroutine PrepareInterpolation_2
  
  !******************************************************************
  ! clc_VelInterpolation_1
  !******************************************************************
  subroutine clc_VelInterpolation_1(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz

    ! locals
    type(real4)::pos
    integer::id,jd,kd,pid,i,j,k,nlocal
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,RatioXc(0:1),RatioYc(0:1),RatioZc(0:1),RatioXp(0:1),RatioYp(0:1),RatioZp(0:1)

    nlocal=GPrtcl_list%nlocal
    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)
      SumXDir=0.0_RK;  SumYDir=0.0_RK;  SumZDir=0.0_RK

      idxc_interp=indxyz(1,pid)
      idxp_interp=indxyz(2,pid)
      idyc_interp=indxyz(3,pid)
      idyp_interp=indxyz(4,pid)
      idzc_interp=indxyz(5,pid)
      idzp_interp=indxyz(6,pid)

      RatioXp(1)=(pos%x-xc(idxp_interp))*rdx+0.5_RK
      RatioXp(0)=1.0_RK-RatioXp(1)
      RatioXc(1)=(pos%x-xc(idxc_interp))*rdx
      RatioXc(0)=1.0_RK-RatioXc(1)
 
      RatioYp(0)= RatioYp_interp(pid)%x
      RatioYp(1)= RatioYp_interp(pid)%y
      RatioYc(0)= RatioYc_interp(pid)%x
      RatioYc(1)= RatioYc_interp(pid)%y

      RatioZp(1)=(pos%z-zc(idzp_interp))*rdz+0.5_RK
      RatioZp(0)=1.0_RK-RatioZp(1)
      RatioZc(1)=(pos%z-zc(idzc_interp))*rdz
      RatioZc(0)=1.0_RK-RatioZc(1)
      
      ! ux grid
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,1
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,1
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,1
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      GPrtcl_Vfluid(pid)  = real3(SumXDir,SumYDir,SumZDir)
    ENDDO
  end subroutine clc_VelInterpolation_1
  
  !******************************************************************
  ! clc_VelInterpolation_2
  !******************************************************************
  subroutine clc_VelInterpolation_2(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz

    ! locals
    type(real4)::pos
    integer::id,jd,kd,pid,i,j,k,nlocal
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,RatioXc(0:2),RatioYc(0:2),RatioZc(0:2),RatioXp(0:2),RatioYp(0:2),RatioZp(0:2)

#ifdef CFDSecondOrder
    call update_halo(ux, mb1, hi_ux_interp)
    call update_halo(uz, mb1, hi_uz_interp)
#endif
    nlocal=GPrtcl_list%nlocal
    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)
      SumXDir=0.0_RK;  SumYDir=0.0_RK;  SumZDir=0.0_RK

      idxc_interp=indxyz(1,pid)
      idxp_interp=indxyz(2,pid)
      idyc_interp=indxyz(3,pid)
      idyp_interp=indxyz(4,pid)
      idzc_interp=indxyz(5,pid)
      idzp_interp=indxyz(6,pid)

      prx=(pos%x-xc(idxp_interp))*rdx+0.5_RK
      RatioXp(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
      RatioXp(1)=      prx*(2.0_RK-prx)
      RatioXp(2)=0.5_RK* prx*(prx-1.0_RK)
      prx=(pos%x-xc(idxc_interp))*rdx
      RatioXc(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
      RatioXc(1)=      prx*(2.0_RK-prx)
      RatioXc(2)=0.5_RK* prx*(prx-1.0_RK)          

      RatioYp(0) = RatioYp_interp(pid)%x
      RatioYp(1) = RatioYp_interp(pid)%y
      RatioYp(2) = RatioYp_interp(pid)%z
      RatioYc(0) = RatioYc_interp(pid)%x
      RatioYc(1) = RatioYc_interp(pid)%y
      RatioYc(2) = RatioYc_interp(pid)%z

      prz=(pos%z-zc(idzp_interp))*rdz+0.5_RK
      RatioZp(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
      RatioZp(1)=      prz*(2.0_RK-prz)
      RatioZp(2)=0.5_RK* prz*(prz-1.0_RK) 
      prz=(pos%z-zc(idzc_interp))*rdz
      RatioZc(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
      RatioZc(1)=      prz*(2.0_RK-prz)
      RatioZc(2)=0.5_RK* prz*(prz-1.0_RK)    
      
      ! ux grid
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      GPrtcl_Vfluid(pid)  = real3(SumXDir,SumYDir,SumZDir)
    ENDDO
  end subroutine clc_VelInterpolation_2

  !******************************************************************
  ! FinalFpForce
  !******************************************************************
  subroutine FinalFpForce()
    implicit none
    
    if(GPrtcl_list%nlocal > 0) then
      deallocate(indxyz,RatioYp_interp,RatioYc_interp)
    endif
  end subroutine FinalFpForce

  !******************************************************************
  ! clc_FpForce
  !******************************************************************
  subroutine clc_FpForce()!clc_FpForce(ux,uy,uz,pressure)
    implicit none
    !real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,pressure

    ! locals
    !type(real4)::pos
    type(real3)::FpVelDiff
    integer::pid,nlocal,itype
    real(RK)::cd,rep,diam,udiff,vdiff,wdiff,veldiff,taup,Mass

    nlocal=GPrtcl_list%nlocal
    DO pid=1,nlocal
      itype= GPrtcl_pType(pid)
      diam = 2.0_RK*GPrtcl_PosR(pid)%w
      taup = LPTProperty%Prtcl_PureProp(itype)%RelaxionTime
      Mass = LPTProperty%Prtcl_PureProp(itype)%Mass

      FpVelDiff = GPrtcl_Vfluid(pid)-GPrtcl_LinVel(1,pid) 
      udiff= FpVelDiff%x; vdiff= FpVelDiff%y; wdiff= FpVelDiff%z
      veldiff=sqrt(udiff*udiff+ vdiff*vdiff+ wdiff*wdiff)

      rep= veldiff *diam/xnu
      cd = 1.0_RK+0.15_RK*(rep**0.687_RK)
      GPrtcl_FpForce(pid)=(Mass*cd/taup)*FpVelDiff
    ENDDO
  end subroutine clc_FpForce

#ifdef CFDLPT_TwoWay
  !******************************************************************
  ! distribute_FpForce_1
  !****************************************************************** 
  subroutine distribute_FpForce_1()
    implicit none

    ! locals
    real(RK)::prx,pry,prz
    type(real3)::Pos,FpForce
    integer::i,j,k,id,jd,kd,pid
    real(RK),dimension(0:1)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    FpForce_x=0.0_RK; FpForce_y=0.0_RK; FpForce_z=0.0_RK
    DO pid=1,GPrtcl_list%nlocal
      Pos     = GPrtcl_PosR(pid)
      FpForce = zero_r3-GPrtcl_FpForce(pid)
      idxc_interp=indxyz(1,pid)
      idxp_interp=indxyz(2,pid)
      idyc_interp=indxyz(3,pid)
      idyp_interp=indxyz(4,pid)
      idzc_interp=indxyz(5,pid)
      idzp_interp=indxyz(6,pid)

      RatioXp(1)=(pos%x-xc(idxp_interp))*rdx+0.5_RK
      RatioXp(0)=1.0_RK-RatioXp(1)
      RatioXc(1)=(pos%x-xc(idxc_interp))*rdx
      RatioXc(0)=1.0_RK-RatioXc(1)
 
      RatioYp(0)= RatioYp_interp(pid)%x
      RatioYp(1)= RatioYp_interp(pid)%y
      RatioYc(0)= RatioYc_interp(pid)%x
      RatioYc(1)= RatioYc_interp(pid)%y

      RatioZp(1)=(pos%z-zc(idzp_interp))*rdz+0.5_RK
      RatioZp(0)=1.0_RK-RatioZp(1)
      RatioZc(1)=(pos%z-zc(idzc_interp))*rdz
      RatioZc(0)=1.0_RK-RatioZc(1)

      ! ux grid
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,1
            id = i+idxp_interp
            prx= RatioXp(i)
            FpForce_x(id,jd,kd)= FpForce_x(id,jd,kd)+ prx*pry* FpForce%x
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd = j+idyp_interp
          pry= RatioYp(j)*prz*iDistRatioYp(jd)
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_y(id,jd,kd)= FpForce_y(id,jd,kd)+ prx*pry* FpForce%y
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,1
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,1
          jd = j+idyc_interp
          pry= RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_z(id,jd,kd)= FpForce_z(id,jd,kd)+ prx*pry* FpForce%z
          enddo
        enddo
      enddo
    ENDDO

    call Gather_Halo_dist_1()
  end subroutine distribute_FpForce_1
  
  !******************************************************************
  ! distribute_FpForce_2
  !****************************************************************** 
  subroutine distribute_FpForce_2()
    implicit none

    ! locals
    real(RK)::prx,pry,prz
    type(real3)::Pos,FpForce
    integer::i,j,k,id,jd,kd,pid
    real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    FpForce_x=0.0_RK; FpForce_y=0.0_RK; FpForce_z=0.0_RK
    DO pid=1,GPrtcl_list%nlocal
      Pos     = GPrtcl_PosR(pid)
      FpForce = zero_r3-GPrtcl_FpForce(pid)
      idxc_interp=indxyz(1,pid)
      idxp_interp=indxyz(2,pid)
      idyc_interp=indxyz(3,pid)
      idyp_interp=indxyz(4,pid)
      idzc_interp=indxyz(5,pid)
      idzp_interp=indxyz(6,pid)
      
      prx=(pos%x-xc(idxp_interp))*rdx+0.5_RK
      RatioXp(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
      RatioXp(1)=      prx*(2.0_RK-prx)
      RatioXp(2)=0.5_RK* prx*(prx-1.0_RK)
      prx=(pos%x-xc(idxc_interp))*rdx
      RatioXc(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
      RatioXc(1)=      prx*(2.0_RK-prx)
      RatioXc(2)=0.5_RK* prx*(prx-1.0_RK)          

      RatioYp(0) = RatioYp_interp(pid)%x
      RatioYp(1) = RatioYp_interp(pid)%y
      RatioYp(2) = RatioYp_interp(pid)%z
      RatioYc(0) = RatioYc_interp(pid)%x
      RatioYc(1) = RatioYc_interp(pid)%y
      RatioYc(2) = RatioYc_interp(pid)%z

      prz=(pos%z-zc(idzp_interp))*rdz+0.5_RK
      RatioZp(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
      RatioZp(1)=      prz*(2.0_RK-prz)
      RatioZp(2)=0.5_RK* prz*(prz-1.0_RK) 
      prz=(pos%z-zc(idzc_interp))*rdz
      RatioZc(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
      RatioZc(1)=      prz*(2.0_RK-prz)
      RatioZc(2)=0.5_RK* prz*(prz-1.0_RK) 

      ! ux grid
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            FpForce_x(id,jd,kd)= FpForce_x(id,jd,kd)+ prx*pry* FpForce%x
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz*iDistRatioYp(jd)
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_y(id,jd,kd)= FpForce_y(id,jd,kd)+ prx*pry* FpForce%y
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_z(id,jd,kd)= FpForce_z(id,jd,kd)+ prx*pry* FpForce%z
          enddo
        enddo
      enddo
    ENDDO
    
    call Gather_Halo_dist_2()
  end subroutine distribute_FpForce_2
#include "LPT_GatherHalo_inc.f90"
#endif
end module LPT_Fpforce

  !******************************************************************
  ! Gather_Halo_dist_1
  !******************************************************************
subroutine Gather_Halo_dist_1()
  implicit none
  type(HaloInfo)::hi_Force
  
  ! Force_x
  hi_Force%pencil = y_pencil
  hi_Force%xmh=0;  hi_Force%xph=1
  hi_Force%ymh=0;  hi_Force%yph=0
  hi_Force%zmh=1;  hi_Force%zph=1
  call sum_halo(FpForce_x,mb1,hi_Force)   
  
  ! Force_y
  hi_Force%pencil = y_pencil
  hi_Force%xmh=1;  hi_Force%xph=1
  hi_Force%ymh=0;  hi_Force%yph=0
  hi_Force%zmh=1;  hi_Force%zph=1
  call sum_halo(FpForce_y,mb1,hi_Force)

  ! Force_z
  hi_Force%pencil = y_pencil
  hi_Force%xmh=1;  hi_Force%xph=1
  hi_Force%ymh=0;  hi_Force%yph=0
  hi_Force%zmh=0;  hi_Force%zph=1
  call sum_halo(FpForce_z,mb1,hi_Force)
end subroutine Gather_Halo_dist_1

!******************************************************************
! Gather_Halo_dist_2
!******************************************************************
subroutine Gather_Halo_dist_2()
  implicit none
  type(HaloInfo)::hi_Force
  
  ! Force_x
  hi_Force%pencil = y_pencil
  hi_Force%xmh=1;  hi_Force%xph=2
  hi_Force%ymh=0;  hi_Force%yph=0
  hi_Force%zmh=1;  hi_Force%zph=1
  call sum_halo(FpForce_x,mb1,hi_Force)   
  
  ! Force_y
  hi_Force%pencil = y_pencil
  hi_Force%xmh=1;  hi_Force%xph=1
  hi_Force%ymh=0;  hi_Force%yph=0
  hi_Force%zmh=1;  hi_Force%zph=1
  call sum_halo(FpForce_y,mb1,hi_Force)

  ! Force_z
  hi_Force%pencil = y_pencil
  hi_Force%xmh=1;  hi_Force%xph=1
  hi_Force%ymh=0;  hi_Force%yph=0
  hi_Force%zmh=1;  hi_Force%zph=2
  call sum_halo(FpForce_z,mb1,hi_Force)
end subroutine Gather_Halo_dist_2

module LPT_Integration
  use m_TypeDef
  use LPT_Property
  use LPT_Variables
  use LPT_Parameters
  implicit none
  private    
  real(RK),parameter,dimension(2):: AB2C = [1.5_RK,-0.5_RK]
  real(RK),parameter,dimension(3):: AB3C = [23.0_RK,-16.0_RK,5.0_RK]/12.0_RK
    
  public::Prtcl_Integrate
contains

  !******************************************************************
  ! Prtcl_Integrate
  !******************************************************************  
  subroutine Prtcl_Integrate(iCountLPT)
    implicit none
    integer,intent(in)::iCountLPT
     
    ! locals
    integer::pid,nlocal,itype
    real(RK)::dt,Mass,TimeIntCoe(3)
    type(real3)::linVel1,linVel2,PosOld,Gravity
    
    dt=LPT_opt%dt
    nlocal = GPrtcl_list%nlocal
    Gravity= LPT_opt%gravity
        
    ! linear position
    if(LPT_Opt%PI_Method==PIM_AB2) then
      if(iCountLPT==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
      else
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
      endif
      DO pid = 1,nlocal 
        itype= GPrtcl_pType(pid)
        Mass = LPTProperty%Prtcl_PureProp(itype)%Mass      
        GPrtcl_linAcc(1,pid)= (1.0_RK/Mass)*GPrtcl_FpForce(pid)+ Gravity
        
        PosOld = GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid) = PosOld

        linVel1=GPrtcl_linVel(1,pid)
        GPrtcl_PosR(pid)=PosOld+(TimeIntCoe(1)*linVel1 + TimeIntCoe(2)*GPrtcl_linVel(2,pid))*dt
        GPrtcl_linVel(1,pid)=linVel1+(TimeIntCoe(1)*GPrtcl_linAcc(1,pid)+TimeIntCoe(2)*GPrtcl_linAcc(2,pid))*dt
        GPrtcl_linVel(2,pid)=linVel1
        GPrtcl_linAcc(2,pid)=GPrtcl_linAcc(1,pid)
      ENDDO

    elseif(LPT_Opt%PI_Method==PIM_AB3 ) then
      if(iCountLPT==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
        TimeIntCoe(3)=0.0_RK
      elseif(iCountLPT==2) then
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
        TimeIntCoe(3)=0.0_RK 
      else
        TimeIntCoe(1)=AB3C(1)
        TimeIntCoe(2)=AB3C(2)
        TimeIntCoe(3)=AB3C(3)    
      endif
      DO pid=1,nlocal
        itype= GPrtcl_pType(pid)
        Mass = LPTProperty%Prtcl_PureProp(itype)%Mass      
        GPrtcl_linAcc(1,pid)= (1.0_RK/Mass)*GPrtcl_FpForce(pid) + Gravity
        
        PosOld = GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid) = PosOld

        linVel1=GPrtcl_linVel(1,pid)
        linVel2=GPrtcl_linVel(2,pid)
                
        GPrtcl_PosR(pid)=PosOld+(TimeIntCoe(1)*linVel1+TimeIntCoe(2)*linVel2+TimeIntCoe(3)*GPrtcl_linVel(3,pid))*dt
        GPrtcl_linVel(1,pid) =linVel1+(TimeIntCoe(1)*GPrtcl_linAcc(1,pid)+TimeIntCoe(2)*GPrtcl_linAcc(2,pid)+ &
                                     TimeIntCoe(3)*GPrtcl_linAcc(3,pid))*dt
        GPrtcl_linVel(3,pid) = linVel2
        GPrtcl_linVel(2,pid) = linVel1
        GPrtcl_linAcc(3,pid) = GPrtcl_linAcc(2,pid)
        GPrtcl_linAcc(2,pid) = GPrtcl_linAcc(1,pid)
      ENDDO
    endif
    
  end subroutine Prtcl_Integrate

end module LPT_Integration

module LPT_IOAndVisu
  use MPI
  use LPT_Comm
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_Variables
  use LPT_Decomp_2d
  use LPT_Parameters
  use m_Decomp2d,only: nrank,nproc
  implicit none
  private

  integer,parameter::IK=4
  integer::Prev_BackUp_itime= 53456791
  logical::saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel, save_LinAcc

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
    procedure:: Write_Restart =>  PIO_Write_Restart
    procedure:: Delete_Prev_Restart =>  PIO_Delete_Prev_Restart
    procedure,private:: Write_XDMF  =>  PIO_Write_XDMF
  end type Prtcl_IO_Visu
  type(Prtcl_IO_Visu),public:: LPT_IO

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_vector,  Prtcl_dump_int_matrix
    module procedure Prtcl_dump_real_vector, Prtcl_dump_real3_vector
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
    integer::nUnitFile,ierror,indent,nflds,ifld
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel, save_LinAcc
    character(128)::XdmfFile
  
    if(iStage==1) then
      open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
      if(ierror/=0)call LPTLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu", "Cannot open file: "//trim(chFile))
      read(nUnitFile, nml=PrtclVisuOption)
      if(nrank==0)write(LPTLogInfo%nUnit, nml=PrtclVisuOption)
      close(nUnitFile,IOSTAT=ierror)
      return
    endif

    ! initialize the XDMF/XDF file
    if(nrank/=0) return
    write(xdmfFile,"(A)") trim(LPT_opt%ResultsDir)//"PartVisuFor"//trim(LPT_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call LPTLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu","Cannot open file: "//trim(XdmfFile))
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! Time series
    indent =  4
    nflds = (LPT_Opt%ilast - LPT_Opt%ifirst +1)/LPT_Opt%SaveVisu  + 1
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    indent = indent + 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Time TimeType="List">'
    indent = indent + 4
    write(nUnitFile,'(A,I6,A)')repeat(' ',indent)//'<DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no') repeat(' ',indent)
    do ifld = 1,nflds
      write(nUnitFile,'(I9)',advance='no') ((ifld-1)*LPT_Opt%SaveVisu + LPT_Opt%ifirst-1)
    enddo
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,fmt='(A)')repeat(' ',indent)//'</Time>'
    close(nUnitFile,IOSTAT=ierror)
    if( .not. saveXDMFOnce) return
    
    do ifld = 1,nflds
      call this%Write_XDMF((ifld-1)*LPT_Opt%SaveVisu + LPT_Opt%ifirst-1)
    enddo

    ! XDMF/XMF Tail
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)    
  end subroutine PIO_Init_visu

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

    !locals
    integer(kind=MPI_OFFSET_KIND)::disp
    integer:: indent,nUnitFile,ierror,np,dims,iprec
    character(128)::XdmfFile

    if(nrank/=0) return 
    np=LPT_Opt%np_InDomain
    write(xdmfFile,"(A)") trim(LPT_opt%ResultsDir)//"PartVisuFor"//trim(LPT_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"PIO_Write_XDMF","Cannot open file: "//trim(XdmfFile))

    indent = 8; disp = 0_MPI_OFFSET_KIND
    write(xdmfFile,"(A)") "PartVisuFor"//trim(LPT_opt%RunName)
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
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(XdmfFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Geometry>'

    IF(save_ID) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"ID","Scalar","Int",disp)
    ENDIF
    IF(save_Diameter) THEN
      dims=1; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"Diameter","Scalar","Float",disp)
    ENDIF
    IF(save_Type) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"Type","Scalar","Int",disp)
    ENDIF
    IF(save_UsrMark) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"UsrMark","Scalar","Int",disp)
    ENDIF
    IF(save_LinVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"LinVel","Vector","Float",disp)
    ENDIF
    IF(save_LinAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"LinAcc","Vector","Float",disp)
    ENDIF
    write(nUnitFile,'(A)')'        </Grid>'
    close(nUnitFile)
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
    !
    write(chFile,"(A,I10.10)") trim(LPT_opt%RestartDir)//"RestartFor"//trim(LPT_opt%RunName),Prev_BackUp_itime 
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    !
    Prev_BackUp_itime = itime
  end subroutine PIO_Delete_Prev_Restart

  !**********************************************************************
  ! PIO_Read_Restart
  !**********************************************************************
  subroutine PIO_Read_Restart(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(128)::chFile
    integer,parameter::NumRead=2000
    real(RK)::xst,xed,yst,yed,zst,zed
    integer,allocatable,dimension(:):: nP_in_bin
    type(real3),allocatable,dimension(:)::real3Vec,PosVec
    integer(kind=8)::disp,disp_pos,disp_int,disp_real3
    integer::itime,nUnit,ierror,nlocal,np,i,k,itype,nlocal_sum,nreal3,tsize,nLeft,nRead,int_t(3)

    itime = LPT_Opt%ifirst - 1
    xst=LPT_decomp%xSt; xed=LPT_decomp%xEd
    yst=LPT_decomp%ySt; yed=LPT_decomp%yEd
    zst=LPT_decomp%zSt; zed=LPT_decomp%zEd

    ! Begin to write Restart file
    write(chFile,"(A,I10.10)") trim(LPT_opt%RestartDir)//"RestartFor"//trim(LPT_opt%RunName),itime
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart","Cannot open file: "//trim(chFile))
    disp =1_8; read(nUnit,pos=disp,IOSTAT=ierror)np; disp=disp+int_byte
    if(np>LPT_Opt%numPrtcl .and. nrank==0) then
      call LPTLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," np_InDomain > numPrtcl " )
    endif
    LPT_Opt%np_InDomain = np

    tsize=GPrtcl_list%tsize
    nreal3 = 2*tsize
    allocate(real3Vec(nreal3),PosVec(NumRead))
    allocate(nP_in_bin(LPT_opt%numPrtcl_Type)); nP_in_bin=0

    nlocal=0; nLeft=np
    disp_pos  = disp
    disp_int  = disp_pos+ real3_byte*np
    disp_real3= disp_int+ int_byte*np*3
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp_pos,IOSTAT=ierror)PosVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      do i=1,nRead
        if(PosVec(i)%x>=xst .and. PosVec(i)%x< xed .and. PosVec(i)%y>=yst .and. &
           PosVec(i)%y< yed .and. PosVec(i)%z>=zst .and. PosVec(i)%z<zed) then
          if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1

          read(nUnit,pos=disp_int,IOSTAT=ierror)int_t(1:3)
          GPrtcl_id(nlocal)=int_t(1)      ! id
          itype=int_t(2)
          GPrtcl_pType(nlocal)=itype      ! pType
          nP_in_bin(itype)= nP_in_bin(itype)+1
          GPrtcl_UsrMark(nlocal)=int_t(3) ! Usr_Mark

          GPrtcl_PosR(nlocal)= PosVec(i)  ! PosR
          GPrtcl_PosR(nlocal)%w=LPTProperty%Prtcl_PureProp(itype)%Radius
          k=0;
          read(nUnit,pos=disp_real3,IOSTAT=ierror)real3Vec(1:nreal3)
          GPrtcl_LinVel(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinVec
          GPrtcl_LinAcc(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinAcc
        endif
        disp_int  = disp_int  + int_byte*3
        disp_real3= disp_real3+ real3_byte*nreal3       
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    deallocate(PosVec,Real3Vec)
    call MPI_ALLREDUCE(nP_in_bin, LPTProperty%nPrtcl_in_Bin,LPT_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    deallocate(nP_in_bin)
    close(nUnit,IOSTAT=ierror)
    GPrtcl_list%nlocal = nlocal
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= np .and. nrank==0) then
      call LPTLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart: "," nlocal_sum/= np_InDomain " )
    endif
  end subroutine PIO_Read_Restart

  !**********************************************************************
  ! PIO_Write_Restart
  !**********************************************************************
  subroutine PIO_Write_Restart(this,itime)
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime
   
    ! locals 
    character(128)::chFile
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    integer,parameter::NumRestart=200
    integer,allocatable,dimension(:,:)::IntMat
    type(real3),allocatable,dimension(:)::real3Vec
    integer(kind=MPI_OFFSET_KIND)::disp,bgn_byte,FileSize
    integer::pid,i,k,nlocal,bgn_ind,prank,pProc,ierror,fh
    integer::tsize,color,key,Prtcl_WORLD,nreal3,nRestart,nLeft

    ! Calculate the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind= clc_bgn_ind(nlocal)

    ! Create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierror)
    if(color==2) return

    ! Begin to write Restart file
    write(chFile,"(A,I10.10)") trim(LPT_opt%RestartDir)//"RestartFor"//trim(LPT_opt%RunName),itime
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_8,ierror)  ! Guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND

    ! Write LPT_Opt%np_InDomain, in the begining of the Restart file
    call MPI_COMM_RANK(Prtcl_WORLD, prank, ierror)
    call MPI_COMM_SIZE(Prtcl_WORLD, pProc, ierror)
    if(prank==0) then
      call MPI_FILE_WRITE_AT(fh,disp,LPT_Opt%np_InDomain,1,int_type,MPI_STATUS_IGNORE,ierror)
    endif
    disp = disp + int_byte
    call MPI_BARRIER(Prtcl_WORLD,ierror)

    ! Begin to write
    pvsize%sizes(1)   = LPT_Opt%np_InDomain
    pvsize%subsizes(1)= nlocal
    pvsize%starts(1)  = bgn_ind
    allocate(real3Vec(nlocal))
    do pid=1,nlocal
      real3Vec(pid)=GPrtcl_PosR(pid)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)
    
    pmsize%sizes(1)   = 3;   pmsize%sizes(2)   = LPT_Opt%np_InDomain
    pmsize%subsizes(1)= 3;   pmsize%subsizes(2)= nlocal
    pmsize%starts(1)  = 0;   pmsize%starts(2)  = bgn_ind
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
    nreal3 = 2*tsize
    call MPI_FILE_OPEN(Prtcl_WORLD,chFile, MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierror)
    call MPI_FILE_GET_SIZE(fh,FileSize,ierror)
    FileSize=FileSize+int(nreal3*real3_byte,8)*int(LPT_Opt%np_InDomain,8)
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    
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
      enddo
      call MPI_FILE_WRITE_AT(fh,bgn_byte,real3Vec,k,real3_type,MPI_STATUS_IGNORE,ierror)
      bgn_byte=bgn_byte+int(nreal3*real3_byte,8)*int(nRestart,8)
      nLeft=nLeft-nRestart
      if(nLeft==0)exit
    ENDDO
    deallocate(real3Vec)
    call MPI_FILE_CLOSE(fh, ierror)
    call MPI_COMM_FREE(Prtcl_WORLD,ierror)
  end subroutine PIO_Write_Restart

  !**********************************************************************
  ! PIO_Init_visu
  !**********************************************************************
  subroutine PIO_Final_visu(this)
    implicit none 
    class(Prtcl_IO_Visu)::this

    ! locals
    integer::nUnitFile,ierror
    character(128)::XdmfFile

    if(nrank/=0 .or. saveXDMFOnce) return
    write(xdmfFile,"(A)") trim(LPT_opt%ResultsDir)//"PartVisuFor"//trim(LPT_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"PIO_Final_visu","Cannot open file:  "//trim(XdmfFile))
    ! XDMF/XMF Tail
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)

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
    integer :: ierror,fh,i,color,key,Prtcl_WORLD,nlocal,bgn_ind

    ! write xdmf file first
    if(.not.saveXDMFOnce) call this%Write_XDMF(itime)

    ! update the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind=clc_bgn_ind(nlocal)

    ! create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierror)
    if(color==2) return
    
    ! begin to dump
    write(chFile,"(A,I10.10)") trim(LPT_opt%ResultsDir)//"PartVisuFor"//trim(LPT_opt%RunName),itime
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_8,ierror)  ! guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND
    pvsize%sizes(1)     = LPT_Opt%np_InDomain
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
    integer:: ierror,newtype
    integer,dimension(1) :: sizes, subsizes, starts

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
    integer:: ierror,newtype
    integer,dimension(1) :: sizes, subsizes, starts

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
    disp = disp + sizes(1) * real3_byte
  end subroutine Prtcl_dump_real3_vector

end module LPT_IOAndVisu

module LPT_Parameters
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d,only:nrank
#ifdef CFDSecondOrder
  use m_Parameters,only: BcOption
#endif
  use m_Parameters,only: dtMax, ifirst, ilast, BackupFreq,xlx,yly,zlz,SaveVisu
  implicit none
  private
  
  ! Log
  type(LogType),public::LPTLogInfo
  
  integer,parameter,public:: x_axis = 1
  integer,parameter,public:: y_axis = 2
  integer,parameter,public:: z_axis = 3 
    
  integer,parameter,public:: PIM_AB2 = 2
  integer,parameter,public:: PIM_AB3 = 3
   
  ! default values  
  type LPT_Options
    logical:: RestartFlag=.false.
    integer:: numPrtcl    = 8000     ! total particle number
    integer:: np_InDomain            ! particle in domain
    integer:: ifirst                 ! first time step
    integer:: ilast                  ! last time step
    real(RK):: dt   =  1.0E-5_RK     ! time step 
    type(real3):: SimDomain_min
    type(real3):: SimDomain_max
    type(real3):: gravity = real3(0.0_RK,-9.81_RK,0.0_RK) ! gravity or other constant body forces if any
    logical,dimension(3):: IsPeriodic = .false.

    integer:: PI_Method = PIM_AB2      ! integration scheme for translational motion   
    integer:: numPrtcl_Type=1          ! number of particle type 
    character(64)::RunName  = "LPTRun" ! run name
    character(64)::ResultsDir  = "."   ! result directory 
    character(64)::RestartDir="."      ! restart directory
    integer:: SaveVisu      = 1000     ! save frequency for visulizing file
    integer:: BackupFreq    = 100000   ! save frequency for restarting file
    integer:: Cmd_LFile_Freq= 500      ! report frequency in the terminal 
    integer:: LF_file_lvl   = 5        ! logfile report level      
    integer:: LF_cmdw_lvl   = 3        ! terminal report level
  contains 
    procedure :: ReadLPTOption => LO_ReadLPTOption
  end type LPT_Options
  type(LPT_Options),public::  LPT_opt
    
contains

  !**********************************************************************
  ! LO_ReadLPTOption
  !**********************************************************************
  subroutine LO_ReadLPTOption(this, chFile)
    implicit none
    class(LPT_Options):: this
    character(*),intent(in)::chFile
           
    ! locals
    real(RK)::gravity(3)
    logical::RestartFlag
    character(64):: RunName, ResultsDir,RestartDir
    integer::numPrtcl,SaveVisuLPT,PI_Method,numPrtcl_Type,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,nUnitFile,ierror
    NAMELIST /LPTOptions/RestartFlag,numPrtcl,gravity,PI_Method,numPrtcl_Type,RunName,ResultsDir, &
                         RestartDir,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl
               
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) then
       print*, "Cannot open file: "//trim(chFile); STOP
    endif
    read(nUnitFile, nml=LPTOptions)
    close(nUnitFile,IOSTAT=ierror)
    
    SaveVisuLPT=SaveVisu
    this%RestartFlag = RestartFlag
    this%numPrtcl    = numPrtcl
    this%dt       = dtMax
    this%ifirst   = ifirst
    this%ilast    = ilast
    this%gravity  = gravity
    this%SimDomain_min = zero_r3
    this%SimDomain_max = real3(xlx,yly,zlz)

#ifdef CFDSecondOrder
    if(BcOption(1)==0)this%IsPeriodic(1)=.true.
    if(BcOption(3)==0)this%IsPeriodic(2)=.true.
    if(BcOption(5)==0)this%IsPeriodic(3)=.true.
#elif CFDFourthOrder
    this%IsPeriodic(1)=.true.
    this%IsPeriodic(2)=.false.
    this%IsPeriodic(3)=.true.
#endif
           
    this%PI_Method = PI_Method
    this%numPrtcl_Type = numPrtcl_Type
           
    write(this%RunName,"(A)") RunName
    write(this%ResultsDir,"(A)") ResultsDir
    write(this%RestartDir,"(A)") RestartDir
    this%SaveVisu = SaveVisuLPT
    this%BackupFreq = BackupFreq
    this%Cmd_LFile_Freq = Cmd_LFile_Freq
    this%LF_file_lvl = LF_file_lvl
    this%LF_cmdw_lvl = LF_cmdw_lvl

  end subroutine LO_ReadLPTOption
end module LPT_Parameters

module LPT_Property
  use m_TypeDef
  use m_LogInfo
  use LPT_Parameters
  use m_Decomp2d,only: nrank
  use m_Parameters,only: xnu,FluidDensity
  implicit none
  private     
    
  type PureProperty
    real(RK):: Radius
    real(RK):: Density
    real(RK):: Mass
    real(RK):: Volume
    real(RK):: RelaxionTime
  end type PureProperty
    
  type PhysicalProperty
    integer,allocatable,dimension(:) :: nPrtcl_in_Bin  
    type(pureProperty),allocatable,dimension(:) :: Prtcl_PureProp     
  contains
    procedure:: InitPrtclProperty
  end type PhysicalProperty
  type(PhysicalProperty),public::LPTProperty

contains

  !*******************************************************
  ! initializing the size distribution with property 
  !*******************************************************
  subroutine InitPrtclProperty( this,  chFile )
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK):: sum_divided,rtemp 
    real(RK),dimension(:),allocatable:: Bin_Divided, Density, Diameter
    namelist/ParticlePhysicalProperty/Bin_Divided, Density, Diameter
    integer:: i,j,iTV(8),nPType,nUnitFile,ierror,sum_prtcl,bin_pnum,prdiff,bin_id
        
    nPType  = LPT_opt%numPrtcl_Type
    allocate( Bin_Divided(nPType))
    allocate( Density(nPType))
    allocate( Diameter(nPType))
    allocate( this%nPrtcl_in_Bin(nPType))
    allocate( this%Prtcl_PureProp(nPType))
          
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror/=0.and.nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=ParticlePhysicalProperty)
    if(nrank==0)write(LPTLogInfo%nUnit, nml=ParticlePhysicalProperty)
    close(nUnitFile, IOSTAT=ierror)
        
    ! calculate this%nPrtcl_in_Bin
    sum_divided=0.0_RK
    do i=1,nPType
      sum_divided = sum_divided + Bin_Divided(i)
    enddo
    sum_prtcl = 0
    do i=1,nPType
      bin_pnum = int(LPT_opt%numPrtcl*Bin_Divided(i)/sum_divided)
      this%nPrtcl_in_Bin(i)= bin_pnum
      sum_prtcl = sum_prtcl + bin_pnum
    enddo
    prdiff = LPT_opt%numPrtcl - sum_prtcl
    if(prdiff>0) then
      call date_and_time(values=iTV); !iTV=0
      call random_seed(size= i)
      call random_seed(put = iTV(7)*iTV(8)+[(j,j=1,i)])
      do i=1, prdiff 
        call random_number(rtemp)
        bin_id = int(rtemp*nPType) + 1
        this%nPrtcl_in_Bin(bin_id) = this%nPrtcl_in_Bin(bin_id)  + 1
      enddo
     endif
        
     ! calculate particle properties
     do i = 1, nPType
       this%Prtcl_PureProp(i)%Density= Density(i)
       this%Prtcl_PureProp(i)%Radius = 0.5_RK*Diameter(i)
       this%Prtcl_PureProp(i)%Volume = 1.3333333333333_RK*Pi*(this%Prtcl_PureProp(i)%Radius)**3
       this%Prtcl_PureProp(i)%Mass   = Density(i)*this%Prtcl_PureProp(i)%Volume
       this%Prtcl_PureProp(i)%RelaxionTime= Density(i)*Diameter(i)*Diameter(i)/(18.00_RK*xnu*FluidDensity)
     enddo       
  end subroutine InitPrtclProperty

end module LPT_Property

module LPT_Statistics
  use MPI
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_Decomp_2d
  use LPT_Variables
  use LPT_Parameters
  use m_Decomp2d,only:nrank
  use m_MeshAndMetries,only:yp
  use m_Parameters,only: ivstats,saveStat,nyc,itime,yly,ilast 
  implicit none
  private

  integer::nslab,nShannon
  integer::npType,npstime
  real(RK),dimension(:),allocatable:: ypForPs  ! y point for particle statistics 
  integer,dimension(:,:),allocatable::npsum
  real(RK),dimension(:,:),allocatable::upsum,vpsum,wpsum,upupsum,vpvpsum,wpwpsum,upvpsum
  real(RK),dimension(:,:),allocatable::ufsum,vfsum,wfsum,ufufsum,vfvfsum,wfwfsum,ufvfsum

  public::InitLPTStatistics,ClcLPTStatistics
contains
  
!********************************************************************************
!   InitLPTStatistics
!********************************************************************************
  subroutine InitLPTStatistics(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    character(len=128)::filename
    integer:: nUnit,ierror,j,k
    integer:: iErr01,iErr02,iErr03,iErr04,iErr05,iErrSum
    NAMELIST/ParticleStatisticOption/nslab,nShannon

    open(newunit=nUnit, file=chFile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics", "Cannot open file: "//trim(chFile))
    read(nUnit, nml=ParticleStatisticOption)
    if(nrank==0)write(LPTLogInfo%nUnit, nml=ParticleStatisticOption)
    close(nUnit,IOSTAT=ierror)
    if(nrank==0 .and. mod(nyc,nslab)/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics","mod(nyc,nslab)/=0")
    if(nrank==0 .and. (nShannon>nyc .or. nShannon<0)) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics","nShannon wrong")    
    
    npType= LPT_opt%numPrtcl_Type
    allocate(ypForPs(nslab+1),      npsum(npType,nslab),   Stat=iErr01)
    allocate(upsum(npType,nslab),   vpsum(npType,nslab),   wpsum(npType,nslab),   Stat=iErr02)
    allocate(ufsum(npType,nslab),   vfsum(npType,nslab),   wfsum(npType,nslab),   Stat=iErr03)
    allocate(upupsum(npType,nslab), vpvpsum(npType,nslab), wpwpsum(npType,nslab), upvpsum(npType,nslab), Stat=iErr04)
    allocate(ufufsum(npType,nslab), vfvfsum(npType,nslab), wfwfsum(npType,nslab), ufvfsum(npType,nslab), Stat=iErr05)
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)
    if(iErrSum/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"LPT_InitStat ","Allocation failed")
    
    k=nyc/nslab
    do j=0,nslab
      ypForPs(j+1)=yp(j*k+1)
    enddo
    call ResetStatVar()
    
    if(nrank/=0) return
    write(filename,'(A,I10.10)') trim(LPT_Opt%ResultsDir)//'shannon',ilast
    open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics","Cannot open file: "//trim(filename))
    close(nUnit,IOSTAT=ierror)    
  end subroutine InitLPTStatistics

  !******************************************************************
  ! ResetStatVar
  !******************************************************************
  subroutine ResetStatVar()
    implicit none

    npstime=0;       npsum=0
    upsum=0.0_RK;      vpsum=0.0_RK;      wpsum=0.0_RK    
    ufsum=0.0_RK;      vfsum=0.0_RK;      wfsum=0.0_RK
    upupsum=0.0_RK;    vpvpsum=0.0_RK;    wpwpsum=0.0_RK;    upvpsum=0.0_RK;
    ufufsum=0.0_RK;    vfvfsum=0.0_RK;    wfwfsum=0.0_RK;    ufvfsum=0.0_RK;
  end subroutine ResetStatVar

  !********************************************************************************
  ! ClcLPTStatistics
  !********************************************************************************  
  subroutine ClcLPTStatistics()
    implicit none

    ! locals
    character(len=128)::filename
    type(real3)::pos,VPrtcl,VFluid
    real(RK)::irnpsum,inpstime,dyShannon,pse
    real(RK),dimension(npType)::shannon_entropy
    integer::pid,nlocal,itype,js,je,jc,ierror,nUnit
    real(RK),dimension(14,npType,nslab)::sumStat
    real(RK),dimension(:,:,:),allocatable::sumStatR
    integer,dimension(:,:),allocatable::npsslot,npsslotR
    
    allocate(npsslot(npType,nShannon),npsslotR(npType,nShannon))
    sumStat=0.0_RK; npsslot=0
    nlocal= GPrtcl_list%nlocal
    dyShannon=yly/real(nShannon,RK)
    do pid=1,nlocal
      pos  = GPrtcl_posR(pid)
      itype= GPrtcl_pType(pid)
      VPrtcl=GPrtcl_linVel(1,pid)
      VFluid=GPrtcl_VFluid(pid)

      ! Shannon entropy
      jc=int(pos%y/dyShannon)+1
      npsslot(itype,jc)=npsslot(itype,jc)+1
      
      ! if pos%y is within [0,yly), jc will be within [1,nslab]
      js=0
      je=nslab+2
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(pos%y< ypForPs(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      npsum(itype,jc)     = npsum(itype,jc) + 1
      sumStat(1,itype,jc) = sumStat(1,itype,jc) + VPrtcl%x
      sumStat(2,itype,jc) = sumStat(2,itype,jc) + VPrtcl%y
      sumStat(3,itype,jc) = sumStat(3,itype,jc) + VPrtcl%z
      sumStat(4,itype,jc) = sumStat(4,itype,jc) + VPrtcl%x *VPrtcl%x
      sumStat(5,itype,jc) = sumStat(5,itype,jc) + VPrtcl%y *VPrtcl%y
      sumStat(6,itype,jc) = sumStat(6,itype,jc) + VPrtcl%z *VPrtcl%z
      sumStat(7,itype,jc) = sumStat(7,itype,jc) + VPrtcl%x *VPrtcl%y
      sumStat(8,itype,jc) = sumStat(8,itype,jc) + VFluid%x
      sumStat(9,itype,jc) = sumStat(9,itype,jc) + VFluid%y
      sumStat(10,itype,jc)= sumStat(10,itype,jc)+ VFluid%z
      sumStat(11,itype,jc)= sumStat(11,itype,jc)+ VFluid%x *VFluid%x
      sumStat(12,itype,jc)= sumStat(12,itype,jc)+ VFluid%y *VFluid%y
      sumStat(13,itype,jc)= sumStat(13,itype,jc)+ VFluid%z *VFluid%z
      sumStat(14,itype,jc)= sumStat(14,itype,jc)+ VFluid%x *VFluid%y
    enddo
    do jc=1,nslab
      do itype=1,npType
        upsum(itype,jc)  = upsum(itype,jc)  + sumStat( 1,itype,jc)
        vpsum(itype,jc)  = vpsum(itype,jc)  + sumStat( 2,itype,jc)
        wpsum(itype,jc)  = wpsum(itype,jc)  + sumStat( 3,itype,jc)
        upupsum(itype,jc)= upupsum(itype,jc)+ sumStat( 4,itype,jc)
        vpvpsum(itype,jc)= vpvpsum(itype,jc)+ sumStat( 5,itype,jc)
        wpwpsum(itype,jc)= wpwpsum(itype,jc)+ sumStat( 6,itype,jc)
        upvpsum(itype,jc)= upvpsum(itype,jc)+ sumStat( 7,itype,jc)
        ufsum(itype,jc)  = ufsum(itype,jc)  + sumStat( 8,itype,jc)
        vfsum(itype,jc)  = vfsum(itype,jc)  + sumStat( 9,itype,jc)
        wfsum(itype,jc)  = wfsum(itype,jc)  + sumStat(10,itype,jc)
        ufufsum(itype,jc)= ufufsum(itype,jc)+ sumStat(11,itype,jc)
        vfvfsum(itype,jc)= vfvfsum(itype,jc)+ sumStat(12,itype,jc)
        wfwfsum(itype,jc)= wfwfsum(itype,jc)+ sumStat(13,itype,jc)
        ufvfsum(itype,jc)= ufvfsum(itype,jc)+ sumStat(14,itype,jc)
      enddo
    enddo
    call MPI_REDUCE(npsslot,npsslotR, npType*nShannon,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
       
    ! Shannon entropy part
    shannon_entropy = 0.0_RK    
    IF(nrank==0) THEN
      DO itype=1,npType
        nlocal=LPTProperty%nPrtcl_in_Bin(itype)
        do jc=1,nShannon
          pse = real(npsslotR(itype,jc),RK)/real(nlocal,RK)
          if(pse > 0.0_RK) then
            shannon_entropy(itype) = shannon_entropy(itype) -pse*log(pse)
          endif
        enddo
        shannon_entropy(itype) = shannon_entropy(itype)/log(real(nShannon,RK))
      ENDDO
      
      write(filename,'(A,I10.10)') trim(LPT_Opt%ResultsDir)//'shannon',ilast
      open(newunit=nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierror)
      if(ierror /= 0) then
        call LPTLogInfo%CheckForError(ErrT_Pass,"ClcLPTStatistics","Cannot open file: "//trim(filename))
      else  
        write(nUnit,*)itime,shannon_entropy
      endif
      close(nUnit,IOSTAT=ierror)
    ENDIF
    deallocate(npsslot,npsslotR)
    
    npstime = npstime + 1
    if(mod(itime,SaveStat)/=0) return
    
    sumStat( 1,:,:)=upsum;    sumStat( 2,:,:)=vpsum;    sumStat( 3,:,:)=wpsum
    sumStat( 4,:,:)=upupsum;  sumStat( 5,:,:)=vpvpsum;  sumStat( 6,:,:)=wpwpsum;  sumStat( 7,:,:)=upvpsum
    sumStat( 8,:,:)=ufsum;    sumStat( 9,:,:)=vfsum;    sumStat(10,:,:)=wfsum
    sumStat(11,:,:)=ufufsum;  sumStat(12,:,:)=vfvfsum;  sumStat(13,:,:)=wfwfsum;  sumStat(14,:,:)=ufvfsum
    allocate(npsslot(npType,nslab),sumStatR(14,npType,nslab))
    call MPI_REDUCE(npsum,  npsslot,    npType*nslab, int_type, MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(sumStat,sumStatR,14*npType*nslab, real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)

    if(nrank==0) then
      inpstime = 1.0_RK/real(npstime,RK)
      do itype=1,npType
        write(filename,"(A,I10.10,A,I4.4)") trim(LPT_Opt%ResultsDir)//'pstats',itime,'_set',itype
        open(newunit=nUnit,file=trim(filename),status='replace',form='formatted',IOSTAT=ierror)
        if(ierror /= 0) then
          call LPTLogInfo%CheckForError(ErrT_Pass,"ClcLPTStatistics","Cannot open file: "//trim(filename))
        else              
          write(nUnit,'(A,I7,A,I7,A,I7)')'    The time step range for this particle statistics is ', &
                                       itime-(npstime-1)*ivstats, ':', ivstats, ':', itime
          write(nUnit,*)
          write(nUnit,'(A)')'  yp, np, up, vp, wp, upup, vpvp, wpwp, upvp, uf, vf, wf, ufuf, vfvf, wfwf, ufvf' 
          do jc=1,nslab
            if(npsslot(itype,jc)==0) then
              write(nUnit,'(20ES24.15)')(ypForPs(jc)+ypForPs(jc+1))*0.5_RK,   & ! 1
                                     0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,   &
                                     0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK
            else
              irnpsum=1.0_RK/real(npsslot(itype,jc),RK)
              write(nUnit,'(20ES24.15)') (ypForPs(jc)+ypForPs(jc+1))*0.5_RK, & ! 1
                                        real(npsslot(itype,jc),RK)*inpstime, & ! 2
                                              sumStatR( 1,itype,jc)*irnpsum, & ! 3
                                              sumStatR( 2,itype,jc)*irnpsum, & ! 4
                                              sumStatR( 3,itype,jc)*irnpsum, & ! 5
                                              sumStatR( 4,itype,jc)*irnpsum, & ! 6
                                              sumStatR( 5,itype,jc)*irnpsum, & ! 7
                                              sumStatR( 6,itype,jc)*irnpsum, & ! 8
                                              sumStatR( 7,itype,jc)*irnpsum, & ! 9
                                              sumStatR( 8,itype,jc)*irnpsum, & ! 10
                                              sumStatR( 9,itype,jc)*irnpsum, & ! 11
                                              sumStatR(10,itype,jc)*irnpsum, & ! 12
                                              sumStatR(11,itype,jc)*irnpsum, & ! 13
                                              sumStatR(12,itype,jc)*irnpsum, & ! 14
                                              sumStatR(13,itype,jc)*irnpsum, & ! 15
                                              sumStatR(14,itype,jc)*irnpsum    ! 16
            endif
          enddo
        endif
        close(nUnit,IOSTAT=ierror)
      enddo
    endif

    call ResetStatVar()
  end subroutine ClcLPTStatistics

end module LPT_Statistics

module LPT_System
  use MPI
  use m_Timer
  use LPT_Comm
  use m_TypeDef
  use LPT_Property
  use LPT_Geometry
  use LPT_decomp_2d
  use LPT_Variables
  use LPT_IOAndVisu
  use LPT_Parameters
  use LPT_Statistics
  use LPT_Integration
  use LPT_ContactSearchPW
  use m_Decomp2d,only:nrank
  use m_Parameters,only:ivstats
  implicit none
  private
    
  !// LPTSystem class 
  type LPTSystem
    integer :: iterNumber   = 0  ! iteration number 
        
    !// timers
    type(timer):: m_total_timer
    type(timer):: m_integration_timer
    type(timer):: m_write_prtcl_timer
    type(timer):: m_comm_exchange_timer
  contains
    procedure:: Initialize => LPT_Initialize
    procedure:: iterate    => LPT_iterate
  end type LPTSystem
  type(LPTSystem),public::LPT
  
  integer::iCountLPT
contains

!********************************************************************************
!   Initializing LPTSystem object with particles which are inserted from a 
!   predefined plane 
!********************************************************************************
  subroutine  LPT_Initialize(this,chLPTPrm)
    implicit none
    class(LPTSystem)::this
    character(*),intent(in)::chLPTPrm
    
    ! locals
    integer::ierror
    character(256)::chStr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    !// Initializing main log info and visu
    iCountLPT=0
    if(LPT_Opt%RestartFlag) iCountLPT=10
    this%IterNumber=LPT_Opt%ifirst-1
    write(chStr,"(A)") 'mkdir -p '//LPT_opt%ResultsDir//' '//LPT_opt%RestartDir//' 2> /dev/null'
    if (nrank==0) call system(trim(adjustl(chStr)))
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call LPTLogInfo%InitLog(LPT_opt%ResultsDir,LPT_opt%RunName,LPT_opt%LF_file_lvl,LPT_opt%LF_cmdw_lvl)
    if(nrank==0) call LPTLogInfo%CreateFile(LPT_opt%RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call LPTLogInfo%OpenFile()
    if(nrank==0)call Write_LPT_Opt_to_Log()
    call InitLPTStatistics(chLPTPrm)

    ! Step1: Physical property
    call LPTProperty%InitPrtclProperty(chLPTPrm)
    if(nrank==0) then
      call LPTLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call LPTLogInfo%OutInfo("Physical properties contains "// trim( num2str(LPT_opt%numPrtcl_Type ) ) //" particle types ",2)
    endif

    ! Step2: set the geometry
    call LPTGeometry%MakeGeometry()
    if(nrank==0) then
      call LPTLogInfo%OutInfo("Step2: Geometry is set", 1 )
      call LPTLogInfo%OutInfo("Geometry Contains "//trim(num2str(LPTGeometry%num_pWall))//" Plane walls.", 2)
    endif

    ! Step3: initilize all the particle variables
    call GPrtcl_list%AllocateAllVar()
    call LPT_IO%Init_visu(chLPTPrm,1)
    t_restart1=MPI_WTIME()
    if(.not.LPT_Opt%RestartFlag) then
      if(LPT_Opt%numPrtcl>0) call GPrtcl_list%MakingAllPrtcl(chLPTPrm)
      if(nrank==0) then
        call LPTLogInfo%OutInfo("Step3: Initial Particle coordinates are MAKING into LPTSystem ...", 1 )
        call LPTLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(LPT_opt%numPrtcl)),2)
      endif
      LPT_opt%np_InDomain = LPT_opt%numPrtcl
    else
      if(LPT_Opt%numPrtcl>0) call LPT_IO%Read_Restart()
      if(nrank==0) then
        call LPTLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call LPTLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(LPT_opt%np_InDomain)),2)
      endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1
    if(nrank==0 .and. LPT_Opt%RestartFlag) call LPTLogInfo%OutInfo("Restart time [sec] :"//trim(num2str(t_res_tot)),2)

    ! Step4: Initializing visu
    call LPT_IO%Init_visu(chLPTPrm,2)
    
    ! Step5: initialize the inter-processors communication
    call LPTComm%InitComm()
    if(nrank==0) call LPTLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    
    ! Step6: timers for recording the execution time of different parts of program
    if(nrank==0) call LPTLogInfo%OutInfo("Step5: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_comm_exchange_timer%reset()
    call this%m_write_prtcl_timer%reset()
  end subroutine LPT_Initialize

  !********************************************************************************
  !   iterating over time 
  !   calls all the required methods to do numIter iterations in the LPT system
  !********************************************************************************
  subroutine LPT_iterate(this,itime)
    implicit none
    class(LPTSystem) this
    integer,intent(in)::itime

    ! locals
    integer::ierror
    
    ! body
    call this%m_total_timer%start()

    ! correcting position and velocities 
    call this%m_integration_timer%start()
    iCountLPT=iCountLPT+1
    call Prtcl_Integrate(iCountLPT)
    call LPTContactSearchPW%FindContactsPW()
    call this%m_integration_timer%finish()

    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call LPTComm%Comm_For_Exchange()
    call this%m_comm_exchange_timer%finish()
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, LPT_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    if( mod(itime,LPT_opt%SaveVisu)== 0)   call LPT_IO%dump_visu(itime)
    if( mod(itime,LPT_opt%BackupFreq)== 0 .or. itime==LPT_opt%ilast) then
      call LPT_IO%Write_Restart(itime)
      call LPT_IO%Delete_Prev_Restart(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()
    
    if(mod(itime,ivstats)==0) call ClcLPTStatistics()
    
    ! output to log file and terminal/command window
    IF((this%IterNumber==1 .or. mod(itime,LPT_opt%Cmd_LFile_Freq)==0) ) THEN
      if(nrank/=0) return
    
      ! command window and log file output
      call LPTLogInfo%OutInfo("LPT performed "//trim(num2str(itime))//" iterations up to here!",1)
      call LPTLogInfo%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(this%m_total_timer%tot_time))//", "// &
      trim(num2str(this%m_total_timer%last_time ))//", "//trim(num2str(this%m_total_timer%average())),2)

      call LPTLogInfo%OutInfo("Integration time [tot, ave]        : "//trim(num2str(this%m_integration_timer%tot_time))//", "// &
      trim(num2str(this%m_integration_timer%average())), 3)

      call LPTLogInfo%OutInfo("Comm_For_Exchange [tot, ave]       : "//trim(num2str(this%m_comm_exchange_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_exchange_timer%average())), 3)

      call LPTLogInfo%OutInfo("Write to file time [tot, ave]      : "//trim(num2str(this%m_write_prtcl_timer%tot_time))//", "// &
      trim(num2str(this%m_write_prtcl_timer%average())), 3)
     
      call LPTLogInfo%OutInfo("Particle number in  domain:  "//trim(num2str(LPT_Opt%np_InDomain)), 2)        
    ENDIF

  end subroutine LPT_iterate

  !**********************************************************************
  ! Write_LPT_Opt_to_Log
  !**********************************************************************
  subroutine Write_LPT_Opt_to_Log()
    implicit none

    ! locals
    real(RK)::dtLPT
    logical::RestartFlag,IsPeriodic(3)
    type(real3):: gravity, minpoint, maxpoint
    character(64):: RunName, ResultsDir,RestartDir
    integer::SaveVisuLPT,BackupFreqLPT, Cmd_LFile_Freq, LF_file_lvl, LF_cmdw_lvl
    integer::numPrtcl,PI_Method,numPrtcl_Type,ifirstLPT,ilastLPT
    NAMELIST /LPTOptions/ RestartFlag,numPrtcl,dtLPT,gravity,minpoint,maxpoint,PI_Method,numPrtcl_Type,RunName,RestartDir,   &
                          ResultsDir,BackupFreqLPT,SaveVisuLPT,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,ifirstLPT,ilastLPT,IsPeriodic

    RestartFlag = LPT_opt%RestartFlag 
    numPrtcl    = LPT_opt%numPrtcl 
    dtLPT       = LPT_opt%dt       
    ifirstLPT   = LPT_opt%ifirst   
    ilastLPT    = LPT_opt%ilast    
    gravity     = LPT_opt%gravity  
    minpoint    = LPT_opt%SimDomain_min 
    maxpoint    = LPT_opt%SimDomain_max 
    IsPeriodic  = LPT_opt%IsPeriodic 
    PI_Method   = LPT_opt%PI_Method
    write(RunName,"(A)")LPT_opt%RunName 
    write(ResultsDir,"(A)")LPT_opt%ResultsDir 

    numPrtcl_Type = LPT_opt%numPrtcl_Type
    write(RestartDir,"(A)") LPT_opt%RestartDir 
    BackupFreqLPT = LPT_opt%BackupFreq 
    SaveVisuLPT   = LPT_opt%SaveVisu 
    Cmd_LFile_Freq= LPT_opt%Cmd_LFile_Freq 
    LF_file_lvl   = LPT_opt%LF_file_lvl 
    LF_cmdw_lvl   = LPT_opt%LF_cmdw_lvl
    write(LPTLogInfo%nUnit, nml=LPTOptions)
  end subroutine Write_LPT_Opt_to_Log

end module LPT_System
