program main_channel3d_2nd
  use MPI
  use m_Decomp2d
  use m_Variables
  use m_IOAndVisu
  use m_Parameters
  use m_ChannelSystem
  use m_Timer,only:time2str  
  use m_Poisson,only:Destory_Poisson_FFT_Plan
  implicit none
  integer::intT,ierror
  character(len=128)::chPrm
  character(len=10)::RowColStr
  
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! read Channel options
  intT=command_argument_count()
  if((intT/=1 .and. intT/=3) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif
  call get_command_argument(1,chPrm)
  call ReadAndInitParameters(chPrm)
  if(intT==3) then
    call get_command_argument(2,RowColStr)
    read(RowColStr,*) p_row
    call get_command_argument(3,RowColStr)
    read(RowColStr,*) p_col
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  ! Initialize Decomp-2d
  call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
 
  call ChannelInitialize(chPrm)    ! Topest level initialing for Channel body
  do itime=ifirst, ilast
    call ChannelIterate()
  enddo
  if(nrank==0)call MainLog%OutInfo("Good job! Channel3d_2nd finished successfully at "//time2str(),1)

  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize()
  call MPI_FINALIZE(ierror)
end program main_channel3d_2nd

module m_BC_and_Halo
  use MPI
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1,hi1,OutFlowInfoX,OutFlowInfoY
  implicit none
  private
  type(HaloInfo)::hi_uxPrSrc,hi_uyPrSrc,hi_uzPrSrc  ! halo info type 1
    
  public:: Init_BC_and_Halo, clcOutFlowVelocity, correctOutFlowFaceVelocity
  public:: SetBC_and_UpdateHalo, SetBC_and_UpdateHaloForPrSrc, SetBC_and_UpdateHalo_pr
contains

  !******************************************************************
  ! Init_BC_and_Halo
  !******************************************************************
  subroutine Init_BC_and_Halo()
    implicit none

    hi1%pencil = y_pencil
    hi1%xmh=1;  hi1%xph=1
    hi1%zmh=1;  hi1%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi1%ymh=1;  hi1%yph=1
    else
      hi1%ymh=0;  hi1%yph=0
    endif
                
    hi_uxPrSrc%pencil = y_pencil
    hi_uxPrSrc%xmh=0;  hi_uxPrSrc%xph=1
    hi_uxPrSrc%ymh=0;  hi_uxPrSrc%yph=0
    hi_uxPrSrc%zmh=0;  hi_uxPrSrc%zph=0

    hi_uyPrSrc%pencil = y_pencil
    hi_uyPrSrc%xmh=0;  hi_uyPrSrc%xph=0;
    hi_uyPrSrc%zmh=0;  hi_uyPrSrc%zph=0;
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_uyPrSrc%ymh=0;  hi_uyPrSrc%yph=1
    else
      hi_uyPrSrc%ymh=0;  hi_uyPrSrc%yph=0
    endif

    hi_uzPrSrc%pencil = y_pencil
    hi_uzPrSrc%xmh=0;  hi_uzPrSrc%xph=0
    hi_uzPrSrc%ymh=0;  hi_uzPrSrc%yph=0
    hi_uzPrSrc%zmh=0;  hi_uzPrSrc%zph=1
  end subroutine Init_BC_and_Halo
  
  !******************************************************************
  ! clcOutFlowVelocity
  !******************************************************************   
  subroutine clcOutFlowVelocity(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,im,jm,ierror
    real(RK)::ConVectVel,rTemp,convEd,rCoeC,rCoeP
    
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyp; jm=jc-1
      ConVectVel=0.0_RK
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ConVectVel=ConVectVel+uy(ic,jc,kc)
        enddo
      enddo
      call MPI_ALLREDUCE(ConVectVel,rTemp,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      ConVectVel=rTemp/(real(nxc,RK)*real(nzc,RK))
      !print*, 'ConVectVel_y=',ConVectVel
      rCoeC=rdyc(jc); rCoeP=rdyp(jm)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          convEd=-ConVectVel*(ux(ic,jc,kc)-ux(ic,jm,kc))*rCoeC
          OutFlowInfoY(4,ic,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoY(1,ic,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc)+OutFlowInfoY(4,ic,kc)
          OutFlowInfoY(1,ic,kc)=convEd

          convEd=-ConVectVel*(uy(ic,jc,kc)-uy(ic,jm,kc))*rCoeP
          OutFlowInfoY(5,ic,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoY(2,ic,kc)
          uy(ic,jc,kc)=uy(ic,jc,kc)+OutFlowInfoY(5,ic,kc)
          OutFlowInfoY(2,ic,kc)=convEd

          convEd=-ConVectVel*(uz(ic,jc,kc)-uz(ic,jm,kc))*rCoeC
          OutFlowInfoY(6,ic,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoY(3,ic,kc)
          uz(ic,jc,kc)=uz(ic,jc,kc)+OutFlowInfoY(6,ic,kc)
          OutFlowInfoY(3,ic,kc)=convEd
        enddo
      enddo    
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then  
      ic=nxp; im=ic-1
      ConVectVel=0.0_RK
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          ConVectVel=ConVectVel+ux(ic,jc,kc)*dyp(jc)
        enddo
      enddo
      call MPI_ALLREDUCE(ConVectVel,rTemp,1,real_type,MPI_SUM,DECOMP_2D_COMM_ROW,ierror)
      ConVectVel=rTemp/(yly*real(nzc,RK))
      !print*, 'ConVectVel_x=',ConVectVel
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          convEd=-ConVectVel*(ux(ic,jc,kc)-ux(im,jc,kc))*rdx
          OutFlowInfoX(4,jc,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoX(1,jc,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc)+OutFlowInfoX(4,jc,kc)
          OutFlowInfoX(1,jc,kc)=convEd        
        
          convEd=-ConVectVel*(uy(ic,jc,kc)-uy(im,jc,kc))*rdx
          OutFlowInfoX(5,jc,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoX(2,jc,kc)
          uy(ic,jc,kc)=uy(ic,jc,kc)+OutFlowInfoX(5,jc,kc)
          OutFlowInfoX(2,jc,kc)=convEd
          
          convEd=-ConVectVel*(uz(ic,jc,kc)-uz(im,jc,kc))*rdx
          OutFlowInfoX(6,jc,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoX(3,jc,kc)
          uz(ic,jc,kc)=uz(ic,jc,kc)+OutFlowInfoX(6,jc,kc)
          OutFlowInfoX(3,jc,kc)=convEd           
        enddo
      enddo
    endif
  end subroutine clcOutFlowVelocity

  !******************************************************************
  ! correctOutFlowFaceVelocity
  !******************************************************************   
  subroutine correctOutFlowFaceVelocity(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,ierror
    real(RK)::CorrectVel,rTemp,OutFlowArea
    
    if(BcOption(xp_dir)/=BC_OutFlow .and. BcOption(yp_dir)/=BC_OutFlow) return
    CorrectVel=0.0_RK
    
    ! yp-dir
    if(BcOption(yp_dir)<0) then
      jc=nyp; rTemp=0.0_RK
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uy(ic,jc,kc)
        enddo
      enddo
      CorrectVel=CorrectVel+rTemp*dx*dz
    endif
    
    ! ym-dir
    if(BcOption(ym_dir)<0) then
      jc=1; rTemp=0.0_RK
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uy(ic,jc,kc)
        enddo
      enddo
      CorrectVel=CorrectVel-rTemp*dx*dz
    endif
    
    ! xp-dir
    if(myProcNghBC(y_pencil,3)<0) then
      ic=nxp; rTemp=0.0_RK
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          rTemp=rTemp+ux(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel+rTemp*dz      
    endif
    
    ! xm-dir
    if(myProcNghBC(y_pencil,4)<0) then
      ic=1; rTemp=0.0_RK
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          rTemp=rTemp+ux(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel-rTemp*dz       
    endif
    
    ! zp-dir
    if(myProcNghBC(y_pencil,1)<0) then
      kc=nzp; rtemp=0.0_RK
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uz(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel+rTemp*dx      
    endif
    
    ! zm-dir
    if(myProcNghBC(y_pencil,2)<0) then
      kc=1; rtemp=0.0_RK
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uz(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel-rTemp*dx      
    endif
    
    OutFlowArea=0.0_RK
    if(BcOption(xp_dir)==BC_outFlow) OutFlowArea=OutFlowArea+yly*zlz
    if(BcOption(yp_dir)==BC_outFlow) OutFlowArea=OutFlowArea+xlx*zlz
    call MPI_ALLREDUCE(CorrectVel,rTemp,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    CorrectVel=-rTemp/OutFlowArea
    
    if(BcOption(yp_dir)==BC_outFlow) then
      jc=nyp
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic,jc,kc)=uy(ic,jc,kc)+CorrectVel
        enddo
      enddo    
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxp
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          ux(ic,jc,kc)=ux(ic,jc,kc)+CorrectVel
        enddo
      enddo       
    endif
  end subroutine correctOutFlowFaceVelocity
        
  !******************************************************************
  ! SetBC_and_UpdateHalo
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    
    ! locals
    integer::ic,kc

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*2.0_RK -ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*2.0_RK -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*2.0_RK -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*dyp(nyp) +ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*2.0_RK   -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*dyp(nyp) +uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = uxBcValue(ym_dir)*2.0_RK-ux(ic, 1, kc)
          uy(ic, 1, kc) = uyBcValue(ym_dir)
          uy(ic, 0, kc) = uyBcValue(ym_dir)*2.0_RK-uy(ic, 2, kc)
          uz(ic, 0, kc) = uzBcValue(ym_dir)*2.0_RK-uz(ic, 1, kc)
        enddo
      enddo
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = ux(ic, 1, kc)-uxBcValue(ym_dir)*dyp(1)
          uy(ic, 1, kc) = uyBcValue(ym_dir)
          uy(ic, 0, kc) = uyBcValue(ym_dir)*2.0_RK-uy(ic, 2, kc)
          uz(ic, 0, kc) = uz(ic, 1, kc)-uzBcValue(ym_dir)*dyp(1)
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NoSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*2.0_RK -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uyBcValue(xp_dir)*2.0_RK -uy(nxc, 0:nyp,y1start(3):y1end(3))
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uzBcValue(xp_dir)*2.0_RK -uz(nxc, 0:nyp,y1start(3):y1end(3))
    CASE(BC_FreeSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*2.0_RK -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uy(nxc, 0:nyp,y1start(3):y1end(3)) +uyBcValue(xp_dir)*dx
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uz(nxc, 0:nyp,y1start(3):y1end(3)) +uzBcValue(xp_dir)*dx        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NoSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*2.0_RK -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uyBcValue(xm_dir)*2.0_RK -uy(1,0:nyp,y1start(3):y1end(3))
      uz(0,0:nyp,y1start(3):y1end(3)) = uzBcValue(xm_dir)*2.0_RK -uz(1,0:nyp,y1start(3):y1end(3))       
    CASE(BC_FreeSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*2.0_RK -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uy(1,0:nyp,y1start(3):y1end(3)) -uyBcValue(xm_dir)*dx
      uz(0,0:nyp,y1start(3):y1end(3)) = uz(1,0:nyp,y1start(3):y1end(3)) -uzBcValue(xm_dir)*dx         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NoSlip) 
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = uxBcValue(zp_dir)*2.0_RK -ux(y1start(1):y1end(1), 0:nyp, nzc) 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uyBcValue(zp_dir)*2.0_RK -uy(y1start(1):y1end(1), 0:nyp, nzc)
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, nzc)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = ux(y1start(1):y1end(1), 0:nyp, nzc) +uxBcValue(zp_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uy(y1start(1):y1end(1), 0:nyp, nzc) +uyBcValue(zp_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, nzc)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NoSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = uxBcValue(zm_dir)*2.0_RK -ux(y1start(1):y1end(1), 0:nyp, 1)
      uy(y1start(1):y1end(1), 0:nyp, 0) = uyBcValue(zm_dir)*2.0_RK -uy(y1start(1):y1end(1), 0:nyp, 1)
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, 2)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = ux(y1start(1):y1end(1), 0:nyp, 1) -uxBcValue(zm_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, 0) = uy(y1start(1):y1end(1), 0:nyp, 1) -uyBcValue(zm_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, 2)     
    END SELECT     

    ! update halo
    call update_halo(ux, mb1, hi1)
    call update_halo(uy, mb1, hi1)
    call update_halo(uz, mb1, hi1)
  end subroutine SetBC_and_UpdateHalo

  !******************************************************************
  ! SetBC_and_UpdateHaloForPrSrc
  !******************************************************************   
  subroutine SetBC_and_UpdateHaloForPrSrc(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    
    ! locals
    integer::ic,kc

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NoSlip,BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic,nyp,kc) = uyBcValue(yp_dir)
        enddo
      enddo
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NoSlip,BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic, 1, kc) = uyBcValue(ym_dir)
        enddo
      enddo
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil,3))
    CASE(BC_NoSlip,BC_FreeSlip)
      ux(nxp,0:nyp,y1start(3):y1end(3))=  uxBcValue(xp_dir)    
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil,4))
    CASE(BC_NoSlip,BC_FreeSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)        
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil,1))
    CASE(BC_NoSlip,BC_FreeSlip)
      uz(y1start(1):y1end(1),0:nyp,nzp) = uzBcValue(zp_dir)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil,2))
    CASE(BC_NoSlip,BC_FreeSlip)
      uz(y1start(1):y1end(1),0:nyp,1) = uzBcValue(zm_dir)  
    END SELECT

    ! update halo
    call update_halo(ux, mb1, hi_uxPrSrc)
    call update_halo(uy, mb1, hi_uyPrSrc)
    call update_halo(uz, mb1, hi_uzPrSrc)
  end subroutine SetBC_and_UpdateHaloForPrSrc
  
  !******************************************************************
  ! SetBC_and_UpdateHalo_pr
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_pr(pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure

    ! locals
    integer::ic,kc

    ! yp-dir
    if(BcOption(yp_dir)<0) then
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          pressure(ic, nyp, kc) = pressure(ic, nyc, kc)
        enddo
      enddo
    endif         
        
    ! ym-dir
    if(BcOption(ym_dir)<0) then
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          pressure(ic, 0, kc) = pressure(ic, 1, kc)
        enddo
      enddo     
    endif

    ! xp-dir
    if(myProcNghBC(y_pencil,3)<0) then
      pressure(nxp,0:nyp,y1start(3):y1end(3))= pressure(nxc,0:nyp,y1start(3):y1end(3))        
    endif
    
    ! xm-dir
    if(myProcNghBC(y_pencil,4)<0) then
      pressure(0,0:nyp,y1start(3):y1end(3))  = pressure(1,0:nyp,y1start(3):y1end(3))       
    endif      
        
    ! zp-dir
    if(myProcNghBC(y_pencil,1)<0) then
      pressure(y1start(1):y1end(1),0:nyp,nzp)= pressure(y1start(1):y1end(1),0:nyp,nzc)         
    endif        
    
    ! zm-dir
    if(myProcNghBC(y_pencil,2)<0) then
      pressure(y1start(1):y1end(1),0:nyp,0)  = pressure(y1start(1):y1end(1),0:nyp,1)       
    endif
    
    call update_halo(pressure, mb1, hi1)
  end subroutine SetBC_and_UpdateHalo_pr
end module m_BC_and_Halo
module m_ChannelSystem
  use MPI
  use m_Timer
  use m_Tools
  use m_Typedef
  use m_LogInfo
  use m_Poisson
  use m_TScheme
  use m_FlowCase
  use m_Decomp2d
  use m_IOAndVisu
  use m_Variables
  use m_DumpPlane
  use m_Parameters
  use m_BC_and_Halo
  use m_MeshAndMetries
  implicit none
  private
    
  !// timers
  type(timer):: total_timer

  public:: ChannelInitialize, ChannelIterate
contains

  !******************************************************************
  ! ChannelInitialize
  !******************************************************************
  subroutine ChannelInitialize(ChannelPrm)
    implicit none
    character(*),intent(in)::ChannelPrm
    
    ! locals
    integer::ierror
    character(256)::chStr
    
    !// Initializing main log info
    write(chStr,"(A)") 'mkdir -p '//ResultsDir//' '//RestartDir//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MainLog%InitLog(ResultsDir,RunName,LF_file_lvl,LF_cmdw_lvl)
    if(nrank==0) call MainLog%CreateFile(RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MainLog%OpenFile()    
    if(nrank==0) call DumpReadedParam()
    
    call InitMeshAndMetries(ChannelPrm)
    call InitVisu(ChannelPrm)
    call Initialize_DumpPlane()

    call AllocateVariables()
    call InitPoissonSolver()
    call Init_BC_and_Halo()
    call InitTimeScheme(ChannelPrm)

    if(.not. RestartFlag) then
      assoDevia: associate(Deviation=>RealArr1)
      call InitVelocity(ux,uy,uz,Deviation,ChannelPrm)
      end associate assoDevia
      if(nrank==0) call MainLog%OutInfo("Initializing all the needed variables into ChannelSystem ...", 1 )
    else
      call read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
      if(nrank==0) call MainLog%OutInfo("Reading all the needed variables into ChannelSystem ...", 1 )
    endif
    call InitStatVar(ChannelPrm) 
    
    SimTime=0.0_RK; dt = dtMax
    call SetBC_and_UpdateHalo(ux,uy,uz)
    call SetBC_and_UpdateHalo_pr( pressure )
    RealArr1=0.0_RK; RealArr2=0.0_RK
       
    ! Timers
    call total_timer%reset()
    call dump_visu(ifirst-1,ux,uy,uz,pressure,RealArr1)
   end subroutine  ChannelInitialize

  !******************************************************************
  ! ChannelIterate
  !******************************************************************
  subroutine ChannelIterate()
    implicit none

    !locals
    integer:: ns
    real(RK)::uddxmax,cflmp,divmax1,divmax2,uxm,vmaxabs(3)

    call total_timer%start()
    call CalcMaxCFL(ux,uy,uz,uddxmax)
    if( icfl==1 ) then
      dt = CFLc/uddxmax
      dt = min(dt, dtMax)
    else
      dt = dtMax  
    endif
    cflmp=uddxmax*dt

    do ns=1, iadvance
      ! step0: Update the Projection Method coefficients.
      call PMcoeUpdate(ns)
      
      ! step1: Calculate the right hand side of the three velocity equations.
      asso_RHS123: associate( RhsX=>RealArr1, RhsY=>RealArr2,  RhsZ=>RealHalo)
      call clcRhsX(ux,uy,uz,RhsX,HistXOld,pressure)
      call clcRhsY(ux,uy,uz,RhsY,HistYOld,pressure)
      call clcRhsZ(ux,uy,uz,RhsZ,HistZOld,pressure)

      ! step2: Calculate the Uhat
      call clcOutFlowVelocity(ux,uy,uz)
      call clcU1Hat(ux,RhsX)
      call clcU2Hat(uy,RhsY)
      call clcU3Hat(uz,RhsZ)
      end associate asso_RHS123

      ! step3: Calculate the source term of the PPE 
      call SetBC_and_UpdateHaloForPrSrc( ux,uy,uz)
      call correctOutFlowFaceVelocity(ux,uy,uz)
      asso_Pr: associate(prsrc =>RealArr1, prphi =>RealArr2, prphiHalo =>RealHalo  )
      call clcPrSrc(ux,uy,uz,prsrc,pressure,divmax1)
      call clcPPE(prsrc,prphiHalo)
      call SetBC_and_UpdateHalo_pr(prphiHalo)
            
      ! step4: Update the velocity field to get the final real velocity.
      call FluidVelUpdate(prphiHalo,ux,uy,uz)
            
      ! step5: Update the real pressure field to get the final pressure.
      call PressureUpdate(pressure, prphiHalo)
      end associate asso_Pr
      call SetBC_and_UpdateHalo( ux,uy,uz)
      call SetBC_and_UpdateHalo_pr( pressure )
    enddo
    if(mod(itime,ivstats)==0) then
      call clcStat(ux,uy,uz,pressure,RealArr1,RealArr2)
      call dump_plane(itime,ux,uy,uz,pressure)
    endif
    if(mod(itime,SaveVisu)== 0)  call dump_visu(itime,ux,uy,uz,pressure,RealArr1)
    if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
      call Write_Restart(itime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
      call Delete_Prev_Restart(itime)
    endif
    call total_timer%finish()

    ! command window and log file output
    IF((itime==ifirst .or. mod(itime, Cmd_LFile_Freq)==0) ) THEN
      call CheckDivergence(ux,uy,uz, divmax2)
      if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelIterate","too big div: "//trim(num2str(divmax2)))
      vmaxabs = CalcVmax(ux,uy,uz)
      if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelIterate","too big velocity: "//trim(num2str(vmaxabs(1)))//", "//trim(num2str(vmaxabs(2)))//", "//trim(num2str(vmaxabs(3))) )

      uxm = CalcUxAver(ux)
      if(nrank==0) then
        call MainLog%OutInfo("Channel3d performed "//trim(num2str(itime))//" iterations up to here!",1)
        call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(total_timer%tot_time))//", "// &
        trim(num2str(total_timer%last_time ))//", "//trim(num2str(total_timer%average())),2)
        call MainLog%OutInfo("SimTime | dt | CFL : "//trim(num2str(SimTime))//' | '//trim(num2str(dt))//' | '//trim(num2str(cflmp)),3)
        call MainLog%OutInfo("Max Abs Div: "//trim(num2str(divmax1))//" | "//trim(num2str(divmax2)) ,3)
        call MainLog%OutInfo("Max Abs Vel: "//trim(num2str(vmaxabs(1)))//" | "//trim(num2str(vmaxabs(2)))//" | "//trim(num2str(vmaxabs(3))), 3)
        call MainLog%OutInfo("Mean Velocity in streamwise: "//trim(num2str(uxm)),3)
      endif
    ENDIF

  end subroutine ChannelIterate

end module m_ChannelSystem
module m_DumpPlane
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_Variables,only: mb1
  implicit none
  private
#define RKP_Dump 4

  logical,parameter::DumpPlaneFlag=.false.
  
  integer,parameter::nyPlane=3
  integer,dimension(nyPlane),parameter::iyPlane=[1, 2, 3]
  
  integer::nTimePlane,iTimePlane
  real(RKP_Dump),dimension(:,:,:,:),allocatable::uxcPlane_y, uycPlane_y, uypPlane_y, uzcPlane_y, prcPlane_y
  
  public::Initialize_DumpPlane, dump_plane
contains

  !******************************************************************
  ! Initialize_DumpPlane
  !******************************************************************
  subroutine Initialize_DumpPlane()
    implicit none
    
    ! locals
    integer::ierrTmp,ierror=0
    
    if(.not. DumpPlaneFlag) return
    if(mod(BackupFreq, ivstats)/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","nTimePlane WRONG")
    endif
    nTimePlane=BackupFreq/ivstats; iTimePlane=0
    allocate(uxcPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(uycPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)    
    allocate(uypPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(uzcPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)    
    allocate(prcPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)          
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","Allocation failed")
  end subroutine Initialize_DumpPlane
  
  !******************************************************************
  ! dump_plane
  !******************************************************************
  subroutine dump_plane(ntime,ux,uy,uz,pr)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr
    
    ! locals
    character(len=128)::filename
    integer,dimension(3)::sizes,subsizes,starts
    integer::ic,jc,kc,jplane1,jplane2,data_type,ierror,newtype,nUnit
        
    if(.not. DumpPlaneFlag) return
    iTimePlane=iTimePlane+1
    
    do kc=y1start(3),y1end(3)
      do jc=1,nyPlane
        jplane1=iyPlane(jc)
        jplane2=jplane1+1
        do ic=y1start(1),y1end(1)
          uxcPlane_y(ic,kc,iTimePlane,jc)=real(ux(ic,jplane1,kc), RKP_Dump)
          uycPlane_y(ic,kc,iTimePlane,jc)=real(uy(ic,jplane1,kc), RKP_Dump)
          uzcPlane_y(ic,kc,iTimePlane,jc)=real(uz(ic,jplane1,kc), RKP_Dump)
          prcPlane_y(ic,kc,iTimePlane,jc)=real(pr(ic,jplane1,kc), RKP_Dump)          
          uypPlane_y(ic,kc,iTimePlane,jc)=real(uy(ic,jplane2,kc), RKP_Dump)
        enddo
      enddo
    enddo
    if(iTimePlane/=nTimePlane) return

    ! Write plane info ========================
    if((RKP_Dump-4)==0) then
      data_type= MPI_REAL
    else    
      data_type= MPI_DOUBLE_PRECISION    
    endif
    sizes(1)= nxc
    sizes(2)= nzc
    sizes(3)= nTimePlane
    subsizes(1)= y1size(1)
    subsizes(2)= y1size(3)
    subsizes(3)= nTimePlane
    starts(1)= y1start(1)-1
    starts(2)= y1start(3)-1
    starts(3)= 0
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
        
    ! uxcPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uxcPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uxcPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo
    
    ! uycPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uycPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uycPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo    

    ! uypPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uypPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uypPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo  

    ! uzcPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uzcPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uzcPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo  

    ! prcPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'prcPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,prcPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo  
    
    call MPI_TYPE_FREE(newtype,ierror)    
    iTimePlane=0
  end subroutine dump_plane
end module m_DumpPlane

#undef RKP_Dump
module m_FlowCase
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use m_FlowType_Channel
  use m_FlowType_TGVortex
  use m_FlowType_AddedNew
  use m_Variables,only: mb1
  implicit none
  private
  
  public:: InitVelocity, InitStatVar, clcStat
contains

  !******************************************************************
  ! SetBC_and_UpdateHalo
  !******************************************************************
  subroutine InitVelocity(ux,uy,uz,Deviation,ChannelPrm)
    implicit none
    character(*),intent(in)::ChannelPrm
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
      
    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitVelocity_CH(ux,uy,uz,Deviation,ChannelPrm)
    case(FT_TG) ! Taylor-Green vortex
      call InitVelocity_TG(ux,uy,uz,Deviation)          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
      call InitVelocity_AN(ux,uy,uz,Deviation)          
    end select    
  end subroutine InitVelocity

  !******************************************************************
  ! InitStatVar
  !******************************************************************
  subroutine InitStatVar(chFile)
    implicit none
    character(*),intent(in)::chFile

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitStatVar_CH(chFile)
    case(FT_TG) ! Taylor-Green vortex
      call InitStatVar_TG()          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
          
    end select

  end subroutine InitStatVar

  !******************************************************************
  ! clcStat
  !****************************************************************** 
  subroutine clcStat(ux,uy,uz,pressure,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp1,ArrTemp2

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call clcStat_CH(ux,uy,uz,pressure,ArrTemp1,ArrTemp2)
    case(FT_TG) ! Taylor-Green vortex
      call clcStat_TG(ux,uy,uz,pressure)          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
          
    end select
  end subroutine clcStat 
    
end module m_FlowCase
module m_FlowType_AddedNew
  use MPI
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
  use m_Tools,only:CalcUxAver
  implicit none
  private    

  public:: InitVelocity_AN, Update_uy_ym_AN
  public:: InitStatVar_AN,  clcStat_AN
contains

  !******************************************************************
  ! InitVelocity_AN
  !******************************************************************
  subroutine InitVelocity_AN(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  

    ! locals
    integer :: ic,jc,kc
    real(RK):: VelRef,Ratiot,yct
    
    VelRef= uxBcValue(ym_dir)
    Ratiot=(uxBcValue(yp_dir)-uxBcValue(ym_dir))/yly
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        yct=yc(jc)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)=  VelRef+Ratiot*yct
        enddo
      enddo
    enddo
    uy=0.0_RK
    uz=0.0_RK
  end subroutine InitVelocity_AN

  !******************************************************************
  ! Update_uy_ym_AN
  !******************************************************************   
  subroutine Update_uy_ym_AN(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(inout):: uy_ym,duy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = uyBcValue(ym_dir)
     
    duy_ym = 0.0_RK!uy_ym - duy_ym
    
  end subroutine Update_uy_ym_AN

  !******************************************************************
  ! InitStatVar_AN
  !******************************************************************
  subroutine InitStatVar_AN()
    implicit none


  end subroutine InitStatVar_AN

  !******************************************************************
  ! clcStat_CH
  !******************************************************************
  subroutine clcStat_AN(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   

  end subroutine clcStat_AN
    
end module m_FlowType_AddedNew
module m_FlowType_Channel
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use iso_c_binding
  use m_MeshAndMetries
  use m_Variables,only:mb1
  use m_Tools,only:CalcUxAver
#ifdef CFDLPT_TwoWay
  use m_Variables,only:FpForce_x,FpForce_y,FpForce_z
#endif
#ifdef CFDACM
  use m_Variables,only:FluidIndicator
#endif
  implicit none
  private
  include "myfftw3.f03"
#define HighOrderGradStat
#define SAVE_SINGLE_Spec2D

  ! statistics variabls
  integer:: nfstime
  real(RK):: PrGradsum
  real(RK),allocatable,dimension(:,:):: SumStat
#ifdef HighOrderGradStat 
  real(RK),allocatable,dimension(:,:)::SumGrad
#endif

  ! SpectraOptions
  integer:: ivSpec,jForLCS(2)
  logical:: clcSpectra1D,clcSpectra2D
  
  ! Spectra variables
  type(C_PTR)::fft_plan_x,fft_plan_z
  integer:: nxh,nxhp,nzh,nzhp,nSpectime
  type(decomp_info),allocatable::decomp_xhzf,decomp_xhzh
  real(RK),allocatable,dimension(:,:,:,:)::EnergySpec2D
  real(RK),allocatable,dimension(:,:,:)::EnergySpecX,EnergySpecZ
    
  public:: InitVelocity_CH,InitStatVar_CH,clcStat_CH

#define iSpec1DUU  1
#define iSpec1DVV  2
#define iSpec1DWW  3
#define iSpec1DPP  4
#define iSpec1DUV  5
#define iSpec1DUV2 6
#define iLCSR1DUU  7
#define iLCSI1DUU  8
#define iLCSR1DUU2 9
#define iLCSR1DVV  10
#define iLCSI1DVV  11
#define iLCSR1DVV2 12
#define iLCSR1DWW  13
#define iLCSI1DWW  14
#define iLCSR1DWW2 15
#define iLCSR1DPP  16
#define iLCSI1DPP  17
#define iLCSR1DPP2 18
#define iSpec1DCC  19
#define iSpec1DUC  20
#define iImag1DUC  21
#define iSpec1DVC  22
#define iImag1DVC  23
#define iLCSR1DCC  24
#define iLCSI1DCC  25
#define iLCSR1DCC2 26

#define iSpec2DUU  1
#define iSpec2DVV  2
#define iSpec2DWW  3
#define iSpec2DPP  4
#define iSpec2DUV  5
#define iLCSR2DUU  6
#define iLCSR2DVV  7
#define iLCSR2DWW  8
#define iLCSR2DPP  9
#define iSpec2DCC  10
#define iLCSR2DCC  11

#ifdef CFDLPT_TwoWay
#define NCHASTAT 45
#elif CFDACM
#define NCHASTAT 51
#else
#define NCHASTAT 35
#endif

#define NEnergySpec1D 18
#define NEnergySpec2D 9

#if defined(HighOrderGradStat)
#define nGradStat     98
#endif

contains
#include "my_FFTW_inc.f90"
#define EnergySpectra_staggered_2nd
#include "EnergySpectra_calcu_fun_inc.f90"
#undef  EnergySpectra_staggered_2nd

  !******************************************************************
  ! InitVelocity_CH
  !******************************************************************
  subroutine InitVelocity_CH(ux,uy,uz,Deviation,chFile)
    implicit none
    character(*),intent(in)::chFile
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    integer::ierror,ic,jc,kc,m1,m2,nUnit,iTV(8)
    real(RK)::uBulkTemp,ybulk1,ybulk2,twopi
    real(RK)::xplus,yplus,zplus,yct,ybar,xpt,zpt,ratiot
    real(RK)::retau_guass,utau_guass,height,rem,wx,wz,xlxPlus,zlzPlus
    NAMELIST/uBulk_Param/ybulk1,ybulk2

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitVelocity_CH","Cannot open file: "//trim(chFile))
    endif
    read(nUnit,nml=uBulk_Param)
    close(nUnit,IOSTAT=ierror)
    height=ybulk2-ybulk1
    uBulkTemp=uBulk*yly/height
    if(FlowType==FT_CH) then
      height=0.5_RK*height
    endif
    rem=uBulkTemp*height/xnu
    ux=0.0_RK; uy=0.0_RK; uz=0.0_RK
    if(abs(uBulk)<1.0E-12)return
    
    retau_guass = 0.1538_RK*rem**0.887741_RK
    utau_guass  = retau_guass*xnu/height
    if(nrank==0)print*,'************** retau_gauss=',retau_guass
    if(nrank==0)print*,'************** utau_gauss= ',utau_guass

    call date_and_time(values=iTV); !iTV=0
    call random_seed(size= ic)
    call random_seed(put = iTV(7)*iTV(8)+[(jc,jc=1,ic)])
    call random_number(Deviation)
    Deviation= 0.2_RK* Deviation + 0.9_RK ! [0.8, 1.2]

    ! modulation of the random noise + initial velocity profile
    twopi=2.0_RK*PI
    wx=twopi/500.0_RK; wz=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   zlzPlus=zlz*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(zlzPlus*wz/twopi)+1; wz=real(m2,RK)*twopi/zlzPlus
    do jc=y1start(2),y1end(2)
      yct = height-abs(height-(yc(jc)-ybulk1))
      if(yc(jc)<ybulk1 .or. yc(jc)>=ybulk2)yct=0.0_RK
      ybar= yct/height; yplus=utau_guass*yct/xnu
      do kc=y1start(3),y1end(3)
        zpt  =real(kc-1,kind=RK)*dz+dz*0.5_RK
        zplus=utau_guass*zpt/xnu
        do ic=y1start(1),y1end(1)
          xpt  =real(ic-1,kind=RK)*dx+dx*0.5_RK
          xplus=utau_guass*xpt/xnu
          !ux(ic,jc,kc) = 0.0052_RK*uBulkTemp*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(ic,jc,kc) ! original expression
          !uz(ic,jc,kc) = 0.0050_RK*uBulkTemp*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(ic,jc,kc) ! original expression
          ux(ic,jc,kc) = uBulkTemp*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(ic,jc,kc)
          uz(ic,jc,kc) = uBulkTemp*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(ic,jc,kc)
          ux(ic,jc,kc) = ux(ic,jc,kc)+ 3.0_RK*uBulkTemp*(ybar-0.5_RK*ybar*ybar)
        enddo
      enddo
    enddo
    ratiot=ubulk/CalcUxAver(ux)
    call MPI_Bcast(ratiot,1,real_type,0,MPI_COMM_WORLD,ierror)
    ux= ux*ratiot
    Deviation=0.0_RK
  end subroutine InitVelocity_CH

  !******************************************************************
  ! InitStatVar_CH
  !******************************************************************
  subroutine InitStatVar_CH(chFile)
    implicit none
    character(*),intent(in)::chFile

    ! locals
    real(RK)::ybulk1,ybulk2
    character(len=128)::filename
    NAMELIST/uBulk_Param/ybulk1,ybulk2
    integer::ierror,nUnit,plan_type,nySpec2D
    real(RK),dimension(:),allocatable::Vec1,Vec2
    NAMELIST/SpectraOptions/clcSpectra1D,clcSpectra2D,ivSpec,jForLCS

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Cannot open file: "//trim(chFile))
    endif
    read(nUnit,nml=uBulk_Param)
    rewind(nUnit)
    read(nUnit, nml=SpectraOptions)
    close(nUnit,IOSTAT=ierror)
        
    if(nrank==0) then
      write(MainLog%nUnit, nml=uBulk_Param)
      write(MainLog%nUnit, nml=SpectraOptions)
      if(mod(saveStat,ivstats)/=0 )  call MainLog%CheckForError(ErrT_Abort,"InitStatVar","ivstats wrong !!!")
      if(clcSpectra1D .and. (mod(saveStat,ivSpec)/=0 .or. mod(ivSpec,ivstats)/=0 )) then
        call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","ivSpec wrong !!!")
      endif
      if(IsUxConst)then
        write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
        open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
        if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","Cannot open file: "//trim(filename))
        close(nUnit,IOSTAT=ierror)
      endif
    endif
    allocate(SumStat(NCHASTAT,nyp),Stat=ierror)
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Allocation failed")  
    nfstime=0; SumStat=0.0_RK; PrGradsum=0.0_RK
#ifdef HighOrderGradStat
    allocate(SumGrad(nGradStat,nyp),Stat=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar: ","Allocation failed 2")  
    SumGrad=0.0_RK
#endif

    if(FFTW_plan_type == 1) then
      plan_type=FFTW_PATIENT
    else
      plan_type=FFTW_ESTIMATE
    endif    
    nxh=nxc/2; nxhp=nxh+1
    nzh=nzc/2; nzhp=nzh+1
    
    ! FFT_plan_x
    allocate(Vec1(x1size(1)),Vec2(x1size(1)))
    fft_plan_x = fftw_plan_r2r_1d(x1size(1),Vec1,Vec2,FFTW_R2HC,plan_type)         
    deallocate(Vec1,Vec2)

    ! FFT_plan_z
    allocate(Vec1(z1size(3)),Vec2(z1size(3)))
    fft_plan_z = fftw_plan_r2r_1d(z1size(3),Vec1,Vec2,FFTW_R2HC,plan_type)         
    deallocate(Vec1,Vec2)
    
    IF(clcSpectra1D) THEN
      allocate(EnergySpecX(nxhp,x1size(2),NEnergySpec1D));EnergySpecX=0.0_RK
      allocate(EnergySpecZ(nzhp,z1size(2),NEnergySpec1D));EnergySpecZ=0.0_RK
    ENDIF
    IF(clcSpectra2D) THEN
      allocate(decomp_xhzf,decomp_xhzh)
      Block
        logical,dimension(6)::initializeIn
        initializeIn=.false.; initializeIn(5:6)=.true.
        call decomp_info_init(nxhp,nyc,nzc, decomp_xhzf,initialize=initializeIn)
        initializeIn=.false.; initializeIn(5)=.true.
        call decomp_info_init(nxhp,nyc,nzhp,decomp_xhzh,initialize=initializeIn)
      end Block
 
      if(FlowType==FT_CH) then
        nySpec2D=nyc/2
      else
        nySpec2D=nyc
      endif
      allocate(EnergySpec2D(decomp_xhzh%y2sz(1),nySpec2D,decomp_xhzh%y2sz(3),NEnergySpec2D),Stat=ierror)
      if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar","Allocation failed For Spectra2D")
      EnergySpec2D=0.0_RK
    ENDIF
  end subroutine InitStatVar_CH

  !******************************************************************
  ! clcStat_CH
  !******************************************************************
  subroutine clcStat_CH(ux,uy,uz,pr,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp1,ArrTemp2
   
    ! locals
    character(len=128)::filename
    integer(kind=8)::disp,disp_inc
    real(RK),allocatable,dimension(:,:,:)::arrx1,arrx2,arrz1,arrz2
    real(RK),allocatable,dimension(:,:,:)::EnergySpecXR,EnergySpecZR
    real(RK)::dudyU,dvdyM,dudzC,dvdzC,dudxC,dvdxU,SumVec(NCHASTAT),rdxh
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,it,jt,kt,ierror,ids,nrankX,nrankZ,nUnit
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(NCHASTAT,nyp)
    real(RK)::uxloc,uyloc,uzloc,prCell,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2
#ifdef HighOrderGradStat
    real(RK)::dpdx,dpdy,dpdz,dpdxx,dpdyy,dpdzz,SumGradR(nGradStat,nyp),SumVec2(nGradStat)
#endif
#ifdef CFDLPT_TwoWay
    real(RK)::uyCellm
#endif

    rdxh=0.5_RK*rdx
    inxz=1.0_RK/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1;
      InterpY1= 0.5_RK*YinterpCoe(jc); InterpY2=0.5_RK-InterpY1
      cac=rdyc(jc);cacU=rdyc(jp); caj=rdyp(jc); SumVec=0.0_RK
#ifdef HighOrderGradStat
      SumVec2=0.0_RK
#endif
      do kc=y1start(3),y1end(3)
        km=kc-1;kp=kc+1
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1
          
          uxloc = ux(ic,jc,kc)
          uyloc = uy(ic,jc,kc)
          uzloc = uz(ic,jc,kc)
          prCell= pr(ic,jc,kc)
          uxCell= 0.5_RK*(ux(ic,jc,kc)+ux(ip,jc,kc))
          uyCell= 0.5_RK*(uy(ic,jc,kc)+uy(ic,jp,kc))
          uzCell= 0.5_RK*(uz(ic,jc,kc)+uz(ic,jc,kp))
          prloc2= InterpY1*(pr(im,jm,kc)+pr(ic,jm,kc))+ InterpY2*(pr(im,jc,kc)+pr(ic,jc,kc))
 
          dudx= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz= (ux(ic,jc,kc)-ux(ic,jc,km))*rdz
          dvdx= (uy(ic,jc,kc)-uy(im,jc,kc))*rdx
          dvdy= (uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          dvdz= (uy(ic,jc,kc)-uy(ic,jc,km))*rdz
          dwdx= (uz(ic,jc,kc)-uz(im,jc,kc))*rdx
          dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
          dwdz= (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz

          dudxx= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          dudyy= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)
          dudzz= (ux(ic,jc,kp)-2.0_RK*ux(ic,jc,kc)+ux(ic,jc,km))*rdz2
          dvdxx= (uy(ip,jc,kc)-2.0_RK*uy(ic,jc,kc)+uy(im,jc,kc))*rdx2
          dvdyy= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          dvdzz= (uy(ic,jc,kp)-2.0_RK*uy(ic,jc,kc)+uy(ic,jc,km))*rdz2
          dwdxx= (uz(ip,jc,kc)-2.0_RK*uz(ic,jc,kc)+uz(im,jc,kc))*rdx2
          dwdyy= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)
          dwdzz= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          vor_x= dwdy-dvdz
          vor_y= dudz-dwdx
          vor_z= dvdx-dudy

          dudxC= (ux(ip,jc,kc)-ux(im,jc,kc))*rdxh
          dvdxU= (uy(ic,jp,kc)-uy(im,jp,kc))*rdx
          dudyU= (ux(ic,jp,kc)-ux(ic,jc,kc))*cacU
          dvdyM= (uy(im,jp,kc)-uy(im,jc,kc))*caj
          dudzC= (InterpY1*(ux(ic,jm,kp)-ux(ic,jm,km)) +InterpY2*(ux(ic,jc,kp)-ux(ic,jc,km)))*rdz
          dvdzC= (uy(im,jc,kp)-uy(im,jc,km)+ uy(ic,jc,kp)-uy(ic,jc,km))*rdz*0.25_RK

          ids=1
          SumVec(ids)=SumVec(ids)+ uxloc;                     ids=ids+1 ! 01
          SumVec(ids)=SumVec(ids)+ uyloc;                     ids=ids+1 ! 02, *yp*
          SumVec(ids)=SumVec(ids)+ uzloc;                     ids=ids+1 ! 03
          SumVec(ids)=SumVec(ids)+ prCell;                    ids=ids+1 ! 04
          SumVec(ids)=SumVec(ids)+ uxloc*uxloc;               ids=ids+1 ! 05
          SumVec(ids)=SumVec(ids)+ uyloc*uyloc;               ids=ids+1 ! 06, *yp*  
          SumVec(ids)=SumVec(ids)+ uzloc*uzloc;               ids=ids+1 ! 07
          SumVec(ids)=SumVec(ids)+ prCell*prCell;             ids=ids+1 ! 08
          SumVec(ids)=SumVec(ids)+ uxCell*uyCell;             ids=ids+1 ! 09 
          SumVec(ids)=SumVec(ids)+ uyCell*uzCell;             ids=ids+1 ! 10
          SumVec(ids)=SumVec(ids)+ uxCell*uzCell;             ids=ids+1 ! 11
          SumVec(ids)=SumVec(ids)+ uxCell*prCell;             ids=ids+1 ! 12
          SumVec(ids)=SumVec(ids)+ uyCell*prCell;             ids=ids+1 ! 13
          SumVec(ids)=SumVec(ids)+ uzCell*prCell;             ids=ids+1 ! 14
          SumVec(ids)=SumVec(ids)+ uxCell*uxCell*uyCell;      ids=ids+1 ! 15
          SumVec(ids)=SumVec(ids)+ uyloc *uyloc *uyloc;       ids=ids+1 ! 16, *yp*
          SumVec(ids)=SumVec(ids)+ uzCell*uzCell*uyCell;      ids=ids+1 ! 17
          SumVec(ids)=SumVec(ids)+ uxCell*uyCell*uyCell;      ids=ids+1 ! 18
          SumVec(ids)=SumVec(ids)+ uxloc*uxloc*uxloc;         ids=ids+1 ! 19 
          SumVec(ids)=SumVec(ids)+ uzloc*uzloc*uzloc;         ids=ids+1 ! 20
          SumVec(ids)=SumVec(ids)+ uxloc*uxloc*uxloc*uxloc;   ids=ids+1 ! 21
          SumVec(ids)=SumVec(ids)+ uyloc*uyloc*uyloc*uyloc;   ids=ids+1 ! 22, *yp*
          SumVec(ids)=SumVec(ids)+ uzloc*uzloc*uzloc*uzloc;   ids=ids+1 ! 23
          SumVec(ids)=SumVec(ids)+ prCell*dudx;               ids=ids+1 ! 24 
          SumVec(ids)=SumVec(ids)+ prCell*dvdy;               ids=ids+1 ! 25 
          SumVec(ids)=SumVec(ids)+ prCell*dwdz;               ids=ids+1 ! 26 
          SumVec(ids)=SumVec(ids)+ prloc2*(dudy+dvdx);        ids=ids+1 ! 27, *yp*
          SumVec(ids)=SumVec(ids)+ uxloc*(dudxx+dudyy+dudzz); ids=ids+1 ! 28 
          SumVec(ids)=SumVec(ids)+ uyloc*(dvdxx+dvdyy+dvdzz); ids=ids+1 ! 29, *yp*
          SumVec(ids)=SumVec(ids)+ uzloc*(dwdxx+dwdyy+dwdzz); ids=ids+1 ! 30
          SumVec(ids)=SumVec(ids)+ dudxC*(dvdxU+dvdx)*0.5_RK  &
                     +(dudyU+dudy)*(dvdyM+dvdy)*0.25_RK;      ids=ids+1 ! 31
          SumVec(ids)=SumVec(ids)+ dudzC*dvdzC;               ids=ids+1 ! 32, *yp*
          SumVec(ids)=SumVec(ids)+ vor_x*vor_x;               ids=ids+1 ! 33, *yp*
          SumVec(ids)=SumVec(ids)+ vor_y*vor_y;               ids=ids+1 ! 34
          SumVec(ids)=SumVec(ids)+ vor_z*vor_z;               ids=ids+1 ! 35, *yp*
#ifdef CFDLPT_TwoWay
          SumVec(ids)=SumVec(ids)+ FpForce_x(ic,jc,kc);       ids=ids+1 ! 36
          SumVec(ids)=SumVec(ids)+ FpForce_y(ic,jc,kc);       ids=ids+1 ! 37, *yp*
          SumVec(ids)=SumVec(ids)+ FpForce_z(ic,jc,kc);       ids=ids+1 ! 38
          SumVec(ids)=SumVec(ids)+ FpForce_x(ic,jc,kc)*FpForce_x(ic,jc,kc); ids=ids+1 ! 39
          SumVec(ids)=SumVec(ids)+ FpForce_y(ic,jc,kc)*FpForce_y(ic,jc,kc); ids=ids+1 ! 40, *yp*
          SumVec(ids)=SumVec(ids)+ FpForce_z(ic,jc,kc)*FpForce_z(ic,jc,kc); ids=ids+1 ! 41          
          SumVec(ids)=SumVec(ids)+ FpForce_x(ic,jc,kc)*uxloc; ids=ids+1 ! 42
          SumVec(ids)=SumVec(ids)+ FpForce_y(ic,jc,kc)*uyloc; ids=ids+1 ! 43, *yp*
          SumVec(ids)=SumVec(ids)+ FpForce_z(ic,jc,kc)*uzloc; ids=ids+1 ! 44
          
          SumVec(ids)=SumVec(ids)+ 0.5_RK*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jp,kc))*uxCell
          uyCellm=(uy(im,jc,kc)+uy(im,jp,kc))*0.5_RK
          SumVec(ids)=SumVec(ids)+FpForce_x(ic,jc,kc)*(uyCellm+uyCell)*0.5_RK; ids=ids+1 ! 45
#endif

#ifdef CFDACM
          SumVec(ids)=SumVec(ids)+ uxCell;                    ids=ids+1 ! 36
          SumVec(ids)=SumVec(ids)+ uyCell;                    ids=ids+1 ! 37
          SumVec(ids)=SumVec(ids)+ uzCell;                    ids=ids+1 ! 38
          SumVec(ids)=SumVec(ids)+ uxCell*uxCell;             ids=ids+1 ! 39
          SumVec(ids)=SumVec(ids)+ uyCell*uyCell;             ids=ids+1 ! 40
          SumVec(ids)=SumVec(ids)+ uzCell*uzCell;             ids=ids+1 ! 41
          if(FluidIndicator(ic,jc,kc)=='P')cycle
          SumVec(ids)=SumVec(ids)+ 1.0_RK;                    ids=ids+1 ! 42          
          SumVec(ids)=SumVec(ids)+ uxCell;                    ids=ids+1 ! 43
          SumVec(ids)=SumVec(ids)+ uyCell;                    ids=ids+1 ! 44
          SumVec(ids)=SumVec(ids)+ uzCell;                    ids=ids+1 ! 45
          SumVec(ids)=SumVec(ids)+ uxCell*uxCell;             ids=ids+1 ! 46
          SumVec(ids)=SumVec(ids)+ uyCell*uyCell;             ids=ids+1 ! 47
          SumVec(ids)=SumVec(ids)+ uzCell*uzCell;             ids=ids+1 ! 48
          SumVec(ids)=SumVec(ids)+ uxCell*uyCell;             ids=ids+1 ! 49
          SumVec(ids)=SumVec(ids)+ prCell;                    ids=ids+1 ! 50
          SumVec(ids)=SumVec(ids)+ prCell*prCell;             ids=ids+1 ! 51
#endif

#ifdef HighOrderGradStat
          dpdx=(pr(ic,jc,kc)-pr(im,jc,kc))*rdx
          dpdy=(pr(ic,jc,kc)-pr(ic,jm,kc))*cac          ! yp ! 
          dpdz=(pr(ic,jc,kc)-pr(ic,jc,km))*rdz
          dpdxx= (pr(ip,jc,kc)-2.0_RK*pr(ic,jc,kc)+pr(im,jc,kc))*rdx2
          dpdyy= ap2c(jc)*pr(ic,jp,kc)+ac2c(jc)*pr(ic,jc,kc)+am2c(jc)*pr(ic,jm,kc)
          dpdzz= (pr(ic,jc,kp)-2.0_RK*pr(ic,jc,kc)+pr(ic,jc,km))*rdz2
          
          ids=1
          SumVec2(ids)=SumVec2(ids)+ prCell*prCell*prCell;        ids=ids+1 ! 01
          SumVec2(ids)=SumVec2(ids)+ prCell*prCell*prCell*prCell; ids=ids+1 ! 02
          
          SumVec2(ids)=SumVec2(ids)+ dudx;                        ids=ids+1 ! 03
          SumVec2(ids)=SumVec2(ids)+ dudx*dudx;                   ids=ids+1 ! 04
          SumVec2(ids)=SumVec2(ids)+ dudx*dudx*dudx;              ids=ids+1 ! 05
          SumVec2(ids)=SumVec2(ids)+ dudx*dudx*dudx*dudx;         ids=ids+1 ! 06
          SumVec2(ids)=SumVec2(ids)+ dudy;                        ids=ids+1 ! 07, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudy*dudy;                   ids=ids+1 ! 08, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudy*dudy*dudy;              ids=ids+1 ! 09, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudy*dudy*dudy*dudy;         ids=ids+1 ! 10, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudz;                        ids=ids+1 ! 11
          SumVec2(ids)=SumVec2(ids)+ dudz*dudz;                   ids=ids+1 ! 12
          SumVec2(ids)=SumVec2(ids)+ dudz*dudz*dudz;              ids=ids+1 ! 13
          SumVec2(ids)=SumVec2(ids)+ dudz*dudz*dudz*dudz;         ids=ids+1 ! 14

          SumVec2(ids)=SumVec2(ids)+ dvdx;                        ids=ids+1 ! 15, *yp*       
          SumVec2(ids)=SumVec2(ids)+ dvdx*dvdx;                   ids=ids+1 ! 16, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdx*dvdx*dvdx;              ids=ids+1 ! 17, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdx*dvdx*dvdx*dvdx;         ids=ids+1 ! 18, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdy;                        ids=ids+1 ! 19                  
          SumVec2(ids)=SumVec2(ids)+ dvdy*dvdy;                   ids=ids+1 ! 20          
          SumVec2(ids)=SumVec2(ids)+ dvdy*dvdy*dvdy;              ids=ids+1 ! 21      
          SumVec2(ids)=SumVec2(ids)+ dvdy*dvdy*dvdy*dvdy;         ids=ids+1 ! 22   
          SumVec2(ids)=SumVec2(ids)+ dvdz;                        ids=ids+1 ! 23, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdz*dvdz;                   ids=ids+1 ! 24, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdz*dvdz*dvdz;              ids=ids+1 ! 25, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdz*dvdz*dvdz*dvdz;         ids=ids+1 ! 26, *yp* 
          
          SumVec2(ids)=SumVec2(ids)+ dwdx;                        ids=ids+1 ! 27
          SumVec2(ids)=SumVec2(ids)+ dwdx*dwdx;                   ids=ids+1 ! 28
          SumVec2(ids)=SumVec2(ids)+ dwdx*dwdx*dwdx;              ids=ids+1 ! 29
          SumVec2(ids)=SumVec2(ids)+ dwdx*dwdx*dwdx*dwdx;         ids=ids+1 ! 30
          SumVec2(ids)=SumVec2(ids)+ dwdy;                        ids=ids+1 ! 31, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdy*dwdy;                   ids=ids+1 ! 32, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdy*dwdy*dwdy;              ids=ids+1 ! 33, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdy*dwdy*dwdy*dwdy;         ids=ids+1 ! 34, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdz;                        ids=ids+1 ! 35
          SumVec2(ids)=SumVec2(ids)+ dwdz*dwdz;                   ids=ids+1 ! 36
          SumVec2(ids)=SumVec2(ids)+ dwdz*dwdz*dwdz;              ids=ids+1 ! 37
          SumVec2(ids)=SumVec2(ids)+ dwdz*dwdz*dwdz*dwdz;         ids=ids+1 ! 38
          
          SumVec2(ids)=SumVec2(ids)+ dpdx;                        ids=ids+1 ! 39
          SumVec2(ids)=SumVec2(ids)+ dpdx*dpdx;                   ids=ids+1 ! 40
          SumVec2(ids)=SumVec2(ids)+ dpdx*dpdx*dpdx;              ids=ids+1 ! 41
          SumVec2(ids)=SumVec2(ids)+ dpdx*dpdx*dpdx*dpdx;         ids=ids+1 ! 42
          SumVec2(ids)=SumVec2(ids)+ dpdy;                        ids=ids+1 ! 43, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdy*dpdy;                   ids=ids+1 ! 44, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdy*dpdy*dpdy;              ids=ids+1 ! 45, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdy*dpdy*dpdy*dpdy;         ids=ids+1 ! 46, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdz;                        ids=ids+1 ! 47
          SumVec2(ids)=SumVec2(ids)+ dpdz*dpdz;                   ids=ids+1 ! 48
          SumVec2(ids)=SumVec2(ids)+ dpdz*dpdz*dpdz;              ids=ids+1 ! 49
          SumVec2(ids)=SumVec2(ids)+ dpdz*dpdz*dpdz*dpdz;         ids=ids+1 ! 50
          
          SumVec2(ids)=SumVec2(ids)+ dudxx;                       ids=ids+1 ! 51
          SumVec2(ids)=SumVec2(ids)+ dudxx*dudxx;                 ids=ids+1 ! 52
          SumVec2(ids)=SumVec2(ids)+ dudxx*dudxx*dudxx;           ids=ids+1 ! 53
          SumVec2(ids)=SumVec2(ids)+ dudxx*dudxx*dudxx*dudxx;     ids=ids+1 ! 54
          SumVec2(ids)=SumVec2(ids)+ dudyy;                       ids=ids+1 ! 55
          SumVec2(ids)=SumVec2(ids)+ dudyy*dudyy;                 ids=ids+1 ! 56
          SumVec2(ids)=SumVec2(ids)+ dudyy*dudyy*dudyy;           ids=ids+1 ! 57
          SumVec2(ids)=SumVec2(ids)+ dudyy*dudyy*dudyy*dudyy;     ids=ids+1 ! 58
          SumVec2(ids)=SumVec2(ids)+ dudzz;                       ids=ids+1 ! 59
          SumVec2(ids)=SumVec2(ids)+ dudzz*dudzz;                 ids=ids+1 ! 60
          SumVec2(ids)=SumVec2(ids)+ dudzz*dudzz*dudzz;           ids=ids+1 ! 61
          SumVec2(ids)=SumVec2(ids)+ dudzz*dudzz*dudzz*dudzz;     ids=ids+1 ! 62
          
          SumVec2(ids)=SumVec2(ids)+ dvdxx;                       ids=ids+1 ! 63, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdxx*dvdxx;                 ids=ids+1 ! 64, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdxx*dvdxx*dvdxx;           ids=ids+1 ! 65, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdxx*dvdxx*dvdxx*dvdxx;     ids=ids+1 ! 66, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy;                       ids=ids+1 ! 67, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy*dvdyy;                 ids=ids+1 ! 68, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy*dvdyy*dvdyy;           ids=ids+1 ! 69, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy*dvdyy*dvdyy*dvdyy;     ids=ids+1 ! 70, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz;                       ids=ids+1 ! 71, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz*dvdzz;                 ids=ids+1 ! 72, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz*dvdzz*dvdzz;           ids=ids+1 ! 73, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz*dvdzz*dvdzz*dvdzz;     ids=ids+1 ! 74, *yp*                                 

          SumVec2(ids)=SumVec2(ids)+ dwdxx;                       ids=ids+1 ! 75
          SumVec2(ids)=SumVec2(ids)+ dwdxx*dwdxx;                 ids=ids+1 ! 76
          SumVec2(ids)=SumVec2(ids)+ dwdxx*dwdxx*dwdxx;           ids=ids+1 ! 77
          SumVec2(ids)=SumVec2(ids)+ dwdxx*dwdxx*dwdxx*dwdxx;     ids=ids+1 ! 78
          SumVec2(ids)=SumVec2(ids)+ dwdyy;                       ids=ids+1 ! 79
          SumVec2(ids)=SumVec2(ids)+ dwdyy*dwdyy;                 ids=ids+1 ! 80
          SumVec2(ids)=SumVec2(ids)+ dwdyy*dwdyy*dwdyy;           ids=ids+1 ! 81
          SumVec2(ids)=SumVec2(ids)+ dwdyy*dwdyy*dwdyy*dwdyy;     ids=ids+1 ! 82
          SumVec2(ids)=SumVec2(ids)+ dwdzz;                       ids=ids+1 ! 83
          SumVec2(ids)=SumVec2(ids)+ dwdzz*dwdzz;                 ids=ids+1 ! 84
          SumVec2(ids)=SumVec2(ids)+ dwdzz*dwdzz*dwdzz;           ids=ids+1 ! 85
          SumVec2(ids)=SumVec2(ids)+ dwdzz*dwdzz*dwdzz*dwdzz;     ids=ids+1 ! 86
          
          SumVec2(ids)=SumVec2(ids)+ dpdxx;                       ids=ids+1 ! 87
          SumVec2(ids)=SumVec2(ids)+ dpdxx*dpdxx;                 ids=ids+1 ! 88
          SumVec2(ids)=SumVec2(ids)+ dpdxx*dpdxx*dpdxx;           ids=ids+1 ! 89
          SumVec2(ids)=SumVec2(ids)+ dpdxx*dpdxx*dpdxx*dpdxx;     ids=ids+1 ! 90
          SumVec2(ids)=SumVec2(ids)+ dpdyy;                       ids=ids+1 ! 91
          SumVec2(ids)=SumVec2(ids)+ dpdyy*dpdyy;                 ids=ids+1 ! 92
          SumVec2(ids)=SumVec2(ids)+ dpdyy*dpdyy*dpdyy;           ids=ids+1 ! 93
          SumVec2(ids)=SumVec2(ids)+ dpdyy*dpdyy*dpdyy*dpdyy;     ids=ids+1 ! 94
          SumVec2(ids)=SumVec2(ids)+ dpdzz;                       ids=ids+1 ! 95
          SumVec2(ids)=SumVec2(ids)+ dpdzz*dpdzz;                 ids=ids+1 ! 96
          SumVec2(ids)=SumVec2(ids)+ dpdzz*dpdzz*dpdzz;           ids=ids+1 ! 97
          SumVec2(ids)=SumVec2(ids)+ dpdzz*dpdzz*dpdzz*dpdzz;     ids=ids+1 ! 98
#endif
        enddo
      enddo
      do kc=1,NCHASTAT
        SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
      enddo
#ifdef HighOrderGradStat
      do kc=1,nGradStat
        SumGrad(kc,jc)=SumGrad(kc,jc)+ SumVec2(kc)*inxz
      enddo
#endif
    ENDDO

    ! nyp only
    jc=nyp; jm=jc-1; cac=rdyc(jc); SumVec=0.0_RK
#ifdef HighOrderGradStat
    SumVec2=0.0_RK
#endif
    InterpY1= 0.5_RK*YinterpCoe(jc); InterpY2=0.5_RK-InterpY1
    do kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do ic=y1start(1),y1end(1)
        im=ic-1;ip=ic+1
        prloc2= InterpY1*(pr(im,jm,kc)+pr(ic,jm,kc))+ InterpY2*(pr(im,jc,kc)+pr(ic,jc,kc))
        dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
        dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
        vor_x=  dwdy
        vor_z= -dudy
        SumVec(27)=SumVec(27)+ prloc2*(dudy+0.0_RK)   ! yp !
        SumVec(33)=SumVec(33)+ vor_x*vor_x            ! yp !
        SumVec(35)=SumVec(35)+ vor_z*vor_z            ! yp !
#if defined(CFDLPT_TwoWay)
        SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)                     ! yp !
        SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*FpForce_y(ic,jc,kc) ! yp !
        SumVec(43)=SumVec(43)+ FpForce_y(ic,jc,kc)*uy(ic,jc,kc)        ! yp !
#endif
#ifdef HighOrderGradStat
        SumVec2( 7)=SumVec2( 7)+ dudy                 ! yp !
        SumVec2( 8)=SumVec2( 8)+ dudy*dudy            ! yp !
        SumVec2( 9)=SumVec2( 9)+ dudy*dudy*dudy       ! yp !
        SumVec2(10)=SumVec2(10)+ dudy*dudy*dudy*dudy  ! yp !
        
        SumVec2(31)=SumVec2(31)+ dwdy                 ! yp !
        SumVec2(32)=SumVec2(32)+ dwdy*dwdy            ! yp !
        SumVec2(33)=SumVec2(33)+ dwdy*dwdy*dwdy       ! yp !
        SumVec2(34)=SumVec2(34)+ dwdy*dwdy*dwdy*dwdy  ! yp !
#endif
      enddo
    enddo
    do kc=1,NCHASTAT
      SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
    enddo
#ifdef HighOrderGradStat
    do kc=1,nGradStat
      SumGrad(kc,jc)=SumGrad(kc,jc)+ SumVec2(kc)*inxz
    enddo
#endif

    ! shear stress and pressure gradient
    PrGradsum   = PrGradsum+ PrGradData(1)
    if(nrank==0 .and. IsUxConst) then
      write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
      open(newunit=nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierror)
      if(ierror/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,2ES24.15)')itime,SimTime,PrGradData(1)
      endif
      close(nUnit,IOSTAT=ierror)
    endif
    nfstime= nfstime + 1
#define EnergySpectra_staggered_2nd
#include "EnergySpectra_staggered_inc.f90"
#undef  EnergySpectra_staggered_2nd
    if(mod(itime,SaveStat)/=0) return

    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,NCHASTAT*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = 1.0_RK/real(nfstime,RK)
      write(filename,"(A,I10.10)") trim(ResultsDir)//'stats',itime
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(nUnit,'(A)')'  '
        if(IsUxConst) then
          write(nUnit,'(A)')'  Constant velocity in x-dir by adding a pressure gradient.'
          write(nUnit,'(A, ES24.15)')'    time averaged pressure gradient is: ',PrGradsum*infstime
        else
          dudy = abs(SumStatR(1,1))*infstime*2.0_RK*rdyc(1)
          if(FlowType==FT_CH) then
            dudyU= abs(SumStatR(1,nyc))*infstime*2.0_RK*rdyc(nyp)
            dudy = dudy+dudyU
          endif
          write(nUnit,'(A, ES24.15)')'  Variable velocity in x-dir while adding a constant body force.',xnu*dudy/yly
        endif
        write(nUnit,'(A)')'  '
        
        Block 
        character(len=128)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',NCHASTAT,'ES24.15)'
        do jc=1,nyp
          write(nUnit,FormatStr)SumStatR(1:NCHASTAT,jc)*infstime
        enddo
        End block
      ENDIF
      close(nUnit,IOSTAT=ierror)
    endif
    
#ifdef HighOrderGradStat
    call MPI_REDUCE(SumGrad,SumGradR,nGradStat*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = 1.0_RK/real(nfstime,RK)
      write(filename,"(A,I10.10)") trim(ResultsDir)//'vGrad',itime
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat","Cannot open file: "//trim(filename))
      ELSE
        Block 
        character(len=128)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',nGradStat,'ES24.15)'
        do jc=1,nyp
          write(nUnit,FormatStr)SumGradR(1:nGradStat,jc)*infstime
        enddo
        End Block
      ENDIF
      close(nUnit,IOSTAT=ierror)  
    endif
    SumGrad=0.0_RK
#endif

#include "EnergySpectra_write_fun_inc.f90"

    nfstime=0; SumStat=0.0_RK; PrGradsum=0.0_RK; nSpectime=0; 
    if(clcSpectra1D) then
      EnergySpecX=0.0_RK; EnergySpecZ=0.0_RK
    endif
    if(clcSpectra2D) EnergySpec2D=0.0_RK
  end subroutine clcStat_CH

#undef iSpec1DUU
#undef iSpec1DVV
#undef iSpec1DWW
#undef iSpec1DPP
#undef iSpec1DUV
#undef iSpec1DUV2
#undef iLCSR1DUU
#undef iLCSI1DUU
#undef iLCSR1DUU2
#undef iLCSR1DVV
#undef iLCSI1DVV
#undef iLCSR1DVV2
#undef iLCSR1DWW
#undef iLCSI1DWW
#undef iLCSR1DWW2
#undef iLCSR1DPP
#undef iLCSI1DPP
#undef iLCSR1DPP2
#undef iSpec1DCC
#undef iSpec1DUC
#undef iImag1DUC
#undef iSpec1DVC
#undef iImag1DVC
#undef iLCSR1DCC
#undef iLCSI1DCC
#undef iLCSR1DCC2

#undef iSpec2DUU
#undef iSpec2DVV
#undef iSpec2DWW
#undef iSpec2DPP
#undef iSpec2DUV
#undef iLCSR2DUU
#undef iLCSR2DVV
#undef iLCSR2DWW
#undef iLCSR2DPP
#undef iSpec2DCC
#undef iLCSR2DCC

#undef NCHASTAT
#undef NEnergySpec1D
#undef NEnergySpec2D
#ifdef SAVE_SINGLE_Spec2D
#undef SAVE_SINGLE_Spec2D
#endif

#ifdef HighOrderGradStat
#undef HighOrderGradStat
#undef nGradStat
#endif
end module m_FlowType_Channel
module m_FlowType_TGVortex
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
  use m_Tools,only:CalcDissipationRate
  implicit none
  private

  public:: InitVelocity_TG, Update_uy_ym_TG
  public:: InitStatVar_TG,  clcStat_TG
contains

  !******************************************************************
  ! InitVelocity_TG
  !******************************************************************
  subroutine InitVelocity_TG(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    integer :: ic,jc,kc
    real(RK):: VelRef,LenRef,xpt,ypt,zpt,xct,yct,zct
      
    VelRef=1.0_RK
    LenRef=1.0_RK
    do kc=y1start(3),y1end(3)
      zpt=real(kc-1,kind=RK)*dz
      zct=zpt+0.5_RK*dz
      do jc=y1start(2),y1end(2)
        ypt=yp(jc)
        yct=yc(jc)
        do ic=y1start(1),y1end(1)
          xpt=real(ic-1,kind=RK)*dx
          xct=xpt+0.5_RK*dx
          ux(ic,jc,kc) =  VelRef*sin(xpt/LenRef)*cos(yct/LenRef)*cos(zct/LenRef)              
          uy(ic,jc,kc) = -VelRef*cos(xct/LenRef)*sin(ypt/LenRef)*cos(zct/LenRef)
          uz(ic,jc,kc) =  0.0_RK !VelRef*cos(xct/LenRef)*cos(yct/LenRef)*sin(zpt/LenRef)
        enddo
      enddo
    enddo
    
  end subroutine InitVelocity_TG

  !******************************************************************
  ! Update_uy_ym_TG
  !******************************************************************   
  subroutine Update_uy_ym_TG(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(inout):: uy_ym,duy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = 0.0_RK
     
    duy_ym = uy_ym - duy_ym
    
  end subroutine Update_uy_ym_TG

  !******************************************************************
  ! InitStatVar_TG
  !******************************************************************
  subroutine InitStatVar_TG()
    implicit none

    ! locals
    integer:: ierror,nUnit
    character(len=128)::filename

    if(nrank/=0) return
    write(filename,'(A,I10.10)')trim(ResultsDir)//"TG_dissp",ilast
    open(newunit=nUnit, file=filename,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_TG","Cannot open file: "//trim(filename))
    close(nUnit,IOSTAT=ierror)

  end subroutine InitStatVar_TG

  !******************************************************************
  ! clcStat_TG
  !******************************************************************
  subroutine clcStat_TG(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   
    ! locals
    real(Rk):: sum_dissp,sumr
    character(len=128)::filename
    integer::ic,jc,kc,nUnit,ierror
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::dissp
 
    sumr= 0.0_RK
    call CalcDissipationRate(ux,uy,uz,dissp)
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          sumr=sumr+ dissp(ic,jc,kc)
        enddo
      enddo
    enddo
    call MPI_REDUCE(sumr,sum_dissp,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)

    if(nrank==0) then
      write(filename,'(A,I10.10)')trim(ResultsDir)//"TG_dissp",ilast
      open(newunit=nUnit, file=filename, status='old',position='append',form='formatted',IOSTAT=ierror )
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_TG","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(2ES24.15)')SimTime,xnu*sum_dissp/real(nxc*nyc*nzc,RK)
      ENDIF
      close(nUnit,IOSTAT=ierror)
    endif
  end subroutine clcStat_TG
    
end module m_FlowType_TGVortex
module m_IOAndVisu
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Tools,only: Clc_Q_vor,Clc_lamda2
  use m_Variables,only:mb1,OutFlowInfoX,OutFlowInfoY
  implicit none
  private

  ! VisuOption
  integer:: iskip,jskip,kskip
  integer:: Prev_BackUp_itime  = 53456791
  logical:: save_ux,save_uy,save_uz,save_wx,save_wy,save_wz,save_wMag
  logical:: save_pr,save_Q_vor,save_lamda2,WriteHistOld,ReadHistOld

  public:: InitVisu, dump_visu, read_restart, write_restart, Delete_Prev_Restart

contains

  !******************************************************************
  ! InitVisu
  !******************************************************************
  subroutine InitVisu(ChannelPrm)
    implicit none
    character(*),intent(in)::ChannelPrm

    ! locals
    character(128)::XdmfFile
    integer::nUnitFile,ierror,nflds,ifld,iprec,i,j,k
    NAMELIST /IO_Options/ save_ux,save_uy,save_uz,save_pr,save_wx,save_wy,save_wz,save_wMag,save_Q_vor, &
                          save_lamda2,WriteHistOld,ReadHistOld,iskip,jskip,kskip
 
    open(newunit=nUnitFile, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitVisu", "Cannot openfile: "//trim(ChannelPrm))
    read(nUnitFile, nml=IO_Options)
    close(nUnitFile,IOSTAT=ierror)
    if(nrank==0) write(MainLog%nUnit,nml=IO_Options)

    ! write XDMF file
    if(nrank/=0) return
    write(xdmfFile,"(A)") trim(ResultsDir)//"VisuFor"//trim(RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitVisu","Cannot open file: "//trim(XdmfFile))
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! grid
    iprec=mytype_save
    write(nUnitFile,'(A,3I7,A)')'    <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',nzc,nyc,nxc,'"/>'
    write(nUnitFile,'(A)')'    <Geometry name="GEO" GeometryType="VXVYVZ">'
    ! x-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nxc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do i=1,nxc
      write(nUnitFile,'(E14.7)',advance='no') (i-1)*dx+dx*0.5_RK
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    ! y-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nyc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do j=1,nyc
      write(nUnitFile,'(E14.7)',advance='no') yc(j)
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    ! z-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do k=1,nzc
      write(nUnitFile,'(E14.7)',advance='no') (k-1)*dz+dz*0.5_RK
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')'    </Geometry>'

    ! Time series
    nflds = (ilast - ifirst +1)/SaveVisu  + 1
    write(nUnitFile,'(A)')'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(nUnitFile,'(A)')'        <Time TimeType="List">'
    write(nUnitFile,'(A,I6,A)')'        <DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no')'        '
    do ifld = ifirst-1,ilast,SaveVisu
      write(nUnitFile,'(I10)',advance='no') ifld
    enddo
    write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'       </Time>'

    ! attribute
    do  ifld=ifirst-1,ilast,SaveVisu
      write(nUnitFile,'(A,I10.10,A)')'        <Grid Name="T',ifld,'" GridType="Uniform">'
      write(nUnitFile,'(A)')'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(nUnitFile,'(A)')'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
      if(save_ux)    call Write_XDMF_One(nUnitFile,ifld,'ux')
      if(save_uy)    call Write_XDMF_One(nUnitFile,ifld,'uy')
      if(save_uz)    call Write_XDMF_One(nUnitFile,ifld,'uz')
      if(save_pr)    call Write_XDMF_One(nUnitFile,ifld,'pr')
      if(save_wx)    call Write_XDMF_One(nUnitFile,ifld,'wx')
      if(save_wy)    call Write_XDMF_One(nUnitFile,ifld,'wy')
      if(save_wz)    call Write_XDMF_One(nUnitFile,ifld,'wz')
      if(save_wMag)  call Write_XDMF_One(nUnitFile,ifld,'wMag')
      if(save_Q_vor) call Write_XDMF_One(nUnitFile,ifld,'Q' )
      if(save_lamda2)call Write_XDMF_One(nUnitFile,ifld,'lamda2')
      write(nUnitFile,'(A)')'        </Grid>'
    enddo

    write(nUnitFile,'(A)')'    </Grid>'
    write(nUnitFile,'(A)')'</Domain>'
    write(nUnitFile,'(A)')'</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)
#ifdef SaveNode
    call MainLog%OutInfo("Choose to save the visualizing file at grid node",2)
#else
    call MainLog%OutInfo("Choose to save the visualizing file at cell center",2)
#endif
  end subroutine InitVisu

  !******************************************************************
  ! Write_XDMF_One
  !******************************************************************
  subroutine Write_XDMF_One(nUnitFile, ifld,chAttribute)
    implicit none
    integer,intent(in)::nUnitFile,ifld
    character(*),intent(in)::chAttribute

    ! locals
    character(128)::chFile
    integer::iprec=mytype_save
    
    write(chFile,'(A,A,I10.10)')"VisuFor"//trim(RunName),"_"//trim(adjustl(chAttribute))//"_",ifld
    write(nUnitFile,'(A)')'            <Attribute Name="'//trim(chAttribute)//'" Center="Node">'
    write(nUnitFile,'(A,I1,A,3I7,A)')'                <DataItem Format="Binary" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzc,nyc,nxc,'">'
    write(nUnitFile,'(A)')'                    '//trim(chFile)
    write(nUnitFile,'(A)')'                </DataItem>'
    write(nUnitFile,'(A)')'            </Attribute>'
  end subroutine Write_XDMF_One

  !******************************************************************
  ! dump_visu
  !******************************************************************
  subroutine dump_visu(ntime,ux,uy,uz,pressure,ArrTemp)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::ArrTemp

    ! locals
    character(128)::chFile
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::dudy,dudz,dvdx,dvdz,dwdx,dwdy
    real(RK)::caj,cac1,cac2,cac12,vor_x,vor_y,vor_z
 
    ! ux
    if(save_ux) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
#ifdef SaveNode
            ArrTemp(ic,jc,kc)=ux(ic,jc,kc)
#else
            ArrTemp(ic,jc,kc)=0.5_RK*(ux(ic+1,jc,kc)+ux(ic,jc,kc))
#endif
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_ux_",ntime
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! uy
    if(save_uy) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
#ifdef SaveNode
            ArrTemp(ic,jc,kc)=uy(ic,jc,kc)
#else
            ArrTemp(ic,jc,kc)=0.5_RK*(uy(ic,jc+1,kc)+uy(ic,jc,kc))
#endif
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_uy_",ntime
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! uz
    if(save_uz) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
#ifdef SaveNode
            ArrTemp(ic,jc,kc)=uz(ic,jc,kc)
#else
            ArrTemp(ic,jc,kc)=0.5_RK*(uz(ic,jc,kc+1)+uz(ic,jc,kc))
#endif
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_uz_",ntime
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! pressure
    if(save_pr) then
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_pr_",ntime
      call decomp_2d_write_every(y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wx
    if(save_wx) then
      do kc=y1start(3),y1end(3)
        kp=kc+1
        km=kc-1
        do jc=y1start(2),y1end(2)
          jp=jc+1
          jm=jc-1
          cac1 = rdyc(jc)
          cac2 = rdyc(jp)    
          cac12= cac1 - cac2
          do ic=y1start(1),y1end(1)
            dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *0.25_RK
            dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                  +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                  -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *0.25_RK
            ArrTemp(ic,jc,kc)= dwdy -dvdz
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wx_",ntime
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wy
    if(save_wy) then
      do kc=y1start(3),y1end(3)
        kp=kc+1
        km=kc-1
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            ip=ic+1
            im=ic-1
            dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *0.25_RK
            dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *0.25_RK
            ArrTemp(ic,jc,kc)= dudz -dwdx
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wy_",ntime
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wz
    if(save_wz) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          jp=jc+1
          jm=jc-1
          caj  = rdyp(jc)
          cac1 = rdyc(jc)
          cac2 = rdyc(jp)    
          cac12= cac1 - cac2
          do ic=y1start(1),y1end(1)
            ip=ic+1
            im=ic-1
            dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2  &
                  +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12 &
                  -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1  )  *0.25_RK
            dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *0.25_RK
            ArrTemp(ic,jc,kc)= dvdx -dudy
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wz_",ntime
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wMag
    if(save_wMag) then
      do kc=y1start(3),y1end(3)
        kp=kc+1
        km=kc-1
        do jc=y1start(2),y1end(2)
          jp=jc+1
          jm=jc-1
          caj  = rdyp(jc)
          cac1 = rdyc(jc)
          cac2 = rdyc(jp)    
          cac12= cac1 - cac2
          do ic=y1start(1),y1end(1)
            ip=ic+1
            im=ic-1

            dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2    &
                  +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12   &
                  -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1    )    *0.25_RK
            dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *0.25_RK
      
            dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *0.25_RK
            dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *0.25_RK

            dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *0.25_RK
            dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                  +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                  -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *0.25_RK

            vor_x= dwdy-dvdz
            vor_y= dudz-dwdx
            vor_z= dvdx-dudy
            ArrTemp(ic,jc,kc)= sqrt(vor_x*vor_x +vor_y*vor_y +vor_z*vor_z)
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wMag_",ntime 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! Q
    if(save_Q_vor) then
      call Clc_Q_vor(ux,uy,uz,ArrTemp)
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_Q_",ntime 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif  
  
    ! lamda2
    if(save_lamda2) then
      call Clc_lamda2(ux,uy,uz,ArrTemp)
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_lamda2_",ntime 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif
  end subroutine dump_visu

  !**********************************************************************
  ! Delete_Prev_Restart
  !**********************************************************************
  subroutine Delete_Prev_Restart(ntime)
    implicit none
    integer,intent(in)::ntime

    ! locals
    integer::nUnit,ierror
    character(128)::chFile

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if(nrank/=0) return
    ! 
    write(chFile,"(A,I10.10)") trim(RestartDir)//"RestartFor"//trim(RunName),Prev_BackUp_itime
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    ! 
    write(chFile,"(A,I10.10)") trim(RestartDir)//"OutFlowFor"//trim(RunName),Prev_BackUp_itime 
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    !
    write(chFile,"(A,I10.10)") trim(RestartDir)//"PrDataFor"//trim(RunName),Prev_BackUp_itime 
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    !
    Prev_BackUp_itime = ntime
  end subroutine Delete_Prev_Restart

  !******************************************************************
  ! write_restart
  !******************************************************************
  subroutine write_restart(ntime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in):: HistXOld,HistYOld,HistZOld

    ! locals
    character(128)::chFile
    integer(MPI_OFFSET_KIND)::disp
    integer::fh,ierror,ic,jc,kc,newtype,sizes(3),subsizes(3),starts(3)

    ! begin to write restart file
    write(chFile,"(A,I10.10)") trim(RestartDir)//"RestartFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_var(fh,disp,y_pencil,      ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,      uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,      uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    if(ischeme==FI_AB2 .and. WriteHistOld) then
      call decomp_2d_write_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistZOld)
    endif
    call MPI_FILE_CLOSE(fh,ierror)

    ! Write PrGradData
    if(IsUxConst .and. nrank==0) then
      write(chFile,"(A,I10.10)") trim(RestartDir)//"PrDataFor"//trim(RunName),ntime
      open(newunit=fh, file=chFile, status='replace', action='write', IOSTAT=ierror)
      if(ierror/=0) then
        call MainLog%CheckForError(ErrT_Abort,"write_restart","Cannot open file: "//trim(chFile))
      else
        write(fh,'(4ES26.17)') PrGradData(1:4)
      endif
      close(fh,IOSTAT=ierror)
    endif
    
    ! Write OutFlow
    if(BcOption(xp_dir)/=BC_OutFlow .and. BcOption(yp_dir)/=BC_OutFlow) return
    write(chFile,"(A,I10.10)") trim(RestartDir)//"OutFlowFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    disp = 0_MPI_OFFSET_KIND 
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyp
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          OutFlowInfoY(4,ic,kc)=ux(ic,jc,kc)
          OutFlowInfoY(5,ic,kc)=uy(ic,jc,kc)
          OutFlowInfoY(6,ic,kc)=uz(ic,jc,kc)         
        enddo
      enddo
      sizes=[3, nxc, nzc]
      subsizes=[3,y1size(1),y1size(3)]
      starts=[0,y1start(1)-1,y1start(3)-1]
      call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(fh,OutFlowInfoY(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_WRITE_ALL(fh,OutFlowInfoY(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
    endif
    if(myProcNghBC(y_pencil,3)==BC_OutFlow) then
      ic=nxp
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          OutFlowInfoX(4,jc,kc)=ux(ic,jc,kc)
          OutFlowInfoX(5,jc,kc)=uy(ic,jc,kc)
          OutFlowInfoX(6,jc,kc)=uz(ic,jc,kc)          
        enddo
      enddo
      sizes=[3, nyc, nzc]
      subsizes=[3,y1size(2),y1size(3)]
      starts=[0,y1start(2)-1,y1start(3)-1]
      call MPI_FILE_OPEN(DECOMP_2D_COMM_ROW, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(fh,OutFlowInfoX(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_WRITE_ALL(fh,OutFlowInfoX(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)    
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
    endif    
  end subroutine write_restart

  !******************************************************************
  ! read_restart
  !******************************************************************
  subroutine read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out):: HistXOld,HistYOld,HistZOld

    ! locals
    character(128)::chFile
    integer(MPI_OFFSET_KIND)::disp,byte_total1,byte_total2,filebyte
    integer::fh,ierror,ntime,ic,jc,kc,newtype,sizes(3),subsizes(3),starts(3)
    
    ! begin to write restart file
    ntime= ifirst - 1
    write(chFile,"(A,I10.10)") trim(RestartDir)//"RestartFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Read_Restart","Cannot open file: "//trim(chFile))
    call MPI_FILE_GET_SIZE(fh,filebyte,ierror)
    byte_total1=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*7_MPI_OFFSET_KIND
    byte_total2=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*4_MPI_OFFSET_KIND
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      if(filebyte /= byte_total1 .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong1")
      endif      
    else
      if((filebyte /= byte_total1 .and. filebyte /= byte_total2) .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong2")
      endif      
    endif
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_var(fh,disp,y_pencil,      ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,      uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,      uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      call decomp_2d_read_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistZOld)
    endif
    call MPI_FILE_CLOSE(fh,ierror)

    ! Read PrGradData
    if(IsUxConst) then
      write(chFile,"(A,I10.10)") trim(RestartDir)//"PrDataFor"//trim(RunName),ntime
      open(newunit=fh, file=chFile, status='old', action='read', IOSTAT=ierror)
      if(ierror/=0) then
        if(nrank==0) then
          call MainLog%OutInfo("read_restart: Cannot open file "//trim(chFile),1)
          call MainLog%OutInfo(" PrGradData=0.0 will be used ! ",2)
        endif
        PrGradData=0.0_RK
      else
        read(fh,*) PrGradData(1:4)
      endif
      close(fh,IOSTAT=ierror)
    endif
    
    ! Read OutFlow
    if(BcOption(xp_dir)/=BC_OutFlow .and. BcOption(yp_dir)/=BC_OutFlow) return
    write(chFile,"(A,I10.10)") trim(RestartDir)//"OutFlowFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Read_Restart","Cannot open file: "//trim(chFile))
    call MPI_FILE_GET_SIZE(fh,filebyte,ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    byte_total1=0; byte_total2=0;
    if(BcOption(xp_dir)==BC_OutFlow) then
      byte_total1=byte_total1+int(mytype_bytes,8)*int(nyc,8)*int(nzc,8)*6_8
      byte_total2=byte_total2+int(mytype_bytes,8)*int(nyc,8)*int(nzc,8)*3_8
    endif
    if(BcOption(yp_dir)==BC_OutFlow) then
      byte_total1=byte_total1+int(mytype_bytes,8)*int(nxc,8)*int(nzc,8)*6_8
      byte_total2=byte_total2+int(mytype_bytes,8)*int(nxc,8)*int(nzc,8)*3_8
    endif
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      if(filebyte /= byte_total1 .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong3")
      endif      
    else
      if((filebyte /= byte_total1 .and. filebyte /= byte_total2) .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong4")
      endif      
    endif
    disp = 0_MPI_OFFSET_KIND 
    if(BcOption(yp_dir)==BC_OutFlow) then
      sizes=[3, nxc, nzc]
      subsizes=[3,y1size(1),y1size(3)]
      starts=[0,y1start(1)-1,y1start(3)-1]
      call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_READ_ALL(fh,OutFlowInfoY(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_READ_ALL(fh,OutFlowInfoY(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
      jc=nyp
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)=OutFlowInfoY(4,ic,kc)
          uy(ic,jc,kc)=OutFlowInfoY(5,ic,kc)
          uz(ic,jc,kc)=OutFlowInfoY(6,ic,kc)       
        enddo
      enddo   
    endif    
    if(myProcNghBC(y_pencil,3)==BC_OutFlow) then
      sizes=[3, nyc, nzc]
      subsizes=[3,y1size(2),y1size(3)]
      starts=[0,y1start(2)-1,y1start(3)-1]
      call MPI_FILE_OPEN(DECOMP_2D_COMM_ROW, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_READ_ALL(fh,OutFlowInfoX(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_READ_ALL(fh,OutFlowInfoX(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
      ic=nxp
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          ux(ic,jc,kc)=OutFlowInfoX(4,jc,kc)
          uy(ic,jc,kc)=OutFlowInfoX(5,jc,kc)
          uz(ic,jc,kc)=OutFlowInfoX(6,jc,kc)         
        enddo
      enddo                      
    endif
  end subroutine read_restart
end module m_IOAndVisu
module m_MeshAndMetries
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  implicit none
  private 
  
  real(RK),public::dx,  dy,  dz   ! average mesh intervals in three directions
  real(RK),public::dx2, dy2, dz2  ! square of average mesh intervals in three directions  
  real(RK),public::rdx, rdy, rdz  ! inverse average mesh intervals in three directions
  real(RK),public::rdx2,rdy2,rdz2 ! square of the inverse average mesh intervals in three directions
#if defined CFDACM
  real(RK),public:: dyUniform,rdyUniform
#endif

  real(RK),public,allocatable,dimension(:):: yp   ! point coordinate in y-dir. Suffix 'v' means 'vector'
  real(RK),public,allocatable,dimension(:):: xc   ! center coordinate in x-dir    
  real(RK),public,allocatable,dimension(:):: yc   ! center coordinate in y-dir 
  real(RK),public,allocatable,dimension(:):: zc   ! center coordinate in z-dir
  real(RK),public,allocatable,dimension(:):: VolCell   ! cell volume
  real(RK),public,allocatable,dimension(:):: DeltaCell ! (dx*dy*dz)^(1/3)

  real(RK),public,allocatable,dimension(:):: dyp  ! point coordinate interval in y-dir
  real(RK),public,allocatable,dimension(:):: dyc  ! center coordinate interval in y-dir
  real(RK),public,allocatable,dimension(:):: rdyp ! inverse point coordinate interval in y-dir
  real(RK),public,allocatable,dimension(:):: rdyc ! inverse center coordinate interval in y-dir

  ! Pressure Laplacian metries in y-dir
  real(RK),public,allocatable,dimension(:):: ap2Pr
  real(RK),public,allocatable,dimension(:):: ac2Pr  
  real(RK),public,allocatable,dimension(:):: am2Pr 

  ! uy/uz Laplacian metries in x-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am1c
  real(RK),public,allocatable,dimension(:):: ac1c
  real(RK),public,allocatable,dimension(:):: ap1c  
  
  ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am2c
  real(RK),public,allocatable,dimension(:):: ac2c
  real(RK),public,allocatable,dimension(:):: ap2c
  
  ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am2p
  real(RK),public,allocatable,dimension(:):: ac2p  
  real(RK),public,allocatable,dimension(:):: ap2p
  
  ! ux/uy Laplacian metries in z-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am3c
  real(RK),public,allocatable,dimension(:):: ac3c
  real(RK),public,allocatable,dimension(:):: ap3c 

  ! for linear interpolation in y-dir
  real(RK),public,allocatable,dimension(:):: YinterpCoe
  
  public:: InitMeshAndMetries
contains

  !******************************************************************
  ! InitMeshAndMetries
  !****************************************************************** 
  subroutine InitMeshAndMetries(ChannelPrm)
    implicit none 
    character(*),intent(in)::ChannelPrm
    
    ! locals
    integer::nSection
    character(128)::chFile
    NAMELIST/MeshSection/nSection
    integer::j,nUnitFile,ierrTmp,ierror=0
    real(RK),allocatable,dimension(:)::SectionLength,SectioncStret,ySectionCoord
    integer,allocatable,dimension(:)::nycSection,StretType,StretOption,nidYSection
    NAMELIST/MeshOptions/SectionLength,SectioncStret,nycSection,StretType,StretOption
    
    dx = xlx/real(nxc,kind=RK)
    dy = yly/real(nyc,kind=RK)
    dz = zlz/real(nzc,kind=RK)
    rdx= real(nxc,kind=RK)/xlx
    rdy= real(nyc,kind=RK)/yly
    rdz= real(nzc,kind=RK)/zlz 
    dx2= dx*dx;  rdx2= rdx*rdx
    dy2= dy*dy;  rdy2= rdy*rdy
    dz2= dz*dz;  rdz2= rdz*rdz
    
    allocate(yp(0:nyp),    Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(yc(0:nyp),    Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(dyp(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(rdyp(0:nyp),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(dyc(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(rdyc(0:nyp),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)    
    
    allocate(ap2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)   
    allocate(am2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap1c(1:nxc),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac1c(1:nxc),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am1c(1:nxc),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)   
    allocate(am2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap3c(1:nzc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac3c(1:nzc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am3c(1:nzc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)

    allocate(xc(0:nxp),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(zc(0:nzp),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(VolCell(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(DeltaCell(0:nyp), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)

    allocate(YinterpCoe(1:nyp), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Allocation failed")

    ! Read mesh options and claculate yp for every section
    open(newunit=nUnitFile, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "Cannot open file: "//trim(ChannelPrm))
    read(nUnitFile, nml=MeshSection)
    if(nSection<0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nSection Wrong !!!")
    allocate(nycSection(nSection),StretType(nSection),StretOption(nSection))
    allocate(SectionLength(nSection),SectioncStret(nSection))
    allocate(nidYSection(nSection+1),ySectionCoord(nSection+1))
    read(nUnitFile, nml=MeshOptions)
    if(nrank==0) then
      do j=1,nSection
        if(SectionLength(j)<0.0_RK)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectionLength Wrong !!!")
        if(SectioncStret(j)<0.0_RK)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectioncStret Wrong !!!")
        if(nycSection(j)<1)      call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 1 !!!")
        if(StretType(j)<0  .or. StretType(j)>3 ) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretType Wrong !!!")
        if(StretOption(j)<0.or. StretOption(j)>1)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretOption Wrong !!!")
      enddo
      if(sum(nycSection) /= nyc) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 2 !!!")
      write(MainLog%nUnit, nml=MeshSection)
      write(MainLog%nUnit, nml=MeshOptions)
    endif
    close(nUnitFile,IOSTAT=ierror)
    nidYSection=1
    ySectionCoord=0.0_RK
    Block 
      real(RK)::SumLength
      SumLength=sum(SectionLength)
      do j=1,nSection
        nidYSection(j+1)  =nidYSection(j)   + nycSection(j)
        ySectionCoord(j+1)=ySectionCoord(j) + SectionLength(j)/SumLength*yly
      enddo
    End Block
    ySectionCoord(nSection+1)=yly
    deallocate(nycSection,SectionLength)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    do j=1,nSection
      call clcCoord(yp(1:nyp),nidYSection(j:j+1),ySectionCoord(j:j+1),StretType(j),SectioncStret(j),StretOption(j))
    enddo
    deallocate(StretType,StretOption,SectioncStret,nidYSection,ySectionCoord)
    write(chFile,"(A)") trim(ResultsDir)//"yMeshFor"//trim(RunName)//".txt"
    if(nrank==0) then
      open(newunit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=ierror)
      if(ierror/=0.and.nrank==0)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file: "//trim(chFile))
      do j=1,nyp
        write(nUnitFile,*)j,yp(j)
      enddo
      close(nUnitFile,IOSTAT=ierror)
    endif

    ! yc, center coordinate interval in y-dir
    do j=1,nyc
      yc(j) = 0.5_RK*(yp(j)+yp(j+1))
    enddo
    yc(0)=2.0_RK*yp(1)-yc(1)
    yc(nyp)=2.0_RK*yp(nyp)-yc(nyc)

    ! xc,zc, center coordinate interval in x-dir and z-dir
    do j=0,nxp
      xc(j)=dx*(real(j,RK)-0.5_RK)
    enddo
    do j=0,nzp
      zc(j)=dz*(real(j,RK)-0.5_RK)
    enddo

    ! dyp, point coordinate interval in y-dir 
    do j=1,nyc
      dyp(j) = yp(j+1) - yp(j)    
    enddo
    dyp(0)  =dyp(1)
    dyp(nyp)=dyp(nyc)
#ifdef CFDACM
    dyUniform=dyp(1)
    rdyUniform=1.0_RK/dyUniform
#endif

    ! VolCell,volume of the cell
    do j=0,nyp
      VolCell(j)=dyp(j)*dx*dz
    enddo
  
    ! DeltaCell, (VolCell)^(1/3)
    do j=0,nyp
      DeltaCell(j)=(dyp(j)*dx*dz)**(0.333333333333333333333333333_RK)
    enddo

    ! dyc, center coordinate interval in y-dir
    do j=2,nyc
      dyc(j) = yc(j)-yc(j-1)
    enddo
    dyc(1)=dyp(1)
    dyc(nyp)=dyp(nyc)
    
    ! rdyp and rdyc, the reverse of the dyp and dyc, respectively.
    do j=0, nyp
      rdyp(j) = 1.0_RK/dyp(j)    
    enddo    
    do j=1,nyp
      rdyc(j) = 1.0_RK/dyc(j)    
    enddo

    ! for linear interpolation in y-dir
    do j=1,nyp
      YinterpCoe(j)= dyp(j)/(dyp(j)+dyp(j-1))
    enddo
    
    ! Pressure Laplacian metries in y-dir
    do j=1,nyc
      am2Pr(j)= rdyp(j)*rdyc(j)
      ap2Pr(j)= rdyp(j)*rdyc(j+1)
    enddo
    if(BcOption(ym_dir)/=BC_Period) then
      am2Pr(1)= 0.0_RK
      ap2Pr(1)= rdyp(1)*rdyc(2) 
    endif
    if(BcOption(yp_dir)/=BC_Period) then
      am2Pr(nyc)= rdyp(nyc)*rdyc(nyc)
      ap2Pr(nyc)= 0.0_RK
    endif    
    ac2Pr = -(am2Pr+ap2Pr)
    
    ! uy/uz Laplacian metries in x-dir (STAGGERED VARIABLE)=====================
    am1c = rdx2;  ap1c=rdx2
    if(BcOption(xm_dir)==BC_NoSlip ) then
      am1c(1)= 4.0_RK/3.0_RK*rdx2
      ap1c(1)= 4.0_RK/3.0_RK*rdx2
    endif    
    if(BcOption(xp_dir)==BC_NoSlip ) then
      am1c(nxc)= 4.0_RK/3.0_RK*rdx2
      ap1c(nxc)= 4.0_RK/3.0_RK*rdx2 
    endif 
    ac1c= -(am1c+ap1c)
    
    ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)=====================
    do j=1,nyc
      am2c(j)= rdyp(j)*rdyc(j)
      ap2c(j)= rdyp(j)*rdyc(j+1)
    enddo
    if(BcOption(ym_dir)==BC_NoSlip ) then
      am2c(1)= 4.0_RK*rdyc(1)/( dyc(1)+2.0_RK*dyc(2) )
      ap2c(1)= 4.0_RK*rdyc(2)/( dyc(1)+2.0_RK*dyc(2) )
    endif
    if(BcOption(yp_dir)==BC_NoSlip ) then
      am2c(nyc)= 4.0_RK*rdyc(nyc)/( dyc(nyp)+2.0_RK*dyc(nyc) )
      ap2c(nyc)= 4.0_RK*rdyc(nyp)/( dyc(nyp)+2.0_RK*dyc(nyc) )
    endif
    ac2c= -(am2c+ap2c)
    
    ! ux/uy Laplacian metries in z-dir (STAGGERED VARIABLE)=====================
    am3c = rdz2;  ap3c=rdz2;
    if(BcOption(zm_dir)==BC_NoSlip ) then
      am3c(1)= 4.0_RK/3.0_RK*rdz2
      ap3c(1)= 4.0_RK/3.0_RK*rdz2
    endif
    if(BcOption(zp_dir)==BC_NoSlip ) then
      am3c(nzc)= 4.0_RK/3.0_RK*rdz2
      ap3c(nzc)= 4.0_RK/3.0_RK*rdz2
    endif
    ac3c= -(am3c+ap3c)

    ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
    do j=1,nyc
      am2p(j)= rdyc(j)*rdyp(j-1)        
      ap2p(j)= rdyc(j)*rdyp(j)        
    enddo
    ac2p= -(am2p+ap2p)    
  end subroutine InitMeshAndMetries

  !******************************************************************
  ! clcCoord
  !******************************************************************
  subroutine clcCoord(coordinate,ncoorId,SectionCoord,StretType,cStret,StretOption)
    implicit none
    real(RK),dimension(:),intent(inout)::coordinate
    real(RK),intent(in)::SectionCoord(2),cStret
    integer,intent(in)::ncoorId(2),StretType,StretOption
    
    ! locals
    integer::j,jt,m
    real(RK)::secLen,tstr,xi,setCrd
    
    m=ncoorId(2)-ncoorId(1)
    secLen=SectionCoord(2)-SectionCoord(1)
    coordinate(ncoorId(1))=SectionCoord(1)
    coordinate(ncoorId(2))=SectionCoord(2)
    SELECT CASE(StretType)
    CASE(0) ! Uniform
      do j=0,m
        jt=ncoorId(1)+j
        coordinate(jt)=SectionCoord(1)+real(j,RK)/real(m,RK)*secLen
      enddo
    CASE(1) ! Tangent hyperbolic function
      tstr= tanh(cStret)
      do j=0,m
        xi= real(j,kind=RK)/real(m)               ! For j: [0,m],  xi: [0, 1]
        setCrd= tanh(cStret*(xi-1.0_RK))/tstr + 1.0_RK  ! For j: [0,m],  setCrd: [0, 1]
        if(StretOption==0) then  ! bottom
          jt=ncoorId(1)+j
          coordinate(jt)=SectionCoord(1) + setCrd*secLen
        else                     ! top
          jt=ncoorId(2)-j
          coordinate(jt)=SectionCoord(2) - setCrd*secLen
        endif        
      enddo
    CASE(2) ! Sine/cosine function
      tstr=sin(0.5_RK*cStret*PI)
      do j=0,m
        xi= real(j,RK)/real(m,RK)-1.0_RK             ! For j: [0,m],  xi: [-1,0]
        setCrd= sin(cStret*xi*PI*0.5_RK)/tstr + 1.0_RK ! For j: [0,m],  setCrd: [0, 1]
        if(StretOption==0) then  ! bottom
          jt=ncoorId(1)+j
          coordinate(jt)=SectionCoord(1) + setCrd*secLen
        else                     ! top
          jt=ncoorId(2)-j
          coordinate(jt)=SectionCoord(2) - setCrd*secLen
        endif
      enddo
    CASE(3) ! Proportional sequence
      if(cStret==1.0_RK) then
        do j=0,m
          jt=j+ncoorId(1)
          coordinate(jt)=SectionCoord(1)+real(j,RK)/real(m,RK)*secLen
        enddo
      else
        tstr=(cStret**m -1.0_RK)/(cStret-1.0_RK)
        do j=0,m-1
          setCrd=(cStret**j)/tstr
          if(StretOption==0) then  ! bottom
            jt=ncoorId(1)+j
            coordinate(jt+1)=coordinate(jt) + setCrd*secLen
          else                     ! top
            jt=ncoorId(2)-j
            coordinate(jt-1)=coordinate(jt) - setCrd*secLen
          endif
        enddo
      endif
    END SELECT
    coordinate(ncoorId(1))=SectionCoord(1)
    coordinate(ncoorId(2))=SectionCoord(2)  
  end subroutine clcCoord

end module m_MeshAndMetries
module m_Parameters
  use m_LogInfo
  use m_TypeDef
  implicit none  
  private

  ! Log 
  type(LogType),public::MainLog
  
  ! Decomp2d options
  integer,public:: p_row,p_col
    
  ! Flow type option
  integer,parameter,public:: FT_CH=1  ! Channel
  integer,parameter,public:: FT_HC=2  ! 0.5_RK Channel
  integer,parameter,public:: FT_TG=3  ! Taylor-Green vortex
  integer,parameter,public:: FT_HI=4  ! Homogenerous isotropic turbulence
  integer,parameter,public:: FT_AN=5  ! Added new
  logical,public::IsUxConst=.false.
  real(RK),public::uBulk=0.0
  integer,public:: FlowType
  real(RK),dimension(4),public:: PrGradData=0.0_RK ! PrGradAve, PrGradAveOld, PrGradNow, ForcedOld
      
  ! mesh options
  real(RK),public::xlx, yly, zlz  ! domain length in three directions
  integer,public::nxp,nyp,nzp     ! grid point number in three directions
  integer,public::nxc,nyc,nzc     ! grid center number in three directions. nxc = nxp-1
  
  integer,parameter,public:: x_pencil=1
  integer,parameter,public:: y_pencil=2
  integer,parameter,public:: z_pencil=3
  integer,parameter,public:: xm_dir=1  ! x- direction
  integer,parameter,public:: xp_dir=2  ! x+ direction
  integer,parameter,public:: ym_dir=3  ! y- direction
  integer,parameter,public:: yp_dir=4  ! y+ direction
  integer,parameter,public:: zm_dir=5  ! z- direction
  integer,parameter,public:: zp_dir=6  ! z+ direction  
  
  ! Physical properties
  real(RK),public:: xnu                    ! Kinematic viscosity
  real(RK),public:: FluidDensity           ! Fluid density 
  real(RK),dimension(3),public:: gravity   ! Gravity or  other constant body forces (if any)
  
  ! Time stepping scheme and Projection method options
  integer,parameter,public:: FI_AB2 =1          !  AB2 for convective term
  integer,parameter,public:: FI_RK2 =2          !  RK2 for convective term
  integer,parameter,public:: FI_RK3 =3          !  RK3 for convective term
  real(RK),public:: dt         ! current time step
  real(RK),public:: dtMax      ! Maxium time step
  real(RK),public:: SimTime=0.0_RK  ! Real simulation time
  integer,public :: iCFL       ! Use CFL condition to change time step dynamically( 1: yes, 2:no ).
  real(RK),public:: CFLc       ! CFL parameter
  integer,public::  itime      ! current time step
  integer,public::  ifirst     ! First iteration
  integer,public::  ilast      ! Last iteration 
  integer,public::  iadvance
  integer,public::  ischeme
  integer,public::  IsImplicit  !(0=full explicit, 1=partial implicit, 2=full implicit )
  real(RK),public::pmGamma
  real(RK),public::pmTheta
  real(RK),public::pmAlpha
  real(RK),public::pmBeta

  ! FFT_option
  integer,public:: FFTW_plan_type
  
  ! Boundary conditions
  !  0: Periodic
  ! -1: NoSlip
  ! -2: Slip
  ! -3: Convective (ONLY AVAILABLE for x+ and y+)
  integer,parameter,public::BC_PERIOD  =  0 
  integer,parameter,public::BC_NoSlip  = -1
  integer,parameter,public::BC_FreeSlip= -2
  integer,parameter,public::BC_OutFlow = -3
  integer, public,dimension(6):: BcOption
  real(RK),public,dimension(6):: uxBcValue
  real(RK),public,dimension(6):: uyBcValue
  real(RK),public,dimension(6):: uzBcValue

  ! I/O, Statistics
  integer,public:: ivstats            ! time step interval for statistics calculation 
  integer,public:: SaveVisu           ! Output visulizing file frequency
  integer,public:: BackupFreq         ! Output Restarting file frequency
  integer,public:: SaveStat           ! Output Statistics file frequency

  logical,public::       RestartFlag  ! restart or not
  character(64),public:: RunName      ! Run name
  character(64),public:: ResultsDir   ! Result directory
  character(64),public:: RestartDir   ! Restart directory
  integer,public:: Cmd_LFile_Freq= 1  ! report frequency in the terminal 
  integer,public:: LF_file_lvl   = 5  ! logfile report level      
  integer,public:: LF_cmdw_lvl   = 3  ! terminal report level

  ! limited velocity and div
  real(RK),public:: vel_limit
  real(RK),public:: div_limit

  public:: ReadAndInitParameters,DumpReadedParam,PMcoeUpdate
contains
    
  !******************************************************************
  ! InitParameters
  !****************************************************************** 
  subroutine ReadAndInitParameters(chFile)
    implicit none 
    character(*),intent(in)::chFile
    
    ! locals
    integer:: nUnitFile,ierror
    NAMELIST/BasicParam/FlowType,IsUxConst,uBulk,xlx,yly,zlz,nxc,nyc,nzc,xnu,dtMax,iCFL,CFLc,ifirst,    &
                        ilast,ischeme,IsImplicit,FFTW_plan_type,BcOption,gravity,uxBcValue,uyBcValue,   &
                        uzBcValue,ivstats,BackupFreq,SaveStat,SaveVisu,RestartFlag,RunName,Cmd_LFile_Freq,   &
                        ResultsDir,RestartDir,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,vel_limit,div_limit,FluidDensity
 
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0) then
      print*,"Cannot open file: "//trim(chFile); STOP
    endif
    read(nUnitFile, nml=BasicParam)
    close(nUnitFile,IOSTAT=ierror)  

    nxp= nxc+1
    nyp= nyc+1
    nzp= nzc+1
    if(ischeme==FI_AB2) then
      iadvance = 1
    elseif(ischeme==FI_RK2) then
      iadvance = 2
    elseif(ischeme==FI_RK3) then
      iadvance = 3
    else
      print*,"Time scheme WRONG!!! ischeme=",ischeme; STOP      
    endif   
  end subroutine ReadAndInitParameters

  !******************************************************************
  ! DumpReadedParam
  !****************************************************************** 
  subroutine DumpReadedParam()
    implicit none

    ! locals
    NAMELIST/BasicParam/FlowType,IsUxConst,uBulk,xlx,yly,zlz,nxc,nyc,nzc,xnu,dtMax,iCFL,CFLc,ifirst,    &
                        ilast,ischeme,IsImplicit,FFTW_plan_type,BcOption,gravity,uxBcValue,uyBcValue,   &
                        uzBcValue,ivstats,BackupFreq,SaveStat,SaveVisu,RestartFlag,RunName,Cmd_LFile_Freq,   &
                        ResultsDir,RestartDir,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,vel_limit,div_limit,FluidDensity
    write(MainLog%nUnit, nml=BasicParam)
  end subroutine DumpReadedParam

  !******************************************************************
  ! PMcoeUpdate
  !******************************************************************
  subroutine PMcoeUpdate(ns)
    implicit none
    integer,intent(in)::ns

    ! locals
    real(RK),dimension(3):: pmGammaConst,pmThetaConst,pmAlphaConst

    if(ischeme==FI_AB2) then
      pmGammaConst= [ 1.5_RK, 0.0_RK, 0.0_RK]
      pmThetaConst= [-0.5_RK, 0.0_RK, 0.0_RK]
      if((.not. RestartFlag) .and. itime==ifirst) then
        pmGammaConst= [ 1.0_RK, 0.0_RK, 0.0_RK]
        pmThetaConst= [ 0.0_RK, 0.0_RK, 0.0_RK]
      endif
    elseif(ischeme==FI_RK2) then
      pmGammaConst= [ 0.5_RK, 1.0_RK, 0.0_RK]
      pmThetaConst= [ 0.0_RK,-0.5_RK, 0.0_RK]
    else
      pmGammaConst= [8.0_RK/15.0_RK, 5.0_RK/12.0_RK, 0.75_RK]
      pmThetaConst= [0.0_RK, -17.0_RK/60.0_RK, -5.0_RK/12.0_RK]   
    endif
    pmAlphaConst = pmGammaConst + pmThetaConst
       
    pmGamma = pmGammaConst(ns) *dt
    pmTheta = pmThetaConst(ns) *dt
    pmAlpha = pmAlphaConst(ns) *dt
    pmBeta  = 0.5_RK*pmAlphaConst(ns) *xnu *dt
    SimTime = SimTime+dt*pmAlphaConst(ns)
    if(ns==1) then
      PrGradData(2)=PrGradData(1)
      PrGradData(1)=0.0_RK
    endif
  end subroutine PMcoeUpdate
end module m_Parameters
module m_Poisson
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters
  use iso_c_binding
  use m_MeshAndMetries
  use m_Variables,only: mb1
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal
  implicit none
  private
  include "myfftw3.f03"
  
  real(RK)::normfft
  type(decomp_info),allocatable::decomp_PPE
  real(RK),allocatable,dimension(:)::WaveNumX,WaveNumZ
  type(C_PTR)::fwd_plan_x,bwd_plan_x,fwd_plan_z,bwd_plan_z
  real(RK),allocatable,dimension(:,:,:)::a_reduce,c_reduce 
  procedure(),pointer::clcPPE=>null(),execute_FFTW_r2r_z=>null()
  
  public:: InitPoissonSolver,clcPPE,Destory_Poisson_FFT_Plan
contains
#define nTime_FFT_Test 5

#define my_FFTW_inc_add_z2
#include "my_FFTW_inc.f90"
#undef  my_FFTW_inc_add_z2

#define my_Poisson_inc_add_Periodic_2d
#include "my_Poisson_inc.f90"
#undef  my_Poisson_inc_add_Periodic_2d
    
  !******************************************************************
  ! InitPoissonSolver
  !******************************************************************     
  subroutine InitPoissonSolver()
    implicit none
    
    ! locals
    real(RK)::WaveCoe,wa1,wa3,best_time,t1(2),t2(2),normTmp
    integer::i,k,iErr01,iErr02,iChoice,iFFTz,IsReduce,IsReduceR
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm)::prphiHalo
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::prsrc
        
    prsrc=0.0_RK
    normfft = 1.0_RK
    
    ! Modified wave number in x-dir
    allocate(WaveNumX(nxc), Stat =iErr01)
    IF(BcOption(xm_dir)==BC_PERIOD) THEN    
      WaveCoe = 2.0_RK*PI
      normTmp = 1.0_RK
    ELSE
      WaveCoe = PI
      normTmp = 2.0_RK
    ENDIF
    do i=1,nxc
      wa1= WaveCoe*real(i-1,RK)/real(nxc,RK)
      WaveNumX(i)=2.0_RK*rdx2*(cos(wa1)-1.0_RK)
    enddo
    normfft = normfft*normTmp*real(nxc,RK)

    ! Modified wave number in z-dir    
    allocate(WaveNumZ(nzc), Stat =iErr02)
    IF(BcOption(zm_dir)==BC_PERIOD ) THEN    
      WaveCoe = 2.0_RK*PI
      normTmp = 1.0_RK
    ELSE
      WaveCoe= PI
      normTmp= 2.0_RK
    ENDIF   
    do k=1,nzc
      wa3= WaveCoe*real(k-1,RK)/real(nzc,RK)
      WaveNumZ(k)=2.0_RK*rdz2*(cos(wa3)-1.0_RK)
    enddo
    normfft = normfft*normTmp*real(nzc,RK); normfft=1.0_RK/normfft
        
    best_time=max(huge(wa1),1.0E+20)
    if(nrank==0) call MainLog%OutInfo("Auto-tuning mode for Poisson Solver......",1)
    IF(BcOption(ym_dir)==BC_PERIOD) THEN
      ! Choice-1
      call Create_Poisson_FFT_Plan(x1size,z2size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1_periodic(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      ! 
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1_periodic(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !       
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-1, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=1; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=1; iFFTz=2
      endif
      call Destory_Poisson_FFT_Plan()
   
      ! Choice-2    
      call Create_Poisson_FFT_Plan(x2size,z1size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1_periodic(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1_periodic(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !       
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-2, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=2; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=2; iFFTz=2
      endif
      call Destory_Poisson_FFT_Plan()    
    ELSE
      ! Choice-1
      call Create_Poisson_FFT_Plan(x1size,z2size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      ! 
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !      
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-1, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=1; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=1; iFFTz=2
      endif
      call Destory_Poisson_FFT_Plan()
   
      ! Choice-2    
      call Create_Poisson_FFT_Plan(x2size,z1size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      ! 
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-2, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=2; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=2; iFFTz=2
      endif
      call Destory_Poisson_FFT_Plan()

      ! Choice-3
      call Create_Poisson_FFT_Plan(x1size,z2size)
      IsReduce=1
      if(z2size(2)<4)IsReduce=0
      call MPI_ALLREDUCE(IsReduce,IsReduceR,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD,iErr01)
      if(IsReduceR==1) then
        allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(z2start,z2end,'z')
        call MPI_BARRIER(MPI_COMM_WORLD,iErr01)        
        !
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
        t1(1)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_x1_reduce(prsrc,prphiHalo)
        enddo
        t2(1)=MPI_WTIME()-t1(1)
        nullify(execute_FFTW_r2r_z)
        ! 
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
        t1(2)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_x1_reduce(prsrc,prphiHalo)
        enddo
        t2(2)=MPI_WTIME()-t1(2)
        nullify(execute_FFTW_r2r_z)
        !        
        call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
        if(nrank==0)call MainLog%OutInfo("Choice-3, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
        if(best_time>t1(1)) then
          best_time=t1(1); iChoice=3; iFFTz=1
        endif
        if(best_time>t1(2)) then
          best_time=t1(2); iChoice=3; iFFTz=2
        endif
        call decomp_info_finalize(decomp_PPE)
        deallocate(decomp_PPE,a_reduce,c_reduce)
      else
        if(nrank==0)call MainLog%OutInfo("Choice-3, z2size<4, ignore",2)
      endif
      call Destory_Poisson_FFT_Plan()
        
      ! Choice-4 
      call Create_Poisson_FFT_Plan(x2size,z1size)
      IsReduce=1
      if(x2size(2)<4)IsReduce=0
      call MPI_ALLREDUCE(IsReduce,IsReduceR,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD,iErr01)
      if(IsReduceR==1) then   
        allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(x2start,x2end,'x')
        call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
        ! 
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
        t1(1)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_z1_reduce(prsrc,prphiHalo)
        enddo
        t2(1)=MPI_WTIME()-t1(1)
        nullify(execute_FFTW_r2r_z)
        ! 
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
        t1(2)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_z1_reduce(prsrc,prphiHalo)
        enddo
        t2(2)=MPI_WTIME()-t1(2)
        nullify(execute_FFTW_r2r_z)
        !
        call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
        if(nrank==0)call MainLog%OutInfo("Choice-4, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
        if(best_time>t1(1)) then
          best_time=t1(1); iChoice=4; iFFTz=1
        endif
        if(best_time>t1(2)) then
          best_time=t1(2); iChoice=4; iFFTz=2
        endif
        call decomp_info_finalize(decomp_PPE)
        deallocate(decomp_PPE,a_reduce,c_reduce)
      else
        if(nrank==0)call MainLog%OutInfo("Choice-4, x2size<4, ignore",2)
      endif
      call Destory_Poisson_FFT_Plan()
    ENDIF
    !
    if(nrank==0) then
      call MainLog%OutInfo("The best Poisson Solver choice is probably Choice-"//num2str(iChoice),2)
      call MainLog%OutInfo("Corresponding Global Data Transpose is:",2)
    endif
    if(iFFTz==1) then
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
    else
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
    endif
    if(iChoice==1) then
      IF(BcOption(ym_dir)==BC_PERIOD ) THEN
        clcPPE => clcPPE_x1_periodic      
      ELSE
        clcPPE => clcPPE_x1
      ENDIF
      call Create_Poisson_FFT_Plan(x1size,z2size)
      if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> y2 -> z2 -> x1 -> y1",3)
    elseif(iChoice==2) then
      IF(BcOption(ym_dir)==BC_PERIOD ) THEN
        clcPPE => clcPPE_z1_periodic      
      ELSE
        clcPPE => clcPPE_z1
      ENDIF
      call Create_Poisson_FFT_Plan(x2size,z1size)
      if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> y2 -> x2 -> z1 -> y1",3)
    elseif(iChoice==3) then
      clcPPE => clcPPE_x1_reduce 
      call Create_Poisson_FFT_Plan(x1size,z2size)
      if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> x1 -> y1",3)
      allocate(decomp_PPE)
      allocate(a_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
      allocate(c_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
      call Initialize_ReduceMatrix(z2start,z2end,'z')
    elseif(iChoice==4) then
      clcPPE => clcPPE_z1_reduce 
      call Create_Poisson_FFT_Plan(x2size,z1size)
      if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> z1 -> y1",3)
      allocate(decomp_PPE)
      allocate(a_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
      allocate(c_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
      call Initialize_ReduceMatrix(x2start,x2end,'x')
    endif
    if(nrank==0) print*," "
    prsrc=0.0_RK; prphiHalo=0.0_RK
  end subroutine InitPoissonSolver

  !******************************************************************
  ! Create_Poisson_FFT_Plan
  !******************************************************************
  subroutine Create_Poisson_FFT_Plan(xsizeIn,zsizeIn)
    implicit none
    integer,dimension(3),intent(in)::xsizeIn,zsizeIn
   
    ! locals
    integer::plan_type
    integer(C_FFTW_R2R_KIND)::kind_fwd,kind_bwd
    real(RK),dimension(:),allocatable::Vec1,Vec2
    
    if(FFTW_plan_type == 1) then
      plan_type=FFTW_PATIENT
    else
      plan_type=FFTW_ESTIMATE
    endif
        
    ! FFT in x
    IF(BcOption(xm_dir)==BC_PERIOD) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
    ELSE
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01    
    ENDIF
    allocate(Vec1(xsizeIn(1)),Vec2(xsizeIn(1)))
    fwd_plan_x= fftw_plan_r2r_1d(xsizeIn(1),Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_x= fftw_plan_r2r_1d(xsizeIn(1),Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
   
    ! FFT in z
    IF(BcOption(zm_dir)==BC_PERIOD ) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
    ELSE
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01    
    ENDIF

    allocate(Vec1(zsizeIn(3)),Vec2(zsizeIn(3)))
    fwd_plan_z= fftw_plan_r2r_1d(zsizeIn(3),Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_z= fftw_plan_r2r_1d(zsizeIn(3),Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
  end subroutine Create_Poisson_FFT_Plan

#undef nTime_FFT_Test
end module m_Poisson
module m_Tools
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
  implicit none
  private
  
  public:: CalcUxAver, Clc_Q_vor, Clc_lamda2
  public:: CalcMaxCFL, CheckDivergence, CalcVmax, CalcDissipationRate
  public:: InverseTridiagonal, InversePeriodicTridiagonal,InversePTriFixedCoe
contains    
   
  !******************************************************************
  ! CalcMaxCFL
  !******************************************************************     
  !  uddxmax = max{ |u1|/dx + |u2|/dy  + |u3|/dz } at cell center
  subroutine CalcMaxCFL(ux,uy,uz,uddxmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),intent(out):: uddxmax
    
    ! locals
    real(RK)::vcf
    integer::ic,jc,kc,ip,jp,kp,ierror
    
    uddxmax=0.0_RK
    do kc=y1start(3),y1end(3)
      kp=kc+1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        do ic=y1start(1),y1end(1)
          ip=ic+1
          vcf= abs(ux(ic,jc,kc)+ux(ip,jc,kc))*rdx+abs(uy(ic,jc,kc)+uy(ic,jp,kc))*rdyp(jc)+ &
               abs(uz(ic,jc,kc)+uz(ic,jc,kp))*rdz
          if(vcf>uddxmax)uddxmax=vcf
        enddo
      enddo
    enddo      
    vcf = uddxmax*0.5_RK

    call MPI_ALLREDUCE(vcf,uddxmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
  end subroutine CalcMaxCFL
    
  !******************************************************************
  ! CheckDivergence
  !******************************************************************     
  subroutine  CheckDivergence(ux,uy,uz, divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),intent(out)::divmax
    
    ! locals
    real(RK)::udiv,divmax1
    integer::ic,jc,kc,ip,jp,kp,ierror

    !  ***** compute the divg(U)
    divmax1=0.0_RK
    do kc=y1start(3),y1end(3)
      kp=kc+1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        do ic=y1start(1),y1end(1)
          ip=ic+1
          udiv= abs((ux(ip,jc,kc)-ux(ic,jc,kc))*rdx +(uy(ic,jp,kc)-uy(ic,jc,kc))*rdyp(jc)+ &
                    (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz)
          if(udiv>divmax1) divmax1=udiv
        enddo
      enddo
    enddo
    call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine CheckDivergence
  
  !******************************************************************
  ! CalcVmax
  !******************************************************************  
  function CalcVmax(ux,uy,uz) result(vmaxabs)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(3)::vmaxabs
    
    ! locals
    integer::ic,jc,kc,ierror
    real(RK)::vfm,vmaxabs1(3)

    vmaxabs=0.0_RK; vmaxabs1=0.0_RK
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          vfm=abs(ux(ic,jc,kc))
          if(vfm>vmaxabs1(1))vmaxabs1(1)=vfm
          vfm=abs(uy(ic,jc,kc))
          if(vfm>vmaxabs1(2))vmaxabs1(2)=vfm
          vfm=abs(uz(ic,jc,kc))
          if(vfm>vmaxabs1(3))vmaxabs1(3)=vfm                
        enddo
      enddo
    enddo
    call MPI_REDUCE(vmaxabs1,vmaxabs,3,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)

  end function CalcVmax

  !******************************************************************
  ! InversePeriodicTridiagonal
  !******************************************************************
  subroutine InversePeriodicTridiagonal(aj,bj,cj,fj,m,n) ! my periodic tridiagonal solvers (variable coefficients )
    implicit none
    integer,intent(in)::m,n
    real(RK),dimension(m,n),intent(in):: aj,bj,cj
    real(RK),dimension(m,n),intent(inout)::fj

    ! locals
    integer::  i,j
    real(RK):: ppj,arrmn1(m,n),arrmn2(m,n),arrmn3(m,n)

    do i=1,m                                                     
      arrmn1(i,1)= -cj(i,1)/bj(i,1)                                        
      arrmn2(i,1)= -aj(i,1)/bj(i,1)                                     
      fj(i,1)= fj(i,1)/bj(i,1)                                         
    enddo                                                          
                                                                      
    ! forward elimination sweep                                                                      
    do j=2,n-1                                                     
      do i=1,m                                                     
        ppj =1.0_RK/(bj(i,j)+ aj(i,j)*arrmn1(i,j-1))                              
        arrmn1(i,j) = -cj(i,j)*ppj                                            
        arrmn2(i,j) = -aj(i,j)*arrmn2(i,j-1)*ppj                                   
        fj(i,j) = (fj(i,j)-aj(i,j)*fj(i,j-1))*ppj                         
      enddo                                                         
    enddo                                         
                                                                      
    ! backward pass                                                                   
    do i=1,m                                                     
      arrmn2(i,n)= 1.0_RK                                                    
      arrmn3(i,n)= 0.0_RK                                                  
    enddo                                                         
    do j=n-1,1,-1                                                                                                           
      do i=1,m                                                     
        arrmn2(i,j)= arrmn2(i,j) + arrmn1(i,j)*arrmn2(i,j+1)                                 
        arrmn3(i,j)= fj(i,j) + arrmn1(i,j)*arrmn3(i,j+1)                               
      enddo                                                          
    enddo                                                          
    do i=1,m                                                     
      fj(i,n)=(fj(i,n)-cj(i,n)*arrmn3(i,1)-aj(i,n)*arrmn3(i,n-1))/(cj(i,n)*arrmn2(i,1)+aj(i,n)*arrmn2(i,n-1)+bj(i,n))           
    enddo
                                                                    
    ! backward elimination pass                                                                     
    do j=n-1,1,-1                                                    
      do i=1,m                                                     
        fj(i,j)= fj(i,n)*arrmn2(i,j)+arrmn3(i,j)                                 
      enddo                                                         
    enddo
  end subroutine InversePeriodicTridiagonal

  !******************************************************************
  ! InversePTriFixedCoe
  !******************************************************************
  subroutine InversePTriFixedCoe(ajf,bjf,cjf,fj,m,n) ! my periodic tridiagonal solvers (fixedcoefficients )
    implicit none
    integer,intent(in)::m,n
    real(RK),intent(in):: ajf,bjf,cjf
    real(RK),dimension(m,n),intent(inout)::fj

    ! locals
    integer::  i,j
    real(RK),dimension(n):: ppj,vecn1,vecn2
    real(RK),dimension(m,n)::arrmn

    vecn1(1)= -cjf/bjf                                     
    vecn2(1)= -ajf/bjf
    vecn2(n)= 1.0_RK
    do j=2,n-1                                                     
      ppj(j)   = 1.0_RK/(bjf+ ajf*vecn1(j-1))                              
      vecn1(j) = -cjf*ppj(j)                                            
      vecn2(j) = -ajf*vecn2(j-1)*ppj(j)                                                       
    enddo 
    do j=n-1,1,-1                                                     
      vecn2(j)= vecn2(j) + vecn1(j)*vecn2(j+1)                                                        
    enddo
         
    ! forward elimination sweep
    do i=1,m                     
      fj(i,1)= fj(i,1)/bjf                                         
    enddo                                                             
    do j=2,n-1                                                     
      do i=1,m                                   
        fj(i,j) = (fj(i,j)-ajf*fj(i,j-1))*ppj(j)                         
      enddo                                                         
    enddo                                         
                                                                      
    ! backward pass                                                 
    do i=1,m                                              
      arrmn(i,n)= 0.0_RK                                                    
    enddo
    do j=n-1,1,-1                                                                                                           
      do i=1,m                                                                      
        arrmn(i,j)= fj(i,j) + vecn1(j)*arrmn(i,j+1)                               
      enddo                                                          
    enddo                                                             
    do i=1,m                                                     
      fj(i,n)=(fj(i,n)-cjf*arrmn(i,1)-ajf*arrmn(i,n-1))/(cjf*vecn2(1)+ajf*vecn2(n-1)+bjf)
    enddo
                                                                    
    ! backward elimination pass                                                                     
    do j=n-1,1,-1                                                    
      do i=1,m                                                     
        fj(i,j)= fj(i,n)*vecn2(j)+arrmn(i,j)                                
      enddo                                                         
    enddo
  end subroutine InversePTriFixedCoe

  !******************************************************************
  ! InverseTridiagonal
  !******************************************************************  
  subroutine InverseTridiagonal(aj,bj,cj,fj,m,n)   ! my tridiagonal solvers (variable coefficients )
    implicit none
    integer,intent(in)::m,n
    real(RK),dimension(m,n),intent(in):: aj,bj,cj
    real(RK),dimension(m,n),intent(inout)::fj
    
    ! locals
    integer:: i,j
    real(RK),dimension(m):: vecm
    real(RK),dimension(m,n)::arrmn
      
    do i=1,m
      vecm(i)=bj(i,1)
      fj(i,1)=fj(i,1)/vecm(i)
    enddo
    do j=2,n
      do i=1,m
        arrmn(i,j)=cj(i,j-1)/vecm(i)
        vecm(i)=bj(i,j)-aj(i,j)*arrmn(i,j)
        fj(i,j)=(fj(i,j)-aj(i,j)*fj(i,j-1))/vecm(i)
      enddo
    enddo
    do j=n-1,1,-1
      do i=1,m
        fj(i,j)=fj(i,j)-arrmn(i,j+1)*fj(i,j+1)
       enddo
    enddo    
  end subroutine InverseTridiagonal

  !******************************************************************
  ! CalcUxAver
  !******************************************************************
  function CalcUxAver(ux) result(res)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux
    real(RK):: res

    ! locals
    integer:: ic,jc,kc,ierror
    real(RK):: res1

    res1=0.0_RK
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          res1=res1+ux(ic,jc,kc)*dyp(jc)
        enddo
      enddo
    enddo
    call MPI_REDUCE(res1,res,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    res=res/(real(nxc*nzc,kind=RK))/yly

  end function CalcUxAver

  !******************************************************************
  ! Clc_Q_vor
  !******************************************************************
  subroutine Clc_Q_vor(ux,uy,uz,Q_vor)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::Q_vor

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::udiv,caj,cac1,cac2,cac12
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    
    ! Q = 0.5*(ui,i *uj,j - ui,j *uj,i)
    ! Similar expression in subroutine "ClcVelStrain"
    do kc=y1start(3),y1end(3)
      kp=kc+1
      km=kc-1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=y1start(1),y1end(1)
          ip=ic+1
          im=ic-1
          dudx=  (ux(ip,jc,kc) -ux(ic,jc,kc))*rdx
          dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2    &
                +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12   &
                -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1    )    *0.25_RK
          dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *0.25_RK
      
          dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *0.25_RK
          dvdy=  (uy(ic,jp,kc) -uy(ic,jc,kc))*rdyp(jc)
          dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *0.25_RK

          dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *0.25_RK
          dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *0.25_RK
          dwdz=  (uz(ic,jc,kp) -uz(ic,jc,kc))*rdz

          udiv= dudx +dvdy +dwdz 
          Q_vor(ic,jc,kc)= 0.5_RK*(udiv*udiv -dudx*dudx -dvdy*dvdy -dwdz*dwdz) -dudy*dvdx -dudz*dwdx -dvdz*dwdy
        enddo
      enddo
    enddo
  end subroutine Clc_Q_vor

  !******************************************************************
  ! Clc_lamda2
  !******************************************************************
  subroutine Clc_lamda2(ux,uy,uz,lamda2_vor)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::lamda2_vor

    ! locals
    logical:: AllRealFlag
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,rb,rc,rd,to
    real(RK)::Ome12,Ome31,Ome23,Eps11,Eps12,Eps13,Eps22,Eps23,Eps33
    real(RK)::caj,cac1,cac2,cac12,Lam11,Lam12,Lam13,Lam22,Lam23,Lam33,Root(3)
    
    ! Similar expression in subroutine "ClcVelStrain"

    ! Ome(i,j)= 0.5*(uj,i - ui,j)
    ! Eps(i,j)= 0.5*(uj,i + ui,j)
    ! Lam(i,j)= Ome(i,k) *Ome(k,j) +Eps(i,k) *Eps(k,j)
    ! lamda2 is the second eigvalue of tensor Lam.

    ! Lam12= Ome(1,k) *Ome(k,2) +Eps(1,k) *Eps(k,2)= Eps11*Eps12 +Eps12*Eps22 +Eps13*Eps23 +Ome31*Ome23
    ! Lam13= Ome(1,k) *Ome(k,3) +Eps(1,k) *Eps(k,3)= Eps11*Eps13 +Eps12*Eps23 +Eps13*Eps33 +Ome12*Ome23
    ! Lam23= Ome(2,k) *Ome(k,3) +Eps(2,k) *Eps(k,3)= Eps12*Eps13 +Eps22*Eps23 +Eps23*Eps33 +Ome12*Ome31
    ! Lam11= Ome(1,k) *Ome(k,1) +Eps(1,k) *Eps(k,1)= Eps11*Eps11 +Eps12*Eps12 +Eps13*Eps13 -Ome12*Ome12 -Ome31*Ome31
    ! Lam22= Ome(2,k) *Ome(k,2) +Eps(2,k) *Eps(k,2)= Eps12*Eps12 +Eps22*Eps22 +Eps23*Eps23 -Ome12*Ome12 -Ome23*Ome23
    ! Lam33= Ome(3,k) *Ome(k,3) +Eps(3,k) *Eps(k,3)= Eps13*Eps13 +Eps23*Eps23 +Eps33*Eps33 -Ome31*Ome31 -Ome23*Ome23

    do kc=y1start(3),y1end(3)
      kp=kc+1
      km=kc-1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=y1start(1),y1end(1)
          ip=ic+1
          im=ic-1
          dudx=  (ux(ip,jc,kc) -ux(ic,jc,kc))*rdx
          dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2    &
                +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12   &
                -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1    )    *0.25_RK
          dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *0.25_RK
      
          dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *0.25_RK
          dvdy=  (uy(ic,jp,kc) -uy(ic,jc,kc))*rdyp(jc)
          dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *0.25_RK

          dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *0.25_RK
          dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *0.25_RK
          dwdz=  (uz(ic,jc,kp) -uz(ic,jc,kc))*rdz

          Ome12= dvdx - dudy
          Ome23= dwdy - dvdz
          Ome31= dudz - dwdx
          Eps11= dudx + dudx
          Eps12= dudy + dvdx
          Eps13= dudz + dwdx
          Eps22= dvdy + dvdy
          Eps23= dvdz + dwdy
          Eps33= dwdz + dwdz
          Lam12= Eps11*Eps12 +Eps12*Eps22 +Eps13*Eps23 +Ome31*Ome23
          Lam13= Eps11*Eps13 +Eps12*Eps23 +Eps13*Eps33 +Ome12*Ome23
          Lam23= Eps12*Eps13 +Eps22*Eps23 +Eps23*Eps33 +Ome12*Ome31
          Lam11= Eps11*Eps11 +Eps12*Eps12 +Eps13*Eps13 -Ome12*Ome12 -Ome31*Ome31
          Lam22= Eps12*Eps12 +Eps22*Eps22 +Eps23*Eps23 -Ome12*Ome12 -Ome23*Ome23
          Lam33= Eps13*Eps13 +Eps23*Eps23 +Eps33*Eps33 -Ome31*Ome31 -Ome23*Ome23

          rb= -(Lam11+Lam22+Lam33)
          rc= Lam11*Lam22+Lam22*Lam33+Lam33*Lam11-Lam12*Lam12-Lam13*Lam13-Lam23*Lam23
          rd= Lam12*Lam12*Lam33+Lam13*Lam13*Lam22+Lam23*Lam23*Lam11-Lam11*Lam22*Lam33 -2.0_RK*Lam12*Lam23*Lam13

          call CubicRoot(rb,rc,rd,Root,AllRealFlag,to)
          if(AllRealFlag) then
            lamda2_vor(ic,jc,kc)=0.25_RK*Root(2)
          else
            call MainLog%CheckForError(ErrT_Abort,"Clc_lamda2","lamda2 wrong")
          endif
        enddo
      enddo
    enddo
  end subroutine Clc_lamda2

  !******************************************************************
  ! CubicRoot
  !******************************************************************
  subroutine CubicRoot(b,c,d,Root,AllRealFlag,to)
    implicit none
    real(RK),intent(in):: b,c,d
    logical,intent(out):: AllRealFlag
    real(RK),dimension(3),intent(out)::Root
    real(RK),intent(out)::to

    ! locals
    real(RK):: SMALL=1.0E-12_RK
    real(RK):: p,q,t,rho,eta,root1,root2,root3
    real(RK):: OneThird=0.3333333333333333333333_RK

    p = c-OneThird*b*b
    q = (2.0_RK*b*b*b-9.0_RK*b*c +27.0_RK*d)/27.0_RK
    t = q*q/4.0_RK+p*p*p/27.0_RK; to=t

    !  as q^2 / 4 + p^3/27 < 0, p < 0, -p and rho > 0
    if(t<0.0_RK) then
      AllRealFlag=.true.
      t  = -OneThird*b
      rho= ((-p)**1.5_RK)/sqrt(27.0_RK)
      eta= acos(0.5_RK*q/rho)
      rho= -2.0_RK*(rho**OneThird)
      root1= cos(OneThird*eta)*rho +t
      root2= cos(OneThird*(2.0_RK*Pi +eta))*rho +t
      root3= cos(OneThird*(4.0_RK*Pi +eta))*rho +t
      if(root1<root2 .and. root1<root3) then
        if(root2<root3) then
          Root=(/root1,root2,root3/);return
        else
          Root=(/root1,root3,root2/);return
        endif
      elseif(root1>root2 .and. root1>root3) then
        if(root2<root3) then
          Root=(/root2,root3,root1/);return
        else
          Root=(/root3,root2,root1/);return
        endif
      else
        if(root2<root3) then
          Root=(/root2,root1,root3/);return
        else
          Root=(/root3,root1,root2/);return
        endif
      endif

    elseif(t<SMALL) then
      AllRealFlag=.true.  
      rho= 0.5_RK*q
      if(rho<0.0_RK) then
        rho= ((-rho)**OneThird)
      else
        rho= -rho**OneThird
      endif
      Root(1)= 2.0_RK*rho-OneThird*b
      Root(2)= -rho-OneThird*b
      Root(3)= Root(2)

    else
      AllRealFlag=.false.
      t= sqrt(t)
      rho= 0.5_RK*q+t
      eta= 0.5_RK*q-t
      if(rho<0.0_RK) then
        rho= ((-rho)**OneThird)
      else
        rho= -rho**OneThird
      endif
      if(eta<0.0_RK) then
        eta= ((-eta)**OneThird)
      else
        eta= -eta**OneThird
      endif
      Root(1)= rho + eta-OneThird*b
      Root(2)= -0.5_RK*(rho+eta)-OneThird*b
      Root(3)= (0.5_RK*sqrt(3.0))*(eta-rho)
    endif
  end subroutine CubicRoot

  !******************************************************************
  ! CalcDissipationRate
  !******************************************************************
  subroutine CalcDissipationRate(ux,uy,uz,dissp)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::dissp

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::st1,st2,st3,st4,st5,st6,caj,cac1,cac2,cac12

    ! epsilon_ij= (2*xnu*<S_ij*S_ij>)/xnu, where S_ij = 0.5*(dui/dxj + duj/dxi)
    ! Similar expression in subroutine "ClcVelStrain"
    DO kc=y1start(3),y1end(3)
      kp=kc+1
      km=kc-1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=y1start(1),y1end(1)
          ip=ic+1
          im=ic-1
          st1=(ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          st2=(uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          st3=(uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
          st4=0.125_RK*((uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx  &
                       +(ux(ip,jp,kc) +ux(ic,jp,kc))*cac2                             &
                       +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12                            &
                       -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1                             )
          st5=0.125_RK*((ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz  &  
                       +(uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx  )
          st6=0.125_RK*((uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz  &
                       +(uz(ic,jp,kp) +uz(ic,jp,kc))*cac2                             &
                       +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12                            &
                       -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1                             )
          dissp(ic,jc,kc)= 2.0_RK*(st1*st1+ st2*st2+ st3*st3)+ 4.0_RK*(st4*st4 + st5*st5 + st6*st6)
        enddo
      enddo
    ENDDO

  end subroutine CalcDissipationRate
      
end module m_Tools
module m_TScheme
  use MPI
  use m_LogInfo
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1,OutFlowInfoX,OutFlowInfoY
#if defined CFDDEM || defined CFDLPT_TwoWay
  use m_Variables,only: FpForce_x,FpForce_y,FpForce_z
#endif
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal,InversePTriFixedCoe
  implicit none
  private
  
  ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am1cForCN,ap1cForCN
  
  ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am2cForCN,ap2cForCN  
  
  ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am3cForCN,ap3cForCN
  
  procedure(),pointer,public::clcRhsX, clcRhsY, clcRhsZ
  procedure(),pointer,public::clcU1Hat,clcU2Hat,clcU3Hat
  procedure(),pointer,public::clcPrSrc,PressureUpdate
  public:: InitTimeScheme,FluidVelUpdate
contains    

#include "m_TSchemeFEXP_inc.f90"
#include "m_TSchemePIMP_inc.f90"
#include "m_TSchemeFIMP_inc.f90"
  !******************************************************************
  ! InitTimeScheme
  !******************************************************************
  subroutine InitTimeScheme(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    integer::nUnit,ierror

    ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
    allocate(ap1cForCN(1:nxc),am1cForCN(1:nxc),Stat=ierror)
    ap1cForCN=ap1c;   am1cForCN=am1c;
    if(BcOption(xm_dir)==BC_NoSlip) then
      am1cForCN(1)= 2.0_RK*am1c(1)
    elseif(BcOption(xm_dir)==BC_FreeSlip) then
      am1cForCN(1)= 0.0_RK      
    endif
    if(BcOption(xp_dir)==BC_NoSlip) then
      ap1cForCN(nxc)= 2.0_RK*ap1c(nxc)
    elseif(BcOption(xp_dir)==BC_FreeSlip) then
      ap1cForCN(nxc)= 0.0_RK    
    endif
        
    ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
    allocate(ap2cForCN(1:nyc),am2cForCN(1:nyc),Stat=ierror)
    ap2cForCN = ap2c;   am2cForCN = am2c;
    if(BcOption(ym_dir)==BC_NoSlip) then
      am2cForCN(1)= 2.0_RK*am2c(1)    
    elseif(BcOption(ym_dir)==BC_FreeSlip) then
      am2cForCN(1)= 0.0_RK    
    endif
    if(BcOption(yp_dir)==BC_NoSlip) then
      ap2cForCN(nyc)= 2.0_RK*ap2c(nyc)
    elseif(BcOption(yp_dir)==BC_FreeSlip) then
      ap2cForCN(nyc)= 0.0_RK    
    endif
    
    ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
    allocate(ap3cForCN(1:nzc),am3cForCN(1:nzc),Stat=ierror)
    ap3cForCN = ap3c;   am3cForCN = am3c;
    if(BcOption(zm_dir)==BC_NoSlip) then
      am3cForCN(1)= 2.0_RK*am3c(1)    
    elseif(BcOption(zm_dir)==BC_FreeSlip) then
      am3cForCN(1)= 0.0_RK
    endif
    if(BcOption(zp_dir)==BC_NoSlip) then
      ap3cForCN(nzc)= 2.0_RK*ap3c(nzc)
    elseif(BcOption(zp_dir)==BC_FreeSlip) then
      ap3cForCN(nzc)= 0.0_RK
    endif
    
    ! FEXP 0: full explicit
    ! PIMP 1: partial implicit, only use C-N in y-dir 
    ! FIMP 2: full implicit, use C-N in all 3 dirs.
    if( (BcOption(xm_dir)==BC_PERIOD .and. BcOption(xp_dir)/=BC_PERIOD) .or.  (BcOption(ym_dir)==BC_PERIOD .and. BcOption(yp_dir)/=BC_PERIOD) .or. &
        (BcOption(zm_dir)==BC_PERIOD .and. BcOption(zp_dir)/=BC_PERIOD)   ) then
      call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Periodic Bc Wrong ")
    endif
    if( BcOption(xm_dir)==BC_OutFlow .or. BcOption(ym_dir)==BC_OutFlow .or. BcOption(zm_dir)==BC_OutFlow &
    .or.BcOption(zp_dir)==BC_OutFlow) call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","OutFlow Is ONLY Supported in xp-dir and yp-dir")

    ! Note here, if use full implicit time scheme, the potential periodic bc should be in x-dir first,a nd then z-dir, then y-dir
    if(IsImplicit==2) then
      if(BcOption(xm_dir) /=BC_PERIOD .and. (BcOption(ym_dir)==BC_PERIOD .or. BcOption(zm_dir)==BC_PERIOD) ) then
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Bc type wrong 1 ")
      endif
      if(BcOption(zm_dir) /=BC_PERIOD .and. BcOption(ym_dir) ==BC_PERIOD ) then
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Bc type wrong 2 ")
      endif
    endif

    SELECT CASE(IsImplicit)
    CASE(0)
      clcRhsX => clcRhsX_FEXP
      clcRhsY => clcRhsY_FEXP
      clcRhsZ => clcRhsZ_FEXP
      clcU1Hat=> clcU1Hat_FEXP
      clcU2Hat=> clcU2Hat_FEXP
      clcU3Hat=> clcU3Hat_FEXP
      clcPrSrc=> clcPrSrcOther
      PressureUpdate => PressureUpdate_FEXP
    CASE(1)
      clcRhsX => clcRhsX_PIMP
      clcRhsY => clcRhsY_PIMP
      clcRhsZ => clcRhsZ_PIMP
      if(BcOption(ym_dir)==BC_PERIOD) then
        clcU1Hat=> clcU1Hat_PIMP_0
        clcU2Hat=> clcU2Hat_PIMP_0
        clcU3Hat=> clcU3Hat_PIMP_0
      else
        clcU1Hat=> clcU1Hat_PIMP
        clcU2Hat=> clcU2Hat_PIMP
        clcU3Hat=> clcU3Hat_PIMP
      endif
      clcPrSrc=> clcPrSrcOther
      PressureUpdate => PressureUpdate_PIMP
    CASE(2)
      clcRhsX => clcRhsX_FIMP
      clcRhsY => clcRhsY_FIMP
      clcRhsZ => clcRhsZ_FIMP
      if(BcOption(xm_dir) ==BC_PERIOD) then ! There can be several periodic Bcs
        if(BcOption(zm_dir) ==BC_PERIOD) then
          if(BcOption(ym_dir) ==BC_PERIOD) then !
            clcU1Hat=> clcU1Hat_FIMP_000
            clcU2Hat=> clcU2Hat_FIMP_000 
            clcU3Hat=> clcU3Hat_FIMP_000    
          else
            clcU1Hat=> clcU1Hat_FIMP_010
            clcU2Hat=> clcU2Hat_FIMP_010
            clcU3Hat=> clcU3Hat_FIMP_010
          endif
        else
          clcU1Hat=> clcU1Hat_FIMP_011
          clcU2Hat=> clcU2Hat_FIMP_011
          clcU3Hat=> clcU3Hat_FIMP_011
        endif
      else                                  ! NO periodic Bc exist
        clcU1Hat=> clcU1Hat_FIMP_111
        clcU2Hat=> clcU2Hat_FIMP_111
        clcU3Hat=> clcU3Hat_FIMP_111
      endif
      clcPrSrc=> clcPrSrc_FIMP
      PressureUpdate => PressureUpdate_FIMP
    END SELECT
  end subroutine InitTimeScheme

  !******************************************************************
  ! clcPrSrcOther
  !******************************************************************  
  subroutine clcPrSrcOther(ux,uy,uz,prsrc,pressure,divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierror
    real(RK)::sudtal,sucaj,rdiv,divmax1

    divmax1=0.0_RK
    sudtal=1.0_RK/pmAlpha
    DO kc=y1start(3),y1end(3)
       kp=kc+1
       do jc=y1start(2),y1end(2)
         jp=jc+1
         sucaj=rdyp(jc)
         do ic=y1start(1),y1end(1)
           ip=ic+1
           rdiv= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx + (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj + &
                 (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
           divmax1=max(abs(rdiv),divmax1)
           prsrc(ic,jc,kc)= sudtal * rdiv
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine clcPrSrcOther    

  !******************************************************************
  ! clcPrSrc_FIMP
  !******************************************************************  
  subroutine clcPrSrc_FIMP(ux,uy,uz,prsrc,pressure,divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierror
    real(RK)::sudtal,sucaj,rdiv,divmax1,xnuhm

    divmax1=0.0_RK
    xnuhm= -0.5_RK*xnu
    sudtal=1.0_RK/pmAlpha
    DO kc=y1start(3),y1end(3)
       kp=kc+1
       do jc=y1start(2),y1end(2)
         jp=jc+1
         sucaj=rdyp(jc)
         do ic=y1start(1),y1end(1)
           ip=ic+1
           rdiv= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx + (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj + &
                 (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
           divmax1=max(abs(rdiv),divmax1)
           prsrc(ic,jc,kc)= sudtal * rdiv
           pressure(ic,jc,kc)= pressure(ic,jc,kc)+ xnuhm*rdiv         ! new added for full implicit scheme
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine clcPrSrc_FIMP
  
  !******************************************************************
  ! FluidVelUpdate
  !******************************************************************
  subroutine FluidVelUpdate(prphiHalo,ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,im,jm,km
    real(RK)::rdxEta2,sucacEta2,rdzEta2,locphi

    rdxEta2=rdx*pmAlpha
    rdzEta2=rdz*pmAlpha
    DO kc=y1start(3),y1end(3)
      km=kc-1
      do jc=y1start(2),y1end(2)
        jm=jc-1
        sucacEta2=rdyc(jc)*pmAlpha
        do ic=y1start(1),y1end(1)
          im=ic-1
          locphi=prphiHalo(ic,jc,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc)-  (locphi-prphiHalo(im,jc,kc))*rdxEta2
          uy(ic,jc,kc)=uy(ic,jc,kc)-  (locphi-prphiHalo(ic,jm,kc))*sucacEta2
          uz(ic,jc,kc)=uz(ic,jc,kc)-  (locphi-prphiHalo(ic,jc,km))*rdzEta2                
        enddo
      enddo
    ENDDO
  end subroutine FluidVelUpdate

end module m_TScheme

  ! This file is included in the module m_TScheme

  !******************************************************************
  ! clcRhsX_FEXP
  !******************************************************************
  subroutine  clcRhsX_FEXP(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistXold
    
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,dp1ns,Forced,Forced_tot,ForcedCoe,ForcedXnu
    real(RK)::d11q1,d22q1,d33q1,dcq123,convEd1,gradp1,InterpY1,InterpY2,InterpY3,InterpY4
   
    Forced=0.0_RK 
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucaj= 0.5_RK*rdyp(jc)
        ForcedCoe=dyp(jc); ForcedXnu=ForcedCoe*xnu 
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1 
          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)) )*sucaj            
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3
                
          d11q1= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          d22q1= ap2c(jc)*ux(ic,jp,kc) +ac2c(jc)*ux(ic,jc,kc) +am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) +ac3c(kc)*ux(ic,jc,kc) +am3c(kc)*ux(ic,jc,km)
          dcq123= d11q1+d22q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
#ifdef CFDDEM
          Forced_tot= 0.5_RK*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
          Forced=Forced +ForcedXnu*dcq123 +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1) +Forced_tot
#elif CFDLPT_TwoWay
          Forced_tot= FpForce_x(ic,jc,kc)
          Forced=Forced +ForcedXnu*dcq123 +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1) +Forced_tot
#else
          Forced=Forced +ForcedXnu*dcq123
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1 +pmTheta*HistXold(ic,jc,kc) -pmAlpha*gradp1
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(Forced,Forced_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      Forced= -Forced_tot/(real(nxc*nzc,kind=RK))/yly
      dp1ns = pmGamma*(Forced-gravity(1)) +pmTheta*PrGradData(4)
      PrGradData(4) =Forced-gravity(1)
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc) +dp1ns
          enddo
        enddo
      enddo
      PrGradData(3)= dp1ns/pmAlpha
      PrGradData(1)= PrGradData(1) +dp1ns/dt
    ENDIF
  end subroutine clcRhsX_FEXP

  !******************************************************************
  ! clcRhsY_FEXP
  !******************************************************************    
  subroutine  clcRhsY_FEXP(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsY
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::hdx1,hdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq123,convEd2,gradp2,InterpY1,InterpY2    
    
    hdx1=0.5_RK*rdx
    hdx3=0.5_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucac = rdyc(jc)
        qsucac= 0.25_RK*sucac
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h21=( (InterpY1*ux(ip,jm,kc)+InterpY2*ux(ip,jc,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(InterpY1*ux(ic,jm,kc)+InterpY2*ux(ic,jc,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*hdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (InterpY1*uz(ic,jm,kp)+InterpY2*uz(ic,jc,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(InterpY1*uz(ic,jm,kc)+InterpY2*uz(ic,jc,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*hdx3
                
          d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
          d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
          dcq123= d11q2+d22q2+d33q2
       
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+xnu*dcq123+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#elif CFDLPT_TwoWay
          convEd2= -h21-h22-h23+xnu*dcq123+ gravity(2) +FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+xnu*dcq123+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_FEXP

  !******************************************************************
  ! clcRhsZ_FEXP
  !****************************************************************** 
  subroutine  clcRhsZ_FEXP(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::d11q3,d22q3,d33q3,dcq123,convEd3,gradp3
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj,InterpY1,InterpY2,InterpY3,InterpY4
    
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucaj =0.5_RK*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (InterpY3*uz(ic,jc,kc) +InterpY4*uz(ic,jp,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (InterpY1*uz(ic,jm,kc) +InterpY2*uz(ic,jc,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
          d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
          d33q3= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          dcq123= d11q3+d22q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz                
#ifdef CFDDEM
          convEd3= -h31-h32-h33+xnu*dcq123+ gravity(3)+0.5_RK*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3= -h31-h32-h33+xnu*dcq123+ gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3= -h31-h32-h33+xnu*dcq123+ gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_FEXP

  !******************************************************************
  ! clcU1Hat_FEXP
  !******************************************************************    
  subroutine clcU1Hat_FEXP(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1) 
          ux(ic,jc,kc)= ux(ic,jc,kc)+ RhsX(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FEXP
  
  !******************************************************************
  ! clcU2Hat_FEXP
  !******************************************************************  
  subroutine clcU2Hat_FEXP(uy,RhsY)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    
    ! locals
    integer::ic,jc,kc

    DO kc=y1start(3),y1end(3) 
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+ RhsY(ic,jc,kc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FEXP  
  
  !******************************************************************
  ! clcU3Hat_FEXP
  !******************************************************************  
  subroutine clcU3Hat_FEXP(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    integer::ic,jc,kc
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+ RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FEXP

  !******************************************************************
  ! PressureUpdate_FEXP
  !******************************************************************
  subroutine PressureUpdate_FEXP(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout):: pressure
    
    ! locals
    integer::ic,jc,kc

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ prphiHalo(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_FEXP

  ! This file is included in the module m_TScheme

  !******************************************************************
  ! clcRhsX_FIMP
  !******************************************************************    
  subroutine  clcRhsX_FIMP(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistXold
    
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,dp1ns,Forced,Forced_tot,ForcedCoe,ForcedXnu
    real(RK)::d11q1,d22q1,d33q1,dcq123,convEd1,gradp1,InterpY1,InterpY2,InterpY3,InterpY4
    
    Forced=0.0_RK
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        sucaj=0.5_RK*rdyp(jc)
        ForcedCoe=dyp(jc); ForcedXnu=ForcedCoe*xnu
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3            
        do ic=y1start(1),y1end(1)
          im=ic-1; ip=ic+1

          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1 
          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)) )*sucaj            
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3
                
          d11q1= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
          dcq123= d11q1+d22q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
#ifdef CFDDEM
          Forced_tot= 0.5_RK*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
          Forced=Forced +ForcedXnu*dcq123 +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + gravity(1) +Forced_tot
#elif CFDLPT_TwoWay
          Forced_tot= FpForce_x(ic,jc,kc)
          Forced=Forced +ForcedXnu*dcq123 +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + gravity(1) +Forced_tot
#else
          Forced=Forced +ForcedXnu*dcq123
          convEd1= -h11-h12-h13 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ 2.0_RK*pmBeta*dcq123
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(Forced,Forced_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      Forced= -Forced_tot/(real(nxc*nzc,kind=RK))/yly
      dp1ns = pmGamma*(Forced-gravity(1)) +pmTheta*PrGradData(4)
      PrGradData(4) =Forced-gravity(1)
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc) +dp1ns
          enddo
        enddo
      enddo
      PrGradData(3) = dp1ns/pmAlpha
      PrGradData(1) = PrGradData(1) +dp1ns/dt
    ENDIF
  end subroutine clcRhsX_FIMP

  !******************************************************************
  ! clcRhsY_FIMP
  !******************************************************************    
  subroutine  clcRhsY_FIMP(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure 

    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsY
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::hdx1,hdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq123,convEd2,gradp2,InterpY1,InterpY2    
    
    hdx1=0.5_RK*rdx
    hdx3=0.5_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        sucac = rdyc(jc)
        qsucac= 0.25_RK*sucac
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1
        do ic=y1start(1),y1end(1)
          im=ic-1; ip=ic+1

          h21=( (InterpY1*ux(ip,jm,kc)+InterpY2*ux(ip,jc,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(InterpY1*ux(ic,jm,kc)+InterpY2*ux(ic,jc,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*hdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (InterpY1*uz(ic,jm,kp)+InterpY2*uz(ic,jc,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(InterpY1*uz(ic,jm,kc)+InterpY2*uz(ic,jc,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*hdx3
                
          d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
          d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
          dcq123= d11q2+d22q2+d33q2
       
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#elif CFDLPT_TwoWay
          convEd2= -h21-h22-h23+ gravity(2) +FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ 2.0_RK*pmBeta*dcq123
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_FIMP

  !******************************************************************
  ! clcRhsZ_FIMP
  !****************************************************************** 
  subroutine  clcRhsZ_FIMP(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::d11q3,d22q3,d33q3,dcq123,convEd3,gradp3
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj,InterpY1,InterpY2,InterpY3,InterpY4
    
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucaj =0.5_RK*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (InterpY3*uz(ic,jc,kc) +InterpY4*uz(ic,jp,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (InterpY1*uz(ic,jm,kc) +InterpY2*uz(ic,jc,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
          d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
          d33q3= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          dcq123= d11q3+d22q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz                
#ifdef CFDDEM
          convEd3= -h31-h32-h33+ gravity(3)+0.5_RK*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3= -h31-h32-h33+ gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3= -h31-h32-h33+ gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ 2.0_RK*pmBeta*dcq123
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_FIMP

  !******************************************************************
  ! clcU1Hat_FIMP_000
  !******************************************************************    
  subroutine clcU1Hat_FIMP_000(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX) 

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsX) 

    ! compute dq2* sweeping in the y direction
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*2.0_RK +1.0_RK
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_000

  !******************************************************************
  ! clcU1Hat_FIMP_010
  !******************************************************************    
  subroutine clcU1Hat_FIMP_010(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX)

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,Rhsx)

    ! compute dq2* sweeping in the y direction
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_010
  
  !******************************************************************
  ! clcU1Hat_FIMP_011
  !******************************************************************    
  subroutine clcU1Hat_FIMP_011(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX) 

    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+1.0_RK
      enddo
    enddo
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsX)

    ! compute dq2* sweeping in the y direction
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_011

  !******************************************************************
  ! clcU1Hat_FIMP_111
  !******************************************************************    
  subroutine clcU1Hat_FIMP_111(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxc; rTemp=pmBeta*rdx2
      DO kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoX(4,jc,kc)
        enddo
      ENDDO      
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    do jc=1,x1size(2)
      tridmi(jc,1) = 0.0_RK
      tridci(jc,1) = 1.0_RK
      tridpi(jc,1) = 0.0_RK
    enddo
    do ic=2,x1size(1)
      do jc=1,x1size(2)  
        tridmi(jc,ic) =  mic
        tridci(jc,ic) =  cic
        tridpi(jc,ic) =  pic
      enddo
    enddo
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)
      do jc=1,x1size(2)
        tridfi(jc,1) = 0.0_RK
        do ic=2,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX)

    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+1.0_RK
      enddo
    enddo
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsX)

    ! compute dq2* sweeping in the y direction
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_111

  !******************************************************************
  ! clcU2Hat_FIMP_000
  !******************************************************************    
  subroutine clcU2Hat_FIMP_000(uy,RhsY)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsY) 

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)

    ! compute dq2* sweeping in the y direction
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*2.0_RK +1.0_RK
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsY(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU2Hat_FIMP_000

  !******************************************************************
  ! clcU2Hat_FIMP_010
  !******************************************************************    
  subroutine clcU2Hat_FIMP_010(uy,RhsY)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy   
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsY)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)
    
    ! compute dq2* sweeping in the y direction
    do ic=y1start(1),y1end(1) 
      tridpj(ic,1)=0.0_RK
      tridcj(ic,1)=1.0_RK
      tridmj(ic,1)=0.0_RK
    enddo
    do jc=2,nyc
      do ic=y1start(1),y1end(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+1.0_RK
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    DO kc=y1start(3),y1end(3) 
      do ic=y1start(1),y1end(1)
        tridfj(ic,1)=0.0_RK
      enddo
      do jc=2,nyc
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=1,nyc
         do ic=y1start(1),y1end(1) 
           uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FIMP_010

  !******************************************************************
  ! clcU2Hat_FIMP_011
  !******************************************************************    
  subroutine clcU2Hat_FIMP_011(uy,RhsY)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsY)
    
    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+1.0_RK
      enddo
    enddo
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)

    ! compute dq2* sweeping in the y direction
    do ic=y1start(1),y1end(1) 
      tridpj(ic,1)=0.0_RK
      tridcj(ic,1)=1.0_RK
      tridmj(ic,1)=0.0_RK
    enddo
    do jc=2,nyc
      do ic=y1start(1),y1end(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+1.0_RK
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    DO kc=y1start(3),y1end(3) 
      do ic=y1start(1),y1end(1)
        tridfj(ic,1)=0.0_RK
      enddo
      do jc=2,nyc
        do ic=y1start(1),y1end(1) 
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=1,nyc
        do ic=y1start(1),y1end(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FIMP_011

  !******************************************************************
  ! clcU2Hat_FIMP_111
  !******************************************************************    
  subroutine clcU2Hat_FIMP_111(uy,RhsY)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy 
    
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(x1size(2),x1size(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxc; rTemp=pmBeta*rdx2
      DO kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoX(5,jc,kc)
        enddo
      ENDDO      
    endif
            
    ! compute dq1* sweeping in the x direction
    do ic=1,x1size(1)
      do jc=1,x1size(2)  
        tridpi(jc,ic)= -pmBeta*ap1cForCN(ic)
        tridmi(jc,ic)= -pmBeta*am1cForCN(ic)
        tridci(jc,ic)= -tridpi(jc,ic)-tridmi(jc,ic)+1.0_RK
      enddo
    enddo
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsY)
    
    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+1.0_RK
      enddo
    enddo
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)

    ! compute dq2* sweeping in the y direction
    do ic=y1start(1),y1end(1) 
      tridpj(ic,1)=0.0_RK
      tridcj(ic,1)=1.0_RK
      tridmj(ic,1)=0.0_RK
    enddo
    do jc=2,nyc
      do ic=y1start(1),y1end(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+1.0_RK
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    DO kc=y1start(3),y1end(3) 
      do ic=y1start(1),y1end(1)
        tridfj(ic,1)=0.0_RK
      enddo
      do jc=2,nyc
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=1,nyc
        do ic=y1start(1),y1end(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FIMP_111

  !******************************************************************
  ! clcU3Hat_FIMP_000
  !******************************************************************    
  subroutine clcU3Hat_FIMP_000(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsZ_temp) 

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 
    call transpose_z1_to_y1(arrz1,RhsZ_temp)
    
    ! compute dq2* sweeping in the y direction
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*2.0_RK +1.0_RK
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsZ_temp(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_000

  !******************************************************************
  ! clcU3Hat_FIMP_010
  !******************************************************************    
  subroutine clcU3Hat_FIMP_010(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
      
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsZ_temp)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsZ_temp)

    ! compute dq2* sweeping in the y direction
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_010

  !******************************************************************
  ! clcU3Hat_FIMP_011
  !******************************************************************    
  subroutine clcU3Hat_FIMP_011(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*2.0_RK +1.0_RK
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsZ_temp)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    do ic=1,z1size(1)
      tridpk(ic,1) =  0.0_RK
      tridck(ic,1) =  1.0_RK
      tridmk(ic,1) =  0.0_RK
    enddo
    dO kc=2,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) =  pkc
        tridck(ic,kc) =  ckc
        tridmk(ic,kc) =  mkc
      enddo
    enddo
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      do ic=1,z1size(1)
        tridfk(ic,1) =  0.0_RK
      enddo
      dO kc=2,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk,tridck,tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsZ_temp)

    ! compute dq2* sweeping in the y direction
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_011

  !******************************************************************
  ! clcU3Hat_FIMP_111 (unfinished !!!)
  !******************************************************************    
  subroutine clcU3Hat_FIMP_111(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxc; rTemp=pmBeta*rdx2
      DO kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoX(6,jc,kc)
        enddo
      ENDDO      
    endif
        
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* sweeping in the x direction
    do jc=1,x1size(2)
      do ic=1,x1size(1)
        tridpi(jc,ic)= -pmBeta*ap1cForCN(ic)
        tridmi(jc,ic)= -pmBeta*am1cForCN(ic)
        tridci(jc,ic)= -tridpi(jc,ic)-tridmi(jc,ic)+1.0_RK
      enddo
    enddo
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsZ_temp)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*2.0_RK +1.0_RK
    do ic=1,z1size(1)
      tridpk(ic,1) =  0.0_RK
      tridck(ic,1) =  1.0_RK
      tridmk(ic,1) =  0.0_RK
    enddo
    dO kc=2,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) =  pkc
        tridck(ic,kc) =  ckc
        tridmk(ic,kc) =  mkc
      enddo
    enddo
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      do ic=1,z1size(1)
        tridfk(ic,1) =  0.0_RK
      enddo
      dO kc=2,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk,tridck,tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsZ_temp)

    ! compute dq2* sweeping in the y direction
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_111

  !******************************************************************
  ! PressureUpdate_FIMP
  !******************************************************************
  subroutine PressureUpdate_FIMP(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout):: pressure
    
    ! locals
    integer::ic,jc,kc

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ prphiHalo(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_FIMP

  ! This file is included in the module m_TScheme

  !******************************************************************
  ! clcRhsX_PIMP
  !******************************************************************    
  subroutine  clcRhsX_PIMP(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistXold
    
    ! locals    
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,dp1ns,Forced,Forced_tot,ForcedCoe,ForcedXnu
    real(RK)::d11q1,d22q1,d33q1,dcq13,convEd1,gradp1,InterpY1,InterpY2,InterpY3,InterpY4
    
    Forced=0.0_RK
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucaj=0.5_RK*rdyp(jc)
        ForcedCoe=dyp(jc); ForcedXnu=ForcedCoe*xnu 
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3            
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1
          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)) )*sucaj           
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3
                
          d11q1= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
          dcq13= d11q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
#ifdef CFDDEM
          Forced_tot= 0.5_RK*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
          Forced=Forced +ForcedXnu*(dcq13+d22q1) +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + xnu*dcq13 + gravity(1) +Forced_tot
#elif CFDLPT_TwoWay
          Forced_tot= FpForce_x(ic,jc,kc)
          Forced=Forced +ForcedXnu*(dcq13+d22q1) +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + xnu*dcq13 + gravity(1) +Forced_tot
#else
          Forced=Forced +ForcedXnu*(dcq13+d22q1)
          convEd1= -h11-h12-h13 + xnu*dcq13 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ 2.0_RK*pmBeta*d22q1
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(Forced,Forced_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      Forced= -Forced_tot/(real(nxc*nzc,kind=RK))/yly
      dp1ns = pmGamma*(Forced-gravity(1)) +pmTheta*PrGradData(4)
      PrGradData(4) =Forced-gravity(1)
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc) +dp1ns
          enddo
        enddo
      enddo
      PrGradData(3) = dp1ns/pmAlpha
      PrGradData(1) = PrGradData(1) +dp1ns/dt
    ENDIF
  end subroutine clcRhsX_PIMP

  !******************************************************************
  ! clcRhsY_PIMP
  !******************************************************************    
  subroutine  clcRhsY_PIMP(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure 
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsY
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::hdx1,hdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq13,convEd2,gradp2,InterpY1,InterpY2     
    
    hdx1=0.5_RK*rdx
    hdx3=0.5_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        sucac = rdyc(jc)
        qsucac= 0.25_RK*sucac
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1
        do ic=y1start(1),y1end(1)
          im=ic-1; ip=ic+1

          h21=( (InterpY1*ux(ip,jm,kc)+InterpY2*ux(ip,jc,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(InterpY1*ux(ic,jm,kc)+InterpY2*ux(ic,jc,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*hdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (InterpY1*uz(ic,jm,kp)+InterpY2*uz(ic,jc,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(InterpY1*uz(ic,jm,kc)+InterpY2*uz(ic,jc,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*hdx3
                
          d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
          d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
          dcq13= d11q2+d33q2
       
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+xnu*dcq13+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#elif CFDLPT_TwoWay
          convEd2= -h21-h22-h23+xnu*dcq13+ gravity(2) +FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+xnu*dcq13+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ 2.0_RK*pmBeta*d22q2
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_PIMP

  !******************************************************************
  ! clcRhsZ_PIMP
  !****************************************************************** 
  subroutine  clcRhsZ_PIMP(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::d11q3,d22q3,d33q3,dcq13,convEd3,gradp3
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj,InterpY1,InterpY2,InterpY3,InterpY4
    
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        sucaj=0.5_RK*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3
        do ic=y1start(1),y1end(1)
          im=ic-1; ip=ic+1

          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (InterpY3*uz(ic,jc,kc) +InterpY4*uz(ic,jp,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (InterpY1*uz(ic,jm,kc) +InterpY2*uz(ic,jc,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
          d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
          d33q3= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          dcq13= d11q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz                
#ifdef CFDDEM
          convEd3= -h31-h32-h33+xnu*dcq13+ gravity(3)+0.5_RK*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3= -h31-h32-h33+xnu*dcq13+ gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3= -h31-h32-h33+xnu*dcq13+ gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ 2.0_RK*pmBeta*d22q3
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_PIMP

  !******************************************************************
  ! clcU1Hat_PIMP
  !******************************************************************    
  subroutine clcU1Hat_PIMP(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
 
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif    
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc) 
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc) 
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo   
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_PIMP

  !******************************************************************
  ! clcU1Hat_PIMP_0
  !******************************************************************    
  subroutine clcU1Hat_PIMP_0(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK)::mjc,cjc,pjc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridfj
    
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*2.0_RK +1.0_RK
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_PIMP_0
  
  !******************************************************************
  ! clcU2Hat_PIMP
  !******************************************************************  
  subroutine clcU2Hat_PIMP(uy,RhsY)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
    do ic=y1start(1),y1end(1) 
      tridpj(ic,1)=0.0_RK
      tridcj(ic,1)=1.0_RK
      tridmj(ic,1)=0.0_RK
    enddo
    do jc=2,nyc
      do ic=y1start(1),y1end(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+1.0_RK
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    DO kc=y1start(3),y1end(3) 
      do ic=y1start(1),y1end(1)
        tridfj(ic,1)=0.0_RK
      enddo
      do jc=2,nyc
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=1,nyc
        do ic=y1start(1),y1end(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_PIMP 

  !******************************************************************
  ! clcU2Hat_PIMP_0
  !******************************************************************  
  subroutine clcU2Hat_PIMP_0(uy,RhsY)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mjc,cjc,pjc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridfj

    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*2.0_RK +1.0_RK
    DO kc=y1start(3),y1end(3) 
      do jc=1,nyc
        do ic=y1start(1),y1end(1) 
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
       enddo
       call InversePTriFixedCoe(mjc,cjc,pjc,tridfj,y1size(1),nyc)
       do jc=1,nyc
         do ic=y1start(1),y1end(1) 
           uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
         enddo
       enddo  
     ENDDO
  end subroutine clcU2Hat_PIMP_0
  
  !******************************************************************
  ! clcU3Hat_PIMP
  !******************************************************************  
  subroutine clcU3Hat_PIMP(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+1.0_RK
      enddo
    enddo  
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsZ(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_PIMP

  !******************************************************************
  ! clcU3Hat_PIMP_0
  !******************************************************************  
  subroutine clcU3Hat_PIMP_0(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    integer::ic,jc,kc
    real(RK):: mjc,cjc,pjc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridfj

    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*2.0_RK +1.0_RK
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsZ(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mjc,cjc,pjc,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_PIMP_0

  !******************************************************************
  ! PressureUpdate_PIMP
  !******************************************************************
  subroutine PressureUpdate_PIMP(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout):: pressure
    
    integer::ic,jc,kc,jp,jm
    real(RK)::pmBetap,pmBetac,pmBetam

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        pmBetap= -pmBeta*ap2Pr(jc)
        pmBetac= -pmBeta*ac2Pr(jc) + 1.0_RK
        pmBetam= -pmBeta*am2Pr(jc)
        do ic=y1start(1),y1end(1)
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ pmBetap*prphiHalo(ic,jp,kc)+ pmBetac*prphiHalo(ic,jc,kc)+ &
                                                  pmBetam*prphiHalo(ic,jm,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_PIMP
module m_Variables
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  implicit none
  private
  
  ! define all major arrays here 
  real(RK), public,allocatable, dimension(:,:,:) :: ux
  real(RK), public,allocatable, dimension(:,:,:) :: uy
  real(RK), public,allocatable, dimension(:,:,:) :: uz
  real(RK), public,allocatable, dimension(:,:,:) :: HistxOld
  real(RK), public,allocatable, dimension(:,:,:) :: HistyOld
  real(RK), public,allocatable, dimension(:,:,:) :: HistzOld
  real(RK), public,allocatable, dimension(:,:,:) :: pressure
  
  real(RK), public,allocatable, dimension(:,:,:) :: RealArr1
  real(RK), public,allocatable, dimension(:,:,:) :: RealArr2
  real(RK), public,allocatable, dimension(:,:,:) :: Realhalo

  real(RK), public,allocatable, dimension(:,:,:)::OutFlowInfoX,OutFlowInfoY

  type(MatBound),public:: mb1  ! matrix bound type 1
  type(HaloInfo),public:: hi1  ! matrix bound type 1
#ifdef CFDDEM
  type(MatBound),public:: mb_dist  ! matrix bound for distribution
  type(HaloInfo),public:: hi_dist  ! halo info for distribution
  real(RK), public,allocatable, dimension(:,:,:)::FpForce_x,FpForce_y,FpForce_z
#endif

#ifdef CFDACM
  type(MatBound),public:: mb_dist  ! matrix bound for distribution
  type(HaloInfo),public:: hi_dist  ! halo info for distribution
  real(RK), public,allocatable, dimension(:,:,:):: IBMArr1,IBMArr2,IBMArr3
  character,public,allocatable, dimension(:,:,:):: FluidIndicator
#endif

#ifdef CFDLPT_TwoWay
  real(RK),public,allocatable, dimension(:,:,:)::FpForce_x,FpForce_y,FpForce_z
#endif
    
  public:: AllocateVariables
contains  

  !******************************************************************
  ! InverseTridiagonal
  !****************************************************************** 
  subroutine AllocateVariables()
    implicit none
      
    ! locals
    integer::ierrTmp,ierror=0

    mb1%pencil = y_pencil  
    mb1%xme=1;  mb1%xpe=2
    mb1%yme=1;  mb1%ype=2
    mb1%zme=1;  mb1%zpe=2
    
    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------
    call myallocate(ux, mb1, opt_global=.true.)
    call myallocate(uy, mb1, opt_global=.true.)
    call myallocate(uz, mb1, opt_global=.true.)
    call myallocate(pressure, mb1, opt_global=.true.)
    call myallocate(RealHalo, mb1, opt_global=.true.)
#ifdef CFDACM
    call myallocate(IBMArr1, mb1, opt_global=.true.)
    call myallocate(IBMArr2, mb1, opt_global=.true.)
    call myallocate(IBMArr3, mb1, opt_global=.true.)
    IBMArr1=0.0_RK; IBMArr2=0.0_RK; IBMArr3=0.0_RK
#endif
#ifdef CFDLPT_TwoWay
    call myallocate(FpForce_x, mb1, opt_global=.true.)
    call myallocate(FpForce_y, mb1, opt_global=.true.)
    call myallocate(FpForce_z, mb1, opt_global=.true.)
    FpForce_x=0.0_RK; FpForce_y=0.0_RK; FpForce_z=0.0_RK
#endif
   
    !-------------------------------------------------
    ! Arrays without ghost cells
    !-------------------------------------------------
    allocate(HistxOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(HistyOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(HistzOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(RealArr1(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(RealArr2(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
#ifdef CFDACM
    allocate(FluidIndicator(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
#endif
    ux=0.0_RK;         uy=0.0_RK;        uz=0.0_RK
    HistxOld=0.0_RK;   HistyOld=0.0_RK;  HistzOld=0.0_RK
    pressure=0.0_RK;   RealArr1=0.0_RK;  RealArr2=0.0_RK;  RealHalo=0.0_RK;
     
    ! xp - outflow
    if(myProcNghBC(y_pencil,3)==BC_OutFlow) then
      allocate(OutFlowInfoX(6,y1start(2):y1end(2),y1start(3):y1end(3)),stat=ierror); OutFlowInfoX=0.0_RK
    endif
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed 1")
    
    ! yp - outflow
    if(BcOption(yp_dir)==BC_OutFlow) then
      allocate(OutFlowInfoY(6,y1start(1):y1end(1),y1start(3):y1end(3)),stat=ierror); OutFlowInfoY=0.0_RK    
    endif
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed 2")
  end subroutine AllocateVariables
    
end module m_Variables
