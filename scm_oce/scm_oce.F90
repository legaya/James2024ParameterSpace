module scm_oce

   implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine obl_stp(  z_r     ,      &            ! Depth at cell centers    [m]
                     z_w     ,      &            ! Depth at cell interfaces [m]
                     Hz      ,      &            ! Cell thickness           [m]
                     unudge  ,      &            ! zonal      geostrophic current [m/s]
                     vnudge  ,      &            ! meridional geostrophic current [m/s]
                     tnudge  ,      &            ! tracer reference profile (for nudging) [Celsius/PSU]
                     u       ,      &            ! zonal velocity [m/s]
                     v       ,      &            ! meridional velocity [m/s]
                     t       ,      &            ! active tracers [Celsius/PSU]
                     turb    ,      &            ! GLS variables TKE + length scale (TKE : m2/s2)
                     lmix    ,      &            ! mixing length scale [m]
                     eps     ,      &            ! TKE dissipation [m2/s3]
                     rho0    ,      &            ! Reference constant density [kg/m3]
                     rho1    ,      &            ! Density perturbation [kg/m3]
                     Akv     ,      &            ! Turbulent viscosity  [m2/s]
                     Akt     ,      &            ! Turbulent diffusivity [m2/s]
                     wx_NL_KPP,     &            ! KPP non local flux ([C m/s] or [psu m/s])
                     c_mu    ,      &            ! Stability function for u,v (GLS), Coeff in the akv equation (TKE)
                     c_mu_prime ,   &            ! Stability function for T,S (GLS), Coeff in the akt equation (TKE)
                     r_D     ,      &            ! bottom drag (r_D = C_D |u1|)
                     sustr   ,      &            ! zonal wind stress [m2/s2 = (N/m2) / (kg/m3) ]
                     svstr   ,      &            ! meridional wind stress  [m2/s2 = (N/m2) / (kg/m3)]
                     srflx   ,      &            ! solar radiation
                     stflx1  ,      &            ! net heat flux
                     stflx2  ,      &            ! net freshwater flux
                     dtdz_bot,      &            ! Vertical derivative of tracers at the bottom (edge maintenance condition)
                     delta   ,      &            ! nudging coefficient
                     f       ,      &            ! Coriolis parameter
                     Ricr    ,      &            ! Critical Richardson number (for KPP)
                     hbls    ,      &            ! Surface boundary layer depth
                     dt      ,      &            ! Time-step
                     dpdx    ,      &            ! forcing-term for u equation
                     trb_scheme ,   &            ! Choice of turbulence scheme
                     sfunc_opt,     &            ! Choice of stability function for GLS
                     lin_eos ,      &            ! Boolean for use of linear equation of state
                     alpha   ,      &            ! Thermal expansion coefficient in linear EOS
                     T0      ,      &            ! Reference temperature in linear EOS
                     beta   ,       &            ! Haline contraction coefficient in linear EOS
                     S0      ,      &            ! Reference salinity in linear EOS
                     Zob     ,      &            ! Bottom roughness length
                     Neu_bot ,      &            ! Bottom boundary condition for GLS prognostic variables
                     EVD     ,      &            ! Activation of EVD for TKE Scheme
                     NonLocalKPP,   &            ! Activation of Non Local Term KPP
                     check_inputs,  &            ! Print values of inputs if True
                     nstp    ,      &            ! time n
                     nnew    ,      &            ! time n+1
                     N       ,      &            ! Number of vertical grid points        )
                     ntra    ,      &
                     ntime   ,      &
                     ngls    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       use scm_gls
       use scm_kpp
       use scm_tke

       implicit none
       integer,                       intent(in   ) :: N,ntra,ntime,ngls
       integer,                       intent(in   ) :: nstp
       integer,                       intent(in   ) :: nnew
       integer,                       intent(in   ) :: trb_scheme,sfunc_opt
       logical,                       intent(in   ) :: lin_eos,Neu_bot,check_inputs,EVD,NonLocalKPP

       real(8),dimension( 1:N, ntime      ), intent(inout) :: u
       real(8),dimension( 1:N, ntime      ), intent(inout) :: v
       real(8),dimension( 1:N, ntime, ntra), intent(inout) :: t
       real(8),dimension( 0:N             ), intent(inout) :: Akv
       real(8),dimension( 0:N, ntra       ), intent(inout) :: Akt
       real(8),dimension( 0:N, ntra       ), intent(inout) :: wx_NL_KPP
       real(8),dimension( 1:N             ), intent(inout) :: rho1
       real(8),dimension( 0:N, ntime, ngls), intent(inout) :: turb
       real(8),dimension( 0:N             ), intent(inout) :: lmix
       real(8),dimension( 0:N             ), intent(inout) :: eps
       real(8),dimension( 0:N             ), intent(inout) :: c_mu
       real(8),dimension( 0:N             ), intent(inout) :: c_mu_prime
       real(8),dimension( 1:N,ntra       ), intent(inout) :: delta
       real(8),dimension( 1:N       ), intent(in   ) :: unudge
       real(8),dimension( 1:N       ), intent(in   ) :: vnudge
       real(8),dimension( 1:N, ntra ), intent(in   ) :: tnudge
       real(8),dimension( ntra      ), intent(in   ) :: dtdz_bot
! Grid variables
       real(8),dimension( 1:N      ), intent(in   ) ::  z_r
       real(8),dimension( 0:N      ), intent(in   ) ::  z_w
       real(8),dimension( 1:N      ), intent(in   ) ::  Hz
       real(8),                       intent(inout) ::  hbls
       real(8),                       intent(in   ) ::  sustr
       real(8),                       intent(in   ) ::  svstr
       real(8),                       intent(in   ) ::  srflx
       real(8),                       intent(in   ) ::  stflx1,stflx2
       real(8),                       intent(in   ) ::  f,Ricr,dt,dpdx,rho0,alpha,T0,beta,S0
       real(8),                       intent(in   ) ::  r_D,zOb
! local variables
       integer                                      ::  k,itrc
       !logical                                      ::  check_inputs = .true.
       real(8)                                      ::  bvf(0:N),ghat(0:N),hbl_avg
       real(8)                                      ::  sig,cff, cff1, cff2, cff3
       real(8)                                      ::  FC(0:N), CF(0:N), Rz(1:N)
       real(8)                                      ::  swr_frac(0:N),stflx(2)
       real(8),parameter                            :: cp = 3985.0d0
       real(8), parameter :: nuwm =1.0e-4     !<-- minimum turbulent viscosity
       real(8), parameter :: nuws =0.1e-4     !<-- minimum turbulent diffusion
       real(8)  :: c1,c2,c3,c4,c5,c6,cb1,cb2,cb3,cb4,cb5,cbb
       real(8),dimension( 0:N             )  :: Akv_tmp
       real(8),dimension( 0:N, ntra       )  :: Akt_tmp

       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       IF( check_inputs ) THEN
          write(*,'(/1x,A/,/1x,A,5x,A/)') 'Vertical coordinate system:',   &
                                          'level','z_r         z_w           Hz'
          DO k=1,N
             write(*,'(I4,F12.4,F12.4,F12.4)') k, z_r(k), z_w(k), Hz(k)
          END DO
          print*,'Boundary layer depth       = ',hbls
          print*,'Time_step                  = ',dt
          print*,'Critical Richardson number = ',Ricr
          print*,'r_D                        = ',r_D
          print*,'Zob                        = ',Zob
          print*,'fcor                       = ',f
          print*,'lin_eos                    = ',lin_eos,T0,alpha
          print*,'-----------------------------------'
          print*,'stress fluxes   = ',sustr,svstr
          print*,'solar radiation = ',srflx
          print*,'net fluxes      = ',stflx1,stflx2
          print*,'ntra      = ',ntra
       END IF
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       stflx(1) = stflx1
       stflx(2) = stflx2
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       call  lmd_swfrac (N, swr_frac, Hz )

       IF(lin_eos) THEN
          call  rho_lin_eos    (N, rho1,bvf, t(:,nstp,1:2),z_r,rho0,alpha,T0,beta,S0)
       ELSE
          call  rho_eos        (N, rho1,bvf, t(:,nstp,1:2),z_r,rho0 )
       ENDIF

       SELECT CASE(trb_scheme)
       CASE(0)
          call  lmd_vmix   (N, Akv ,Akt, u(:,nstp),v(:,nstp),rho1,bvf,z_r)
          call  lmd_kpp    (N, Akv, Akt, hbls, u(:,nstp),v(:,nstp), t(:,nstp,:), bvf, &
                                 z_r,z_w,Hz,Ricr,f,sustr,svstr,srflx,stflx,rho0,swr_frac,ghat )
          !call  lmd_bkpp   (N, Akv, Akt, hbls, u(:,nstp),v(:,nstp),bvf,z_r,z_w,Hz,Ricr,f,r_D,Zob)
       CASE(4)
          call  lmd_vmix   (N, Akv ,Akt, u(:,nstp),v(:,nstp),rho1,bvf,z_r)
          call  lmd_kpp_LMD94    (N, Akv, Akt, hbls, u(:,nstp),v(:,nstp), t(:,nstp,:), bvf, rho1, &
                                 z_r,z_w,Hz,Ricr,f,sustr,svstr,srflx,stflx,rho0,swr_frac,ghat )
          !call  lmd_bkpp_LMD94(N,Akv,Akt,hbls,u(:,nstp),v(:,nstp),bvf,rho0,rho1,z_r,z_w,Hz,Ricr,f,r_D,Zob)
       CASE(5)
          call  tke_stp(Hz,z_r,u,v,bvf,turb(:,:,1),turb(:,:,2),lmix,eps,Akv,Akt,c_mu,c_mu_prime,r_D,sustr,svstr,      &
                                        dt,Zob,Neu_bot,EVD,nstp,nnew,N,ntra,ngls,ntime)
          !remark: turb(:,:,2) = k/leps = eps/ceps/k**(1/2)
          ghat(0:N) = 0.d0
       CASE DEFAULT
          call  gls_stp (Hz,u,v,bvf,turb,lmix,eps,Akv,Akt,c_mu,c_mu_prime,r_D,sustr,svstr,      &
                    trb_scheme,sfunc_opt,dt,Zob,Neu_bot,nstp,nnew,N,ntra,ngls,ntime)
          ghat(0:N) = 0.d0                                                      ! no non-local term for GLS
          !remark: turb(:,:,2) = psi
       END SELECT

       IF( .NOT. NonLocalKPP ) THEN
         ghat(0:N) = 0.d0     ! Deactivate non local term if asked
       END IF
!
! Tracer equations: vertical viscosity and solar penetration
!--------- ---------- -------- ----- --- -------- ---------
!
        do itrc=1,ntra
!-----------
          if (itrc.eq.1) then
             FC(N)=stflx(itrc)     !<-- surface heat flux (including
                                    !    latent and solar components)
             do k=N-1,1,-1
               !FC(k)=srflx*swr_frac(k)-(stflx(1)-srflx)*ghat(k)    !<-- penetration of solar heat flux, signe initial
               FC(k)=srflx*swr_frac(k)+(stflx(1)-srflx)*ghat(k)    !<-- penetration of solar heat flux, signe que je pense etre le bon
               wx_NL_KPP(k,1)=-(stflx(1)-srflx)*ghat(k)
             enddo
          else
             FC(N)=stflx(itrc)      !<-- salinity (fresh water flux)
             do k=N-1,1,-1
                 !FC(k)=-stflx(itrc)*ghat(k)     ! signe initial
                 FC(k)=+stflx(itrc)*ghat(k)    ! signe que je pense etre le bon
                 wx_NL_KPP(k,2)=-stflx(itrc)*ghat(k)
                 !FC(k)=0.D0
             enddo
          endif
!
! Bottom flux
!------
          FC(0) = nuws*dTdz_bot(itrc)  !<-- Neumann BC at the bottom
!
!  Lateral and vertical heat flux convergence due to 3D processes
!-------
          Do k=1,N
               Rz(k) = 0.
          Enddo

!
! Implicit integration for vertical diffusion
!------
          do k=1,N
            t(k,nnew,itrc)=Hz(k)*t(k,nstp,itrc)       &
                             +dt*(FC(k)-FC(k-1))      &     !++ dQs/dz
                             +dt*Hz(k)*Rz(k)                !++ R(z)
          enddo

          FC(1)=2.D0*dt*Akt(1,itrc)/(Hz(2)+Hz(1))       !--> resolve
          cff=1./(Hz(1)+FC(1))                          ! tri-diagonal
          CF(1)=cff*FC(1)                               ! system...
          t(1,nnew,itrc)=cff*t(1,nnew,itrc)

          do k=2,N-1
            FC(k)=2.D0*dt*Akt(k,itrc)/(Hz(k+1)+Hz(k))
            cff=1./(Hz(k)+FC(k)+FC(k-1)*(1.-CF(k-1)))
            CF(k)=cff*FC(k)
            t(k,nnew,itrc)=cff*( t(k,nnew,itrc) +FC(k-1)             &
                                      *t(k-1,nnew,itrc))
          enddo
          t(N,nnew,itrc)=(t(N,nnew,itrc)+FC(N-1)*t(N-1,nnew,itrc))   &
                                  /(Hz(N)+FC(N-1)*(1.-CF(N-1)))
          do k=N-1,1,-1
             t(k,nnew,itrc)=t(k,nnew,itrc)+CF(k)*t(k+1,nnew,itrc)  !<-- tracer value of implicit diffusion
          enddo
!
! Nudging (toward tnudge)
!------
          do k=N-1,1,-1
              t(k,nnew,itrc)=t(k,nnew,itrc) - dt*delta(k,itrc)*(           &
                            t(k,nnew,itrc)-tnudge(k,itrc) )
          enddo
!-----------
         enddo   ! <-- itrc
!-----------


!
! Momentum equations: Coriolis terms and vertical viscosity
!--------- ---------- -------- ----- --- -------- ---------
!


!
! Coriolis term  (forward-backward)
!------

        IF(nstp==1) THEN
           do k=1,N
              cff       = f * (   v(k,nstp) - vnudge(k) )
              u(k,nnew) =         u(k,nstp) + dt * cff
              cff       = f * (   u(k,nnew) - unudge(k) )
              v(k,nnew) = Hz(k)*( v(k,nstp) - dt * cff  )
              u(k,nnew) = Hz(k)*( u(k,nnew)             )
           enddo
        ELSE
           do k=1,N
              cff       = f * (   u(k,nstp) - unudge(k) )
              v(k,nnew) =         v(k,nstp) - dt * cff
              cff       = f * (   v(k,nnew) - vnudge(k) )
              u(k,nnew) = Hz(k)*( u(k,nstp) + dt * cff  )
              v(k,nnew) = Hz(k)*( v(k,nnew)             )
           enddo
        ENDIF

!
! Apply surface forcing
!------
        u(N,nnew)=u(N,nnew) + dt*sustr    !<-- sustr is in m2/s2 here
        v(N,nnew)=v(N,nnew) + dt*svstr
!
! Resolve tri-diagonal system
!------
        FC(1)=2.D0*dt*Akv(1)/(Hz(2)+Hz(1)) !<--     c(1)     ! system
        cff=1./(Hz(1)+FC(1)+dt*r_D)        !<-- 1 / b(1) implicit bottom drag appears here
        CF(1)=cff*FC(1)                    !<-- q(1)
        u(1,nnew)=cff*u(1,nnew)
        v(1,nnew)=cff*v(1,nnew)
        do k=2,N-1
          FC(k)=2.D0*dt*Akv(k)/(Hz(k+1)+Hz(k))
          cff=1.D0/(Hz(k)+FC(k)+FC(k-1)*(1.D0-CF(k-1)))
          CF(k)=cff*FC(k)
          u(k,nnew)=cff*(u(k,nnew)+FC(k-1)*u(k-1,nnew))
          v(k,nnew)=cff*(v(k,nnew)+FC(k-1)*v(k-1,nnew))
        enddo
        cff=1./( Hz(N) +FC(N-1)*(1.-CF(N-1)) )
        u(N,nnew)=cff*(u(N,nnew)+FC(N-1)*u(N-1,nnew))
        v(N,nnew)=cff*(v(N,nnew)+FC(N-1)*v(N-1,nnew))
!
! Finalize and apply damping term
!------

        do k=N-1,1,-1
          u(k,nnew)=u(k,nnew)+CF(k)*u(k+1,nnew) - dt*dpdx
          v(k,nnew)=v(k,nnew)+CF(k)*v(k+1,nnew)
        enddo
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Compute hbls for diagnostics
        Akv_tmp(:  ) = nuwm
        Akt_tmp(:,1) = nuws
        Akt_tmp(:,2) = nuws
        if (ntra>2) then
           Akt_tmp(:,2:ntra) = nuws
        endif
        call  lmd_kpp_LMD94    (N, Akv_tmp, Akt_tmp, hbls, u(:,nstp),v(:,nstp), t(:,nstp,:), bvf, rho1, &
                                 z_r,z_w,Hz,0.25,f,sustr,svstr,srflx,stflx,rho0,swr_frac,ghat )
!





       return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine obl_stp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











!===================================================================================================
       SUBROUTINE rho_eos (N,rho1,bvf,t,z_r,rho0)
!---------------------------------------------------------------------------------------------------
      implicit none
!
!-- Equation Of State variables to compute oceanic density ------------------------------------
!
      integer, intent(in   )                   :: N
      real(8), intent(  out)                   :: bvf (0:N)
      real(8), intent(  out)                   :: rho1(1:N)
      real(8), intent(in   )                   :: t  (1:N,2)
      real(8), intent(in   )                   :: z_r(1:N  )
      real(8), intent(in   )                   :: rho0
! local variables
      real(8)                                  :: Ts,Tt,sqrtTs
      integer                                  :: k
      real(8)                                  :: K0(N),K1(N),K2(N)
      real(8)                                  :: r0, cff
      real(8) A00, A01, A02, A03, A04, A10, A11, A12, A13
      real(8) AS0, AS1, AS2, B00, B01, B02, B03, B10, B11
      real(8) B12, BS1, E00, E01, E02, E10, E11, E12
      real(8) QR , Q01, Q02, Q03, Q04, Q05, Q10, Q11
      real(8) Q12, Q13, Q14, QS0, QS1, QS2, Q20
      real(8), parameter                       :: g=9.81
! parameter values
      parameter(A00=+19092.56 ,  A01=+209.8925   , A02=-3.041638,     &
               A03=-1.852732e-3, A04=-1.361629e-5, A10=104.4077  ,    &
               A11=-6.500517   , A12=+0.1553190  , A13=2.326469e-4,   &
               AS0=-5.587545   , AS1=+0.7390729  , AS2=-1.909078e-2,  &
               B00=+4.721788e-1, B01=+1.028859e-2, B02=-2.512549e-4,  &
               B03=-5.939910e-7, B10=-1.571896e-2, B11=-2.598241e-4,  &
               B12=+7.267926e-6, BS1=+2.042967e-3,                    &
               E00=+1.045941e-5, E01=-5.782165e-10,E02=+1.296821e-7,  &
               E10=-2.595994e-7, E11=-1.248266e-9, E12=-3.508914e-9)
      parameter(QR=+999.842594 , Q01=+6.793952e-2, Q02=-9.095290e-3,  &
               Q03=+1.001685e-4, Q04=-1.120083e-6, Q05=+6.536332e-9,  &
               Q10=+0.824493   , Q11=-4.08990e-3 , Q12=+7.64380e-5,   &
               Q13=-8.24670e-7 , Q14=+5.38750e-9 , QS0=-5.72466e-3,   &
               QS1=+1.02270e-4 , QS2=-1.65460e-6 , Q20=+4.8314e-4)
!---------------------------------------------------------------------------------------------------
      r0 = QR-1000.d0
! Compute density anomaly via Equation Of State (EOS) for seawater
!-------
      do k=1,N
!----
        Tt       = t(k,1)
        Ts       = t(k,2)
        sqrtTs   = sqrt(Ts)

        K0(k)    =    A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))            &
                + Ts*( A10+Tt*(A11+Tt*(A12+Tt*A13))+sqrtTs*(          &
                                    AS0+Tt*(AS1+Tt*AS2) ))

        K1(k)   =    B00+Tt*(B01+Tt*(B02+Tt*B03))+Ts*( B10            &
                +    Tt*(B11+Tt*B12)+sqrtTs*BS1)

        K2(k)   =    E00+Tt*(E01+Tt*E02)+Ts*(E10+Tt*(E11+Tt*E12))

        rho1(k) =   r0+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))      &
                       +Ts*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))      &
                            +sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20)
!----
    enddo
!----
    do k=1,N-1
        cff    = 1./(z_r(k+1)-z_r(k))
        bvf(k) = -cff*(g/rho0)*(rho1(k+1)-rho1(k))  ! Brunt-Vaisala frequency
    enddo
    bvf(0) = bvf(1  )
    bvf(N) = bvf(N-1)
!---------------------------------------------------------------------------------------------------
END SUBROUTINE rho_eos
!===================================================================================================


!===================================================================================================
       SUBROUTINE rho_lin_eos (N,rho1,bvf,t,z_r,rho0,Tcoef,T0,Scoef,S0)
!---------------------------------------------------------------------------------------------------
      implicit none
!
!-- Equation Of State variables to compute oceanic density ------------------------------------
!
      integer, intent(in   )                   :: N
      real(8), intent(  out)                   :: bvf (0:N)
      real(8), intent(  out)                   :: rho1(1:N)
      real(8), intent(in   )                   :: t  (1:N,2)
      real(8), intent(in   )                   :: z_r(1:N  )
      real(8), intent(in   )                   :: rho0, Tcoef, T0, Scoef, S0
! local variables
      real(8)                                  :: Ts,Tt,sqrtTs
      integer                                  :: k
      real(8)                                  :: K0(N),K1(N),K2(N)
      real(8)                                  :: r0, cff
      real(8), parameter                       :: g=9.81
!---------------------------------------------------------------------------------------------------
! Compute density anomaly via linear Equation Of State (EOS)
!-------
      do k=1,N
!----
         rho1(k)= rho0*( 1. - Tcoef*( t(k,1) - T0 ) + Scoef*( t(k,2) - S0 ))
!----
    enddo
!----
    do k=1,N-1
        cff    = 1./(z_r(k+1)-z_r(k))
        bvf(k) = -cff*(g/rho0)*(rho1(k+1)-rho1(k))  ! Brunt-Vaisala frequency
    enddo
    bvf(0) = bvf(1  )
    bvf(N) = bvf(N-1)

!---------------------------------------------------------------------------------------------------
END SUBROUTINE rho_lin_eos
!===================================================================================================








!
!===================================================================================================
subroutine lmd_swfrac(N,swr_frac,Hz)
!---------------------------------------------------------------------------------------------------
      implicit none
      integer,intent(in   )   :: N
      real(8),intent(  out)   :: swr_frac(0:N)
      real(8),intent(in   )   :: Hz      (1:N)
! local variables
      integer k, Jwt
      real(8) swdk1,swdk2,xi1,xi2
      real(8) mu1(5),mu2(5), r1(5), attn1, attn2
!
! Compute fraction of solar shortwave flux penetrating to specified
! depth due to exponential decay in Jerlov water type.
!
! output: swr_frac     shortwave (radiation) fractional decay.
!
! Reference:
! Paulson, C.A., and J.J. Simpson, 1977: Irradiance measurements
! in the upper ocean, J. Phys. Oceanogr., 7, 952-956.
!----------

      mu1(1)=0.35    ! reciprocal of the absorption coefficient
      mu1(2)=0.6     ! for each of the two solar wavelength bands
      mu1(3)=1.0     ! as a function of Jerlov water type (Paulson
      mu1(4)=1.5     ! and Simpson, 1977) [dimensioned as length,
      mu1(5)=1.4     ! meters];

      mu2(1)=23.0
      mu2(2)=20.0
      mu2(3)=17.0
      mu2(4)=14.0
      mu2(5)=7.9

      r1(1)=0.58     ! fraction of the total radiance for
      r1(2)=0.62     ! wavelength band 1 as a function of Jerlov
      r1(3)=0.67     ! water type (fraction for band 2 is always
      r1(4)=0.77     ! r2=1-r1);
      r1(5)=0.78
                     ! set Jerlov water type to assign everywhere
      Jwt=1          ! (an integer from 1 to 5).

      attn1=-1./mu1(Jwt)
      attn2=-1./mu2(Jwt)

      swdk1=r1(Jwt)               ! surface, then attenuate
      swdk2=1.-swdk1              ! them separately throughout
      swr_frac(N)=1.              ! the water column.

      do k=N,1,-1
         xi1=attn1*Hz(k)
         if (xi1 .gt. -20.) then        ! this logic to avoid
            swdk1=swdk1*exp(xi1)        ! computing exponent for
         else                           ! a very large argument
            swdk1=0.
         endif
         xi2=attn2*Hz(k)
         if (xi2 .gt. -20.) then
            swdk2=swdk2*exp(xi2)
         else
            swdk2=0.
         endif
         swr_frac(k-1)=swdk1+swdk2
      enddo
!---------------------------------------------------------------------------------------------------
end subroutine lmd_swfrac
!===================================================================================================
!





!===================================================================================================
subroutine set_forces(tsec,flxtime,sustr_m,svstr_m,srflx_m,stflx_m,ssflx_m,dqdt_m,sst_m,sss,S_sfc,rho0,nfrc, &
                      sustr,svstr,stflx1,stflx2,srflx)   !<-- Time interpolation of surface fluxes at the current time-step
!---------------------------------------------------------------------------------------------------
      implicit none
      real(8),intent(in)  :: tsec,sss,S_sfc,rho0
      integer,intent(in)  :: nfrc
      real(8),intent(in)  :: flxtime(nfrc)
      real(8),intent(in)  :: sustr_m(nfrc)
      real(8),intent(in)  :: svstr_m(nfrc)
      real(8),intent(in)  :: srflx_m(nfrc)
      real(8),intent(in)  :: stflx_m(nfrc)
      real(8),intent(in)  :: ssflx_m(nfrc)
      real(8),intent(in)  ::  dqdt_m(nfrc)
      real(8),intent(in)  ::   sst_m(nfrc)
      real(8),intent(out) :: sustr,svstr,srflx
      real(8),intent(out) :: stflx1,stflx2
      real(8)             :: td,cff2,cff,dqdt,tdays
      real(8),parameter   :: cp = 3985.0d0
      integer             :: kt,k1,k2


      ! convert tsec in tdays
      tdays = tsec / (3600.d0 * 24.d0)
      td    = mod(tdays,360.d0)

      if( td.le.flxtime(1) .or. td.ge.flxtime(12)) then
         k2 = 1
         k1 = 12
         cff  = 1./ (360.+flxtime(k2)-flxtime(k1))
         if(td.le.flxtime(k2)) cff2 = cff*(     flxtime(k2) - td)
         if(td.ge.flxtime(k1)) cff2 = cff*(360.+flxtime(k2) - td)
      else
        do kt=1,11
          if( td.ge.flxtime(kt) .and. td.lt.flxtime(kt+1) ) then
            k1   = kt
            k2   = kt+1
            cff  = 1./ (flxtime(k2)-flxtime(k1))
            cff2 = cff*(flxtime(k2) - td)
          endif
        enddo
      endif

      cff=0.01/86400.
!
! Linear interpolation
!-----
      sustr        =      (cff2*sustr_m(k1)  + (1.-cff2)*sustr_m(k2))/ rho0
      svstr        =      (cff2*svstr_m(k1)  + (1.-cff2)*svstr_m(k2))/ rho0
      srflx        =      (cff2*srflx_m(k1)  + (1.-cff2)*srflx_m(k2))/(rho0*Cp)
      stflx1     =      (cff2*stflx_m(k1)  + (1.-cff2)*stflx_m(k2))/(rho0*Cp)
      stflx2     =  cff*(cff2*ssflx_m(k1)  + (1.-cff2)*ssflx_m(k2))*S_sfc
      dqdt         =       cff2*dqdt_m (k1)  + (1.-cff2)*dqdt_m (k2)
!      sst          =       cff2*sst_m  (k1)  + (1.-cff2)*sst_m  (k2)
!
! Qcorrection for uncoupled cases
!-----

!       stflx(itemp)=stflx(itemp)+dqdt*(t(N,nstp,itemp)-sst)/(rho0*Cp)
       stflx2=stflx2+dqdt*(S_sfc-sss)*cff
!
!---------------------------------------------------------------------------------------------------
end subroutine set_forces
!===================================================================================================







end module scm_oce
