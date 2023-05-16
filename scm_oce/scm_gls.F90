module scm_gls

   implicit none

contains

!----------------------------------------------------------------------------------------
subroutine gls_stp( Hz      ,      &            ! Cell thickness           [m]
                    u       ,      &            ! zonal velocity [m/s]
                    v       ,      &            ! meridional velocity [m/s]
                    bvf     ,      &            ! Brunt-Vaisala frequency
                    trb     ,      &            ! GLS variables TKE + length scale
                    lmix    ,      &            ! mixing length scale [m]
                    eps     ,      &            ! TKE dissipation [m2/s3]
                    Akv     ,      &            ! Turbulent viscosity  [m2/s]
                    Akt     ,      &            ! Turbulent diffusivity [m2/s]
                    c_mu    ,      &            ! Stability function for u,v
                    c_mu_prime ,   &            ! Stability function for T,S
                    r_D     ,      &
                    sustr   ,      &            ! zonal wind stress [m2/s2 = (N/m2) / (kg/m3) ]
                    svstr   ,      &            ! meridional wind stress  [m2/s2 = (N/m2) / (kg/m3)]
                    gls_scheme ,   &
                    sfunc_opt  ,   &
                    dt      ,      &            ! Time-step
                    Zob     ,      &            ! bottom roughness length
                    Neu_bot ,      &            ! Nature of bottom boundary condition
                    nstp    ,      &            ! time n
                    nnew    ,      &            ! time n+1
                    N       ,      &            ! Number of vertical grid points        )
                    ntra    ,      &
                    ngls    ,      &
                    ntime          )
       !------------------------------------------------------------------------
       integer,                              intent(in   ) :: N,ntra,ntime,ngls
       integer,                              intent(in   ) :: nstp
       integer,                              intent(in   ) :: nnew
       integer,                              intent(in   ) :: gls_scheme
       integer,                              intent(in   ) :: sfunc_opt
       real(8),dimension( 1:N, ntime      ), intent(in   ) :: u
       real(8),dimension( 1:N, ntime      ), intent(in   ) :: v
       real(8),dimension( 0:N             ), intent(in   ) :: bvf
       real(8),dimension( 0:N             ), intent(inout) :: Akv
       real(8),dimension( 0:N, ntra       ), intent(inout) :: Akt
       real(8),dimension( 0:N, ntime, ngls), intent(inout) :: trb
       real(8),dimension( 0:N             ), intent(inout) :: lmix
       real(8),dimension( 0:N             ), intent(inout) :: eps
       real(8),dimension( 0:N             ), intent(inout) :: c_mu
       real(8),dimension( 0:N             ), intent(inout) :: c_mu_prime
! Grid variables
       real(8),dimension( 1:N      ),        intent(in   ) ::  Hz
       real(8),                              intent(in   ) ::  sustr
       real(8),                              intent(in   ) ::  svstr
       real(8),                              intent(in   ) ::  dt
       real(8),                              intent(in   ) ::  r_D,Zob
       logical,                              intent(in   ) ::  Neu_bot
       !------------------------------------------------------------------------
! local variables
       integer   ::  k,ig,ig1,ig2,igls,itke,tind
       real(8)   :: diss  (1:N-1)
       real(8)   :: shear2(1:N-1)
       real(8)   :: FC    (1:N  )
       real(8)   :: DC    (1:N-1)
       real(8)   :: CF(1:N-1)
       real(8)   :: RH(1:N-1)
       real(8)   :: rp,    rm,    rn                !<-- n,m and p exponents
       real(8)   :: beta1, beta2, beta3m, beta3p    !<-- beta terms for the psi equation
       real(8)   :: OneOverSig(2)                   !<-- inverse of Schmidt number for tke and psi
       real(8)   :: e1,e2,e3
       real(8)   :: c1   ,c2    ,c3    ,c4    ,c5    , c6
       real(8)   :: cb1   ,cb2   ,cb3   ,cb4   ,cb5   ,cbb
       real(8)   :: a1    ,a2    ,a3    ,a5    ,nn
       real(8)   :: ab1   ,ab2   ,ab3   ,ab5   ,nb
       real(8)   :: sf_d0 ,sf_d1 ,sf_d2 ,sf_d3 ,sf_d4 , sf_d5
       real(8)   :: sf_n0 ,sf_n1 ,sf_n2
       real(8)   :: sf_nb0,sf_nb1,sf_nb2
       real(8)   :: lim_am0,lim_am1,lim_am2,lim_am3,lim_am4,lim_am5,lim_am6
       real(8)   :: z0_s,ustar_sfc_sq,ustar_bot_sq,L_lim,z0_b,trb_min(2)
       real(8)   :: cff,cff1,cff2,cff3m,cff3p,lgthsc,flux_top,flux_bot,trb_sfc,trb_bot
       real(8)   :: invk       , invG       , Bprod , Sprod    , epsilon
       real(8)   :: alpha_n    , alpha_m    , c_mu_k2eps, c_mu_prime_k2eps, Denom, gls_min
       real(8)   :: alpha_n_min, alpha_m_max, cm0   , cm0inv2  , gls, du, dv
       real(8),parameter :: vonKar = 0.4
       real(8),parameter :: nuwm   =   1e-04
       real(8),parameter :: nuws   = 0.1e-04
       real(8),parameter :: eps_bvf   =   1e-14
       real(8),parameter :: eps_min   = 1.0e-12
       real(8),parameter :: tke_min   = 1.0e-06
       real(8),parameter :: galp   =  0.53
       real(8),parameter ::  chk   =  1400./9.81
       !--------------------------------------------------
       igls = 2
       itke = 1
       !--------------------------------------------------
       SELECT CASE( gls_scheme )
       CASE(1)     ! k-omega
          rp    = -1.0 ; rm    = 0.5  ; rn     = -1.0
          beta1 = 0.555; beta2 = 0.833; beta3m = -0.6; beta3p = 1.0
          OneOverSig = (/ 0.5, 0.5 /)
       CASE(2)     ! k-epsilon
          rp    = 3.0 ; rm    = 1.5 ; rn     = -1.0
          beta1 = 1.44; beta2 = 1.92; beta3m = -0.4; beta3p = 1.0
          OneOverSig = (/ 1.0, 0.7692 /)
       CASE(3) ! gen-model
          rp    = 0.0; rm    = 1.0 ; rn     = -0.67
          beta1 = 1.0; beta2 = 1.22; beta3m =  0.05; beta3p = 1.0
          OneOverSig = (/ 1.25, 0.9345 /)
       CASE DEFAULT
          print*,'Error in the definition of the closure scheme'
          stop
       END SELECT
       e1 =  3.0 + 1.*rp / rn
       e2 =  1.5 + 1.*rm / rn
       e3 = -1.0 / rn
       !--------------------------------------------------
       Call stab_func(sfunc_opt,c1,c2,c3,c4,c5,c6,cb1,cb2,cb3,cb4,cb5,cbb)
       !--------------------------------------------------
       a1 = 0.66666666667-0.5*c2;a2 = 1.-0.5*c3;a3 =1.-0.5*c4;a5 = 0.5-0.5*c6
       ab1 = 1.-cb2;ab2 = 1.-cb3;ab3 = 2.*(1.-cb4);ab5 = 2.*cbb*(1.-cb5)
       nn  = 0.5*c1;nb  = cb1
       sf_d0 = 36.0*nn*nn*nn*nb*nb
       sf_d1 = 84.0*a5*ab3*nn*nn*nb+36.0*ab5*nn*nn*nn*nb
       sf_d2 = 9.0*(ab2*ab2-ab1*ab1)*nn*nn*nn-12.0*(a2*a2-3.*a3*a3)*nn*nb*nb
       sf_d3 = 12.0*a5*ab3*(a2*ab1-3.0*a3*ab2)* nn       &
                         + 12.0*a5*ab3*(    a3*a3-a2*a2)* nb   &
                         + 12.0*   ab5*(3.0*a3*a3-a2*a2)*nn*nb
       sf_d4 = 48.0*a5*a5*ab3*ab3*nn + 36.0*a5*ab3*ab5*nn*nn
       sf_d5 = 3.0*(a2*a2-3.0*a3*a3)*(ab1*ab1-ab2*ab2)*nn
       sf_n0  = 36.0*a1*nn*nn*nb*nb
       sf_n1  = - 12.0*a5*ab3*(ab1+ab2)*nn*nn                    &
                         + 8.0*a5*ab3*(6.0*a1-a2-3.0*a3)*nn*nb   &
                         + 36.0*a1*ab5*nn*nn*Nb
       sf_n2  = 9.0*a1*(ab2*ab2-ab1*ab1)*nn*nn
       sf_nb0 = 12.0*ab3*nn*nn*nn*nb
       sf_nb1 = 12.0*a5*ab3*ab3*nn*nn
       sf_nb2 = 9.0*a1*ab3*(ab1-ab2)*nn*nn + ( 6.0*a1*(a2-3.0*a3)   &
                - 4.0*(a2*a2-3.0*a3*a3) )*ab3 * nn * nb


       lim_am0 = sf_d0*sf_n0
       lim_am1 = sf_d0*sf_n1 + sf_d1*sf_n0
       lim_am2 = sf_d1*sf_n1 + sf_d4*sf_n0
       lim_am3 = sf_d4*sf_n1
       lim_am4 = sf_d2*sf_n0
       lim_am5 = sf_d2*sf_n1+sf_d3*sf_n0
       lim_am6 = sf_d3*sf_n1
       !--------------------------------------------------
       ! Initialization of various constants
       cm0     =  ( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn) )**0.25  ! Compute cmu0
       cm0inv2 = 1./cm0**2                                               ! inverse of cmu0 squared
       ! minimum value of alpha_n to ensure that alpha_m is positive
       alpha_n_min = 0.5*( - ( sf_d1 + sf_nb0 )  + sqrt(  ( sf_d1 + sf_nb0 )**2     &
                - 4. * sf_d0 *( sf_d4 + sf_nb1 ) ) ) / ( sf_d4 + sf_nb1 )
       cff     = (cm0**3 )*(tke_min**1.5) / eps_min                      ! Compute gls_min consistently
       gls_min = (cm0**rp)*(tke_min**rm ) * ( cff**rn )                  !  with eps_min/tke_min

       trb_min(itke) = tke_min
       trb_min(igls) = gls_min

       !--------------------------------------------------
       ! Compute the vertical shear
       tind = nstp
       DO k=1,N-1
          cff = 2. / ( Hz(k ) + Hz(k+1 ) )
          du  = cff*( u(k+1, tind)-u(k, tind) )
          dv  = cff*( v(k+1, tind)-v(k, tind) )
          shear2(k) = du*du + dv*dv
       ENDDO

       !--------------------------------------------------
       ! Compute ustar squared at the surface and at the bottom
       ustar_sfc_sq = sqrt( sustr**2+svstr**2 )
       ustar_bot_sq = r_D * sqrt( u(1,tind)**2 + v(1,tind)**2  )
       ! Compute the dissipation rate
       DO k=1,N-1
          cff       = (cm0**e1) * ( trb( k,nstp,itke )**e2 )  &
                                * ( trb( k,nstp,igls )**e3 )
          diss(k)   = MAX( cff , eps_min )
       ENDDO

       !--------------------------------------------------
       DO ig = 1,ngls     ! ig = 2 for gls and = 1 for tke
       !--------------------------------------------------
          ! Off-diagonal terms for the tridiagonal problem
          cff=-0.5*dt
          DO k=2,N-1
             FC(k) = cff*OneOverSig(ig)*( Akv(k)+Akv(k-1) ) / Hz(k)
          ENDDO

          IF(Neu_bot) THEN
             FC(1) = 0.
          ELSE
             FC(1) = cff*OneOverSig(ig)*( Akv(1)+Akv(0) ) / Hz(1)
          END IF

          FC(N)=0.
          ! Production/Dissipation terms and diagonal term
          DO k=1,N-1
             ig1   = (igls-ig); ig2 = (ig-itke)  ! tke: (ig1 = 1,ig2 = 0) ; gls: (ig1 = 0,ig2 = 1)
             invk  =     1. / trb( k,nstp,itke ) ! 1/tke
             gls   =          trb( k,nstp,igls )
             ! invG = 1 for tke invg=1/psi for gls
             invG  =  ig1+ig2*(1./gls)
             cff1  =  ig1+ig2*beta1   * invk*gls
             cff2  = (ig1+ig2*beta2 ) * invk
             cff3m =  ig1+ig2*beta3m  * invk*gls
             cff3p =  ig1+ig2*beta3p  * invk*gls
             ! Shear and buoyancy production
             Sprod =  cff1*Akv(k) * shear2(k)
             Bprod = -Akt(k,1)*( cff3m*MAX(bvf(k),0.) + cff3p*MIN(bvf(k),0.) )
             ! Patankar trick to ensure non-negative solutions
             cff   =       0.5*(Hz(k)+Hz(k+1))
             IF( (Bprod + Sprod) .gt. 0.) THEN
                RH(k) = cff*( trb(k,nnew,ig) + dt*(Bprod+Sprod) )
                DC(k) = cff*(1.+dt*cff2*diss(k))-FC(k)-FC(k+1)
             ELSE
                RH(k) = cff*( trb(k,nnew,ig) + dt*       Sprod  )
                DC(k) = cff*(1.+dt*(cff2*diss(k)                    &
                                  -invG*Bprod)) - FC(k) - FC(k+1)
             ENDIF
          ENDDO

          ! Boundary conditions
          IF( ig == itke ) THEN
             ! surface
             trb_sfc      = MAX( tke_min, cm0inv2*ustar_sfc_sq )
             flux_top     = 0.
             ! bottom
             trb_bot      = MAX( tke_min, cm0inv2*ustar_bot_sq )
             flux_bot     = 0.
             ! finalize
             IF(Neu_bot) THEN
                RH(1   ) = RH(  1) + dt*flux_bot
             ELSE
                RH(1   ) = RH(  1) - FC(1)*trb_bot
             ENDIF
             RH(N-1 ) = RH(N-1) + dt*flux_top
             trb(N,nnew,ig ) = trb_sfc
             trb(0,nnew,ig ) = trb_bot
          ELSE
             ! surface
             z0_s = MAX( 1.e-2 , chk*ustar_sfc_sq )   !<-- Charnock
!                  cff     = 30.*tanh( 0.6 / (28.*sqrt( ustar_sfc_sq(i,j) )) )
!                  z0_s    = MAX( 1.e-2 ,
!     &                       1.3*( 782.353/g )*ustar_sfc_sq(i,j)*(cff**1.5) )
              cff = 0.5*( trb(N-1,nnew,itke )+trb( N  ,nnew,itke ) )
              lgthsc      = vonKar*(0.5*Hz(N)+z0_s)
              trb_sfc     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)*(cff**rm))
              flux_top    = -rn*cm0**(rp+1.)*vonKar*OneOverSig(igls)  &
                                       *(cff**(rm+0.5))*(lgthsc**rn)
              ! bottom
              z0_b        = MAX( Zob , 1.E-04 )
              cff         = 0.5*( trb(1,nnew,itke ) + trb(0,nnew,itke ) )
              lgthsc      = vonKar*(0.5*Hz(1)+z0_b)
              trb_bot     = MAX(gls_min,(cm0**rp)*(lgthsc**rn)*(cff**rm))
              flux_bot    =-rn*cm0**(rp+1.)                &
                                *vonKar*OneOverSig(igls)   &
                                *(cff**(rm+0.5))*(lgthsc**rn)

              IF( ustar_bot_sq == 0. ) THEN
                 flux_bot = 0.
                 trb_bot  = gls_min
              ENDIF
              ! finalize
              IF(Neu_bot) THEN
                 RH(  1  ) = RH(  1) + dt*flux_bot
              ELSE
                 RH(  1  ) = RH(  1) - FC(1)*trb_bot
              END IF
              RH( N-1 ) = RH(N-1) + dt*flux_top
              trb( N,nnew,ig ) = trb_sfc
              trb( 0,nnew,ig ) = trb_bot
            ENDIF

            ! tridiagonal resolution
            cff       =  1./DC(N-1)
            CF(N-1) = cff*FC(N-1)
            RH(N-1) = cff*RH(N-1)

            DO k=N-2,1,-1
               cff   =   1./(DC(k)-CF(k+1)*FC(k+1))
               CF(k) = cff*FC(k)
               RH(k) = cff*( RH(k)-FC(k+1)*RH(k+1))
            ENDDO

            trb(1,nnew,ig ) = MAX( RH(1), trb_min(ig) )

            DO k=2,N-1
               RH(k) = RH(k)-CF(k)*RH(k-1)
               trb(k,nnew,ig ) = MAX( RH(k), trb_min(ig) )
            ENDDO
         !--------------------------------------------------
         ENDDO     ! ig loop
         !--------------------------------------------------

         DO k=1,N-1
            !
            ! Galperin limitation : l <= l_lim
            L_lim = galp * sqrt( 2.* trb(k,nnew,itke)) /        &
                                  ( sqrt(max(eps_bvf, bvf(k)))  )
            !
            ! Limitation on psi (use MAX because rn is negative)
            cff = (cm0**rp) * (L_lim**rn) * (trb(k,nnew,itke)**rm)
            trb( k,nnew,igls ) = MAX( trb( k,nnew,igls ),cff )
            !
            ! Dissipation rate
            epsilon = (cm0**e1) * ( trb(k,nnew,itke )**e2 )   &
                                * ( trb(k,nnew,igls )**e3 )
            epsilon = MAX(epsilon,eps_min)
            eps(k) = epsilon
            !
            ! Compute alpha_n and alpha_m
            cff     = ( trb(k,nnew,itke)/epsilon )**2
            alpha_m     = cff*  shear2(k)
            alpha_n     = cff*     bvf(k)
            !
            ! Limitation of alpha_n and alpha_m
            alpha_n     = MIN(  MAX( 0.73*alpha_n_min , alpha_n ) , 1.0e10 )



            alpha_m_max = ( lim_am0 + lim_am1 * alpha_n                         &
                           +  lim_am2 * alpha_n**2 + lim_am3 * alpha_n**3) /   &
                      ( lim_am4 + lim_am5 * alpha_n + lim_am6 * alpha_n**2)
            alpha_m = MIN(alpha_m , alpha_m_max)
            !
            ! Compute stability functions
            Denom = sf_d0  +  sf_d1*alpha_n +  sf_d2*alpha_m   &
                 + sf_d3*alpha_n*alpha_m + sf_d4*alpha_n**2 + sf_d5*alpha_m**2
            c_mu_k2eps      = (sf_n0  +  sf_n1*alpha_n +  sf_n2*alpha_m)/Denom
            c_mu_prime_k2eps = (sf_nb0 + sf_nb1*alpha_n + sf_nb2*alpha_m)/Denom
            !
            ! Finalize the computation of Akv and Akt
            cff = trb( k,nnew,itke )**2 / epsilon
            Akv(k  )= MAX( cff*c_mu_k2eps     ,nuwm )
            Akt(k,1)= MAX( cff*c_mu_prime_k2eps,nuws )

            Akt(k,2:ntra)= Akt(k,1)
            !sh Akt(k,2)= Akt(k,1)
            lmix( k ) =  cm0**3 * cff / sqrt( trb( k,nnew,itke ) )
            c_mu(k) = c_mu_k2eps / cm0**3
            c_mu_prime(k) = c_mu_prime_k2eps / cm0**3
         ENDDO

         epsilon = (cm0**e1) * ( trb(N,nnew,itke )**e2 )   &
                             * ( trb(N,nnew,igls )**e3 )
         epsilon = MAX(epsilon,eps_min)
         eps(N)  = epsilon
         cff     = ( trb(N,nnew,itke)/epsilon )**2
         lmix(N) =  cm0**3 * cff / sqrt( trb( N,nnew,itke ) )

         epsilon = (cm0**e1) * ( trb(1,nnew,itke )**e2 )   &
                             * ( trb(1,nnew,igls )**e3 )
         epsilon = MAX(epsilon,eps_min)
         eps(0)  = epsilon
         cff     = ( trb(1,nnew,itke)/epsilon )**2
         lmix(0) =  cm0**3 * cff / sqrt( trb( 0,nnew,itke ) )

         Akv(N)   = MAX(  1.5*Akv(N-1  )-0.5*Akv(N-2  ), nuwm )
         Akv(0)   = MAX(  1.5*Akv(  1  )-0.5*Akv(  2  ), nuwm )
         Akt(N,1) = MAX(  1.5*Akt(N-1,1)-0.5*Akt(N-2,1), nuws )
         Akt(0,1) = MAX(  1.5*Akt(  1,1)-0.5*Akt(  2,1), nuws )
         Akt(N,2:ntra) = Akt(N,1)
         Akt(0,2:ntra) = Akt(0,1)

         c_mu(0) = 0.         ! Not used at 0 and N so we put them at 0
         c_mu_prime(0) = 0.
         c_mu(N) = 0.
         c_mu_prime(N) = 0.

         ! print*,'sf_d0       = ',sf_d0/sf_d0
         ! print*,'sf_d1       = ',sf_d1/sf_d0
         ! print*,'sf_d2       = ',sf_d2/sf_d0
         ! print*,'sf_d3       = ',sf_d3/sf_d0
         ! print*,'sf_d4       = ',sf_d4/sf_d0
         ! print*,'sf_d5       = ',sf_d5/sf_d0
         ! print*,'sf_n0       = ',sf_n0/sf_d0
         ! print*,'sf_n1       = ',sf_n1/sf_d0
         ! print*,'sf_n2       = ',sf_n2/sf_d0
         ! print*,'sf_nb0      = ',sf_nb0/sf_d0
         ! print*,'sf_nb1      = ',sf_nb1/sf_d0
         ! print*,'sf_nb2      = ',sf_nb2/sf_d0
         ! print*,'cm0      = ',cm0
         ! print*,'alpha_n_min      = ',alpha_n_min
         ! print*,'alpha_m_max      = ',alpha_m_max
         ! print*,'keps_old'

         return

end subroutine gls_stp


subroutine stab_func(sfunc_opt,c1,c2,c3,c4,c5,c6,cb1,cb2,cb3,cb4,cb5,cbb)
implicit none
integer,intent(in)    :: sfunc_opt
real(8),intent(out)   :: c1,c2,c3,c4,c5,c6,cb1,cb2,cb3,cb4,cb5,cbb
SELECT CASE(sfunc_opt)

CASE(1)   ! Gibson Launder 1978
c1=3.6;c2=0.8;c3=1.2;c4=1.2;c5=0.;c6=0.5
cb1=3.0;cb2=0.3333;cb3=0.333;cb4=0.;cb5=0.3333;cbb=0.8
CASE(2)   ! Mellor Yamada  1982
c1=6.;c2=0.32;c3=0.;c4=0.;c5=0.;c6=0.
cb1=3.728;cb2=0.;cb3=0.;cb4=0.;cb5=0.;cbb=0.6102
CASE(3)   ! Kantha Clayson 1994
c1=6.;c2=0.32;c3=0.;c4=0.;c5=0.;c6=0.
cb1=3.728;cb2=0.7;cb3=0.7;cb4=0.;cb5=0.2;cbb=0.6102
CASE(4)   ! Luyten & al.   1996
c1=3.;c2=0.8;c3=2.;c4=1.118;c5=0.;c6=0.5
cb1=3.;cb2=0.3333;cb3=0.3333;cb4=0.;cb5=0.3333;cbb=0.8
CASE(5)   ! Canuto & al. B 2001
c1=5.;c2=0.6983;c3=1.9664;c4=1.094;c5=0.;c6=0.495
cb1=5.6;cb2=0.6;cb3=1.;cb4=0.;cb5=0.3333;cbb=0.477
CASE(6)   ! Cheng 2002
c1=5.;c2=0.7983;c3=1.968;c4=1.136;c5=0.;c6=0.5
cb1=5.52;cb2=0.2134;cb3=0.3570;cb4=0.;cb5=0.3333;cbb=0.82
CASE DEFAULT ! Canuto A
c1=5.;c2=0.8;c3=1.968;c4=1.136;c5=0.;c6=0.4
cb1=5.95;cb2=0.6;cb3=1.;cb4=0.;cb5=0.3333;cbb=0.72
END SELECT
end subroutine stab_func

end module scm_gls
