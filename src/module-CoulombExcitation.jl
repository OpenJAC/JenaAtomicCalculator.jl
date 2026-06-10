
#== September 2025
    The code is running through (together with example-Ds.jl) to compute partial and total Coulomb-excitation cross sections.
    Until now, however, the reduced matrix elements of <alphaf Jf || H_{tL} )(q) || alphai Ji> are still missing ... and, hence,
    cross sections are wrong.
    
    Below, the procedures from RATIP has been collected for stimulating the further development; in practice, however,
    the present code and data structures are rather different to keep it close to the standards in JAC.
    
    
    
   ! Define an internal structure type(coulex_transition) which stores all necessary information for an excitation 
   ! transition
   !
   type :: coulex_mline
      integer          :: L, M
      complex(kind=dp) :: me  
   end type coulex_mline
   !
   type :: coulex_qamp
      integer          :: number_of_mlines
      real(kind=dp)    :: q, w   
      complex(kind=dp) :: amplitude   
      type(coulex_mline), dimension(:), pointer :: mline  
   end type coulex_qamp
   !
   type :: coulex_subtransition
      integer          :: totalM_i, totalM_f
      integer          :: number_of_qs
      real(kind=dp)    :: cs  
      type(coulex_qamp), dimension(:), pointer :: qamp
   end type coulex_subtransition
   !
   type :: coulex_transition
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: number_of_subtransitions
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: energy, cs, alignment2, alignment4
      real(kind=dp)    :: beta, q0
      type(coulex_subtransition), dimension(:), pointer :: subtransition
   end type coulex_transition
   integer          :: number_of_transitions                 = 0, &
                       number_of_GLq_points                  = 0, &
		       number_of_beta_values
   !
   !
   ! Define global logical flags for the control of the COULEX program; the default values for these flags may be 
   ! overwritten interactively during input time
   logical, public :: coulex_apply_exp_energies          = .false.,  &
		      coulex_print_csf_scheme            = .false.,  &
		      coulex_print_cs_in_hartree         = .false.,  &
		      coulex_print_each_line             = .false.,  &
		      coulex_print_selected_trans        = .false.,  &
		      coulex_sort_transition_energy      = .false.,  &
		      coulex_write_transition_file       = .false.
   !
   subroutine coulex_multipole_amplitude(i,j,k,iq)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the excitation amplitude for the transition i, subtransition j and multipole component k by the 
   ! summation over the 'pure' multipole matrix and by using the proper weights and integrals.
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: i, j, k, iq
      !
      integer	          :: asfi, asff, l, r, rr, s, ss
      complex(kind=dp)    :: value
      !
      asfi  = transition(i)%asfi;  asff = transition(i)%asff
      value = zero
      do  r = 1,coulex%no_f
	 rr = coulex%ndx_f(r)
	 do  s = 1,coulex%no_i
	    ss = coulex%ndx_i(s)
	    value = value + asf_final%asf(asff)%eigenvector(rr) * &
		            coulex%matrix(r,s) * asf_initial%asf(asfi)%eigenvector(ss)
	 end do
      end do
      !
      transition(i)%subtransition(j)%qamp(iq)%mline(k)%me = value
      !
      print 1,    i,j,k,iq,transition(i)%subtransition(j)%qamp(iq)%mline(k)%me
    1 format("*** i,j,k,iq,transition(i)%subtransition(j)%qamp(iq)%mline(k)%me", &
             4i5,3x,2(es16.8,2x))
      !
   end subroutine coulex_multipole_amplitude
   !
   !
   subroutine coulex_print_results(stream)
   !--------------------------------------------------------------------------------------------------------------------
   ! Writes the Coulomb excitation cross sections and alignment parameters in a neat summary to stream.
   !
   ! Calls:  angular_momentum_string().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      integer             :: i, j, k
      real(kind=dp)       :: energy, wa
      character(len=5)    :: string_Mi, string_Mf
      !
    1 format(//32x,"============================================================", &
              /32x,"|  Summary of all cross sections and alignment parameters  |", &
              /32x,"============================================================"  )
    2 format(//,"     I-level-F     I--J^P--F      Transition Energy    ",                          &
                "  beta     M_I   M_F    No_Glq   cross section     Alignment_2     Alignment_4 ",  &
              /,"                                     (in ",a4,")       ",                          &  
                "                                     (in ",a4,")    ",                             &
             /,4x,132("-") )
    3 format(/,4x,i4," -",i4,3x,a4,a1,3x,a4,a1,5x,es14.7,6x, e8.2,41x,2(es14.7,2x))
    4 format(66x,a5,1x,a5,1x,i5,4x,es14.7,2x)
    5 format(66x,"Total cs:            ",es14.7)
    6 format(4x,132("-") )
      !
      write(stream,1)
      write(stream,2) trim(energy_unit), trim(coulex_cs_unit)
      !	      
      do  i = 1,number_of_transitions
         energy = convert_energy("from a.u.",transition(i)%energy)
         !!x if (energy_inverse) then
         !!x    energy = energy_factor / transition(i)%energy
         !!x else
         !!x    energy = energy_factor * transition(i)%energy
         !!x end if
         write(stream,3) transition(i)%level_i,transition(i)%level_f,                                     &
                         trim(angular_momentum_string(transition(i)%totalJ_i,4)), transition(i)%parity_i, &
                         trim(angular_momentum_string(transition(i)%totalJ_f,4)), transition(i)%parity_f, &
                         energy,transition(i)%beta, transition(i)%alignment2, transition(i)%alignment4
	 !
	 wa = zero
	 !
	 do  j = 1,transition(i)%number_of_subtransitions
	    string_Mi = adjustl(angular_momentum_string(abs(transition(i)%subtransition(j)%totalM_i)))
	    string_Mf = adjustl(angular_momentum_string(abs(transition(i)%subtransition(j)%totalM_f)))
	    if (transition(i)%subtransition(j)%totalM_i < 0) then
	       string_Mi = "-"//string_Mi
	    else
	       string_Mi = " "//string_Mi
	    end if
	    if (transition(i)%subtransition(j)%totalM_f < 0) then
	       string_Mf = "-"//string_Mf
	    else
	       string_Mf = " "//string_Mf
	    end if
	    !
	    write(stream,4) string_Mi, string_Mf, transition(i)%subtransition(j)%number_of_qs,           &
			                          transition(i)%subtransition(j)%cs * coulex_cs_factor
	    wa = wa + transition(i)%subtransition(j)%cs * coulex_cs_factor	 
	 end do
         write(stream,5) wa
	 ! 
      end do
      write(stream,6) 
      !
   end subroutine coulex_print_results
   !
   !
   subroutine coulex_print_transitions(stream)
   !--------------------------------------------------------------------------------------------------------------------
   ! Prints a neat table of all selected transitions on stream before the actual computation starts; only the quantum 
   ! numbers of the atomic states and the transition energies are displayed.
   !
   ! Calls: angular_momentum_string().
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      integer		  :: i, j, k
      real(kind=dp)	  :: energy
      character(len=5)    :: string_Mi, string_Mf
      !
      write(stream,1) number_of_transitions,trim(energy_unit)
      do  i = 1,number_of_transitions
         energy = convert_energy("from a.u.",transition(i)%energy)
	 !!x if (energy_inverse) then
         !!x    energy = energy_factor / transition(i)%energy
         !!x else
         !!x    energy = energy_factor * transition(i)%energy
         !!x end if
         write(stream,2) transition(i)%level_i,transition(i)%level_f,                                                &
                         trim(angular_momentum_string(transition(i)%totalJ_i,4)), transition(i)%parity_i,            &
                         trim(angular_momentum_string(transition(i)%totalJ_f,4)), transition(i)%parity_f,            &
                         energy,transition(i)%beta, number_of_GLq_points, transition(i)%subtransition(1)%qamp(1)%q,  &
			 transition(i)%subtransition(1)%qamp(number_of_GLq_points)%q
	 do  j = 1,transition(i)%number_of_subtransitions
	    string_Mi = adjustl(angular_momentum_string(abs(transition(i)%subtransition(j)%totalM_i)))
	    string_Mf = adjustl(angular_momentum_string(abs(transition(i)%subtransition(j)%totalM_f)))
	    if (transition(i)%subtransition(j)%totalM_i < 0) then
	       string_Mi = "-"//string_Mi
	    else
	       string_Mi = " "//string_Mi
	    end if
	    if (transition(i)%subtransition(j)%totalM_f < 0) then
	       string_Mf = "-"//string_Mf
	    else
	       string_Mf = " "//string_Mf
	    end if
	    !
	    write(stream,3) string_Mi, string_Mf,                            &
	           (transition(i)%subtransition(j)%qamp(1)%mline(k)%L,       &
		    transition(i)%subtransition(j)%qamp(1)%mline(k)%M,       &
	        k=1,transition(i)%subtransition(j)%qamp(1)%number_of_mlines)
	 end do
      end do
      write(stream,4) 
      !
    1 format( "The following ",i5," transitions are selected:",              &
           //,"     I-level-F     I--J^P--F      Transition Energy    ",        &
              " beta   No_Glq    q_0 ... q_max                  .. first line", &
            /,"                                     (in ",a4,")       ",        &  
              "   M_I   M_F        (L, M) multipole components ",               &
              "   .. subsequent lines",                                         &
            /,4x,134("-") )
    2 format(/,4x,i4," -",i4,3x,a4,a1,3x,a4,a1,5x,es14.7,4x, e8.2,i5,6x,2(es14.7,2x))
    3 format(55x,a5,1x,a5,7x,10( "(",i1,",",i2,") " ), 75x,10( "(",i1,",",i2,") " ))
    4 format(4x,134("-") )
      !
   end subroutine coulex_print_transitions
   !
   !
   subroutine coulex_pure_matrix(tr,str,ch,iq,csf_set)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates the 'pure' multipole matrix for the given configuration scheme csf_set. This pure matrix contains 
   ! already the contributions from all the individual (total) ranks as associated with the multipole pair L,M. The 
   ! first no_f CSF belong to the final-state representation and the following no_f+1,...,no_f+no_i to the 
   ! initial states. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      integer, intent(in)	  :: tr, str, ch, iq
      type(csf_basis), intent(in) :: csf_set
      !
      integer			  :: i, ia, ib, no_T_coeff, L, M, nu, parity, r, s, ss, tt, totalJ_i, totalJ_f, &
                                     totalM_i, totalM_f	   
      type(nkappa)		  :: aa, bb, cc, dd
      real(kind=dp)		  :: aweight, bweight
      type(nkappa)		  :: Ta, Tb
      !
      if (rabs_use_stop   .and.  csf_set%nocsf /= coulex%no_f+coulex%no_i) then
	 stop "coulex_pure_matrix(): program stop A."
      end if
      !
      coulex%matrix = cmplx(zero,zero)
      !
      totalJ_i = transition(tr)%totalJ_i
      totalJ_f = transition(tr)%totalJ_f
      totalM_i = transition(tr)%subtransition(str)%totalM_i
      totalM_f = transition(tr)%subtransition(str)%totalM_f
      L        = transition(tr)%subtransition(str)%qamp(iq)%mline(ch)%L
      M        = transition(tr)%subtransition(str)%qamp(iq)%mline(ch)%M
    1 format(a,6i6)
      !
      do  r = 1,coulex%No_f
	 do  s = coulex%No_f+1,coulex%No_f+coulex%No_i
	    ss = s - coulex%No_f
            !
            ! Calculate the individual contributions to the matrix in turn
	    ! (i) Direct term with total rank L
	    !
	    if (abs(Clebsch_Gordan(totalJ_i,totalM_i,L+L,M+M,totalJ_f,totalM_f)) > eps10) then
	       bweight = Clebsch_Gordan(totalJ_i,totalM_i,L+L,M+M,totalJ_f,totalM_f)   !!x / sqrt(totalJ_f+one)
	       nu      = L
	       parity  = 1;   if (mod(L+32,2) == 1) parity = - 1
	       !! call anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
	       !
	       number_of_mct_coefficients = 0
	       call mct_generate_coefficients(r,r,s,s,csf_set,parity,nu)
	       no_T_coeff = number_of_mct_coefficients
	       !
	       do  tt = 1,no_T_coeff
	          ia	= mct_list(tt)%a
	          ib       = mct_list(tt)%b
	          aweight  = mct_list(tt)%T
	          Ta	= wave_final%rwf(ia)%orbital
	          Tb	= wave_initial%rwf(ib)%orbital
	          !
	          coulex%matrix(r,ss) = coulex%matrix(r,ss) + aweight*bweight &
		                      * multipole_reduced_F1_integral(transition(tr)%subtransition(str)%qamp(iq)%q, &
		                                                      L,parity,wave_final%rwf(ia),wave_initial%rwf(ib))
	       end do
	    end if
	    !
	    ! (ii) Total rank t = L - 1
	    !
	    if (abs(Clebsch_Gordan(L+L,M+M,2,0,L+L-2,M+M)  &
                  * Clebsch_Gordan(totalJ_i,totalM_i,L+L-2,M+M,totalJ_f,totalM_f)) > eps10)   then
	       bweight = Clebsch_Gordan(L+L,M+M,2,0,L+L-2,M+M) *              &
	                 Clebsch_Gordan(totalJ_i,totalM_i,L+L-2,M+M,totalJ_f,totalM_f) &  !!x / sqrt(totalJ_f+one) 
			 * (-transition(tr)%beta)
	       nu      = L-1
	       parity  = 1;   if (mod(L+1+32,2) == 1) parity = - 1
	       !! call anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
	       !
	       number_of_mct_coefficients = 0
	       call mct_generate_coefficients(r,r,s,s,csf_set,parity,nu)
	       no_T_coeff = number_of_mct_coefficients
	       !
	       do  tt = 1,no_T_coeff
	          ia	= mct_list(tt)%a
	          ib       = mct_list(tt)%b
	          aweight  = mct_list(tt)%T
	          Ta	= wave_final%rwf(ia)%orbital
	          Tb	= wave_initial%rwf(ib)%orbital
	          !
	          coulex%matrix(r,ss) = coulex%matrix(r,ss) + cmplx(zero,-one)* aweight*bweight                    &
		              * multipole_reduced_F23_integral(transition(tr)%subtransition(str)%qamp(iq)%q,       &
		                                               L,L-1,parity,wave_final%rwf(ia),wave_initial%rwf(ib))
	       end do
	    end if
	    !
	    ! (iii) Total rank t = L
	    !
	    if (abs(Clebsch_Gordan(L+L,M+M,2,0,L+L,M+M) *                                    &
	            Clebsch_Gordan(totalJ_i,totalM_i,L+L,M+M,totalJ_f,totalM_f)) > eps10) then
	       bweight = Clebsch_Gordan(L+L,M+M,2,0,L+L,M+M) *                &
	                 Clebsch_Gordan(totalJ_i,totalM_i,L+L,M+M, totalJ_f,totalM_f)  &  !!x / sqrt(totalJ_f+one) 
			 * (-transition(tr)%beta)
	       nu      = L
	       parity  = 1;   if (mod(L+1+32,2) == 1) parity = - 1
	       !! call anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
	       !
	       number_of_mct_coefficients = 0
	       call mct_generate_coefficients(r,r,s,s,csf_set,parity,nu)
	       no_T_coeff = number_of_mct_coefficients
	       !
	       do  tt = 1,no_T_coeff
	          ia	= mct_list(tt)%a
	          ib       = mct_list(tt)%b
	          aweight  = mct_list(tt)%T
	          Ta	= wave_final%rwf(ia)%orbital
	          Tb	= wave_initial%rwf(ib)%orbital
	          !
	          coulex%matrix(r,ss) = coulex%matrix(r,ss) + cmplx(zero,-one)* aweight*bweight                  &
		              * multipole_reduced_F23_integral(transition(tr)%subtransition(str)%qamp(iq)%q,     &
		                                               L,L,parity,wave_final%rwf(ia),wave_initial%rwf(ib))
	       end do
	    end if
	    !
	    ! (iv) Total rank t = L + 1
	    !
	    if (abs(Clebsch_Gordan(L+L,M+M,2,0,L+L+2,M+M) *                                    &
	            Clebsch_Gordan(totalJ_i,totalM_i,L+L+2,M+M,totalJ_f,totalM_f)) > eps10) then
	       bweight = Clebsch_Gordan(L+L,M+M,2,0,L+L+2,M+M) *                               &
	                 Clebsch_Gordan(totalJ_i,totalM_i,L+L+2,M+M,totalJ_f,totalM_f)         & !!x / sqrt(totalJ_f+one) 
			 * (-transition(tr)%beta)
	       nu      = L+1
	       parity  = 1;   if (mod(L+1+32,2) == 1) parity = - 1
	       !! call anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
	       !
	       number_of_mct_coefficients = 0
	       call mct_generate_coefficients(r,r,s,s,csf_set,parity,nu)
	       no_T_coeff = number_of_mct_coefficients
	       !
	       do  tt = 1,no_T_coeff
	          ia	= mct_list(tt)%a
	          ib       = mct_list(tt)%b
	          aweight  = mct_list(tt)%T
	          Ta	= wave_final%rwf(ia)%orbital
	          Tb	= wave_initial%rwf(ib)%orbital
	          !
	          coulex%matrix(r,ss) = coulex%matrix(r,ss) + cmplx(zero,-one)* aweight*bweight                    &
		              * multipole_reduced_F23_integral(transition(tr)%subtransition(str)%qamp(iq)%q,       &
		                                               L,L+1,parity,wave_final%rwf(ia),wave_initial%rwf(ib))
	       end do
	    end if
	    !
	 end do
      end do
      !
   end subroutine coulex_pure_matrix
   !
   !
   subroutine coulex_transition_properties(tline)
   !--------------------------------------------------------------------------------------------------------------------
   ! Calculates all selected transition properties of transition i and its subtransitions from the amplitudes of the 
   ! individual multipole lines. 
   !--------------------------------------------------------------------------------------------------------------------
      !
      type(coulex_transition), intent(inout) :: tline
      !
      integer	       :: iq, j, k, L, M, stream, J_f, M_f
      real(kind=dp)    :: energy, acq, q, sigma, cs, wa2, wa4
      complex(kind=dp) :: M_fi_q, me
      !
      ! Calculate the q-dependent and independent amplitudes M_fi(q; ...)
      !
      cs = zero
      do  j = 1,tline%number_of_subtransitions
         sigma = zero
         do  iq = 1,tline%subtransition(j)%number_of_qs
	    M_fi_q = zero
	    do  k = 1,tline%subtransition(j)%qamp(iq)%number_of_mlines
	       L  = tline%subtransition(j)%qamp(iq)%mline(k)%L
	       M  = tline%subtransition(j)%qamp(iq)%mline(k)%M
	       me = tline%subtransition(j)%qamp(iq)%mline(k)%me
	       !
	       acq    = acos(tline%q0 / tline%subtransition(j)%qamp(iq)%q)
	       M_fi_q = M_fi_q + (cmplx(zero,one)**L) * conjg(spherical_Ylm(L,M,acq,zero)) * me
	    end do
	    q     = tline%subtransition(j)%qamp(iq)%q
	    tline%subtransition(j)%qamp(iq)%amplitude = M_fi_q
	    !
	    sigma = sigma + tline%subtransition(j)%qamp(iq)%w * abs( M_fi_q ) * abs( M_fi_q ) * q   &
	            / ((q*q - tline%beta*tline%beta * tline%q0*tline%q0)**two)
	 end do
	 tline%subtransition(j)%cs = sigma * two*pi * ((two*two*two*pi / (c*tline%beta))**two)
	 cs = cs + tline%subtransition(j)%cs			     
      end do
      !
      tline%cs = cs
      wa2 = zero;   wa4 = zero
      do  j = 1,tline%number_of_subtransitions
         J_f = tline%totalJ_f
	 M_f = tline%subtransition(j)%totalM_f
         wa2 = wa2 + (-1)**((J_f-M_f)/two) * Clebsch_Gordan(J_f, M_f, J_f, -M_f,2+2,0) * tline%subtransition(j)%cs
         wa4 = wa4 + (-1)**((J_f-M_f)/two) * Clebsch_Gordan(J_f, M_f, J_f, -M_f,4+4,0) * tline%subtransition(j)%cs
      end do
      !
      tline%alignment2 = sqrt(J_f+one) * wa2 / tline%cs
      tline%alignment4 = sqrt(J_f+one) * wa4 / tline%cs
      !
      ! Print a short summary of the line if required
      if (coulex_print_each_line) then
	 stream = 6
      13 write(stream,*) " "
	 write(stream,1) " Results for the transition ", tline%level_i," - ",tline%level_f,": ",	    &
			 trim(angular_momentum_string(tline%totalJ_i))," ",  tline%parity_i,"	----> ",    &
			 trim(angular_momentum_string(tline%totalJ_f))," ",  tline%parity_f		  
	 write(stream,*) " --------------------------------------------------------------------"
	 write(stream,*) " "
         energy = convert_energy("from a.u.", tline%energy)
	 !!x if (energy_inverse) then
	 !!x    energy = energy_factor / tline%energy
	 !!x else
	 !!x    energy = energy_factor * tline%energy
	 !!x end if
	 !
         write(stream,*) "   Transition energy             = ",energy,trim(energy_unit)
         write(stream,*) "   Ion velocity v/c              = ",tline%beta
         write(stream,*) "   Minimum momentum transfer q_0 = ",tline%q0
         write(stream,*) "   Total cross section [a.u.]    = ",tline%cs
         write(stream,*) "   Total cross section [barn]    = ",tline%cs * coulex_cs_factor
         write(stream,*) "   Alignment A_2                 = ",tline%alignment2
         write(stream,*) "   Alignment A_4                 = ",tline%alignment4
	 !
       1 format(1x,a,i5,a,i5,a,a,a,a,a,a,a,a,a,a)
       2 format(1x,a,i3,a,i3,/,3x,45("-"))
       3 format(4x,i5,")  ",es14.7,5x,2(es14.7,2x))
	 !
         do  j = 1,tline%number_of_subtransitions  
	    write(stream,*) " "
            write(stream,2) "   Subtransition:  M_i = ", tline%subtransition(j)%totalM_i,"   --->   M_f = ",  &
	                                                 tline%subtransition(j)%totalM_f
	    write(stream,*) " "
            write(stream,*) "   Cross section [a.u.]          = ", tline%subtransition(j)%cs
            write(stream,*) "                 [barn]          = ", tline%subtransition(j)%cs * coulex_cs_factor
	    !
	    write(stream,*) " "
	    write(stream,*) "       i)      q [a.u.]                amplitude"
	    write(stream,*) "   ----------------------------------------------------------"
            do  iq = 1,tline%subtransition(j)%number_of_qs
	       write(stream,3) iq, tline%subtransition(j)%qamp(iq)%q, tline%subtransition(j)%qamp(iq)%amplitude
	    end do
	 end do
	 !
	 ! Re-cycle once on stream 24
	 if (stream == 6) then
	    stream = 24;   goto 13
	 end if
      end if
      !
   end subroutine coulex_transition_properties
   !  

==# 


#######################################################################################################################
#######################################################################################################################

"""
`module  JAC.CoulombExcitation`  
... a submodel of JAC that contains all methods for computing Coulomb excitation cross sections and alignment parameters 
    for the excitation of target or projectile electrons by fast ion impact.
"""
module CoulombExcitation


using  Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings

"""
`struct  CoulombExcitation.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing Coulomb-excitation amplitudes as well as 
        partial and total cross of selected lines. All cross sections are always calculated in momentum space and by
        assuming the excitation by a proton with given ion energy. For other targets, these cross sections need to be 
        multipled with (Z_target)^2.

    + ionEnergies             ::Array{Float64,1}    ... List of ion energies [MeV/u].
    + calcAlignment           ::Bool                ... True, if alignment parameters to be calculated and false otherwise.
    + printBefore             ::Bool                ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection       ... Specifies the selected levels, if any.
    + zerosGL                 ::Int64               ... Number of Gauss-Legendre zeros in the integration over the momentum transfer.
"""
struct Settings  <:  AbstractProcessSettings 
    ionEnergies               ::Array{Float64,1} 
    calcAlignment             ::Bool 
    printBefore               ::Bool 
    lineSelection             ::LineSelection
    zerosGL                   ::Int64
end 


"""
`CoulombExcitation.Settings()`  ... constructor for the default values of Coulomb-excitation line computations.
"""
function Settings()
    Settings(Float64[], false, false, LineSelection(), 0)
end


"""
`CoulombExcitation.Settings(set::CoulombExcitation.Settings;`

        ionEnergies=..,         calcAlignment=..,           printBefore=..,         lineSelection=..,
        zerosGL=..)
                    
    ... constructor for modifying the given CoulombExcitation.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::CoulombExcitation.Settings;    
    ionEnergies::Union{Nothing,Array{Float64,1}}=nothing,           calcAlignment::Union{Nothing,Bool}=nothing,    
    printBefore::Union{Nothing,Bool}=nothing,                       lineSelection::Union{Nothing,LineSelection}=nothing,
    zerosGL::Union{Nothing,Int64}=nothing)  
    
    if  isnothing(ionEnergies)   ionEnergies          = set.ionEnergies             else  ionEnergiesx         = ionEnergies       end 
    if  isnothing(calcAlignment)   calcAlignmentx       = set.calcAlignment           else  calcAlignmentx       = calcAlignment     end 
    if  isnothing(printBefore)   printBeforex         = set.printBefore             else  printBeforex         = printBefore       end 
    if  isnothing(lineSelection)   lineSelectionx       = set.lineSelection           else  lineSelectionx       = lineSelection     end 
    if  isnothing(zerosGL)   zerosGLx             = set.zerosGL                 else  zerosGLx             = zerosGL           end 
    
    Settings( ionEnergiesx, calcAlignmentx, printBeforex, lineSelectionx, zerosGLx)
end


# `Base.show(io::IO, settings::CoulombExcitation.Settings)`  ... prepares a proper printout of the variable settings::CoulombExcitation.Settings.
function Base.show(io::IO, settings::CoulombExcitation.Settings) 
    println(io, "ionEnergies:              $(settings.ionEnergies)  ")
    println(io, "calcAlignment:            $(settings.calcAlignment)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "zerosGL:                  $(settings.zerosGL)  ")
end


"""
`struct  CoulombExcitation.Channel`  
    ... defines a type for a Coulomb-excitation channel to characterize a single q- and w-dependent (complex) amplitude.
        Each of these amplitudes are built on a number of (L, M)-dependent matrix elements which are all computed "on fly".
        
    + q              ::Float64              ... q-value of the magnetic subline.
    + w              ::Float64              ... Gauss-Legendre weight of this amplitude in the integration.
    + amplitude      ::Complex{Float64}     
        ... Coulomb-excitation amplitude K^(Coulex) (vec{q}; alpha_i J_i M_i --> alpha_f J_f M_f) associated with the given channel.
"""
struct  Channel
    q                ::Float64  
    w                ::Float64 
    amplitude        ::Complex{Float64} 
end


# `Base.show(io::IO, channel::CoulombExcitation.Channel)`  ... prepares a proper printout of the variable channel::CoulombExcitation.Channel.
function Base.show(io::IO, channel::CoulombExcitation.Channel) 
    println(io, "q:             $(channel.q)  ")
    println(io, "w:             $(channel.w)  ")
    println(io, "amplitude:     $(channel.amplitude)  ")
end


"""
`struct  CoulombExcitation.MagneticLine`  
    ... defines a type for a Coulomb-excitation magnetic lines to characterize a pair of Mi and Mf values, along with an
        (sub-) cross section and a list of channels
        Each of these amplitudes are built on a number of (L, M)-dependent matrix elements which are all computed "on fly".
        
    + Mi             ::AngularM64           ... magnetic quantum number of initial level.
    + Mf             ::AngularM64           ... magnetic quantum number of final level.
    + partialCs      ::Float64              ... partial cross section (alpha_i J_i M_i --> alpha_f J_f M_f) of this magnetic line
    + channels       ::Array{CoulombExcitation.Channel,1}  ... channels of the magnetic line.
"""
struct  MagneticLine
    Mi               ::AngularM64
    Mf               ::AngularM64
    partialCs        ::Float64
    channels         ::Array{CoulombExcitation.Channel,1}
end


# `Base.show(io::IO, channel::CoulombExcitation.MagneticLine)`  ... prepares a proper printout of the variable channel::CoulombExcitation.MagneticLine.
function Base.show(io::IO, mLine::CoulombExcitation.MagneticLine) 
    println(io, "Mi:            $(mLine.Mi)  ")
    println(io, "Mf:            $(mLine.Mf)  ")
    println(io, "partialCs:     $(mLine.partialCs)  ")
    println(io, "channels:      $(mLine.channels)  ")
end


"""
`struct  CoulombExcitation.Line`  ... defines a type for a Coulomb-excitation line that may include the definition of channels.

    + initialLevel   ::Level                  ... initial-(state) level
    + finalLevel     ::Level                  ... final-(state) level
    + ionEnergy      ::Float64                ... ion energy [MeV/u].
    + q0             ::Float64                ... minimum momentum transfer q0 that, for a given transition, is equivalent to ionEnergy.
    + totalCs        ::Float64                ... total cross section (alpha_i J_i --> alpha_f J_f) for this line.
    + alignmentA2    ::Float64                ... Alignment A_2 parameter.
    + alignmentA4    ::Float64                ... Alignment A_4 parameter.
    + mLines         ::Array{MagneticLine,1}  ... List of CoulombExcitation.MagneticLine's of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    ionEnergy        ::Float64
    q0               ::Float64
    totalCs          ::Float64
    alignmentA2      ::Float64
    alignmentA4      ::Float64
    mLines           ::Array{CoulombExcitation.MagneticLine,1}
end


# `Base.show(io::IO, line::CoulombExcitation.Line)`  ... prepares a proper printout of the variable line::CoulombExcitation.Line.
function Base.show(io::IO, line::CoulombExcitation.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "ionEnergy:         $(line.ionEnergy)  ")
    println(io, "q0:                $(line.q0)  ")
    println(io, "totalCs:           $(line.totalCs)  ")
    println(io, "alignmentA2:       $(line.alignmentA2)  ")
    println(io, "alignmentA4:       $(line.alignmentA4)  ")
    println(io, "mLines:            $(line.mLines)  ")
end

#######################################################################################################################
#######################################################################################################################


"""
`CoulombExcitation.betaProjectile(ionEnergy)`  
    ... returns the (projectile) beta = v/c for ions of given ionEnergy [MeV/u]; a beta::Float64 is returned.
"""
function  betaProjectile(ionEnergy)
    gamma  = 1.0 + ionEnergy / 938.272;    beta = sqrt(1.0 - 1.0/gamma^2)
    return( beta )
end


"""
`CoulombExcitation.computeKjYme(finalLevel::Level, L::AngularJ64, initialLevel::Level)`  
    ... computes the (many-electron) matrix elements <alpha_f J_f || K^(L, Coulex:jY) || alpha_f J_f>; a me::ComplexF64 is returned.
"""
function  computeKjYme(finalLevel::Level, L::AngularJ64, initialLevel::Level)
    @warn("Not yet implemented: computeKjYme")
    me = ComplexF64(1.0)
    
    return( me )
end


"""
`CoulombExcitation.computeKjTme(finalLevel::Level, t::AngularJ64, initialLevel::Level)`  
    ... computes the (many-electron) matrix elements <alpha_f J_f || K^(t, Coulex:jT) || alpha_f J_f>; a me::ComplexF64 is returned.
"""
function  computeKjTme(finalLevel::Level, t::AngularJ64, initialLevel::Level)
    @warn("Not yet implemented: computeKjTme")
    me = ComplexF64(2.0)
    
    return( me )
end


"""
`CoulombExcitation.computeAmplitude(channel::CoulombExcitation.Channel, Mi::AngularM64, Mf::AngularM64, 
                                    line::CoulombExcitation.Line, grid::Radial.Grid, 
                                    settings::CoulombExcitation.Settings; printout::Bool=true)`  
    ... to compute the amplitudes K^(Coulex) (vec{q}; alpha_i J_i M_i --> alpha_f J_f M_f);
        an amplitude::ComplexF64 is returned.
"""
function  computeAmplitude(channel::CoulombExcitation.Channel, Mi::AngularM64, Mf::AngularM64, 
                                    line::CoulombExcitation.Line, grid::Radial.Grid, 
                                    settings::CoulombExcitation.Settings; printout::Bool=true)
    Ji  = line.initialLevel.J;        Jf  = line.finalLevel.J;    amplitude = ComplexF64(0.)
    Jix = AngularMomentum.oneJ(Ji);   Jfx = AngularMomentum.oneJ(Jf)
    Mix = AngularMomentum.oneM(Mi);   Mfx = AngularMomentum.oneM(Mf)
    
    beta = CoulombExcitation.betaProjectile(line.ionEnergy)
    
    # Simply sum over all summations; determine the range of t and L, when it becomes relevant
    ts = AngularMomentum.j_values(Ji, Jf)
    for  t  in  ts
        tx = AngularMomentum.oneJ(t)
        wa = AngularMomentum.ClebschGordan(Jix, Mix, tx, Mfx-Mix, Jfx, Mfx) / sqrt(2*Jfx + 1.0)
        Ls = AngularMomentum.j_values(t, AngularJ64(1))
        wb = ComplexF64(0.)
        for  L  in Ls
            Lx = AngularMomentum.oneJ(L);   Lint = Int64(Lx);   Mint = Int64(Mfx-Mix);   M = Int64( min(Lx, tx) )
            wc = ComplexF64(0.)
            if  t == L   wc = wc + CoulombExcitation.computeKjYme(line.finalLevel, L, line.initialLevel)   end 
            wc = wc - beta * AngularMomentum.ClebschGordan(Lx, 0., 1., 0., tx, 0.) *
                      CoulombExcitation.computeKjTme(line.finalLevel, t, line.initialLevel)
            if  abs(Mint) > Lint   continue   end
            wc = im^Lx * conj( AngularMomentum.sphericalYlm(Lint, Mint, acos(line.q0/channel.q), 0.) ) * wc
            wb = wb + wc
        end
        amplitude = amplitude + wa * wb
    end
    
    return( amplitude )
end


"""
`CoulombExcitation.computeAmplitudesProperties(mLine::CoulombExcitation.MagneticLine, line::CoulombExcitation.Line,
                                               grid::Radial.Grid, settings::CoulombExcitation.Settings; printout::Bool=true)`  
    ... to compute all amplitudes and properties of the given magnetic line by using various parameters from line;
        a mline::CoulombExcitation.MagneticLine is returned for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(mLine::CoulombExcitation.MagneticLine, line::CoulombExcitation.Line,
                                      grid::Radial.Grid, settings::CoulombExcitation.Settings; printout::Bool=true)
    # Compute the amplitudes K^(Coulex) for the given magnetic line
    newChannels = CoulombExcitation.Channel[]
    for  channel in mLine.channels
        amplitude = CoulombExcitation.computeAmplitude(channel, mLine.Mi, mLine.Mf, line, grid, settings; printout=printout)
        push!(newChannels, CoulombExcitation.Channel(channel.q, channel.w, amplitude))
    end
    
    # Compute the partial cross section; collect parameters for line
    Ji2 = Basics.twice(line.initialLevel.J);    beta = CoulombExcitation.betaProjectile(line.ionEnergy);    partialCs = 0.;    
    for ch in newChannels
        wa = (ch.q^2 - (line.q0^2 * beta^2) )
        wa = ch.q * ch.w / (wa^2) * ch.amplitude * conj(ch.amplitude) 
        partialCs = partialCs + wa    
    end 
    partialCs = 2pi * (8pi * Defaults.getDefaults("alpha") / beta)^2 / (Ji2 + 1)  
    
    newmLine = CoulombExcitation.MagneticLine(mLine.Mi, mLine.Mf, partialCs, newChannels)
    return( newmLine )
end


"""
`CoulombExcitation.computeAmplitudesProperties(line::CoulombExcitation.Line, grid::Radial.Grid, 
                                               settings::CoulombExcitation.Settings; printout::Bool=true)`  
    ... to compute the amplitudes, cross sections and alignment parameters of the given line; a line::CoulombExcitation.Line is 
        returned for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(line::CoulombExcitation.Line, grid::Radial.Grid, 
                                      settings::CoulombExcitation.Settings; printout::Bool=true)
    newmLines = CoulombExcitation.MagneticLine[]
    for  mLine in line.mLines
        newmLine = CoulombExcitation.computeAmplitudesProperties(mLine, line, grid, settings; printout=printout)
        push!(newmLines, newmLine)
    end
    
    # Compute the total cross section and alignment parameters
    totalCs     = 0.;    for  mLine in newmLines    totalCs = totalCs + mLine.partialCs   end 
    alignmentA2 = 0.
    alignmentA4 = 0.
    
    
    newLine = CoulombExcitation.Line(line.initialLevel, line.finalLevel, line.ionEnergy, line.q0, totalCs, alignmentA2, alignmentA2, newmLines)
    return( newLine )
end


"""
`CoulombExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                settings::CoulombExcitation.Settings; output=true)`  
    ... to compute the Coulomb excitation amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{CoulombExcitation.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::CoulombExcitation.Settings; output=true)
    println("")
    printstyled("CoulombExcitation.computeLines(): The computation of Coulomb excitation cross sections starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
        
    sa =    "\n* Coulomb excitation cross sections for many-electron atoms and ions can be computed within different representations " *
            "and methods: " *
            "\n  they are all rather tricky and only approximate. The following assumptions and approximations are presently made: \n" *
            "\n    + All (projectile) ion energies are given in [MeV/u] which are converted into relative velocities beta = v/c. " *
            "\n    + All cross sections are computed for a (target) proton Z_t=1; multiply with Z_t^2 to find the correct cross sections. \n"
    println(sa)
    
    lines = CoulombExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    CoulombExcitation.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = CoulombExcitation.Line[]
    for  line in lines
        newLine = CoulombExcitation.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    CoulombExcitation.displayCrossSections(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    CoulombExcitation.displayCrossSections(iostream, newLines, settings)   end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`CoulombExcitation.determineChannels(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)`  
    ... to determine a list of CoulombExcitation.Channel for a magnetic line (mLine) from the initial to final level 
        and by taking into account the particular settings of for this computation; an Array{CoulombExcitation.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)
    channels = CoulombExcitation.Channel[];  
    # Compute the q's and associated weights in the interval [q0, 10*q0]
    gaussLegendre = Radial.GridGL("Finite", q0, 10*q0, settings.zerosGL);     qs = gaussLegendre.t;     ws = gaussLegendre.wt 
    for  (iq, q)  in  enumerate(qs)
        push!(channels, CoulombExcitation.Channel(q, ws[iq], ComplexF64(0.)) )
    end

    return( channels )  
end


"""
`CoulombExcitation.determineMagneticLines(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)`  
    ... to determine a list of CoulombExcitation.MagneticLine's for line from the initial to final level 
        and by taking into account the particular settings of for this computation; an Array{CoulombExcitation.Channel,1} 
        is returned.
"""
function determineMagneticLines(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)
    mLines = CoulombExcitation.MagneticLine[]
    
    for  Mi in AngularMomentum.m_values(initialLevel.J)
        for  Mf in AngularMomentum.m_values(finalLevel.J)
            channels = CoulombExcitation.determineChannels(finalLevel, initialLevel, q0, settings)
            if   length(channels) == 0   continue   end
            push!( mLines, CoulombExcitation.MagneticLine(Mi, Mf, 0., channels) )
        end
    end

    return( mLines )  
end


"""
`CoulombExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::CoulombExcitation.Settings)`  
    ... to determine a list of Coulomb-excitation Line's for transitions between the levels from the given initial- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{CoulombExcitation.Line,1} is returned. Apart from the level specification, all physical properties are set to 
        zero during the initialization process.  
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::CoulombExcitation.Settings)
    lines = CoulombExcitation.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                if   fLevel.energy - iLevel.energy == 0.   continue   end  
                for  ionEnergy  in  settings.ionEnergies
                    # Determine q0 for the given transition and ion energy     
                    beta = CoulombExcitation.betaProjectile(ionEnergy)
                    q0   = (fLevel.energy - iLevel.energy) / (beta * Defaults.getDefaults("speed of light: c"))
                    @show beta, q0
                    mLines = CoulombExcitation.determineMagneticLines(fLevel, iLevel, q0, settings)
                    if   length(mLines) == 0   continue   end
                    push!( lines, CoulombExcitation.Line(iLevel, fLevel, ionEnergy, q0, 0., 0., 0., mLines) )
                end 
            end
        end
    end
    
    return( lines )
end


"""
`CoulombExcitation.displayCrossSections(stream::IO, lines::Array{CoulombExcitation.Line,1}, settings::CoulombExcitation.Settings)`  
    ... to display a list of lines, magnetic lines and channels that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayCrossSections(stream::IO, lines::Array{CoulombExcitation.Line,1}, settings::CoulombExcitation.Settings)
    nx = 120
    println(stream, " ")
    println(stream, "  Partial and total Coulomb-excitation lines, magnetic lines and channel parameters:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Ion energy"; na=4);              
    sb = sb * TableStrings.center(14, " [MeV/u]"; na=4)
    sa = sa * TableStrings.center(16, "Total CS"; na=2);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("cross section"); na=4)
    sa = sa * TableStrings.center(16, "Mi    Mf"; na=2);                          sb = sb * TableStrings.hBlank(16)   
    sa = sa * TableStrings.center(16, "Partial CS"; na=2);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("cross section"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                     fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", line.ionEnergy)  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.totalCs)) * "    "
        println(stream, sa )        
        for  mL in line.mLines 
            sc = "      " * string(mL.Mi) * "    " * string(mL.Mf) * "           "
            sc = sc[1:22]
            sc = sc * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", mL.partialCs))
            println(stream, TableStrings.hBlank(length(sa)) * sc )   
        end
     end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`CoulombExcitation.displayLines(stream::IO, lines::Array{CoulombExcitation.Line,1})`  
    ... to display a list of lines, magnetic lines and channels that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{CoulombExcitation.Line,1})
    nx = 160
    println(stream, " ")
    println(stream, "  Selected Coulomb-excitation lines, magnetic lines and channel parameters:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Ion energy"; na=4);              
    sb = sb * TableStrings.center(14, " [MeV/u]"; na=4)
    sa = sa * TableStrings.flushleft(50, "List of (Mi, Mf, q1, w1) for mLines + channel_1"; na=4);       
    sb = sb * TableStrings.hBlank(34)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                     fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", line.ionEnergy)  * "    "
        MMtuples = Tuple{AngularM64, AngularM64, Float64, Float64}[]
        for  mLine  in  line.mLines
            push!( MMtuples, (mLine.Mi, mLine.Mf, mLine.channels[1].q, mLine.channels[1].w) )
        end
        wa = TableStrings.MMffTupels(100, MMtuples, "mL")
        println(stream, sa * wa[1] )        
        for  ia = 2:length(wa)   println(stream, TableStrings.hBlank(length(sa)) * wa[ia] )   end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
