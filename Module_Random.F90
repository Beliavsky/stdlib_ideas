module Random
! ########################################################################
!
! Copyright 2020 IRD
!
! This file is part of statpack.
!
! statpack is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of 
! the License, or (at your option) any later version.
!
! statpack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You can find a copy of the GNU Lesser General Public License
! in the statpack/doc directory.
!
! ########################################################################
!                                                                        *
! ************************************************************************
! THIS MODULE REPLACES THE FORTRAN 90 INTRINSICS                         *
! random_number AND random_seed BY SEVERAL IMPLEMENTATIONS               *
! OF THE KISS (Keep It Simple Stupid), L'Ecuyer's LFSR113,               *
! MERSENNE TWISTER MT19937 AND MEMT19937-II RANDOM                       *
! NUMBER GENERATORS.                                                     *
!                                                                        *
! IN ADDITION TO 10 DIFFERENT UNIFORM RANDOM GENERATORS, GAUSSIAN        *
! RANDOM GENERATORS, SHUFFLING AND SAMPLING SUBROUTINES ARE ALSO         *
! PROVIDED, AS WELL AS SUBROUTINES FOR GENERATING PSEUDO-RANDOM          *
! ORTHOGONAL MATRICES FOLLOWING THE HAAR DISTRIBUTION OVER THE GROUP OF  *
! ORTHOGONAL MATRICES, PSEUDO_RANDOM SYMMETRIC MATRICES WITH A           *
! PRESCRIBED SPECTRUM OR PSEUDO_RANDOM MATRICES WITH A PRESCRIBED        *
! SINGULAR VALUE DISTRIBUTION.                                           *
!                                                                        *
! MANY PARTS OF THIS MODULE ARE ADAPTED FROM :                           *
!                                                                        *
! Hennecke, M., 1995: A Fortran90 interface to random                    *
!                    number generation.                                  *
!                    Computer Physics Communications.                    *
!                    Volume 90, Number 1, 117-120                        *
!                                                                        *
! LATEST REVISION : 20/10/2020                                           *
!                                                                        *
! ########################################################################
!                                                                        *
! ************************************************************************
!                                                                        *
! THE C PROCESSOR MACROS USED IN THIS MODULE ARE:                        *
!                                                                        *
!  _OPENMP           FOR ACTIVATING OPENMP PARALLELIZATION               *
!  _ALLOC            FOR ALLOCATING LOCAL VARIABLES INSTEAD OF PLACING   *
!                    THEM ON THE STACK IN SOME SUBROUTINES AND FUNCTIONS *
!  _BLAS             FOR USING BLAS WHEN POSSIBLE                        *
!  _MATMUL           FOR REPLACING THE matmul INTRINSIC FUNCTION WITH    *
!                    STATPACK matmul2 FUNCTION                           *
!  _TRANSPOSE        FOR REPLACING THE transpose INTRINSIC FUNCTION WITH *
!                    STATPACK transpose2 FUNCTION                        *
!  _RANDOM_WITH0     FOR GENERATING REAL FLOATING POINT NUMBERS          *
!                    IN THE [0,1[ INTERVAL INSTEAD OF ]0,1[ INTERVAL     *
!  _RANDOM_NOINT32   FOR SIGNALING THAT 32 BIT INTEGERS ARE NOT          *
!                    AVAILABLE                                           *
!  _RANDOM_NOUNIX    FOR SIGNALING THAT THE OPERATING SYSTEM IS NOT UNIX *
!  _RANDOM_GFORTRAN  FOR SIGNALING THAT THE INTEGER Unix FUNCTION getpid *
!                    IS CONSIDERED AS AN INTRINSIC RATHER AS AN EXTERNAL *
!                    PROCEDURE AS FOR THE gfortran COMPILER              *
!  _RANDOM_NAGWARE   FOR SIGNALING THAT THE INTEGER Unix FUNCTION getpid *
!                    IS LOCATED IN THE f90_unix_env MODULE FOR THE NAG   *
!                    FORTRAN COMPILER                                    *
!  _USE_GNU          FOR SIGNALING THAT THE GNU GFORTRAN COMPILER IS     *
!                    USED. HERE, THIS INCLUDES ONLY THE ACTIVATION OF    *
!                    THE _RANDOM_GFORTRAN MACRO                          *
!  _USE_NAGWARE      FOR SIGNALING THAT THE NAG FORTRAN COMPILER IS      *
!                    USED. HERE, THIS INCLUDES ONLY THE ACTIVATION OF    *
!                    THE _RANDOM_NAGWARE MACRO                           *
!  _USE_PGI          FOR DESACTIVATING SOME OPENMP CONSTRUCTS FOR THE    *
!                    PORTLAND PGFORTRAN COMPILER                         *
!                                                                        *
! ************************************************************************
!
#ifdef _USE_GNU
#define _RANDOM_GFORTRAN
#endif
!
#ifdef _USE_NAGWARE
#define _RANDOM_NAGWARE
#endif
!
!
! USED MODULES
! ============
!
    use Select_Parameters,  only : lgl, stnd, extd, i4b, urandom_file
    use Utilities,          only : merror, arth
    use Reals_Constants,    only : zero, half, one, two, twopi
    use Logical_Constants
#ifdef _MATMUL
    use Utilities,          only : matmul=>matmul2
#endif
#ifdef _TRANSPOSE
    use Utilities,          only : transpose=>transpose2
#endif
#ifdef _OPENMP
    use Select_Parameters,  only : omp_limit
    use omp_lib,            only : omp_get_max_threads, omp_in_parallel, omp_get_nested
#endif
#ifdef _RANDOM_NAGWARE
    use f90_unix_env,       only : getpid
#endif
#ifdef _BLAS
    use BLAS_interfaces,    only : gemm
#endif
!
! OTHER MODULES USED IN SPECIFIC SUBROUTINES
! ==========================================
!
!   use Char_Constants
!   use Prob_Procedures
!   use Hous_Procedures
!   use QR_Procedures
!   use BLAS_interfaces
!
! STRONG TYPING IMPOSED
! =====================
!    
    implicit none
!
! PUBLIC ENTITIES 
! ===============
!    
! ALL SUBROUTINES, FUNCTIONS, VARIABLES AND PARAMETERS ARE PRIVATE BY DEFAULT.
!
    private
    public ::                                                        &
              rand_number,    rand_integer32,    rand_integer31,     &
              random_number_, random_integer32_, random_integer31_,  &
              random_seed_,   init_mt19937,      init_memt19937,     &
              normal_rand_number,  normal_random_number_,            &
              normal_rand_number2, normal_random_number2_,           &
              normal_rand_number3, normal_random_number3_,           &
              random_qr_cmp, ortho_gen_random_qr,                    &
              gen_random_sym_mat, gen_random_mat,                    &
              simple_shuffle, drawsample, drawbootsample
!
! GENERIC INTERFACES FOR ROUTINES WITH OVERLOADED VERSIONS
! ========================================================
!
    interface random_number_
        module procedure     random_r0,                        &
                             random_r1,                        &
                             random_r2,                        &
                             random_r3,                        &
                             random_r4,                        &
                             random_r5,                        &
                             random_r6,                        &
                             random_r7
    end interface
!
    interface random_integer32_
        module procedure     random_i32_0,                        &
                             random_i32_1,                        &
                             random_i32_2,                        &
                             random_i32_3,                        &
                             random_i32_4,                        &
                             random_i32_5,                        &
                             random_i32_6,                        &
                             random_i32_7
    end interface
!
    interface random_integer31_
        module procedure     random_i31_0,                        &
                             random_i31_1,                        &
                             random_i31_2,                        &
                             random_i31_3,                        &
                             random_i31_4,                        &
                             random_i31_5,                        &
                             random_i31_6,                        &
                             random_i31_7
    end interface
!
    interface init_mt19937
        module procedure     init_mt19937_r0,                     &
                             init_mt19937_r1
    end interface
!
    interface init_memt19937
        module procedure     init_memt19937_r0,                   &
                             init_memt19937_r1
    end interface
!
    interface normal_random_number_
        module procedure     normal_random_r0,                    &
                             normal_random_r1,                    &
                             normal_random_r2
    end interface
!
    interface normal_random_number2_
        module procedure     normal_random2_r0,                   &
                             normal_random2_r1,                   &
                             normal_random2_r2
    end interface
!
    interface normal_random_number3_
        module procedure     normal_random3_r0,                   &
                             normal_random3_r1,                   &
                             normal_random3_r2
    end interface
!
    interface simple_shuffle
        module procedure     simple_shuffle_rv, simple_shuffle_cv, simple_shuffle_iv
    end interface
!
!
! MODULE PARAMETERS
! =================
!
!   ki_sel IS THE KIND PARAMETER FOR 32-BIT INTEGER.
!
    integer, parameter :: ki0        = kind(0)
    integer, parameter :: ki9        = selected_int_kind(9)
    integer, parameter :: ki_sel     = max(ki9,sign(ki0,-ki9))
!
    integer, parameter :: ratio_i    = (int(bit_size(0_ki_sel)) + bit_size(0) - 1)/bit_size(0)
!
!   bit_size(0_ki_sel) MUST BE EQUAL TO fullbitsize IF 32_BIT INTEGERS ARE AVAILABLE.
!
    integer(ki_sel), parameter :: fullbitsize = 32_ki_sel
!
!   THE CONSTANTS fbs, hbs, qbs AND tbs ARE USED FOR EMULATING UNSIGNED 32-BIT INTEGER ARITHMETIC
!   WITHOUT OVERFLOWS IN INTERNAL ROUTINES. 
!
    integer, parameter :: fbs = int(fullbitsize), hbs = fbs/2, qbs = hbs/2, tbs = 3*qbs, topbit = fbs - 1
!
!   SPECIFY HERE THE LOCATION OF THE URANDOM DEVICE ON YOUR SYSTEM IF IT EXISTS.
!   FOR UNIX SYSTEM, NORMALLY IT IS /dev/urandom .
!
!    character(len=*),  parameter :: urandom_file='/dev/urandom'
!
#ifdef _RANDOM_WITH0
!  THE CONSTANT m_ran_32 IS 2^32.
!  THE CONSTANT m_ran_invm32 IS THE RECIPROCAL OF 2^32.
!
    real(stnd), parameter :: m_ran_32     = 4294967296.0_stnd,        &
                             m_ran_invm32 = one/m_ran_32
#else
!  THE CONSTANT m_ran_32 IS THE REAL NUMBER FOLLOWING 2^32.
!  THE CONSTANT m_ran_invm32 IS THE RECIPROCAL OF THIS REAL NUMBER.
!
    real(stnd), parameter :: m_ran_32     = 4294967296.0_stnd*(one + epsilon(half) ),        &
                             m_ran_invm32 = one/m_ran_32
#endif
!
!   num_bits_stnd IS THE NUMBER OF BITS IN THE MANTISSA OF REAL NUMBERS OF PRECISION stnd. WE
!   ASSUME THAT THE MANTISSA IS NORMALIZED AS IN THE IEEE STANDARD.
!
    integer,  parameter :: num_bits_stnd = min(digits(half),64)-1,       &
                           num_shift=max(min(num_bits_stnd-64,-1),-topbit)
!
!   THE FOLLOWING INTEGER PARAMETERS ARE USED BY THE 32-BIT XORSHIFT RNG WHICH IS ONE OF THE
!   COMPONENT RNGS IN THE KISS RANDOM NUMBER GENERATORS.
!
    integer,  parameter ::     kiss_shift1   = 13,  &
                               kiss_shift2   = -17, &
                               kiss_shift3   = 5
!
!   THE FOLLOWING INTEGER PARAMETERS ARE USED BY THE MERSENNE TWISTER RANDOM NUMBER GENERATORS.
!
!   PERIOD PARAMETERS FOR THE MERSENNE TWISTER MT19937 AND MEMT19937-II RNGS.
!
    integer(ki_sel),  parameter ::  mt_n    = 624_ki_sel,   &
                                    mt_m    = 397_ki_sel
!
#ifdef _RANDOM_NOINT32
    integer(ki_sel),  parameter ::  mt_upper_mask  = ishft( 1_ki_sel, topbit),                &
                                    mt_lower_mask  = ibits( not( mt_upper_mask ), 0, fbs ),   &
                                    mt_matrix_a    = ior( 419999967_ki_sel, mt_upper_mask ),  &
                                    mt_mag01(0:1)  = (/ 0_ki_sel, mt_matrix_a /)
#else
    integer(ki_sel),  parameter ::  mt_upper_mask  = ishft( 1_ki_sel, topbit),                &
                                    mt_lower_mask  = not( mt_upper_mask ),                    &
                                    mt_matrix_a    = ior( 419999967_ki_sel, mt_upper_mask ),  &
                                    mt_mag01(0:1)  = (/ 0_ki_sel, mt_matrix_a /)
#endif
!
!   THESE PARAMETERS MAY BE ALSO SPECIFIED AS FOLLOW ASSUMING 32-BIT INTEGER AND INTEGER REPRESENTATION
!   IS TWO'S COMPLEMENT. HOWEVER, THE INITIALIZATION OF mt_upper_mask TO THE NEGATIVE 32-BIT INTEGER
!   OF THE LARGEST AMPLITUDE AT COMPILATION TIME MAY CAUSE PROBLEMS ON SOME SYSTEMS (E.G. SOME VERSIONS
!   OF GFORTRAN FOR EXAMPLE). THIS IS DUE TO THE FACT THAT THE COMPILER WILL TRY TO FORM
!   THE POSITIVE INTEGER 2147483648 (WHICH IS NOT REPRESENTABLE IF TWO'S COMPLEMENT NOTATION
!   IS USED), BEFORE NEGATING IT. THIS WILL RESULT IN INTEGER OVERFLOW, WHICH IS USUALLY
!   CHECKED AT COMPILATION TIME FOR CONSTANT PARAMETERS.
!
!    integer(ki_sel),  parameter ::     mt_upper_mask  = -2147483648_ki_sel,        &
!                                       mt_lower_mask  =  2147483647_ki_sel,        &
!                                       mt_matrix_a    = -1727483681_ki_sel,        &
!                                       mt_mag01(0:1)  = (/ 0_ki_sel, mt_matrix_a /)
!
!   THESE PARAMETERS MAY BE ALSO SPECIFIED AS FOLLOW, WITH HEXADECIMAL CONSTANTS. HOWEVER, THESE HEXADECIMAL
!   CONSTANTS MAY CAUSE PROBLEMS ON SOME COMPILERS (E.G. GFORTRAN) BECAUSE THE "NEGATIVE" 32-BIT INTEGERS ARE
!   CONSIDERED AS UNSIGNED 32-BIT INTEGERS INSTEAD OF SIGNED 32-BIT INTEGERS DURING COMPILATION. THIS WILL
!   RESULT IN INTEGER OVERFLOW, WHICH IS USUALLY CHECKED AT COMPILATION TIME FOR CONSTANT PARAMETERS.
!
!    integer(ki_sel),  parameter ::     mt_upper_mask  = int( z'80000000', ki_sel),      &
!                                       mt_lower_mask  = int( z'7fffffff', ki_sel),      &
!                                       mt_matrix_a    = int( z'9908b0df', ki_sel),      &
!                                       mt_mag01(0:1)  = (/ 0_ki_sel, mt_matrix_a /)
!
!   DEFAULT SEEDS FOR THE MERSENNE TWISTER MT19937 AND MEMT19937-II RNGS.
!
    integer(ki_sel),  parameter :: mt_default_seed = 5489_ki_sel
!    integer(ki_sel),  parameter :: mt_default_seed = 21641_ki_sel
!
    integer(ki_sel),  dimension(4), parameter ::       &
           mt_default_seed_array = (/ 291_ki_sel, 564_ki_sel, 837_ki_sel, 1110_ki_sel /)
!
!   TEMPERING PARAMETERS FOR THE MERSENNE TWISTER MT19937 RNG.
!
    integer,  parameter ::     mt_t1_shift = 7,   &
                               mt_t2_shift = 15,  &
                               mt_shift0   = -11, &
                               mt_shift1   = -18
!
    integer(ki_sel),  parameter ::  mt_t1_mask     = ior( 489444992_ki_sel , mt_upper_mask),  &
                                    mt_t2_mask     = ior( 1875247104_ki_sel, mt_upper_mask)
!
!   THESE PARAMETERS MAY BE ALSO SPECIFIED AS FOLLOW ASSUMING 32-BIT INTEGER AND INTEGER REPRESENTATION
!   IS TWO'S COMPLEMENT.
!
!    integer(ki_sel),  parameter ::     mt_t1_mask     = -1658038656_ki_sel,      &
!                                       mt_t2_mask     = -272236544_ki_sel
!
!   THESE PARAMETERS MAY BE ALSO SPECIFIED AS FOLLOW, WITH HEXADECIMAL CONSTANTS. HOWEVER, THESE HEXADECIMAL
!   CONSTANTS MAY CAUSE PROBLEMS ON SOME COMPILERS (E.G. GFORTRAN) BECAUSE THE "NEGATIVE" 32-BIT INTEGERS ARE
!   CONSIDERED AS UNSIGNED 32-BIT INTEGERS INSTEAD OF SIGNED 32-BIT INTEGERS DURING COMPILATION.
!
!    integer(ki_sel),  parameter ::     mt_t1_mask     = int( z'9d2c5680', ki_sel),      &
!                                       mt_t2_mask     = int( z'efc60000', ki_sel)
!
!   TEMPERING PARAMETERS FOR THE MERSENNE TWISTER MEMT19937-II RNG.
!
    integer,  parameter ::     memt_lag1     = 151, &
                               memt_lag2     = 36,  &
                               memt_lag1over = 473, &
                               memt_lag2over = 588, &
                               memt_shift1   = 8,   &
                               memt_shift2   = 14
!
    integer(ki_sel), parameter ::  memt_mask1 = ior( 840548011_ki_sel, mt_upper_mask ),  &
                                   memt_mask2 = 1455285546_ki_sel
!
!   THESE PARAMETERS MAY BE ALSO SPECIFIED AS FOLLOW ASSUMING 32-BIT INTEGER AND INTEGER REPRESENTATION
!   IS TWO'S COMPLEMENT.
!
!    integer(ki_sel), parameter ::  memt_mask1 = -1306935637_ki_sel,  &
!                                   memt_mask2 = 1455285546_ki_sel
!
!   THESE PARAMETERS MAY BE ALSO SPECIFIED AS FOLLOW, WITH HEXADECIMAL CONSTANTS. HOWEVER, THESE HEXADECIMAL
!   CONSTANTS MAY CAUSE PROBLEMS ON SOME COMPILERS (E.G. GFORTRAN) BECAUSE THE "NEGATIVE" 32-BIT INTEGERS ARE
!   CONSIDERED AS UNSIGNED 32-BIT INTEGERS INSTEAD OF SIGNED 32-BIT INTEGERS DURING COMPILATION.
!
!    integer(ki_sel),  parameter ::     memt_mask1     = int( z'b219beab', ki_sel),      &
!                                       memt_mask2     = int( z'56bde52a', ki_sel)
!
#ifdef _RANDOM_NOINT32
!   THE FOLLOWING PARAMETER IS USED TO SET THE LAST 32 BITS OF INTEGER OF KIND ki_sel TO 1,
!   IF 32 BIT-INTEGERS ARE NOT AVAILABLE.
!
    integer(ki_sel),  parameter ::  noint32_mask  = not( ior( mt_lower_mask, mt_upper_mask ) )
!
#endif
!
!
! MODULE VARIABLES
! ================
!
!   first     = true MEANS THAT THE random_seed_ SUBROUTINE HAS NEVER BEEN CALLED BEFORE.
!   test_kiss = true MEANS THAT THE KISS RNGS HAVE NOT BEEN YET TESTED FOR OVERFLOWS.
!
    logical(lgl), save :: first = true, test_kiss = true
!
!   THE PRIVATE INTEGER VARIABLE method SETS THE RANDOM NUMBER GENERATOR USED FOR SUBSEQUENT
!   CALLS OF THE rand_number, rand_integer32, rand_integer31, random_number_, random_integer32_,
!   random_integer31_ ROUTINES:
!
!     method = 1  : THE Marsaglia's KISS RNG IS USED.
!     method = 2  : THE FAST Marsaglia's KISS RNG, WHICH USES ONLY
!                   ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS,
!                   IS USED (SAME RESULTS AS THE C VERSION).
!     method = 3  : THE L'Ecuyer's LFSR113 RNG IS USED.
!     method = 4  : THE MERSENNE TWISTER MT19937 RNG IS USED.
!     method = 5  : THE MERSENNE TWISTER MEMT19937-II RNG IS USED.
!     method = 6  : THE EXTENDED VERSION OF THE Marsaglia's KISS RNG IS USED.
!     method = 7  : THE EXTENDED VERSION OF THE FAST Marsaglia's KISS RNG IS USED.
!     method = 8  : THE EXTENDED VERSION OF THE L'Ecuyer's LFSR113 RNG IS USED.
!     method = 9  : THE EXTENDED VERSION OF THE MERSENNE TWISTER MT19937 RNG IS USED.
!     method = 10 : THE EXTENDED VERSION OF THE MERSENNE TWISTER MEMT19937-II RNG IS USED.
!
!   IN ORDER TO CHANGE THE RNG, YOU MUST USE A CALL TO THE random_seed_ SUBROUTINE WITH
!   THE OPTIONAL alg ARGUMENT SET TO 1, 2, 3, 4, 5, 6, 7, 8, 9 OR 10. NOTE THAT FOR RANDOM INTEGER
!   GENERATION, METHODS 6, 7, 8, 9 AND 10 ARE EXACTLY IDENTICAL TO METHODS 1, 2, 3, 4 AND 5,
!   RESPECTIVELY. THEY DIFFER ONLY FOR RANDOM REAL GENERATION.
!
!   THE PRIVATE INTEGER VARIABLE method STORED THE METHOD CURRENTLY IN USE, num_seeds
!   STORED THE NUMBER OF SEEDS USED BY THIS METHOD AND seed_size IS THE CORRESPONDING SIZE
!   OF THE SEED VECTOR AS RETURNED BY A CALL TO random_seed_ SUBROUTINE WITH ARGUMENT size.
!
!   THE DEFAULT METHOD CAN BE SAFELY CHANGED HERE, BETWEEN 1, 2 AND 3.
!   NOTE THAT THE OTHER METHODS CANNOT BE USED SAFELY AS THE DEFAULT METHOD.
!   FURTHERMORE IF 32-BIT INTEGERS ARE NOT AVAILABLE THE DEFAULT METHOD MUST
!   BE 3, BECAUSE THE KISS RNGS WORK PROPERLY ONLY WITH 32-BIT INTEGERS.
!
!    integer,  save :: method  =  1
!    integer,  save :: method  =  2
    integer,  save :: method  =  3
!
!   num_seeds IS THE SIZE OF THE STATE VECTOR OF THE CURRENT RNG.
!   seed_size IS THE SIZE OF THE SEED ARRAY OF THE RANDOM GENERATOR AS
!   RETURNED BY SUBROUTINE random_seed_. NOTE THAT THE seed_size IS NOT
!   NECESSARILY EQUAL TO num_seeds DEPENDING ON THE RNG AND IF 32-BIT INTEGERS
!   ARE OR NOT THE DEFAULT INTEGERS ON YOUR SYSTEM.
!
    integer,  save :: num_seeds, seed_size
!
!   THE NEXT VARIABLES STORE THE STATE VARIABLES FOR EACH RANDOM NUMBER GENERATORS.
!
!   PRIVATE INTEGER STATE VARIABLES USED BY THE Marsaglia's KISS RNG (method = 1 or 6).
!    
    integer(ki_sel), save  :: x = 123456789_ki_sel, y = 362436069_ki_sel,             &
                              z = 521288629_ki_sel, w = 916191069_ki_sel
!
!   PRIVATE INTEGER STATE VARIABLES USED BY THE FAST Marsaglia's KISS RNG, WHICH USES ONLY
!   ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS (method = 2 or 7).
!    
    integer(ki_sel), save  :: x2 = 123456789_ki_sel, y2 = 362436069_ki_sel,            &
                              z2 = 21288629_ki_sel,  w2 = 14921776_ki_sel, c=0_ki_sel
!
!   PRIVATE INTEGER STATE VARIABLES USED BY THE L'Ecuyer's LFSR113 RNG (method = 3 or 8).
!    
    integer(ki_sel), save  :: s1 = 153587801_ki_sel,  s2 = -759022222_ki_sel,        &
                              s3 = 1288503317_ki_sel, s4 = -1718083407_ki_sel
!
!   PRIVATE INTEGER STATE VARIABLES USED BY THE THE MERSENNE TWISTER MT19937 RNG (method = 4 or 9).
!   mt_initialized = false MEANS THAT THE STATE OF MT19937 IS NOT INITIALIZED.
!   mti = mt_n             MEANS THAT THE ARRAY mt(:) IS NOT INITIALIZED.
!   THE ARRAY mt(:) IS THE STATE VECTOR OF THE MERSENNE TWISTER MT19937 RNG.
!    
    logical(lgl),    save  :: mt_initialized = false
    integer(ki_sel), save  :: mti = mt_n, mt(0_ki_sel:mt_n-1_ki_sel) = 0_ki_sel
!
!   PRIVATE INTEGER STATE VARIABLES USED BY THE THE MERSENNE TWISTER MEMT19937-II RNG (method = 5 or 10).
!   memt_initialized = false MEANS THAT THE STATE OF MEMT19937-II IS NOT INITIALIZED.
!   memti = -1_ki_sel        MEANS THAT THE ARRAY memt(:) IS NOT INITIALIZED.
!   THE ARRAY memt(:) IS THE STATE VECTOR OF THE MERSENNE TWISTER MEMT19937-II RNG.
!    
    logical(lgl),    save  :: memt_initialized = false
    integer(ki_sel), save  :: memti = -1_ki_sel, memt(0_ki_sel:mt_n-1_ki_sel) = 0_ki_sel
!
!   THE FOLLOWING PRIVATE VARIABLE m_ran_invmbits IS USED TO TRANSFORM SIGNED 32-BIT INTEGERS
!   INTO REAL FLOATING POINT NUMBERS FOR THE EXTENDED VERSIONS OF THE RNGS.
!
    real(stnd), save :: m_ran_invmbits = -one
!
#ifndef _RANDOM_NOUNIX
!    
! INTEGER Unix FUNCTION USED FOR INITIALIZING THE SEEDS OF THE RNGS ON UNIX SYSTEMS.
!
#ifdef _RANDOM_GFORTRAN
    intrinsic getpid
#else
!
#ifndef _RANDOM_NAGWARE
!
    integer :: getpid
!
    external getpid
#endif
!
#endif
!
#endif
!
!
! =========================================================================================
!
                                  contains
!                                 ========
!
! =========================================================================================
!                            RANDOM FLOATING POINT NUMBER SUBROUTINES
! =========================================================================================
!
    function rand_number() result( harvest )
!
! purpose
! _______
!
!   This function returns a uniformly distributed random number between 0 and 1,
!   exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   None
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation,
!   this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _RANDOM_NOINT32
    use Char_Constants,   only : random_error7, random_error9
#else
    use Char_Constants,   only : random_error1, random_error2,    &
                                 random_error7, random_error8
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel)               :: k, j, tmp
    integer(ki_sel), dimension(2) :: tmpvec
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='rand_number'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( first ) then
#ifdef _RANDOM_NOINT32
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method/=3 ) then
            call merror( name_proc//random_error9 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ibclr(-1_ki_sel, bit_size(tmp)-1_ki_sel)
!
        if ( radix(tmp)/=2 .or. tmp/=huge(tmp) ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
#else
!
!       TEST IF 32-BIT INTEGERS ARE AVAILABLE.
!
        if ( bit_size( tmp )/=fullbitsize ) then
            call merror( name_proc//random_error1 )
        endif
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method<1 .or. method>3 ) then
            call merror( name_proc//random_error2 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ishftc( -8_ki_sel, -3)
!
        if ( radix(tmp)/=2 .or. tmp/=536870911_ki_sel ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
        if ( method/=3 ) then
!
!           TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS.
!
            tmp   = huge( tmp )
!
            if ( integer_not_safe( tmp ) ) then
!
                call merror( name_proc//random_error8 )
!
            endif
!
        endif
#endif
!
    endif
!
!$OMP CRITICAL (ran_num_)
    select case (method)
!
        case (1)
!
!           Marsaglia's KISS RNG.
!
            x = 69069_ki_sel * x + 1327217885_ki_sel
            y = ieor( y, ishft(y, kiss_shift1) )
            y = ieor( y, ishft(y, kiss_shift2) )
            y = ieor( y, ishft(y, kiss_shift3) )
!            y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
            z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
            w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
            tmp = x + y + ishft( z, 16) + w
!
        case (2)
!
!           FAST Marsaglia's KISS RNG WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            x2 = x2 + 545925293_ki_sel
!            x2 = uiadd( x2, 545925293_ki_sel )
            y2 = ieor( y2, ishft(y2, kiss_shift1) )
            y2 = ieor( y2, ishft(y2, kiss_shift2) )
            y2 = ieor( y2, ishft(y2, kiss_shift3) )
!            y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
            tmp = z2 + w2 + c
            z2 = w2
            c  = ishft( tmp, -topbit)
            w2 = iand( tmp, 2147483647_ki_sel)
!
            tmp = x2 + y2 + w2
!
        case (3)
!
!           L'Ecuyer's LFSR113 RNG.
!
#ifdef _RANDOM_NOINT32
            tmp = ieor( ishft(s1,6), s1)
            tmp = ishft( ibits( tmp, 0, fbs ), -13)
            tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
            s1  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s2,2), s2)
            tmp = ishft( ibits( tmp, 0, fbs ), -27)
            tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
            s2  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s3,13), s3)
            tmp = ishft( ibits( tmp, 0, fbs ), -21)
            tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
            s3  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s4,3), s4)
            tmp = ishft( ibits( tmp, 0, fbs ), -12)
            tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
            s4  = ibits( tmp, 0, fbs )
#else
            tmp  = ishft( ieor( ishft(s1,6), s1), -13)
            s1 = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
            tmp  = ishft( ieor( ishft(s2,2), s2), -27)
            s2 = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
            tmp  = ishft( ieor( ishft(s3,13), s3), -21)
            s3 = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
            tmp  = ishft( ieor( ishft(s4,3), s4), -12)
            s4 = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
            tmp = ieor( ieor( ieor(s1,s2), s3), s4)
!
        case (4)
!
!           MERSENNE TWISTER MT19937 RNG.
!
            if ( mti>=mt_n ) then
!
!               GENERATE mt_n WORDS AT ONE TIME.
!
                do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),            &
                                 iand(mt(j+1_ki_sel),mt_lower_mask)    )
!
                    mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)),      &
                                  mt_mag01(iand(tmp,1_ki_sel))         )
!
                end do
!
                do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),            &
                                 iand(mt(j+1_ki_sel),mt_lower_mask)    )
!
                    mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)), ishft(tmp,-1)), &
                                  mt_mag01(iand(tmp,1_ki_sel))            )
!
                end do
!
                tmp   = ior( iand(mt(j),mt_upper_mask),         &
                             iand(mt(0_ki_sel),mt_lower_mask)   )
!
                mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                              mt_mag01(iand(tmp,1_ki_sel))             )
!
                mti = 0_ki_sel
!
            endif
!
            tmp = mt(mti)
            mti = mti + 1_ki_sel
!
!           TEMPERING PHASE.
!
            tmp = ieor( tmp, ishft(tmp,mt_shift0) )
            tmp = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
            tmp = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
            tmp = ieor( tmp, ishft(tmp,mt_shift1) )
!
        case (5)
!
!           MERSENNE TWISTER MT19937-II RNG.
!
            select case (memti)
!
                case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))          )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag2over:mt_n-2_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-1_ki_sel)
!
                    tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                               iand(memt(0_ki_sel),mt_lower_mask)        )
!
                    memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                mt_mag01(iand(tmp,1_ki_sel))               )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = 0_ki_sel
!
            end select
!
#ifdef _RANDOM_NOINT32
            tmp = ibits( tmp, 0, fbs )
#endif
!
        case (6)
!
!           EXTENDED PRECISION VERSION OF Marsaglia's KISS RNG.
!
            do k = 1_ki_sel, 2_ki_sel
!
                x = 69069_ki_sel * x + 1327217885_ki_sel
                y = ieor( y, ishft(y, kiss_shift1) )
                y = ieor( y, ishft(y, kiss_shift2) )
                y = ieor( y, ishft(y, kiss_shift3) )
!                y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
                z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
                w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
                tmpvec(k) = x + y + ishft( z, 16) + w
!
            end do
!
        case (7)
!
!           EXTENDED PRECISION VERSION OF THE FAST Marsaglia's KISS RNG,
!           WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            do k = 1_ki_sel, 2_ki_sel
!
                x2 = x2 + 545925293_ki_sel
!                x2 = uiadd( x2, 545925293_ki_sel )
                y2 = ieor( y2, ishft(y2, kiss_shift1) )
                y2 = ieor( y2, ishft(y2, kiss_shift2) )
                y2 = ieor( y2, ishft(y2, kiss_shift3) )
!                y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
                tmp = z2 + w2 + c
                z2 = w2
                c  = ishft( tmp, -topbit)
                w2 = iand( tmp, 2147483647_ki_sel)
!
                tmpvec(k) = x2 + y2 + w2
!
            end do
!
        case (8)
!
!           EXTENDED PRECISION VERSION OF L'Ecuyer's LFSR113 RNG.
!
            do k = 1_ki_sel, 2_ki_sel
!
#ifdef _RANDOM_NOINT32
                tmp = ieor( ishft(s1,6), s1)
                tmp = ishft( ibits( tmp, 0, fbs ), -13)
                tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
                s1  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s2,2), s2)
                tmp = ishft( ibits( tmp, 0, fbs ), -27)
                tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
                s2  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s3,13), s3)
                tmp = ishft( ibits( tmp, 0, fbs ), -21)
                tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
                s3  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s4,3), s4)
                tmp = ishft( ibits( tmp, 0, fbs ), -12)
                tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
                s4  = ibits( tmp, 0, fbs )
#else
                tmp = ishft( ieor( ishft(s1,6), s1), -13)
                s1  = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
                tmp = ishft( ieor( ishft(s2,2), s2), -27)
                s2  = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
                tmp = ishft( ieor( ishft(s3,13), s3), -21)
                s3  = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
                tmp = ishft( ieor( ishft(s4,3), s4), -12)
                s4  = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
                tmpvec(k) = ieor( ieor( ieor(s1,s2), s3), s4)
!
            end do
!
        case (9)
!
!           EXTENDED PRECISION VERSION OF THE MERSENNE TWISTER MT19937 RNG.
!
            do k = 1_ki_sel, 2_ki_sel
!
                if ( mti>=mt_n ) then
!
!                   GENERATE mt_n WORDS AT ONE TIME.
!
                    do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)),   &
                                      mt_mag01(iand(tmp,1_ki_sel))      )
!
                    end do
!
                    do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)), ishft(tmp,-1)),  &
                                      mt_mag01(iand(tmp,1_ki_sel))             )
!
                    end do
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(0_ki_sel),mt_lower_mask)   )
!
                    mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)),  &
                                  mt_mag01(iand(tmp,1_ki_sel))            )
!
                    mti = 0_ki_sel
!
                endif
!
                tmp = mt(mti)
                mti = mti + 1_ki_sel
!
!               TEMPERING PHASE.
!
                tmp = ieor( tmp, ishft(tmp,mt_shift0) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
!
                tmpvec(k) = ieor( tmp, ishft(tmp,mt_shift1) )
!
            end do
!
        case (10)
!
!           EXTENDED PRECISION VERSION OF THE MERSENNE TWISTER MEMT19937-II RNG.
!
            do k = 1_ki_sel, 2_ki_sel
!
              select case (memti)
!
                case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))          )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmpvec(k) = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmpvec(k) = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmpvec(k) = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag2over:mt_n-2_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmpvec(k) = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-1_ki_sel)
!
                    tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                               iand(memt(0_ki_sel),mt_lower_mask)        )
!
                    memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                mt_mag01(iand(tmp,1_ki_sel))               )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmpvec(k) = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = 0_ki_sel
!
              end select
!
#ifdef _RANDOM_NOINT32
              tmpvec(k) = ibits( tmpvec(k), 0, fbs )
#endif
!
            end do
!
    end select
!$OMP END CRITICAL (ran_num_)
!
    if ( method<=5 ) then
!
#ifdef _RANDOM_NOINT32
            if ( btest(tmp,topbit) ) then
                tmp = ior( tmp, noint32_mask )
            end if
#endif
!
!           CONVERT THE 32-BIT INTEGER INTO A RANDOM FLOATING POINT NUMBER IN THE
!           INTERVAL BETWEEN 0 AND 1.
!
            harvest = half + m_ran_invm32*real( tmp, kind=stnd)
!
    else
!
#ifdef _RANDOM_NOINT32
            if ( btest(tmpvec(1_ki_sel),topbit) ) then
                tmpvec(1_ki_sel) = ior( tmpvec(1_ki_sel), noint32_mask )
            end if
#endif
!
!           CONVERT THE TWO 32-BIT INTEGERS INTO A RANDOM FLOATING POINT NUMBER IN THE
!           INTERVAL BETWEEN 0 AND 1.
!
            harvest = half + m_ran_invm32*real( tmpvec(1_ki_sel),                      kind=stnd) &
                           + m_ran_invmbits*real( ishft( tmpvec(2_ki_sel), num_shift), kind=stnd)
!
    end if
!
!    harvest = max( harvest, zero)
!    if ( harvest>=one )  harvest = nearest( one, -one )
!
#ifdef _RANDOM_WITH0
    harvest = min( harvest, nearest( one, -one ) )
#endif
!
!
! END OF FUNCTION rand_number
! ___________________________
!
    end function rand_number
!
! =========================================================================================
!
    subroutine random_r0( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random number HARVEST between 0 and 1,
!   exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd)
!            A uniformly distributed random real number between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation, this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), dimension(1) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    call random_r1( harvest1 )
!
    harvest = harvest1(1)
!
!
! END OF SUBROUTINE random_r0
! ___________________________
!
    end subroutine random_r0
!
! =========================================================================================
!
    subroutine random_r1( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random vector HARVEST between 0 and 1,
!   exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:)
!            A uniformly distributed random real vector between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation,
!   this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _RANDOM_NOINT32
    use Char_Constants,   only : random_error7, random_error9
#else
    use Char_Constants,   only : random_error1, random_error2,    &
                                 random_error7, random_error8
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)                  :: i, nharvest
    integer(ki_sel)               :: j, k, tmp
    integer(ki_sel), dimension(2) :: tmpvec
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='random_number_'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   QUICK RETURN IF POSSIBLE.
!
    nharvest = size( harvest )
!
    if ( nharvest<=0_i4b ) return
!
    if ( first ) then
#ifdef _RANDOM_NOINT32
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method/=3 ) then
            call merror( name_proc//random_error9 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ibclr(-1_ki_sel, bit_size(tmp)-1_ki_sel)
!
        if ( radix(tmp)/=2 .or. tmp/=huge(tmp) ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
#else
!
!       TEST IF 32-BIT INTEGERS ARE AVAILABLE.
!
        if ( bit_size( tmp )/=fullbitsize ) then
            call merror( name_proc//random_error1 )
        endif
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method<1 .or. method>3 ) then
            call merror( name_proc//random_error2 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ishftc( -8_ki_sel, -3)
!
        if ( radix(tmp)/=2 .or. tmp/=536870911_ki_sel ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
        if ( method/=3 ) then
!
!           TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS.
!
            tmp   = huge( tmp )
!
            if ( integer_not_safe( tmp ) ) then
!
                call merror( name_proc//random_error8 )
!
            endif
!
        endif
#endif
!
    endif
!
!$OMP CRITICAL (ran_num_)
    select case (method)
!
        case (1)
!
!           Marsaglia's KISS RNG.
!
            do i = 1_i4b, nharvest
!
                x = 69069_ki_sel * x + 1327217885_ki_sel
                y = ieor( y, ishft(y, kiss_shift1) )
                y = ieor( y, ishft(y, kiss_shift2) )
                y = ieor( y, ishft(y, kiss_shift3) )
!                y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
                z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
                w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
                tmp = x + y + ishft( z, 16) + w
!
!               CONVERT THE 32-BIT INTEGER INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmp, kind=stnd)
!
            end do
!
        case (2)
!
!           FAST Marsaglia's KISS RNG WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            do i = 1_i4b, nharvest
!
!BUG: For STATPACK with some versions of the GFORTRAN compiler
!     and the optimization level -O3 the following line of code
!     must be replaced with the following commented line below
!     due to a bug in the GFORTRAN compiler.
!
                x2 = x2 + 545925293_ki_sel
!
!                x2 = uiadd( x2, 545925293_ki_sel )
!
                y2 = ieor( y2, ishft(y2, kiss_shift1) )
                y2 = ieor( y2, ishft(y2, kiss_shift2) )
                y2 = ieor( y2, ishft(y2, kiss_shift3) )
!                y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
                tmp = z2 + w2 + c
                z2 = w2
                c  = ishft( tmp, -topbit)
                w2 = iand( tmp, 2147483647_ki_sel)
!
                tmp = x2 + y2 + w2
!
!               CONVERT THE 32-BIT INTEGER INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmp, kind=stnd)
!
            end do
!
        case (3)
!
!           L'Ecuyer's LFSR113 RNG.
!
            do i = 1_i4b, nharvest
!
#ifdef _RANDOM_NOINT32
                tmp = ieor( ishft(s1,6), s1)
                tmp = ishft( ibits( tmp, 0, fbs ), -13)
                tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
                s1  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s2,2), s2)
                tmp = ishft( ibits( tmp, 0, fbs ), -27)
                tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
                s2  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s3,13), s3)
                tmp = ishft( ibits( tmp, 0, fbs ), -21)
                tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
                s3  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s4,3), s4)
                tmp = ishft( ibits( tmp, 0, fbs ), -12)
                tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
                s4  = ibits( tmp, 0, fbs )
#else
                tmp = ishft( ieor( ishft(s1,6), s1), -13)
                s1  = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
                tmp = ishft( ieor( ishft(s2,2), s2), -27)
                s2  = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
                tmp = ishft( ieor( ishft(s3,13), s3), -21)
                s3  = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
                tmp = ishft( ieor( ishft(s4,3), s4), -12)
                s4  = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
                tmp = ieor( ieor( ieor(s1,s2), s3), s4)
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmp,topbit) ) then
                    tmp = ior( tmp, noint32_mask )
                end if
#endif
!
!               CONVERT THE 32-BIT INTEGER INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmp, kind=stnd)
!
            end do
!
        case (4)
!
!           MERSENNE TWISTER MT19937 RNG.
!
            do i = 1_i4b, nharvest
!
                if ( mti>=mt_n ) then
!
!                   GENERATE mt_n WORDS AT ONE TIME.
!
                    do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)),   &
                                      mt_mag01(iand(tmp,1_ki_sel))      )
!
                    end do
!
                    do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)), ishft(tmp,-1)),  &
                                      mt_mag01(iand(tmp,1_ki_sel))             )
!
                    end do
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(0_ki_sel),mt_lower_mask)   )
!
                    mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)),  &
                                  mt_mag01(iand(tmp,1_ki_sel))            )
!
                    mti = 0_ki_sel
!
                endif
!
                tmp = mt(mti)
                mti = mti + 1_ki_sel
!
!               TEMPERING PHASE.
!
                tmp = ieor( tmp, ishft(tmp,mt_shift0) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
                tmp = ieor( tmp, ishft(tmp,mt_shift1) )
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmp,topbit) ) then
                    tmp = ior( tmp, noint32_mask )
                end if
#endif
!
!               CONVERT THE 32-BIT INTEGER INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmp, kind=stnd)
!
            end do
!
        case (5)
!
!           MERSENNE TWISTER MEMT19937-II RNG.
!
            do i = 1_i4b, nharvest
!
                select case (memti)
!
                    case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))          )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag2over:mt_n-2_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-1_ki_sel)
!
                        tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                                   iand(memt(0_ki_sel),mt_lower_mask)        )
!
                        memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                    mt_mag01(iand(tmp,1_ki_sel))               )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = 0_ki_sel
!
                end select
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmp,topbit) ) then
                    tmp = ior( tmp, noint32_mask )
                else
                    tmp = ibits( tmp, 0, fbs )
                end if
#endif
!
!               CONVERT THE 32-BIT INTEGER INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmp, kind=stnd)
!
            end do
!
        case (6)
!
!           EXTENDED PRECISION VERSION OF Marsaglia's KISS RNG.
!
            do i = 1_i4b, nharvest
!
                do k = 1_ki_sel, 2_ki_sel
!
                    x = 69069_ki_sel * x + 1327217885_ki_sel
                    y = ieor( y, ishft(y, kiss_shift1) )
                    y = ieor( y, ishft(y, kiss_shift2) )
                    y = ieor( y, ishft(y, kiss_shift3) )
!                    y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
                    z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
                    w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
                    tmpvec(k) = x + y + ishft( z, 16) + w
!
                end do
!
!               CONVERT THE TWO 32-BIT INTEGERS INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmpvec(1_ki_sel),                      kind=stnd) &
                                  + m_ran_invmbits*real( ishft( tmpvec(2_ki_sel), num_shift), kind=stnd)
!
            end do
!
        case (7)
!
!           EXTENDED PRECISION VERSION OF THE FAST Marsaglia's KISS RNG,
!           WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            do i = 1_i4b, nharvest
!
                do k = 1_ki_sel, 2_ki_sel
!
!BUG: For STATPACK with some versions of the GFORTRAN compiler
!     and the optimization level -O3 the following line of code
!     must be replaced with the following commented line below
!     due to a bug in the GFORTRAN compiler.
!
                    x2 = x2 + 545925293_ki_sel
!
!                    x2 = uiadd( x2, 545925293_ki_sel )
!
                    y2 = ieor( y2, ishft(y2, kiss_shift1) )
                    y2 = ieor( y2, ishft(y2, kiss_shift2) )
                    y2 = ieor( y2, ishft(y2, kiss_shift3) )
!                    y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
                    tmp = z2 + w2 + c
                    z2 = w2
                    c  = ishft( tmp, -topbit)
                    w2 = iand( tmp, 2147483647_ki_sel)
!
                    tmpvec(k) = x2 + y2 + w2
!
                end do
!
!               CONVERT THE TWO 32-BIT INTEGERS INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmpvec(1_ki_sel),                      kind=stnd) &
                                  + m_ran_invmbits*real( ishft( tmpvec(2_ki_sel), num_shift), kind=stnd)
!
            end do
!
        case (8)
!
!           EXTENDED PRECISION VERSION OF L'Ecuyer's LFSR113 RNG.
!
            do i = 1_i4b, nharvest
!
                do k = 1_ki_sel, 2_ki_sel
!
#ifdef _RANDOM_NOINT32
                    tmp = ieor( ishft(s1,6), s1)
                    tmp = ishft( ibits( tmp, 0, fbs ), -13)
                    tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
                    s1  = ibits( tmp, 0, fbs )
!
                    tmp = ieor( ishft(s2,2), s2)
                    tmp = ishft( ibits( tmp, 0, fbs ), -27)
                    tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
                    s2  = ibits( tmp, 0, fbs )
!
                    tmp = ieor( ishft(s3,13), s3)
                    tmp = ishft( ibits( tmp, 0, fbs ), -21)
                    tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
                    s3  = ibits( tmp, 0, fbs )
!
                    tmp = ieor( ishft(s4,3), s4)
                    tmp = ishft( ibits( tmp, 0, fbs ), -12)
                    tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
                    s4  = ibits( tmp, 0, fbs )
#else
                    tmp = ishft( ieor( ishft(s1,6), s1), -13)
                    s1  = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
                    tmp = ishft( ieor( ishft(s2,2), s2), -27)
                    s2  = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
                    tmp = ishft( ieor( ishft(s3,13), s3), -21)
                    s3  = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
                    tmp = ishft( ieor( ishft(s4,3), s4), -12)
                    s4  = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
                    tmpvec(k) = ieor( ieor( ieor(s1,s2), s3), s4)
!
                end do
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmpvec(1_ki_sel),topbit) ) then
                    tmpvec(1_ki_sel) = ior( tmpvec(1_ki_sel), noint32_mask )
                end if
#endif
!
!               CONVERT THE TWO 32-BIT INTEGERS INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmpvec(1_ki_sel),                      kind=stnd) &
                                  + m_ran_invmbits*real( ishft( tmpvec(2_ki_sel), num_shift), kind=stnd)
!
            end do
!
        case (9)
!
!           EXTENDED PRECISION VERSION OF THE MERSENNE TWISTER MT19937 RNG.
!
            do i = 1_i4b, nharvest
!
                do k = 1_ki_sel, 2_ki_sel
!
                    if ( mti>=mt_n ) then
!
!                       GENERATE mt_n WORDS AT ONE TIME.
!
                        do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
                            tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                         iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                            mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)), &
                                          mt_mag01(iand(tmp,1_ki_sel))    )
!
                        end do
!
                        do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                            tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                         iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                            mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)),                     &
                                          ishft(tmp,-1)),mt_mag01(iand(tmp,1_ki_sel)) )
!
                        end do
!
                        tmp   = ior( iand(mt(mt_n-1_ki_sel),mt_upper_mask),    &
                                     iand(mt(0_ki_sel),mt_lower_mask)          )
!
                        mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)), &
                                      mt_mag01(iand(tmp,1_ki_sel))           )
!
                        mti = 0_ki_sel
!
                    endif
!
                    tmp = mt(mti)
                    mti = mti + 1_ki_sel
!
!                   TEMPERING PHASE.
!
                    tmp       = ieor( tmp, ishft(tmp,mt_shift0) )
                    tmp       = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
                    tmp       = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
                    tmpvec(k) = ieor( tmp, ishft(tmp,mt_shift1) )
!
                end do
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmpvec(1_ki_sel),topbit) ) then
                    tmpvec(1_ki_sel) = ior( tmpvec(1_ki_sel), noint32_mask )
                end if
#endif
!
!               CONVERT THE TWO 32-BIT INTEGERS INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmpvec(1_ki_sel),                      kind=stnd) &
                                  + m_ran_invmbits*real( ishft( tmpvec(2_ki_sel), num_shift), kind=stnd)
!
            end do
!
!
        case (10)
!
!           EXTENDED PRECISION VERSION OF THE MERSENNE TWISTER MEMT19937-II RNG.
!
            do i = 1_i4b, nharvest
!
                do k = 1_ki_sel, 2_ki_sel
!
                  select case (memti)
!
                    case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))          )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmpvec(k) = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmpvec(k) = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmpvec(k) = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag2over:mt_n-2_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmpvec(k) = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-1_ki_sel)
!
                        tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                                   iand(memt(0_ki_sel),mt_lower_mask)        )
!
                        memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                    mt_mag01(iand(tmp,1_ki_sel))               )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmpvec(k) = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = 0_ki_sel
!
                  end select
!
#ifdef _RANDOM_NOINT32
                  tmpvec(k) = ibits( tmpvec(k), 0, fbs )
#endif
!
                end do
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmpvec(1_ki_sel),topbit) ) then
                    tmpvec(1_ki_sel) = ior( tmpvec(1_ki_sel), noint32_mask )
                end if
#endif
!
!               CONVERT THE TWO 32-BIT INTEGERS INTO A RANDOM FLOATING POINT NUMBER IN THE
!               INTERVAL BETWEEN 0 AND 1.
!
                harvest(i) = half + m_ran_invm32*real( tmpvec(1_ki_sel),                      kind=stnd) &
                                  + m_ran_invmbits*real( ishft( tmpvec(2_ki_sel), num_shift), kind=stnd)
!
            end do
!
    end select
!$OMP END CRITICAL (ran_num_)
!
!    where( harvest(:)<zero )  harvest(:) = zero
!    where( harvest(:)>=one )  harvest(:) = nearest( one, -one )
!
#ifdef _RANDOM_WITH0
    harvest(:) = min( harvest(:), nearest( one, -one ) )
#endif
!
!
! END OF SUBROUTINE random_r1
! ___________________________
!
    end subroutine random_r1
!
! =========================================================================================
!
    subroutine random_r2( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random matrix HARVEST between 0 and 1,
!   exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:)
!            A uniformly distributed random real matrix between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation, this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 2 )
!
    do i = 1_i4b, nharvest
!
        call random_r1( harvest(:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_r2
! ___________________________
!
    end subroutine random_r2
!
! =========================================================================================
!
    subroutine random_r3( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random array of dimension 3 HARVEST
!   between 0 and 1, exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:,:)
!            A uniformly distributed random real array of dimension 3 between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation, this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 3 )
!
    do i = 1_i4b, nharvest
!
        call random_r2( harvest(:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_r3
! ___________________________
!
    end subroutine random_r3
!
! =========================================================================================
!
    subroutine random_r4( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random array of dimension 4 HARVEST
!   between 0 and 1, exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:,:,:)
!            A uniformly distributed random real array of dimension 4 between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation, this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:,:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 4 )
!
    do i = 1_i4b, nharvest
!
        call random_r3( harvest(:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_r4
! ___________________________
!
    end subroutine random_r4
!
! =========================================================================================
!
    subroutine random_r5( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random array of dimension 5 HARVEST
!   between 0 and 1, exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:,:,:,:)
!            A uniformly distributed random real array of dimension 5 between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation, this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:,:,:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 5 )
!
    do i = 1_i4b, nharvest
!
        call random_r4( harvest(:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_r5
! ___________________________
!
    end subroutine random_r5
!
! =========================================================================================
!
    subroutine random_r6( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random array of dimension 6 HARVEST
!   between 0 and 1, exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:,:,:,:,:)
!            A uniformly distributed random real array of dimension 6 between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation, this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:,:,:,:,:,:)    :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 6 )
!
    do i = 1_i4b, nharvest
!
        call random_r5( harvest(:,:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_r6
! ___________________________
!
    end subroutine random_r6
!
! =========================================================================================
!
    subroutine random_r7( harvest )
!
! purpose
! _______
!
!   This subroutine returns a uniformly distributed random array of dimension 7 HARVEST
!   between 0 and 1, exclusive of the two endpoints 0 and 1.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:,:,:,:,:,:)
!            A uniformly distributed random real array of dimension 7 between 0 and 1.
!
!
! Further Details
! _______________
!
!   If the CPP macro _RANDOM_WITH0 is used during compilation, this routine may return 0 value.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),  intent(out), dimension(:,:,:,:,:,:,:)  :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 7 )
!
    do i = 1_i4b, nharvest
!
        call random_r6( harvest(:,:,:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_r7
! ___________________________
!
    end subroutine random_r7
!
!
! =========================================================================================
!                            RANDOM 32-BIT INTEGER ROUTINES
! =========================================================================================
!
!
    function rand_integer32() result( harvest )
!
! purpose
! _______
!
!   This function returns a random integer in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integer is equivalent to
!   a signed 32-bit integer.
!
!
! Arguments
! _________
!
!   None
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _RANDOM_NOINT32
    use Char_Constants,   only : random_error7, random_error9
#else
    use Char_Constants,   only : random_error1, random_error2,    &
                                 random_error7, random_error8
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel)  :: j, tmp
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='rand_integer32'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( first ) then
#ifdef _RANDOM_NOINT32
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method/=3 ) then
            call merror( name_proc//random_error9 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ibclr(-1_ki_sel, bit_size(tmp)-1_ki_sel)
!
        if ( radix(tmp)/=2 .or. tmp/=huge(tmp) ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
#else
!
!       TEST IF 32-BIT INTEGERS ARE AVAILABLE.
!
        if ( bit_size( tmp )/=fullbitsize ) then
            call merror( name_proc//random_error1 )
        endif
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method<1 .or. method>3 ) then
            call merror( name_proc//random_error2 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ishftc( -8_ki_sel, -3)
!
        if ( radix(tmp)/=2 .or. tmp/=536870911_ki_sel ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
        if ( method/=3 ) then
!
!           TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS.
!
            tmp   = huge( tmp )
!
            if ( integer_not_safe( tmp ) ) then
!
                call merror( name_proc//random_error8 )
!
            endif
!
        endif
#endif
!
    endif
!
!$OMP CRITICAL (ran_num_)
    select case (method)
!
        case (1,6)
!
!           Marsaglia's KISS RNG.
!
            x = 69069_ki_sel * x + 1327217885_ki_sel
            y = ieor( y, ishft(y, kiss_shift1) )
            y = ieor( y, ishft(y, kiss_shift2) )
            y = ieor( y, ishft(y, kiss_shift3) )
!            y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
            z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
            w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
            tmp = x + y + ishft( z, 16) + w
!
        case (2,7)
!
!           FAST Marsaglia's KISS RNG WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            x2 = x2 + 545925293_ki_sel
!            x2 = uiadd( x2, 545925293_ki_sel )
            y2 = ieor( y2, ishft(y2, kiss_shift1) )
            y2 = ieor( y2, ishft(y2, kiss_shift2) )
            y2 = ieor( y2, ishft(y2, kiss_shift3) )
!            y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
            tmp = z2 + w2 + c
            z2 = w2
            c  = ishft( tmp, -topbit)
            w2 = iand( tmp, 2147483647_ki_sel)
!
            tmp = x2 + y2 + w2
!
        case (3,8)
!
!           L'Ecuyer's LFSR113 RNG.
!
#ifdef _RANDOM_NOINT32
            tmp = ieor( ishft(s1,6), s1)
            tmp = ishft( ibits( tmp, 0, fbs ), -13)
            tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
            s1  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s2,2), s2)
            tmp = ishft( ibits( tmp, 0, fbs ), -27)
            tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
            s2  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s3,13), s3)
            tmp = ishft( ibits( tmp, 0, fbs ), -21)
            tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
            s3  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s4,3), s4)
            tmp = ishft( ibits( tmp, 0, fbs ), -12)
            tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
            s4  = ibits( tmp, 0, fbs )
#else
            tmp  = ishft( ieor( ishft(s1,6), s1), -13)
            s1 = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
            tmp  = ishft( ieor( ishft(s2,2), s2), -27)
            s2 = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
            tmp  = ishft( ieor( ishft(s3,13), s3), -21)
            s3 = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
            tmp  = ishft( ieor( ishft(s4,3), s4), -12)
            s4 = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
            tmp = ieor( ieor( ieor(s1,s2), s3), s4)
!
        case (4,9)
!
!           MERSENNE TWISTER MT19937 RNG.
!
            if ( mti>=mt_n ) then
!
                do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
!                   GENERATE mt_n WORDS AT ONE TIME.
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                    mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)),   &
                                  mt_mag01(iand(tmp,1_ki_sel))      )
!
                end do
!
                do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                    mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)), ishft(tmp,-1)),  &
                                  mt_mag01(iand(tmp,1_ki_sel))             )
!
                end do
!
                tmp   = ior( iand(mt(j),mt_upper_mask),        &
                             iand(mt(0_ki_sel),mt_lower_mask)  )
!
                mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)), &
                              mt_mag01(iand(tmp,1_ki_sel))           )
!
                mti = 0_ki_sel
!
            endif
!
            tmp = mt(mti)
            mti = mti + 1_ki_sel
!
!           TEMPERING PHASE.
!
            tmp = ieor( tmp, ishft(tmp,mt_shift0) )
            tmp = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
            tmp = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
            tmp = ieor( tmp, ishft(tmp,mt_shift1) )
!
        case (5,10)
!
!           MERSENNE TWISTER MEMT19937-II RNG.
!
            select case (memti)
!
                case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))          )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag2over:mt_n-2_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-1_ki_sel)
!
                    tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                               iand(memt(0_ki_sel),mt_lower_mask)        )
!
                    memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                mt_mag01(iand(tmp,1_ki_sel))               )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = 0_ki_sel
!
            end select
!
#ifdef _RANDOM_NOINT32
            tmp = ibits( tmp, 0, fbs )
#endif
!
    end select
!$OMP END CRITICAL (ran_num_)
!
#ifdef _RANDOM_NOINT32
    if ( btest(tmp,topbit) ) then
        tmp = ior( tmp, noint32_mask )
    end if
#endif
!
    harvest = tmp
!
! END OF FUNCTION rand_integer32
! ______________________________
!
    end function rand_integer32
!
! =========================================================================================
!
    subroutine random_i32_0( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random integer in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integer is equivalent to
!   a signed 32-bit integer.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b)
!            A random integer in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b), dimension(1) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    call random_i32_1( harvest1 )
!
    harvest = harvest1(1)
!
!
! END OF SUBROUTINE random_i32_0
! ______________________________
!
    end subroutine random_i32_0
!
! =========================================================================================
!
    subroutine random_i32_1( harvest )
!
! purpose
! _______
!
!   This subroutine returns a vector of random integers in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   signed 32-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:)
!            A vector of random integers in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _RANDOM_NOINT32
    use Char_Constants,   only : random_error7, random_error9
#else
    use Char_Constants,   only : random_error1, random_error2,    &
                                 random_error7, random_error8
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)    :: i, nharvest
    integer(ki_sel) :: j, tmp
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='random_integer32_'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   QUICK RETURN IF POSSIBLE.
!
    nharvest = size( harvest )
!
    if ( nharvest<=0_i4b ) return
!
    if ( first ) then
#ifdef _RANDOM_NOINT32
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method/=3 ) then
            call merror( name_proc//random_error9 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ibclr(-1_ki_sel, bit_size(tmp)-1_ki_sel)
!
        if ( radix(tmp)/=2 .or. tmp/=huge(tmp) ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
#else
!
!       TEST IF 32-BIT INTEGERS ARE AVAILABLE.
!
        if ( bit_size( tmp )/=fullbitsize ) then
            call merror( name_proc//random_error1 )
        endif
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method<1 .or. method>3 ) then
            call merror( name_proc//random_error2 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ishftc( -8_ki_sel, -3)
!
        if ( radix(tmp)/=2 .or. tmp/=536870911_ki_sel ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
        if ( method/=3 ) then
!
!           TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS.
!
            tmp   = huge( tmp )
!
            if ( integer_not_safe( tmp ) ) then
!
                call merror( name_proc//random_error8 )
!
            endif
!
        endif
#endif
!
    endif
!
!$OMP CRITICAL (ran_num_)
    select case (method)
!
        case (1,6)
!
!           Marsaglia's KISS RNG.
!
            do i = 1_i4b, nharvest
!
                x = 69069_ki_sel * x + 1327217885_ki_sel
                y = ieor( y, ishft(y, kiss_shift1) )
                y = ieor( y, ishft(y, kiss_shift2) )
                y = ieor( y, ishft(y, kiss_shift3) )
!                y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
                z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
                w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
!               GENERATE A RANDOM SIGNED 32-BIT INTEGER.
!
                harvest(i) = x + y + ishft( z, 16) + w
!
            end do
!
        case (2,7)
!
!           FAST Marsaglia's KISS RNG WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            do i = 1_i4b, nharvest
!
!BUG: For STATPACK with some versions of the GFORTRAN compiler
!     and the optimization level -O3 the following line of code
!     must be replaced with the following commented line below
!     due to a bug in the GFORTRAN compiler.
!
                x2 = x2 + 545925293_ki_sel
!
!                x2 = uiadd( x2, 545925293_ki_sel )
!
                y2 = ieor( y2, ishft(y2, kiss_shift1) )
                y2 = ieor( y2, ishft(y2, kiss_shift2) )
                y2 = ieor( y2, ishft(y2, kiss_shift3) )
!                y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
                tmp = z2 + w2 + c
                z2 = w2
                c  = ishft( tmp, -topbit)
                w2 = iand( tmp, 2147483647_ki_sel)
!
!               GENERATE A RANDOM SIGNED 32-BIT INTEGER.
!
                harvest(i) = x2 + y2 + w2
!
            end do
!
        case (3,8)
!
!           L'Ecuyer's LFSR113 RNG.
!
            do i = 1_i4b, nharvest
!
#ifdef _RANDOM_NOINT32
                tmp = ieor( ishft(s1,6), s1)
                tmp = ishft( ibits( tmp, 0, fbs ), -13)
                tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
                s1  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s2,2), s2)
                tmp = ishft( ibits( tmp, 0, fbs ), -27)
                tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
                s2  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s3,13), s3)
                tmp = ishft( ibits( tmp, 0, fbs ), -21)
                tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
                s3  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s4,3), s4)
                tmp = ishft( ibits( tmp, 0, fbs ), -12)
                tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
                s4  = ibits( tmp, 0, fbs )
#else
                tmp = ishft( ieor( ishft(s1,6), s1), -13)
                s1  = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
                tmp = ishft( ieor( ishft(s2,2), s2), -27)
                s2  = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
                tmp = ishft( ieor( ishft(s3,13), s3), -21)
                s3  = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
                tmp = ishft( ieor( ishft(s4,3), s4), -12)
                s4  = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
                tmp = ieor( ieor( ieor(s1,s2), s3), s4)
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmp,topbit) ) then
                    tmp = ior( tmp, noint32_mask )
                end if
#endif
!
!               GENERATE A RANDOM SIGNED 32-BIT INTEGER.
!
                harvest(i) = tmp
!
            end do
!
        case (4,9)
!
!           MERSENNE TWISTER MT19937 RNG.
!
            do i = 1_i4b, nharvest
!
                if ( mti>=mt_n ) then
!
!                   GENERATE mt_n WORDS AT ONE TIME.
!
                    do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)),   &
                                      mt_mag01(iand(tmp,1_ki_sel))      )
!
                    end do
!
                    do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)), ishft(tmp,-1)),  &
                                      mt_mag01(iand(tmp,1_ki_sel))             )
!
                    end do
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(0_ki_sel),mt_lower_mask)   )
!
                    mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)),  &
                                  mt_mag01(iand(tmp,1_ki_sel))            )
!
                    mti = 0_ki_sel
!
                endif
!
                tmp = mt(mti)
                mti = mti + 1_ki_sel
!
!               TEMPERING PHASE.
!
                tmp = ieor( tmp, ishft(tmp,mt_shift0) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
                tmp = ieor( tmp, ishft(tmp,mt_shift1) )
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmp,topbit) ) then
                    tmp = ior( tmp, noint32_mask )
                end if
#endif
!
!               GENERATE A RANDOM SIGNED 32-BIT INTEGER.
!
                harvest(i) = tmp
!
            end do
!
        case (5,10)
!
!           MERSENNE TWISTER MT19937-II RNG.
!
            do i = 1_i4b, nharvest
!
                select case (memti)
!
                    case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))          )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag2over:mt_n-2_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-1_ki_sel)
!
                        tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                                   iand(memt(0_ki_sel),mt_lower_mask)        )
!
                        memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                    mt_mag01(iand(tmp,1_ki_sel))               )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = 0_ki_sel
!
                end select
!
!               GENERATE A RANDOM SIGNED 32-BIT INTEGER.
!
#ifdef _RANDOM_NOINT32
                if ( btest(tmp,topbit) ) then
                    harvest(i) = ior( tmp, noint32_mask )
                else
                    harvest(i) = ibits( tmp, 0, fbs )
                end if
#else
                harvest(i) = tmp
#endif
!
            end do
!
    end select
!$OMP END CRITICAL (ran_num_)
!
!
! END OF SUBROUTINE random_i32_1
! ______________________________
!
    end subroutine random_i32_1
!
! =========================================================================================
!
    subroutine random_i32_2( harvest )
!
! purpose
! _______
!
!   This subroutine returns a matrix of random integers in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   signed 32-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:)
!            A matrix of random integers in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 2 )
!
    do i = 1_i4b, nharvest
!
        call random_i32_1( harvest(:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i32_2
! ______________________________
!
    end subroutine random_i32_2
!
! =========================================================================================
!
    subroutine random_i32_3( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   signed 32-bit integers.
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:)
!            An array of random integers in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 3 )
!
    do i = 1_i4b, nharvest
!
        call random_i32_2( harvest(:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i32_3
! ______________________________
!
    end subroutine random_i32_3
!
! =========================================================================================
!
    subroutine random_i32_4( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   signed 32-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:)
!            An array of random integers in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 4 )
!
    do i = 1_i4b, nharvest
!
        call random_i32_3( harvest(:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i32_4
! ______________________________
!
    end subroutine random_i32_4
!
! =========================================================================================
!
    subroutine random_i32_5( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   signed 32-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:,:)
!            An array of random integers in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 5 )
!
    do i = 1_i4b, nharvest
!
        call random_i32_4( harvest(:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i32_5
! ______________________________
!
    end subroutine random_i32_5
!
! =========================================================================================
!
    subroutine random_i32_6( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   signed 32-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:,:,:)
!            An array of random integers in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:,:,:,:)    :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 6 )
!
    do i = 1_i4b, nharvest
!
        call random_i32_5( harvest(:,:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i32_6
! ______________________________
!
    end subroutine random_i32_6
!
! =========================================================================================
!
    subroutine random_i32_7( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (-2147483648,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   signed 32-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:,:,:,:)
!            An array of random integers in the interval (-2147483648,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b),  intent(out), dimension(:,:,:,:,:,:,:)  :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 7 )
!
    do i = 1_i4b, nharvest
!
        call random_i32_6( harvest(:,:,:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i32_7
! ______________________________
!
    end subroutine random_i32_7
!
!
! =========================================================================================
!                            RANDOM UNSIGNED 31-BIT INTEGER ROUTINES
! =========================================================================================
!
!
    function rand_integer31() result( harvest )
!
! purpose
! _______
!
!   This function returns a random integer in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integer is equivalent to
!   an unsigned 31-bit integer.
!
!
! Arguments
! _________
!
!   None
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _RANDOM_NOINT32
    use Char_Constants,   only : random_error7, random_error9
#else
    use Char_Constants,   only : random_error1, random_error2,    &
                                 random_error7, random_error8
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel)  :: j, tmp
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='rand_integer31'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( first ) then
#ifdef _RANDOM_NOINT32
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method/=3 ) then
            call merror( name_proc//random_error9 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ibclr(-1_ki_sel, bit_size(tmp)-1_ki_sel)
!
        if ( radix(tmp)/=2 .or. tmp/=huge(tmp) ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
#else
!
!       TEST IF 32-BIT INTEGERS ARE AVAILABLE.
!
        if ( bit_size( tmp )/=fullbitsize ) then
            call merror( name_proc//random_error1 )
        endif
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method<1 .or. method>3 ) then
            call merror( name_proc//random_error2 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ishftc( -8_ki_sel, -3)
!
        if ( radix(tmp)/=2 .or. tmp/=536870911_ki_sel ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
        if ( method/=3 ) then
!
!           TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS.
!
            tmp   = huge( tmp )
!
            if ( integer_not_safe( tmp ) ) then
!
                call merror( name_proc//random_error8 )
!
            endif
!
        endif
#endif
!
    endif
!
!$OMP CRITICAL (ran_num_)
    select case (method)
!
        case (1,6)
!
!           Marsaglia's KISS RNG.
!
            x = 69069_ki_sel * x + 1327217885_ki_sel
            y = ieor( y, ishft(y, kiss_shift1) )
            y = ieor( y, ishft(y, kiss_shift2) )
            y = ieor( y, ishft(y, kiss_shift3) )
!            y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
            z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
            w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
            tmp = x + y + ishft( z, 16) + w
!
        case (2,7)
!
!           FAST Marsaglia's KISS RNG WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            x2 = x2 + 545925293_ki_sel
!            x2 = uiadd( x2, 545925293_ki_sel )
            y2 = ieor( y2, ishft(y2, kiss_shift1) )
            y2 = ieor( y2, ishft(y2, kiss_shift2) )
            y2 = ieor( y2, ishft(y2, kiss_shift3) )
!            y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
            tmp = z2 + w2 + c
            z2 = w2
            c  = ishft( tmp, -topbit)
            w2 = iand( tmp, 2147483647_ki_sel)
!
            tmp = x2 + y2 + w2
!
        case (3,8)
!
!           L'Ecuyer's LFSR113 RNG.
!
#ifdef _RANDOM_NOINT32
            tmp = ieor( ishft(s1,6), s1)
            tmp = ishft( ibits( tmp, 0, fbs ), -13)
            tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
            s1  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s2,2), s2)
            tmp = ishft( ibits( tmp, 0, fbs ), -27)
            tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
            s2  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s3,13), s3)
            tmp = ishft( ibits( tmp, 0, fbs ), -21)
            tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
            s3  = ibits( tmp, 0, fbs )
!
            tmp = ieor( ishft(s4,3), s4)
            tmp = ishft( ibits( tmp, 0, fbs ), -12)
            tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
            s4  = ibits( tmp, 0, fbs )
#else
            tmp = ishft( ieor( ishft(s1,6), s1), -13)
            s1  = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
            tmp = ishft( ieor( ishft(s2,2), s2), -27)
            s2  = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
            tmp = ishft( ieor( ishft(s3,13), s3), -21)
            s3  = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
            tmp = ishft( ieor( ishft(s4,3), s4), -12)
            s4  = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
            tmp = ieor( ieor( ieor(s1,s2), s3), s4)
!
        case (4,9)
!
!           MERSENNE TWISTER MT19937 RNG.
!
            if ( mti>=mt_n ) then
!
                do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
!                   GENERATE mt_n WORDS AT ONE TIME.
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                    mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)),   &
                                  mt_mag01(iand(tmp,1_ki_sel))      )
!
                end do
!
                do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                    mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)), ishft(tmp,-1)),  &
                                  mt_mag01(iand(tmp,1_ki_sel))             )
!
                end do
!
                tmp   = ior( iand(mt(j),mt_upper_mask),        &
                             iand(mt(0_ki_sel),mt_lower_mask)  )
!
                mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)), &
                              mt_mag01(iand(tmp,1_ki_sel))           )
!
                mti = 0_ki_sel
!
            endif
!
            tmp = mt(mti)
            mti = mti + 1_ki_sel
!
!           TEMPERING PHASE.
!
            tmp = ieor( tmp, ishft(tmp,mt_shift0) )
            tmp = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
            tmp = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
            tmp = ieor( tmp, ishft(tmp,mt_shift1) )
!
        case (5,10)
!
!           MERSENNE TWISTER MEMT19937-II RNG.
!
            select case (memti)
!
                case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))          )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (memt_lag2over:mt_n-2_ki_sel)
!
                    tmp = ior( iand(memt(memti),mt_upper_mask),            &
                               iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                    memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                        mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = memti + 1_ki_sel
!
                case (mt_n-1_ki_sel)
!
                    tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                               iand(memt(0_ki_sel),mt_lower_mask)        )
!
                    memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                mt_mag01(iand(tmp,1_ki_sel))               )
!
!                   TEMPERING PHASE.
!
                    tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                    tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                    tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                    memti = 0_ki_sel
!
            end select
!
#ifdef _RANDOM_NOINT32
            tmp = ibits( tmp, 0, fbs )
#endif
!
    end select
!$OMP END CRITICAL (ran_num_)
!
!   GENERATE A RANDOM UNSIGNED 31-BIT INTEGER.
!
    harvest = ibclr( tmp, topbit )
!
! END OF FUNCTION rand_integer31
! ______________________________
!
    end function rand_integer31
!
! =========================================================================================
!
    subroutine random_i31_0( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random integer in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integer is equivalent to
!   an unsigned 31-bit integer.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b)
!            A random integer in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b), dimension(1) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    call random_i31_1( harvest1 )
!
    harvest = harvest1(1)
!
!
! END OF SUBROUTINE random_i31_0
! ______________________________
!
    end subroutine random_i31_0
!
! =========================================================================================
!
    subroutine random_i31_1( harvest )
!
! purpose
! _______
!
!   This subroutine returns a vector of random integers in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   unsigned 31-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:)
!            A vector of random integers in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _RANDOM_NOINT32
    use Char_Constants,   only : random_error7, random_error9
#else
    use Char_Constants,   only : random_error1, random_error2,    &
                                 random_error7, random_error8
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b)    :: i, nharvest
    integer(ki_sel) :: j, tmp
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='random_integer31_'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   QUICK RETURN IF POSSIBLE.
!
    nharvest = size( harvest )
!
    if ( nharvest<=0_i4b ) return
!
    if ( first ) then
#ifdef _RANDOM_NOINT32
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method/=3 ) then
            call merror( name_proc//random_error9 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ibclr(-1_ki_sel, bit_size(tmp)-1_ki_sel)
!
        if ( radix(tmp)/=2 .or. tmp/=huge(tmp) ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
#else
!
!       TEST IF 32-BIT INTEGERS ARE AVAILABLE.
!
        if ( bit_size( tmp )/=fullbitsize ) then
            call merror( name_proc//random_error1 )
        endif
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method<1 .or. method>3 ) then
            call merror( name_proc//random_error2 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ishftc( -8_ki_sel, -3)
!
        if ( radix(tmp)/=2 .or. tmp/=536870911_ki_sel ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
        if ( method/=3 ) then
!
!           TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS.
!
            tmp   = huge( tmp )
!
            if ( integer_not_safe( tmp ) ) then
!
                call merror( name_proc//random_error8 )
!
            endif
!
        endif
#endif
!
    endif
!
!$OMP CRITICAL (ran_num_)
    select case (method)
!
        case (1,6)
!
!           Marsaglia's KISS RNG.
!
            do i = 1_i4b, nharvest
!
                x = 69069_ki_sel * x + 1327217885_ki_sel
                y = ieor( y, ishft(y, kiss_shift1) )
                y = ieor( y, ishft(y, kiss_shift2) )
                y = ieor( y, ishft(y, kiss_shift3) )
!                y = m( m( m( y, kiss_shift1), kiss_shift2), kiss_shift3)
                z = 18000_ki_sel * iand( z, 65535_ki_sel) + ishft( z, -16)
                w = 30903_ki_sel * iand( w, 65535_ki_sel) + ishft( w, -16)
!
                tmp = x + y + ishft( z, 16) + w
!
!               GENERATE A RANDOM UNSIGNED 31-BIT INTEGER.
!
                harvest(i) = ibclr( tmp, topbit)
!
            end do
!
        case (2,7)
!
!           FAST Marsaglia's KISS RNG WHICH USES ONLY ADD, SHIFT, EXCLUSIVE-OR AND "AND" OPERATIONS.
!
            do i = 1_i4b, nharvest
!
!BUG: For STATPACK with some versions of the GFORTRAN compiler
!     and the optimization level -O3 the following line of code
!     must be replaced with the following commented line below
!     due to a bug in the GFORTRAN compiler.
!
                x2 = x2 + 545925293_ki_sel
!
!                x2 = uiadd( x2, 545925293_ki_sel )
!
                y2 = ieor( y2, ishft(y2, kiss_shift1) )
                y2 = ieor( y2, ishft(y2, kiss_shift2) )
                y2 = ieor( y2, ishft(y2, kiss_shift3) )
!                y2 = m( m( m(y2, kiss_shift1), kiss_shift2), kiss_shift3)
                tmp = z2 + w2 + c
                z2 = w2
                c  = ishft( tmp, -topbit)
                w2 = iand( tmp, 2147483647_ki_sel)
!
                tmp = x2 + y2 + w2
!
!               GENERATE A RANDOM UNSIGNED 31-BIT INTEGER.
!
                harvest(i) = ibclr( tmp, topbit)
!
            end do
!
        case (3,8)
!
!           L'Ecuyer's LFSR113 RNG.
!
            do i = 1_i4b, nharvest
!
#ifdef _RANDOM_NOINT32
                tmp = ieor( ishft(s1,6), s1)
                tmp = ishft( ibits( tmp, 0, fbs ), -13)
                tmp = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
                s1  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s2,2), s2)
                tmp = ishft( ibits( tmp, 0, fbs ), -27)
                tmp = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
                s2  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s3,13), s3)
                tmp = ishft( ibits( tmp, 0, fbs ), -21)
                tmp = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
                s3  = ibits( tmp, 0, fbs )
!
                tmp = ieor( ishft(s4,3), s4)
                tmp = ishft( ibits( tmp, 0, fbs ), -12)
                tmp = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
                s4  = ibits( tmp, 0, fbs )
#else
                tmp = ishft( ieor( ishft(s1,6), s1), -13)
                s1  = ieor( ishft( iand(s1,-2_ki_sel), 18), tmp)
!
                tmp = ishft( ieor( ishft(s2,2), s2), -27)
                s2  = ieor( ishft( iand(s2,-8_ki_sel), 2), tmp)
!
                tmp = ishft( ieor( ishft(s3,13), s3), -21)
                s3  = ieor( ishft( iand(s3,-16_ki_sel), 7), tmp)
!
                tmp = ishft( ieor( ishft(s4,3), s4), -12)
                s4  = ieor( ishft( iand(s4,-128_ki_sel), 13), tmp)
#endif
!
                tmp = ieor( ieor( ieor(s1,s2), s3), s4)
!
!               GENERATE A RANDOM UNSIGNED 31-BIT INTEGER.
!
                harvest(i) = ibclr( tmp, topbit)
!
            end do
!
        case (4,9)
!
!           MERSENNE TWISTER MT19937 RNG.
!
            do i = 1_i4b, nharvest
!
                if ( mti>=mt_n ) then
!
!                   GENERATE mt_n WORDS AT ONE TIME.
!
                    do j = 0_ki_sel, mt_n - mt_m - 1_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+mt_m),ishft(tmp,-1)),   &
                                      mt_mag01(iand(tmp,1_ki_sel))      )
!
                    end do
!
                    do j = mt_n - mt_m, mt_n - 2_ki_sel
!
                        tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                     iand(mt(j+1_ki_sel),mt_lower_mask) )
!
                        mt(j) = ieor( ieor(mt(j+(mt_m-mt_n)), ishft(tmp,-1)),  &
                                      mt_mag01(iand(tmp,1_ki_sel))             )
!
                    end do
!
                    tmp   = ior( iand(mt(j),mt_upper_mask),         &
                                 iand(mt(0_ki_sel),mt_lower_mask)   )
!
                    mt(j) = ieor( ieor(mt(mt_m-1_ki_sel),ishft(tmp,-1)),  &
                                  mt_mag01(iand(tmp,1_ki_sel))            )
!
                    mti = 0_ki_sel
!
                endif
!
                tmp = mt(mti)
                mti = mti + 1_ki_sel
!
!               TEMPERING PHASE.
!
                tmp = ieor( tmp, ishft(tmp,mt_shift0) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t1_shift),mt_t1_mask) )
                tmp = ieor( tmp, iand(ishft(tmp,mt_t2_shift),mt_t2_mask) )
                tmp = ieor( tmp, ishft(tmp,mt_shift1) )
!
!               GENERATE A RANDOM UNSIGNED 31-BIT INTEGER.
!
                harvest(i) = ibclr( tmp, topbit)
!
            end do
!
        case (5,10)
!
!           MERSENNE TWISTER MT19937-II RNG.
!
            do i = 1_i4b, nharvest
!
                select case (memti)
!
                    case (0_ki_sel:mt_n-mt_m-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+mt_m),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))          )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-mt_m:memt_lag1over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti+memt_lag1),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag1over:memt_lag2over-1_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti+memt_lag2),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (memt_lag2over:mt_n-2_ki_sel)
!
                        tmp = ior( iand(memt(memti),mt_upper_mask),            &
                                   iand(memt(memti+1_ki_sel),mt_lower_mask)    )
!
                        memt(memti) = ieor( ieor(memt(memti+(mt_m-mt_n)),ishft(tmp,-1)), &
                                            mt_mag01(iand(tmp,1_ki_sel))                 )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = memti + 1_ki_sel
!
                    case (mt_n-1_ki_sel)
!
                        tmp = ior( iand(memt(mt_n-1_ki_sel),mt_upper_mask),  &
                                   iand(memt(0_ki_sel),mt_lower_mask)        )
!
                        memt(mt_n-1_ki_sel) = ieor( ieor(memt(mt_m-1_ki_sel),ishft(tmp,-1)),   &
                                                    mt_mag01(iand(tmp,1_ki_sel))               )
!
!                       TEMPERING PHASE.
!
                        tmp = ieor( memt(memti), iand(memt(memti-memt_lag1over),memt_mask1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift1) )
                        tmp = ieor( tmp, ishft(tmp,memt_shift2) )
                        tmp = ieor( tmp, iand(memt(memti-memt_lag2over),memt_mask2) )
!
                        memti = 0_ki_sel
!
                end select
!
!               GENERATE A RANDOM UNSIGNED 31-BIT INTEGER.
!
#ifdef _RANDOM_NOINT32
                harvest(i) = ibits( tmp, 0, topbit )
#else
                harvest(i) = ibclr( tmp, topbit)
#endif
!
            end do
!
    end select
!$OMP END CRITICAL (ran_num_)
!
!
! END OF SUBROUTINE random_i31_1
! ______________________________
!
    end subroutine random_i31_1
!
! =========================================================================================
!
    subroutine random_i31_2( harvest )
!
! purpose
! _______
!
!   This subroutine returns a matrix of random integers in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   unsigned 31-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:)
!            A matrix of random integers in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 2 )
!
    do i = 1_i4b, nharvest
!
        call random_i31_1( harvest(:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i31_2
! ______________________________
!
    end subroutine random_i31_2
!
! =========================================================================================
!
    subroutine random_i31_3( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   unsigned 31-bit integers.
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:)
!            An array of random integers in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 3 )
!
    do i = 1_i4b, nharvest
!
        call random_i31_2( harvest(:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i31_3
! ______________________________
!
    end subroutine random_i31_3
!
! =========================================================================================
!
    subroutine random_i31_4( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   unsigned 31-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:)
!            An array of random integers in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 4 )
!
    do i = 1_i4b, nharvest
!
        call random_i31_3( harvest(:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i31_4
! ______________________________
!
    end subroutine random_i31_4
!
! =========================================================================================
!
    subroutine random_i31_5( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   unsigned 31-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:,:)
!            An array of random integers in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:,:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 5 )
!
    do i = 1_i4b, nharvest
!
        call random_i31_4( harvest(:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i31_5
! ______________________________
!
    end subroutine random_i31_5
!
! =========================================================================================
!
    subroutine random_i31_6( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   unsigned 31-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:,:,:)
!            An array of random integers in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(out), dimension(:,:,:,:,:,:)    :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 6 )
!
    do i = 1_i4b, nharvest
!
        call random_i31_5( harvest(:,:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i31_6
! ______________________________
!
    end subroutine random_i31_6
!
! =========================================================================================
!
    subroutine random_i31_7( harvest )
!
! purpose
! _______
!
!   This subroutine returns an array of random integers in the interval (0,2147483647)
!   inclusive of the two endpoints. The returned integers are equivalent to
!   unsigned 31-bit integers.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) integer(i4b), dimension(:,:,:,:,:,:,:)
!            An array of random integers in the interval (0,2147483647).
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b),  intent(out), dimension(:,:,:,:,:,:,:)  :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: i, nharvest
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nharvest = size( harvest, 7 )
!
    do i = 1_i4b, nharvest
!
        call random_i31_6( harvest(:,:,:,:,:,:,i) )
!
    end do
!
!
! END OF SUBROUTINE random_i31_7
! ______________________________
!
    end subroutine random_i31_7
!
!
! =========================================================================================
!             SUBROUTINES FOR SEEDING AND MANIPULATING THE STATE OF THE RNGS
! =========================================================================================
!
!
    subroutine init_mt19937_r0( seed )
!
! purpose
! _______
!
!   User interface subroutine for initializing the state of the MT19937 Random Number Generator (RNG)
!   with a scalar seed, directly, without using the subroutine RANDOM_SEED\ _ and its interface.
!
!
! Arguments
! _________
!
!   SEED (INPUT) integer(i4b)
!        On entry, a scalar integer that will be used to initialize the MT19937 RNG.
!
!
! Further Details
! _______________
!
!   Only the first 32 bits of the scalar SEED will be used.
!
!   For more informations on the MT19937 RNG, see:
!
!   (1)  Matsumoto, M., and Nishimura, T., 1998:
!            Mersenne Twister: A 623-dimensionally
!            equidistributed uniform pseudorandom number generator.
!            ACM Trans. on Modeling and Computer Simulation, Vol. 8, No. 1, January pp.3-30
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in)  :: seed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: tmp
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   COPY THE FIRST 32 BITS OF seed IN tmp WITHOUT INTEGER OVERFLOW.
!
    tmp = ibits( seed, 0, topbit )
!
    if ( btest( seed,  topbit ) ) then
        tmp = ibset( tmp, topbit )
    end if
!
    call init_mt_by_scalar( tmp, mt )
!
!   ENSURE THAT INDEX mti IS INITIALIZED CORRECTLY.
!
    mti = mt_n
!
    mt_initialized = true
!
!
! END OF SUBROUTINE init_mt19937_r0
! _________________________________
!
    end subroutine init_mt19937_r0
!
! =========================================================================================
!
    subroutine init_mt19937_r1( seed )
!
! purpose
! _______
!
!   User interface subroutine for initializing the state of the MT19937 Random Number Generator (RNG)
!   with an array of seeds, directly, without using the subroutine RANDOM_SEED\ _ and its interface.
!
!
! Arguments
! _________
!
!   SEED (INPUT) integer(i4b), dimension(:)
!        On entry, a vector of integers that will be used to initialize the MT19937 RNG.
!        If size(SEED) is zero, a default scalar seed will be used instead.
!
!
! Further Details
! _______________
!
!   Only the first 32 bits of each element of the array SEED will be used.
!
!   For more informations on the MT19937 RNG, see:
!
!   (1)  Matsumoto, M., and Nishimura, T., 1998:
!            Mersenne Twister: A 623-dimensionally
!            equidistributed uniform pseudorandom number generator.
!            ACM Trans. on Modeling and Computer Simulation, Vol. 8, No. 1, January pp.3-30
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), dimension(:), intent(in)  :: seed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nseed
    integer(ki_sel), dimension(size(seed)) :: tmp
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nseed = size(seed)
!
    if ( nseed>=1_i4b ) then
!
!       COPY THE FIRST 32 BITS OF EACH ELEMENT OF seed IN THE CORRESPONDING ELEMENT
!       OF tmp WITHOUT INTEGER OVERFLOW.
!
        tmp(:nseed) = ibits( seed(:nseed), 0, topbit )
!
        where(  btest( seed(:nseed),  topbit ) ) tmp(:nseed) = ibset( tmp(:nseed), topbit )
!
        call init_mt_by_array( tmp(:nseed), mt )
!
    else
!
!       IF THE ARRAY seed IS OF ZERO SIZE, USE
!       A DEFAULT SEED.
!
        call init_mt_by_scalar( mt_default_seed, mt )
!
    end if
!
!   ENSURE THAT INDEX mti IS INITIALIZED CORRECTLY.
!
    mti = mt_n
!
    mt_initialized = true
!
!
! END OF SUBROUTINE init_mt19937_r1
! _________________________________
!
    end subroutine init_mt19937_r1
!
! =========================================================================================
!
    subroutine init_memt19937_r0( seed )
!
! purpose
! _______
!
!   User interface subroutine for initializing the state of the MEMT19937-II Random Number Generator (RNG)
!   with a seed, directly, without using the subroutine RANDOM_SEED\ _ and its interface.
!
!
! Arguments
! _________
!
!   SEED (INPUT) integer(i4b)
!        On entry, a scalar integer that will be used to initialize the MEMT19937-II RNG.
!
!
! Further Details
! _______________
!
!   Only the first 32 bits of the scalar SEED will be used.
!
!   For more informations on the MEMT19937-II RNG, see:
!
!   (1)  Harase, S., 2014:
!            On the F2-linear relations of Mersenne Twister pseudorandom number generators.
!            Mathematics and Computers in Simulation, Volume 100, Pages 103-113.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in)  :: seed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: tmp
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   COPY THE FIRST 32 BITS OF seed IN tmp WITHOUT INTEGER OVERFLOW.
!
    tmp = ibits( seed, 0, topbit )
!
    if ( btest( seed,  topbit ) ) then
        tmp = ibset( tmp, topbit )
    end if
!
    call init_mt_by_scalar( tmp, memt )
!
!   ENSURE THAT INDEX memti IS INITIALIZED CORRECTLY.
!
    memti = 0_ki_sel
!
    memt_initialized = true
!
!
! END OF SUBROUTINE init_memt19937_r0
! ___________________________________
!
    end subroutine init_memt19937_r0
!
! =========================================================================================
!
    subroutine init_memt19937_r1( seed )
!
! purpose
! _______
!
!   User interface subroutine for initializing the state of the MEMT19937-II Random Number Generator (RNG)
!   with an array of seeds, directly, without using the subroutine RANDOM_SEED\ _ and its interface.
!
!
! Arguments
! _________
!
!   SEED (INPUT) integer(i4b), dimension(:)
!        On entry, a vector of integers that will be used to initialize the MT19937 RNG.
!        If size(SEED) is zero, a default scalar seed will be used instead.
!
!
! Further Details
! _______________
!
!   Only the first 32 bits of each element of the array SEED will be used.
!
!   For more informations on the MEMT19937-II RNG, see:
!
!   (1)  Harase, S., 2014:
!            On the F2-linear relations of Mersenne Twister pseudorandom number generators.
!            Mathematics and Computers in Simulation, Volume 100, Pages 103-113.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), dimension(:), intent(in)  :: seed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nseed
    integer(ki_sel), dimension(size(seed)) :: tmp
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nseed = size(seed)
!
    if ( nseed>=1_i4b ) then
!
!       COPY THE FIRST 32 BITS OF EACH ELEMENT OF seed IN THE CORRESPONDING ELEMENT
!       OF tmp WITHOUT INTEGER OVERFLOW.
!
        tmp(:nseed) = ibits( seed(:nseed), 0, topbit )
!
        where(  btest( seed(:nseed),  topbit ) ) tmp(:nseed) = ibset( tmp(:nseed), topbit )
!
        call init_mt_by_array( tmp(:nseed), memt )
!
    else
!
!       IF THE ARRAY seed IS OF ZERO SIZE, USE
!       A DEFAULT SEED.
!
        call init_mt_by_scalar( mt_default_seed, memt )
!
    end if
!
!   ENSURE THAT INDEX memti IS INITIALIZED CORRECTLY.
!
    memti = 0_ki_sel
!
    memt_initialized = true
!
!
! END OF SUBROUTINE init_memt19937_r1
! ___________________________________
!
    end subroutine init_memt19937_r1
!
! =========================================================================================
!
    subroutine random_seed_( alg, size, put, get )
!
! purpose
! _______
!
!   User interface for seeding the random number routines in module RANDOM.
!
!   Syntax is like RANDOM_SEED intrinsic subroutine and a call to RANDOM_SEED\ _()
!   without arguments initiates a non-repeatable reset of the seeds used by the
!   random number subroutines in module RANDOM.
!
!   As for RANDOM_SEED intrinsic subroutine, no more than one argument may be specified
!   in a call to RANDOM_SEED\ _ .
!
!
! Arguments
! _________
!
!   ALG  (INPUT, OPTIONAL) integer
!        On entry, a scalar default integer to select the random number generator used in
!        subsequent calls to subroutines RANDOM_NUMBER\ _ , RANDOM_INTEGER32\ _ , RANDOM_INTEGER31\ _
!        and functions RAND_NUMBER, RAND_INTEGER32 and RAND_INTEGER31. The possible values are:
!
!        - ALG=1  : selects the Marsaglia's KISS random number generator;
!        - ALG=2  : selects the fast Marsaglia's KISS random number generator;
!        - ALG=3  : selects the L'Ecuyer's LFSR113 random number generator;
!        - ALG=4  : selects the Mersenne Twister random number generator;
!        - ALG=5  : selects the maximally equidistributed Mersenne Twister random number generator;
!        - ALG=6  : selects the extended precision of the Marsaglia's KISS random number generator;
!        - ALG=7  : selects the extended precision of the fast Marsaglia's KISS random number generator;
!        - ALG=8  : selects the extended precision of the L'Ecuyer's LFSR113 random number generator.
!        - ALG=9  : selects the extended precision of Mersenne Twister random number generator;
!        - ALG=10 : selects the extended precision of maximally equidistributed Mersenne Twister
!          random number generator;
!
!        For other values the random number generator is not changed.
!        The default value is the L'Ecuyer's LFSR113 random number generator (e.g. ALG=3).
!
!   SIZE (OUTPUT, OPTIONAL) integer
!        On exit, the size of the seed array used by the random number generators.
!
!   PUT  (INPUT, OPTIONAL) integer, dimension(:)
!        On entry, a scalar default integer vector of size at least equal to the size of
!        the seed array (as returned by a call to RANDOM_SEED\ _ with argument SIZE) that
!        will be used to reset the seed for subsequent calls to subroutine RANDOM_NUMBER\ _.
!
!   GET  (OUTPUT, OPTIONAL) integer, dimension(:)
!        On exit, a scalar default integer vector, which is the current value of the seed array.
!        Argument GET must be of size at least equal to the size of the seed array (as
!        returned by a call to RANDOM_SEED\ _ with argument SIZE).
!
!
! Further Details
! _______________
!
!   This subroutine is not thread-safe and must not be called in parallel when OPENMP
!   is used. On the other hand, the associated routines RAND_NUMBER, RANDOM_NUMBER\ _,
!   RAND_INTEGER32, RANDOM_INTEGER32\ _, RAND_INTEGER31 and RANDOM_INTEGER31\ _ are
!   thread-safe if used with OpenMP directives and their states will be consistent
!   while called from multiple OpenMP threads.
!
!   The Marsaglia's KISS (Keep It Simple Stupid) random number generator combines:
!
!     1) The congruential generator x(n) = 69069 \cdot x(n-1) + 1327217885 with
!        a period of 2^32;
!     2) A 3-shift shift-register generator with a period of 2^32 - 1;
!     3) Two 16-bit multiply-with-carry generators with a period of
!        597273182964842497 > 2^59.
!
!   The overall period of this KISS random number generator exceeds 2^123. More details
!   on this Marsaglia's KISS random number generator are available in the references (3) and (4).
!   This generator is also the one used by the intrinsic subroutine RANDOM_NUMBER as implemented
!   in the GNU gfortran compiler.
!
!   The "fast" version of the Marsaglia's KISS random number generator uses only
!   add, shift, exclusive-or and 'and' operations to produces exactly the same 32-bit integer
!   output, which C views as unsigned and Fortran views as signed integers. This version avoids
!   multiplication and is probably faster. More details are available in the reference (5).
!
!   The LFSR113 random number generator is described in the reference (2). This
!   random number generator has a period length of about 2^113.
!
!   The MT19937 Mersenne Twister random number generator is described in the reference (7).
!   This random number generator has a period length of about 2^19937-1, and 623-dimensional
!   equidistribution property is assured.
!   
!   The MEMT19937-II Mersenne Twister random number generator is described in the reference (8).
!   This random number generator has also a period length of about 2^19937-1, and a new set of
!   parameters is introduced in the tempering phase of MT19937, which gives a maximally
!   equidistributed Mersenne Twister random number generator.
!
!   Note that the size of the seed array varies according to the selected random
!   number generator.
!
!   For all the random number generators described above, extended precision versions are
!   also available to generate full precision random real numbers of kind STND (up to 63-bit precision),
!   using the method described in the reference (6).
!
!   The FORTRAN versions of these random number generators as implemented here require that
!   32-bits integer type is available on your computer and that 32-bits integers are represented
!   in base 2 with two's complement notation.
!
!   However, the LFSR113, MT19937 and MEMT19937-II Mersenne Twister random number generators will
!   also work if only 64-bits integer type is available on your system, but in that case you must
!   specify the CPP macro _RANDOM_NOINT32 at compilation. Note, however that the other random number
!   generators will not work properly with 64-bits integer type so they cannot be used on such system.
!
!   The KISS random number generators also assumed that integer overflows do not cause crashes.
!   These assumptions are checked before using these random number generators.
!
!   On the other hand, the LFSR113, MT19937 and MEMT19937-II random number generators do not
!   use integer arithmetic and are free of such assumptions.
!
!   This subroutine is adapted from:
!
!   (1)  Hennecke, M., 1995:
!            A Fortran90 interface to random number generation.
!            Computer Physics Communications, Volume 90, Number 1, 117-120
!
!   (2)  L'Ecuyer, P., 1999:
!            Tables of Maximally-Equidistributed Combined LFSR Generators.
!            Mathematics of Computation, 68, 225, 261-269.
!
!   (3)  Marsaglia, G., 1999:
!            Random number generators for Fortran.
!            Posted to the computer-programming-forum. See:
!            http://computer-programming-forum.com/49-fortran/b89977aa62f72ee8.htm
!
!   (4)  Marsaglia, G., 2005:
!            Double precision RNGs.
!            Posted to the electronic billboard sci.math.num-analysis. See:
!            http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
!   (5)  Marsaglia, G., 2007:
!            Fortran and C: United with a KISS.
!            Posted to the Google comp.lang.forum . See:
!            http://groups.google.co.uk/group/comp.lang.fortran/msg/6edb8ad6ec5421a5
!
!   (6)  Doornik, J.A, 2007:
!            Conversion of high-period random number to floating point.
!            ACM Transactions on Modeling and Computer Simulation, 
!            Volume 17, Issue 1, Article No. 3.
!
!   (7)  Matsumoto, M., and Nishimura, T., 1998:
!            Mersenne Twister: A 623-dimensionally
!            equidistributed uniform pseudorandom number generator.
!            ACM Trans. on Modeling and Computer Simulation, Vol. 8, No. 1, January pp.3-30
!
!   (8)  Harase, S., 2014:
!            On the F2-linear relations of Mersenne Twister pseudorandom number generators.
!            Mathematics and Computers in Simulation, Volume 100, Pages 103-113.
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _RANDOM_NOINT32
    use Char_Constants,   only : random_error3, random_error4, random_error7,   &
                                 random_error9, random_error10
#else
    use Char_Constants,   only : random_error1, random_error2, random_error3,   &
                                 random_error4, random_error7, random_error8
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer, optional, intent(in)                :: alg
    integer, optional, intent(out)               :: size
    integer, optional, intent(in),  dimension(:) :: put
    integer, optional, intent(out), dimension(:) :: get
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: tmp
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='random_seed_'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( ( present(size) .and. present(put) ) .or.                    &
         ( present(size) .and. present(get) ) .or.                    &
         ( present(put)  .and. present(get) ) .or.                    &
         ( present(alg)  .and. present(get) ) .or.                    &
         ( present(alg)  .and. present(put) ) .or.                    &
         ( present(alg ) .and. present(size) )      ) then
!
         call merror( name_proc//random_error3 )
!
    endif
!
    if ( first ) then
#ifdef _RANDOM_NOINT32
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method/=3 ) then
            call merror( name_proc//random_error9 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ibclr(-1_ki_sel, bit_size(tmp)-1_ki_sel)
!
        if ( radix(tmp)/=2 .or. tmp/=huge(tmp) ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
#else
!
!       TEST IF 32-BIT INTEGERS ARE AVAILABLE.
!
        if ( bit_size( tmp )/=fullbitsize ) then
            call merror( name_proc//random_error1 )
        endif
!
!       TEST IF THE DEFAULT RNG IS SET CORRECTLY.
!
        if ( method<1 .or. method>3 ) then
            call merror( name_proc//random_error2 )
        endif
!
!       TEST IF THE BASE OF THE INTEGER SYSTEM IS 2 AND
!       IF INTEGER REPRESENTATION IS TWO'S COMPLEMENT NOTATION.
!
        tmp = ishftc( -8_ki_sel, -3)
!
        if ( radix(tmp)/=2 .or. tmp/=536870911_ki_sel ) then
!
            call merror( name_proc//random_error7 )
!
        endif
!
        if ( method/=3 ) then
!
!           TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS.
!
            tmp   = huge( tmp )
!
            if ( integer_not_safe( tmp ) ) then
!
                call merror( name_proc//random_error8 )
!
            endif
!
            test_kiss = false
!
        endif
#endif
!
!       SAVE THE NUMBER OF SEEDS AND SIZE OF THE SEED ARRAY OF
!       THE DEFAULT RNG FOR LATER USE.
!
        if ( method==2 ) then
            num_seeds = 5
        else
            num_seeds = 4
        end if
!
        seed_size = num_seeds*ratio_i
!
        first = false
!
    endif
!
    if ( present(alg ) ) then
!
        if ( alg>0 .and. alg<=10 ) then
!
            method = alg
!
!           REVERT TO THE STANDARD VERSION OF THE RNGS IF THE NUMBER OF
!           BITS IN THE MANTISSA IS LESS THAN 32.
!
            if ( method>5 .and. num_bits_stnd<=int(fullbitsize) ) method = method - 5
!
            if ( method>5 .and. m_ran_invmbits==-one ) then
!
!               COMPUTE THE CONSTANTS NEEDED TO CONVERT RANDOM INTEGERS
!               INTO RANDOM FLOATING POINT NUMBERS IN THE EXTENDED
!               VERSIONS OF THE RNGS.
!
                m_ran_invmbits = half**(num_bits_stnd)
!
            end if
!
            if ( test_kiss ) then
!
!               TEST IF INTEGER OVERFLOWS ARE SAFE FOR THE MARSAGLIA'S KISS RNGS
!               IF THEY ARE SELECTED.
!
                select case (method)
!
                    case(1,2,6,7)
!
#ifdef _RANDOM_NOINT32
                        call merror( name_proc//random_error10 )
#else
                        tmp  = huge( tmp )
!
                        if ( integer_not_safe( tmp ) ) then
!
                            call merror( name_proc//random_error8 )
!
                        end if
!
                        test_kiss = false
#endif
!
                end select
!
            endif
!
!           SAVE THE NUMBER OF SEEDS AND SIZE OF THE SEED ARRAY OF
!           THE CURRENT RNG FOR LATER USE.
!
            select case (method)
!
                case(1,3,6,8)
!
                    num_seeds = 4
                    seed_size = num_seeds*ratio_i
!
                case(2,7)
!
                    num_seeds = 5
                    seed_size = num_seeds*ratio_i
!
                case(4,5,9,10)
!
                    num_seeds = int( mt_n )
                    seed_size = (num_seeds+1)*ratio_i
!
!                   INITIALIZE THE MERSENNE TWISTER RNGS.
!
                    if ( method==4 .or. method==9 ) then
!
                        if ( .not.mt_initialized ) then
!
                            call init_mt_by_array( mt_default_seed_array, mt )
!                            call init_mt_by_scalar( mt_default_seed, mt )
!
!                           ENSURE THAT INDEX mti IS INITIALIZED CORRECTLY.
!
                            mti = mt_n
!
                            mt_initialized = true
!
                        end if
!
                    else
!
                        if ( .not.memt_initialized ) then
!
                            call init_mt_by_array( mt_default_seed_array, memt )
!                            call init_mt_by_scalar( mt_default_seed, memt )
!
!                           ENSURE THAT INDEX memti IS INITIALIZED CORRECTLY.
!
                            memti = 0_ki_sel
!
                            memt_initialized = true
!
                        end if
!
                    end if
!
            end select
!
        end if
!
    else if ( present(size) ) then
!
        size = seed_size
!
    else if ( present(get ) ) then
!
        call seed_get( get )
!
    else if ( present(put ) ) then
!
        call seed_put( put )
!
    else
!
        call seed_set ( )
!
    endif
!
!----------------------------------------------------------------------
                         contains
!----------------------------------------------------------------------
!
        subroutine seed_set ( )
!
! purpose
! _______
!
!   Internal routine used for initializing the state of the random number generators.
!
! ___________________________________________________________________________________
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
        integer               :: count, rate, i, unit, istat
        integer, dimension(8) :: dt
!
        integer(ki_sel)                       :: jsr, jsr_old
#ifdef _USE_PGI
        integer(ki_sel), dimension(mt_n) :: iseed
#else
        integer(ki_sel), dimension(num_seeds) :: iseed
#endif
!        integer(ki_sel), dimension(5) :: iseed
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!       FIRST TRY IF THE OS PROVIDES A RNG.
!
        open( unit=unit, file=urandom_file, access="stream",                &
              form="unformatted", action="read", status="old", iostat=istat )
!
        if ( istat==0 ) then
!
            read(unit) iseed(:num_seeds)
            close(unit)
!
        else
!
!           OTHERWISE USE THE CURRENT TIME AND PID. THE PID IS USEFUL
!           IN CASE ONE LAUNCHES MULTIPLE INSTANCES OF THE SAME PROGRAM
!           IN PARALLEL.
!
            call system_clock( count=count, count_rate=rate )
!
            if ( rate==0 ) then
!
                call merror( name_proc//random_error4 )
!
            else
!
                if ( count==0 ) then
!
                    call date_and_time( values=dt )
!
                    count = ( dt(1) - 2010 )*365*24*60*60*1000  &
                            + dt(2)         *31*24*60*60*1000   &
                            + dt(3)         *24*60*60*1000      &
                            + dt(5)         *60*60*1000         &
                            + dt(6)         *60*1000            &
                            + dt(7)         *1000               &
                            + dt(8)
                end if
!
#ifndef _RANDOM_NOUNIX
                i = getpid()
                count = ieor( count, i )
#endif
!
                jsr = int( count, ki_sel )
!
!                jsr = ieor( 777755555_ki_sel, jsr )
!
!               NOW USE A FAST BASIC 32-BIT XORSHIFT RNG TO GENERATE
!               THE SEEDS FOR THE SELECTED RNG.
!
#ifdef _RANDOM_NOINT32
                jsr = ibits( jsr, 0, fbs )
!
                do i = 1, num_seeds
!
                    jsr_old  = jsr
!
                    jsr      = ibits( ieor( jsr, ishft(jsr, 5) ), 0, fbs )
                    jsr      = ieor( jsr, ishft(jsr, -7) )
                    jsr      = ibits( ieor( jsr, ishft(jsr, 22) ), 0, fbs )
                    iseed(i) = uiadd( jsr, jsr_old )
!
                end do
#else
                do i = 1, num_seeds
!
                    jsr_old  = jsr
!
                    jsr      = ieor( jsr, ishft(jsr, 5) )
                    jsr      = ieor( jsr, ishft(jsr, -7) )
                    jsr      = ieor( jsr, ishft(jsr, 22) )
!                    jsr      = m( m( m(jsr, 5), -7), 22)
                    iseed(i) = uiadd( jsr, jsr_old )
!                    iseed(i) = jsr + jsr_old
!
                end do
#endif
!
            end if
!
        end if
!
        select case (method)
!
            case(1,6)
!
!               SET THE SEED FOR THE Marsaglia's KISS RNG.
!
                x = iseed(1)
                y = iseed(2)
                z = iseed(3)
                w = iseed(4)
!
!               y MUST NOT BE ZERO.
!
                y = ieor( 777755555_ki_sel, y )
!
            case(2,7)
!
!               SET THE SEEDS FOR THE FAST VERSION OF Marsaglia's KISS RNG.
!
                x2 = iseed(1)
                y2 = iseed(2)
                z2 = ibclr( iseed(3), topbit )
                w2 = ibclr( iseed(4), topbit )
                c  = ishft( iseed(5), -topbit)
!
!               THE FOLLOWING INSTRUCTIONS ENSURE THAT THE RNG WILL HAVE
!               A FULL PERIOD.
!
                if ( y2==0_ki_sel) y2 = 777755555_ki_sel
                if ( mod(z2,7559_ki_sel)==0_ki_sel ) z2 = z2 + 1_ki_sel
                if ( mod(w2,7559_ki_sel)==0_ki_sel ) w2 = w2 + 7558_ki_sel
!
            case(3,8)
!
!               SET THE SEEDS FOR THE L'Ecuyer's LFSR113 RNG.
!
#ifdef _RANDOM_NOINT32
                s1 = ibits( iseed(1), 0, fbs )
                s2 = ibits( iseed(2), 0, fbs )
                s3 = ibits( iseed(3), 0, fbs )
                s4 = ibits( iseed(4), 0, fbs )
!
!               THE INITIAL SEEDS s1, s2, s3 ANS s4 MUST BE LARGER THAN 1, 7, 15 AND 127,
!               RESPECTIVELY, WHEN CONSIDERED AS UNSIGNED 32-BIT INTEGERS.
!               THE FOLLOWING INSTRUCTIONS CHECK THESE CONDITIONS ASSUMING THAT
!               THE BASE OF THE INTEGER SYSTEM IS 2 AND INTEGER REPRESENTATION IS
!               TWO'S COMPLEMENT.
!
                if ( iand(s1,  -2_ki_sel)==0_ki_sel ) s1 = uisub( s1, 1023_ki_sel )
                if ( iand(s2,  -8_ki_sel)==0_ki_sel ) s2 = uisub( s2, 1023_ki_sel )
                if ( iand(s3, -16_ki_sel)==0_ki_sel ) s3 = uisub( s3, 1023_ki_sel )
                if ( iand(s4,-128_ki_sel)==0_ki_sel ) s4 = uisub( s4, 1023_ki_sel )
#else
                s1 = iseed(1)
                s2 = iseed(2)
                s3 = iseed(3)
                s4 = iseed(4)
!
!               THE INITIAL SEEDS s1, s2, s3 ANS s4 MUST BE LARGER THAN 1, 7, 15 AND 127,
!               RESPECTIVELY, WHEN CONSIDERED AS UNSIGNED 32-BIT INTEGERS.
!               THE FOLLOWING INSTRUCTIONS CHECK THESE CONDITIONS ASSUMING THAT
!               THE BASE OF THE INTEGER SYSTEM IS 2 AND INTEGER REPRESENTATION IS
!               TWO'S COMPLEMENT.
!
                if ( iand(s1,  -2_ki_sel)==0_ki_sel ) s1 = s1 - 1023_ki_sel
                if ( iand(s2,  -8_ki_sel)==0_ki_sel ) s2 = s2 - 1023_ki_sel
                if ( iand(s3, -16_ki_sel)==0_ki_sel ) s3 = s3 - 1023_ki_sel
                if ( iand(s4,-128_ki_sel)==0_ki_sel ) s4 = s4 - 1023_ki_sel
#endif
!
            case(4,9)
!
!               SET THE SEEDS FOR THE MT19937 Mersenne Twister RNG.
!
                call init_mt_by_array( iseed, mt )
!
!               ENSURE THAT INDEX mti IS REINITIALIZED CORRECTLY.
!
                mti = mt_n
!
            case(5,10)
!
!               SET THE SEEDS FOR THE MEMT19937-II Mersenne Twister RNG.
!
                call init_mt_by_array( iseed, memt )
!
!               ENSURE THAT INDEX memti IS REINITIALIZED CORRECTLY.
!
                memti = 0_ki_sel
!
        end select
!
!
!    END OF subroutine seed_set
! _____________________________
!
        end subroutine seed_set
!
!----------------------------------------------------------------------
!
!
! END OF SUBROUTINE random_seed_
! ______________________________
!
    end subroutine random_seed_
!
!
! =========================================================================================
!                                   INTERNAL PRIVATE SUBROUTINES
! =========================================================================================
!
!
!    function m(k, n)
!
! purpose
! _______
!
!   Internal routine used by the Marsaglia Shift Register RNGs.
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
!    integer(ki_sel) :: m, k
!
!    integer         :: n
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!    m = ieor( k, ishft(k, n) )
!
!
! END OF function m
! _________________
!
!    end function m
!
! =========================================================================================
!
    function integer_not_safe( ival )
!
! purpose
! _______
!
!   Internal routine used for testing integer arithmetic. This function must be called
!   with IVAL set to huge(1_ki_sel) to perform properly.
!   This function returns .true. if integer overflows don't work as expected
!   and are unsafe. The function returns .false. otherwise.
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(ki_sel), intent(in) :: ival
!
    logical :: integer_not_safe
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: tmp
!
    logical :: test1, test2, test3
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    tmp = ival + 1_ki_sel
!
    test1 = tmp /= -ival - 1_ki_sel
    test2 = ival - tmp + 1_ki_sel /= 0_ki_sel
    test3 = ival + tmp + 1_ki_sel /= 0_ki_sel
!
    integer_not_safe = test1 .or. test2 .or. test3
!
!
! END OF function integer_not_safe
! ________________________________
!
    end function integer_not_safe
!
! =========================================================================================
!
    subroutine init_mt_by_array( seed, state )
!
! purpose
! _______
!
!   Internal routine used for initializing the states of the MT19937 or MEMT19937-II 
!   RNGs with an array of seeds.
!
!   This is done without integer overflows.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(ki_sel), intent(in)    :: seed(0_ki_sel:)
    integer(ki_sel), intent(inout) :: state(0_ki_sel:mt_n-1_ki_sel)
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: nseed, i, im1, j, k, tmp
!
!
! PARAMETERS
! __________
!
    integer(ki_sel), parameter :: seed_d =    19650218_ki_sel !z'12BD6AA'
    integer(ki_sel), parameter :: mult_a =     1664525_ki_sel !z'19660D'
    integer(ki_sel), parameter :: mult_b =  1566083941_ki_sel !z'5D588B65'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nseed = size( seed )
!
    call init_mt_by_scalar( seed_d, state )
!
    i = 1_ki_sel
    j = 0_ki_sel
!
    do k = max(mt_n,nseed), 1_ki_sel, -1_ki_sel
!
        im1 = i - 1_ki_sel
!
        tmp      = ieor( state(im1), ishft(state(im1),-30) )
        tmp      = uimlt( tmp, mult_a )
        state(i) = ieor( state(i), tmp )
        state(i) = uiadd( state(i), uiadd( seed(j), j ) )
         
!        state(i) = ieor( state(i),ieor(state(im1),ishft(state(im1),-30))*mult_a ) + seed(j) + j
!
#ifdef _RANDOM_NOINT32
        state(i) = ibits( state(i), 0, fbs )
#endif
!
        i = i + 1_ki_sel
        j = j + 1_ki_sel
!
        if ( i>=mt_n ) then
          state(0_ki_sel) = state(mt_n-1_ki_sel)
          i = 1_ki_sel
        endif
!
        if ( j>=nseed ) then
          j = 0_ki_sel
        endif
!
      end do
!
      do k = mt_n-1_ki_sel, 1_ki_sel, -1_ki_sel
!
        im1 = i - 1_ki_sel
!
        tmp = ieor( state(im1), ishft(state(im1),-30) ) 
        tmp = uimlt( tmp, mult_b )
        state(i) = ieor( state(i), tmp )
        state(i) = uisub( state(i), i )
         
!        state(i) = ieor( state(i), ieor(state(im1),ishft(state(im1),-30))*mult_b ) - i
!
#ifdef _RANDOM_NOINT32
        state(i) = ibits( state(i), 0, fbs )
#endif
!
        i = i + 1_ki_sel
!
        if ( i>=mt_n ) then
          state(0_ki_sel) = state(mt_n-1_ki_sel)
          i = 1_ki_sel
        endif
!
      end do
!
      state(0_ki_sel) = mt_upper_mask
!
!
! END OF SUBROUTINE init_mt_by_array
! __________________________________
!
    end subroutine init_mt_by_array
!
! =========================================================================================
!
    subroutine init_mt_by_scalar( seed, state )
!
! purpose
! _______
!
!   Internal routine used for initializing the states of the MT19937 or MEMT19937-II 
!   RNGs with a seed.
!
!   This is done without integer overflows.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(ki_sel), intent(in)    :: seed
    integer(ki_sel), intent(inout) :: state(0_ki_sel:mt_n-1_ki_sel)
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: i, tmp
!
!
! PARAMETERS
! __________
!
    integer(ki_sel), parameter :: mult_a = 1812433253_ki_sel  !z'6C078965'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    state(0_ki_sel) = seed
#ifdef _RANDOM_NOINT32
    state(0_ki_sel) = ibits( seed, 0, fbs )
#endif
!
    do  i = 1_ki_sel, mt_n-1_ki_sel
!
        tmp      = ieor( state(i-1_ki_sel), ishft(state(i-1_ki_sel),-30) )
        tmp      = uimlt( tmp, mult_a )
        state(i) = uiadd( tmp, i )
!
!        state(i) = mult_a*ieor(state(i-1_ki_sel),ishft(state(i-1_ki_sel),-30)) + i
!
#ifdef _RANDOM_NOINT32
        state(i) = ibits( state(i), 0, fbs )
#endif
!
     end do
!
!
! END OF SUBROUTINE init_mt_by_scalar
! ___________________________________
!
    end subroutine init_mt_by_scalar
!
! =========================================================================================
!
    elemental function uiadd( a, b ) result( c )
!
! purpose
! _______
!
!   Internal routine used to simulate addition for unsigned 32-bit integers.
!   This is done without overflows
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(ki_sel), intent(in) :: a, b
!
    integer(ki_sel) :: c
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: a1, a2, b1, b2, s1, s2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
!
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
!
    s1 = a1 + b1
    s2 = a2 + b2 + ibits( s1, hbs, hbs )
!
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
!
!
! END OF FUNCTION uiadd
! _____________________
!
    end function uiadd
!
! =========================================================================================
!
    elemental function uisub( a, b ) result( c )
!
! purpose
! _______
!
!   Internal routine used to simulate substraction for unsigned 32-bit integers.
!   This is done without underflows or overflows.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(ki_sel), intent(in) :: a, b
!
    integer(ki_sel) :: c
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: a1, a2, b1, b2, s1, s2
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
!
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
!
    s1 = a1 - b1
    s2 = a2 - b2 + ibits( s1, hbs, hbs )
!
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
!
!
! END OF FUNCTION uisub
! _____________________
!
    end function uisub
!
! =========================================================================================
!
    elemental function uimlt( a, b ) result( c )
!
! purpose
! _______
!
!   Internal routine used to simulate multiplication for unsigned 32-bit integers.
!   This is done without underflows or overflows.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(ki_sel), intent(in) :: a, b
!
    integer(ki_sel) :: c
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(ki_sel) :: a0, a1, a2, a3,  &
                       b0, b1, b2, b3,  &
                       p0, p1, p2, p3
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    a0 = ibits( a, 0, qbs )
    a1 = ibits( a, qbs, qbs )
    a2 = ibits( a, hbs, qbs )
    a3 = ibits( a, tbs, qbs )
!
    b0 = ibits( b, 0, qbs )
    b1 = ibits( b, qbs, qbs )
    b2 = ibits( b, hbs, qbs )
    b3 = ibits( b, tbs, qbs )
!
    p0 = a0 * b0
    p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
    p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
    p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
!
    c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
    c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
    c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )
!
!
! END OF FUNCTION uimlt
! _____________________
!
    end function uimlt
!
! =========================================================================================
!
    subroutine seed_put( put )
!
! purpose
! _______
!
!   Internal subroutine for the STATPACK RANDOM_SEED_ subroutine.
!
! Further Details
! _______________
!
!   This subroutine is adapted from:
!
!   (1)  Hennecke, M., 1995:
!             A Fortran90 interface to random number generation.
!             Computer Physics Communications, Volume 90, Number 1, 117-120
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Char_Constants,   only : random_error5
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer, intent(in), dimension(:) :: put
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer                               :: i, i1, i2
    integer,         dimension(ratio_i)   :: put2
    integer(ki_sel), dimension(num_seeds) :: iseed
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='random_seed_'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( size(put) >= seed_size ) then
!
        do i = 1, num_seeds
!
            i1 = 1+((i-1)*ratio_i)
            i2 = i*ratio_i
            put2(:ratio_i) = put(i1:i2)
            iseed(i)       = transfer(put2(:ratio_i), s1)
!
!BUG: For STATPACK on MACOS10 machines with the G95 compiler,
!     replace the following line with the four lines of codes
!     above due to a bug in the TRANSFER intrinsic function.
!
!            iseed(i) = transfer(put(1+((i-1)*ratio_i):i*ratio_i), s1)
!
#ifdef _RANDOM_NOINT32
            iseed(i) = ibits( iseed(i), 0, fbs )
#endif
!
        end do
!
        select case (method)
!
            case(1,6)
!
!               SET THE SEED FOR THE Marsaglia's KISS RANDOM NUMBER GENERATOR.
!
                x = iseed(1)
                y = iseed(2)
                z = iseed(3)
                w = iseed(4)
!
            case(2,7)
!
!               SET THE SEED FOR THE Marsaglia's KISS RANDOM NUMBER GENERATOR.
!
                x2 = iseed(1)
                y2 = iseed(2)
                z2 = iseed(3)
                w2 = iseed(4)
                c  = iseed(5)
!
            case(3,8)
!
!               SET THE SEED FOR THE L'Ecuyer's LFSR113 RANDOM NUMBER GENERATOR.
!
                s1 = iseed(1)
                s2 = iseed(2)
                s3 = iseed(3)
                s4 = iseed(4)
!
            case(4,9)
!
!               SET THE SEEDS FOR THE MT19937 Mersenne Twister RNG.
!
                mt(0:num_seeds-1) = iseed(:num_seeds)
!
                i1 = 1+((i-1)*ratio_i)
                i2 = i*ratio_i
!
                put2(:ratio_i) = put(i1:i2)
                mti = transfer(put2(:ratio_i), s1)
#ifdef _RANDOM_NOINT32
                mti = ibits( mti, 0, fbs )
#endif
!
            case(5,10)
!
!               SET THE SEEDS FOR THE MEMT19937-II Mersenne Twister RNG.
!
                memt(0:num_seeds-1) = iseed(:num_seeds)
!
                i1 = 1+((i-1)*ratio_i)
                i2 = i*ratio_i
!
                put2(:ratio_i) = put(i1:i2)
                memti = transfer(put2(:ratio_i), s1)
#ifdef _RANDOM_NOINT32
                memti = ibits( memti, 0, fbs )
#endif
!
        end select
!
    else
!
        call merror( name_proc//random_error5 )
!
    endif
!
!
! END OF SUBROUTINE seed_put
! __________________________
!
    end subroutine seed_put
!
! =========================================================================================
!
    subroutine seed_get ( get )
!
! purpose
! _______
!
!   Internal subroutine for the STATPACK RANDOM_SEED_ subroutine.
!
! Further Details
! _______________
!
!   This subroutine is adapted from:
!
!   (1)  Hennecke, M., 1995:
!             A Fortran90 interface to random number generation.
!             Computer Physics Communications, Volume 90, Number 1, 117-120
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Char_Constants,   only : random_error6
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer, intent(out), dimension(:) :: get
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer                                 :: i, nseeds
    integer(ki_sel), dimension(num_seeds+1) :: iseed
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='random_seed_'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    if ( size(get) >= seed_size ) then
!
        nseeds = num_seeds
!
        select case (method)
!
            case(1,6)
!
!               GET THE SEED FOR THE Marsaglia's KISS RANDOM NUMBER GENERATOR.
!
                iseed(:4) = (/ x, y, z, w /)
!
            case(2,7)
!
!               GET THE SEED FOR THE Marsaglia's KISS RANDOM NUMBER GENERATOR.
!
                iseed(:5) = (/ x2, y2, z2, w2, c /)
!
            case(3,8)
!
!               GET THE SEED FOR THE L'Ecuyer's LFSR113 RANDOM NUMBER GENERATOR.
!
                iseed(:4) = (/ s1, s2, s3, s4 /)
!
            case(4,9)
!
!               GET THE SEEDS FOR THE MT19937 Mersenne Twister RNG.
!
                iseed(:nseeds) = mt(0:nseeds-1)
!
                nseeds = nseeds + 1
                iseed(nseeds) = mti

!
            case(5,10)
!
!               GET THE SEEDS FOR THE MEMT19937-II Mersenne Twister RNG.
!
                iseed(:nseeds) = memt(0:nseeds-1)
!
                nseeds = nseeds + 1
                iseed(nseeds) = memti
!
        end select
!
        do i = 1, nseeds
!
            get(1+((i-1)*ratio_i):i*ratio_i) = transfer(iseed(i), get, ratio_i)
!
        end do
!
    else
!
        call merror( name_proc//random_error6 )
!
    endif
!
!
! END OF SUBROUTINE seed_get
! __________________________
!
    end subroutine seed_get
!
!
! =========================================================================================
!                            GAUSSIAN RANDOM FLOATING POINT NUMBER SUBROUTINES
! =========================================================================================
!
!
    function normal_rand_number() result( harvest )
!
! purpose
! _______
!
!   This function returns a Gaussian distributed random real number.
!
!
! Arguments
! _________
!
!   None
!
!
! Further Details
! _______________
!
!   This function uses the Cumulative Density Function (CDF) inversion method to generate
!   a Gaussian random real number. Starting with a random number produced by the STATPACK
!   uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce a
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random real number.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 2 and 3) by the subroutine PPND7 described in the reference (1).
!   This method gives 7 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see:
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Prob_Procedures, only : pinvn
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    harvest1 = rand_number()
!
#ifdef _RANDOM_WITH0
    if ( harvest1<=zero )  harvest1 = nearest( zero, one )
#endif
!
    harvest = pinvn( harvest1 )
!
!
! END OF FUNCTION normal_rand_number
! __________________________________
!
    end function normal_rand_number
!
! =========================================================================================
!
    subroutine normal_random_r0( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random real number HARVEST following the standard
!   Gaussian distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd)
!            A Gaussian distributed random real number.
!
!
! Further Details
! _______________
!
!   This subroutine uses the Cumulative Density Function (CDF) inversion method to generate
!   a Gaussian random real number. Starting with a random number produced by the STATPACK
!   uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce a
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random real number.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 2 and 3) by the subroutine PPND7 described in the reference (1).
!   This method gives 7 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see:
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Prob_Procedures, only : pinvn
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), dimension(1) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    call random_r1( harvest1 )
!
#ifdef _RANDOM_WITH0
    if ( harvest1(1)<=zero )  harvest1(1) = nearest( zero, one )
#endif
!
    harvest = pinvn( harvest1(1) )
!
!
! END OF SUBROUTINE normal_random_r0
! __________________________________
!
    end subroutine normal_random_r0
!
! =========================================================================================
!
    subroutine normal_random_r1( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random vector HARVEST following the standard
!   normal (Gaussian) distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:)
!            A Gaussian distributed random real vector.
!
!
! Further Details
! _______________
!
!   This subroutines uses the Cumulative Density Function (CDF) inversion method to generate
!   Gaussian random numbers. Starting with random numbers produced by the STATPACK
!   uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random numbers.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 2 and 3) by the subroutine PPND7 described in the reference (1).
!   This method gives 7 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see:
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _ALLOC
    use Char_Constants,  only : allocate_error
#endif
    use Prob_Procedures, only : pinvn
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
#ifdef _ALLOC
    real(stnd), dimension(:), allocatable :: harvest1
!
    integer      :: iok
!
#else
    real(stnd), dimension(size(harvest)) :: harvest1
#endif
!
    integer(i4b) :: m
!
!
! PARAMETERS
! __________
!
#ifdef _ALLOC
    character(len=*),  parameter :: name_proc='normal_random_number_'
#endif
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    m = size( harvest )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b )  return
!
#ifdef _ALLOC
!
!   ALLOCATE WORK VARIABLE.
!
    allocate( harvest1(m), stat = iok )
!
    if ( iok/=0 ) then
!
        call merror( name_proc//allocate_error )
!
    end if
#endif
!
    call random_r1( harvest1(:m) )
!
#ifdef _RANDOM_WITH0
    where ( harvest1(:m)<=zero )  harvest1(:m) = nearest( zero, one )
#endif
!
    harvest(:m) = pinvn( harvest1(:m) )
!
#ifdef _ALLOC
!
!   DEALLOCATE WORK VARIABLE.
!
    deallocate( harvest1 )
#endif
!
!
! END OF SUBROUTINE normal_random_r1
! __________________________________
!
    end subroutine normal_random_r1
!
! =========================================================================================
!
    subroutine normal_random_r2( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random matrix HARVEST following the standard
!   normal (Gaussian) distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:)
!            A Gaussian distributed random real matrix.
!
!
! Further Details
! _______________
!
!   This subroutines uses the Cumulative Density Function (CDF) inversion method to generate
!   Gaussian random real numbers. Starting with random numbers produced by the STATPACK
!   uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random real numbers.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 2 and 3) by the subroutine PPND7 described in the reference (1).
!   This method gives 7 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _ALLOC
    use Char_Constants,  only : allocate_error
#endif
    use Prob_Procedures, only : pinvn
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
#ifdef _ALLOC
    real(stnd), dimension(:), allocatable :: harvest1
#else
    real(stnd), dimension(size(harvest,1)) :: harvest1
#endif
!
    integer(i4b) :: m, n, i
!
#ifdef _OPENMP
    logical :: test_par, ompparallel, ompnested
#endif
!
!
! PARAMETERS
! __________
!
#ifdef _ALLOC
    character(len=*),  parameter :: name_proc='normal_random_number_'
#endif
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    m = size( harvest, 1 )
    n = size( harvest, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b .or. n<=0_i4b )  return
!
#ifdef _OPENMP
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    ompnested   = omp_get_nested() .and. ompparallel
    test_par    = ( .not.ompparallel .or. ompnested ) .and.      &
                  m*n>=omp_limit                      .and.      &
                  i>1_i4b
#endif
!
!$OMP PARALLEL IF(test_par),PRIVATE(i,harvest1)
!
#ifdef _ALLOC
!
!   ALLOCATE WORK VARIABLE.
!
    allocate( harvest1(m), stat = i )
!
    if ( i/=0_i4b ) then
!
        call merror( name_proc//allocate_error )
!
    end if
#endif
!
!$OMP DO SCHEDULE(STATIC)
!
    do i = 1_i4b, n
!
        call random_r1( harvest1(:m) )
!
#ifdef _RANDOM_WITH0
        where ( harvest1(:m)<=zero )  harvest1(:m) = nearest( zero, one )
#endif
!
        harvest(:m,i) = pinvn( harvest1(:m) )
!
    end do
!
!$OMP END DO
!
#ifdef _ALLOC
!
!   DEALLOCATE WORK VARIABLE.
!
    deallocate( harvest1 )
#endif
!
!$OMP END PARALLEL
!
!
! END OF SUBROUTINE normal_random_r2
! __________________________________
!
    end subroutine normal_random_r2
!
! =========================================================================================
!
    function normal_rand_number2() result( harvest )
!
! purpose
! _______
!
!   This function returns a Gaussian distributed real random number of kind extd.
!
!
! Arguments
! _________
!
!   None
!
!
! Further Details
! _______________
!
!   This function uses the Cumulative Density Function (CDF) inversion method to generate
!   a Gaussian random real number of kind extd. Starting with a random number produced by
!   the STATPACK uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce a
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random real number of kind extd.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 7) by the subroutine PPND16 described in the reference (1).
!   This method gives about 16 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   This function is more accurate than NORMAL_RAND_NUMBER function, but it is slower.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Prob_Procedures, only : pinvn2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(extd) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    harvest1 = rand_number()
!
#ifdef _RANDOM_WITH0
    if ( harvest1<=zero )  harvest1 = nearest( zero, one )
#endif
!
    harvest = pinvn2( real( harvest1, extd ) )
!
!
! END OF FUNCTION normal_rand_number2
! ___________________________________
!
    end function normal_rand_number2
!
! =========================================================================================
!
    subroutine normal_random2_r0( harvest )
!
! purpose
! _______
!
!   This subroutine  returns a Gaussian distributed real random number of kind extd.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(extd)
!            A Gaussian distributed random real number of kind extd.
!
!
! Further Details
! _______________
!
!   This subroutine uses the Cumulative Density Function (CDF) inversion method to generate
!   a Gaussian random ral number of kind extd. Starting with a random number produced by
!   the STATPACK uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce a
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random real number of kind extd.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 7) by the subroutine PPND16 described in the reference (1).
!   This method gives about 16 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Prob_Procedures, only : pinvn2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(extd), intent(out) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), dimension(1) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    call random_r1( harvest1 )
!
#ifdef _RANDOM_WITH0
    if ( harvest1(1)<=zero )  harvest1(1) = nearest( zero, one )
#endif
!
    harvest = pinvn2( real( harvest1(1), extd) )
!
!
! END OF SUBROUTINE normal_random2_r0
! ___________________________________
!
    end subroutine normal_random2_r0
!
! =========================================================================================
!
    subroutine normal_random2_r1( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random real vector of kind extd following the standard
!   normal (Gaussian) distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(extd), dimension(:)
!            A Gaussian distributed random real vector of kind extd.
!
!
! Further Details
! _______________
!
!   This subroutines uses the Cumulative Density Function (CDF) inversion method to generate
!   Gaussian random real numbers of kind extd. Starting with random numbers produced by the
!   STATPACK uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random real numbers of kind extd.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 7) by the subroutine PPND16 described in the reference (1).
!   This method gives about 16 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _ALLOC
    use Char_Constants,  only : allocate_error
#endif
    use Prob_Procedures, only : pinvn2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(extd), intent(out), dimension(:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
#ifdef _ALLOC
    real(stnd), dimension(:), allocatable :: harvest1
!
    integer      :: iok
!
#else
    real(stnd), dimension(size(harvest)) :: harvest1
#endif
!
    integer(i4b) :: m
!
!
! PARAMETERS
! __________
!
#ifdef _ALLOC
    character(len=*),  parameter :: name_proc='normal_random_number2_'
#endif
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    m = size( harvest )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b )  return
!
#ifdef _ALLOC
!
!   ALLOCATE WORK VARIABLE.
!
    allocate( harvest1(m), stat = iok )
!
    if ( iok/=0 ) then
!
        call merror( name_proc//allocate_error )
!
    end if
#endif
!
    call random_r1( harvest1(:m) )
!
#ifdef _RANDOM_WITH0
    where ( harvest1(:m)<=zero )  harvest1(:m) = nearest( zero, one )
#endif
!
    harvest(:m) = pinvn2( real( harvest1(:m), extd ) )
!
#ifdef _ALLOC
!
!   DEALLOCATE WORK VARIABLE.
!
    deallocate( harvest1 )
#endif
!
!
! END OF SUBROUTINE normal_random2_r1
! __________________________________
!
    end subroutine normal_random2_r1
!
! =========================================================================================
!
    subroutine normal_random2_r2( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random real matrix of kind extd following the standard
!   normal (Gaussian) distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(extd), dimension(:,:)
!            A Gaussian distributed random real matrix of kind extd.
!
!
! Further Details
! _______________
!
!   This subroutines uses the Cumulative Density Function (CDF) inversion method to generate
!   Gaussian random real numbers of kind extd. Starting with random numbers produced by the
!   STATPACK uniform random number generator that can produce random numbers with the
!   uniform distribution over the continuous range (0, 1) (denoted U(0, 1)), the CDF
!   method simply inverts the CDF of a standard Gaussian distribution to produce
!   standard Gaussian  (e.g. a Gaussian distribution with mean zero and standard
!   deviation one) random real numbers of kind extd.
!
!   The inverse Gaussian CDF is approximated to high precision using rational approximations
!   (polynomials with degree 7) by the subroutine PPND16 described in the reference (1).
!   This method gives about 16 decimal digits of accuracy in the range [10**(-316), 1-10**(-316)]
!   if computations are done in double or higher precision.
!
!   For more details on Uniform and Gaussian random number generators or the approximation
!   of the inverse Gaussian CDF used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Wichura, M.J., 1988:
!           Algorithm AS 241: The percentage points of the normal distribution.
!           Appl. Statis. 37, 3, 477-484.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _ALLOC
    use Char_Constants,  only : allocate_error
#endif
    use Prob_Procedures, only : pinvn2
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(extd), intent(out), dimension(:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
#ifdef _ALLOC
    real(stnd), dimension(:), allocatable :: harvest1
#else
    real(stnd), dimension(size(harvest,1)) :: harvest1
#endif
!
    integer(i4b) :: m, n, i
!
#ifdef _OPENMP
    logical :: test_par, ompparallel, ompnested
#endif
!
!
! PARAMETERS
! __________
!
#ifdef _ALLOC
    character(len=*),  parameter :: name_proc='normal_random_number2_'
#endif
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    m = size( harvest, 1 )
    n = size( harvest, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b .or. n<=0_i4b )  return
!
#ifdef _OPENMP
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    ompnested   = omp_get_nested() .and. ompparallel
    test_par    = ( .not.ompparallel .or. ompnested ) .and.      &
                  m*n>=omp_limit                      .and.      &
                  i>1_i4b
#endif
!
!$OMP PARALLEL IF(test_par),PRIVATE(i,harvest1)
!
#ifdef _ALLOC
!
!   ALLOCATE WORK VARIABLE.
!
    allocate( harvest1(m), stat = i )
!
    if ( i/=0_i4b ) then
!
        call merror( name_proc//allocate_error )
!
    end if
#endif
!
!$OMP DO SCHEDULE(STATIC)
!
    do i = 1_i4b, n
!
        call random_r1( harvest1(:m) )
!
#ifdef _RANDOM_WITH0
        where ( harvest1(:m)<=zero )  harvest1(:m) = nearest( zero, one )
#endif
!
        harvest(:m,i) = pinvn2( real( harvest1(:m), extd ) )
!
    end do
!
!$OMP END DO
!
#ifdef _ALLOC
!
!   DEALLOCATE WORK VARIABLE.
!
    deallocate( harvest1 )
#endif
!
!$OMP END PARALLEL
!
!
! END OF SUBROUTINE normal_random2_r2
! ___________________________________
!
    end subroutine normal_random2_r2
!
! =========================================================================================
!
    function normal_rand_number3() result( harvest )
!
! purpose
! _______
!
!   This function returns a Gaussian distributed random real number.
!
!
! Arguments
! _________
!
!   None
!
!
! Further Details
! _______________
!
!   This function uses the classical Box-Muller method to generate a Gaussian random real number.
!
!   For more details on Uniform and Gaussian random number generators or the Box-Muller method
!   used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Brent, R.P., 1993:
!           Fast Normal Random Generators for vector processors.
!           Report TR-CS-93-04, Computer Sciences Laboratory, Australian National University.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: r, theta
    real(stnd), dimension(2) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    call random_r1( harvest1(:2) )
!
#ifdef _RANDOM_WITH0
    r = sqrt( -two*log( one-harvest1(1) ) )
#else
    r = sqrt( -two*log( harvest1(1) ) )
#endif
    theta   = twopi*harvest1(2)
    harvest = r*sin( theta )
!
!
! END OF FUNCTION normal_rand_number3
! ___________________________________
!
    end function normal_rand_number3
!
! =========================================================================================
!
    subroutine normal_random3_r0( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random real number HARVEST following the standard
!   Gaussian distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd)
!            A Gaussian distributed random real number.
!
!
! Further Details
! _______________
!
!   This subroutine uses the classical Box-Muller method to generate a Gaussian random real number.
!
!   For more details on Uniform and Gaussian random number generators or the Box-Muller method
!   used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Brent, R.P., 1993:
!           Fast Normal Random Generators for vector processors.
!           Report TR-CS-93-04, Computer Sciences Laboratory, Australian National University.
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)               :: r, theta
    real(stnd), dimension(2) :: harvest1
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    call random_r1( harvest1(:2) )
!
#ifdef _RANDOM_WITH0
    r = sqrt( -two*log( one-harvest1(1) ) )
#else
    r = sqrt( -two*log( harvest1(1) ) )
#endif
    theta   = twopi*harvest1(2)
    harvest = r*sin( theta )
!
!
! END OF SUBROUTINE normal_random3_r0
! ___________________________________
!
    end subroutine normal_random3_r0
!
! =========================================================================================
!
    subroutine normal_random3_r1( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random real vector HARVEST following the standard
!   normal (Gaussian) distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:)
!            A Gaussian distributed random real vector.
!
!
! Further Details
! _______________
!
!   This subroutine uses the classical Box-Muller method to generate Gaussian random real numbers.
!   The computations are parallelized if OPENMP is used.
!
!   For more details on Uniform and Gaussian random number generators or the Box-Muller method
!   used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Brent, R.P., 1993:
!           Fast Normal Random Generators for vector processors.
!           Report TR-CS-93-04, Computer Sciences Laboratory, Australian National University.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _ALLOC
    use Char_Constants,  only : allocate_error
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: r, theta
!
#ifdef _ALLOC
    real(stnd), dimension(:), allocatable :: harvest1
!
#else
    real(stnd), dimension(size(harvest)+1) :: harvest1
#endif
!
    integer(i4b) :: m, j
!
#ifdef _OPENMP
    logical :: test_par, ompparallel, ompnested
#endif
!
!
! PARAMETERS
! __________
!
#ifdef _ALLOC
    character(len=*),  parameter :: name_proc='normal_random_number3_'
#endif
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    m = size( harvest )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b )  return
!
#ifdef _ALLOC
!
!   ALLOCATE WORK VARIABLE.
!
    allocate( harvest1(m+1_i4b), stat = j )
!
    if ( j/=0_i4b ) then
!
        call merror( name_proc//allocate_error )
!
    end if
#endif
!
    call random_r1( harvest1(:m+1_i4b) )
!
#ifdef _OPENMP
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    ompnested   = omp_get_nested() .and. ompparallel
    test_par    = ( .not.ompparallel .or. ompnested ) .and.      &
                  m>=omp_limit                        .and.      &
                  j>1_i4b
#endif
!
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC)     &
!$OMP            ,PRIVATE(j,r,theta)
!
    do j = 1_i4b, m-1_i4b, 2_i4b
!
#ifdef _RANDOM_WITH0
        r = sqrt( -two*log( one-harvest1(j) ) )
#else
        r = sqrt( -two*log( harvest1(j) ) )
#endif
!
        theta            = twopi*harvest1(j+1_i4b)
        harvest(j)       = r*sin( theta )
        harvest(j+1_i4b) = r*cos( theta )
!
    end do
!
!$OMP END PARALLEL DO
!
    if ( m/=2_i4b*(m/2_i4b) ) then
!
#ifdef _RANDOM_WITH0
        r = sqrt( -two*log( one-harvest1(m) ) )
#else
        r = sqrt( -two*log( harvest1(m) ) )
#endif
!
        theta      = twopi*harvest1(m+1_i4b)
        harvest(m) = r*sin( theta )
!
    end if
!
#ifdef _ALLOC
!
!   DEALLOCATE WORK VARIABLE.
!
    deallocate( harvest1 )
#endif
!
!
! END OF SUBROUTINE normal_random3_r1
! ___________________________________
!
    end subroutine normal_random3_r1
!
! =========================================================================================
!
    subroutine normal_random3_r2( harvest )
!
! purpose
! _______
!
!   This subroutine returns a random matrix HARVEST following the standard
!   normal (Gaussian) distribution.
!
!
! Arguments
! _________
!
!   HARVEST  (OUTPUT) real(stnd), dimension(:,:)
!            A Gaussian distributed random real matrix.
!
!
! Further Details
! _______________
!
!   This subroutine uses the classical Box-Muller method to generate Gaussian random real numbers.
!   The computations are parallelized if OPENMP is used.
!
!   For more details on Uniform and Gaussian random number generators or the Box-Muller method
!   used here, see :
!
!   (1) Devroye, L., 1986:
!           Non-Uniform Random Variate Generation. Springer-Verlag,
!           http://cg.scs.carleton.ca/~luc/rnbookindex.html, New York.
!
!   (2) Thomas, D.B., Luk, W., Leong, P.H.W., and Villasenor, J.D., 2007:
!           Gaussian random number generators.
!           ACM Comput. Surv. 39, 4, Article 11 (October 2007), 38 pages,
!           DOI = 10.1145/1287620.1287622 (http://doi.acm.org/10.1145/1287620.1287622)
!
!   (3) Brent, R.P., 1993:
!           Fast Normal Random Generators for vector processors.
!           Report TR-CS-93-04, Computer Sciences Laboratory, Australian National University.
!
!
! __________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
#ifdef _ALLOC
    use Char_Constants,  only : allocate_error
#endif
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(out), dimension(:,:) :: harvest
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd) :: r, theta
!
#ifdef _ALLOC
    real(stnd), dimension(:), allocatable :: harvest1
!
#else
    real(stnd), dimension(size(harvest,1)+1) :: harvest1
#endif
!
    integer(i4b) :: m, n, j, i
!
#ifdef _OPENMP
    logical :: test_par, ompparallel, ompnested
#endif
!
!
! PARAMETERS
! __________
!
#ifdef _ALLOC
    character(len=*),  parameter :: name_proc='normal_random_number3_'
#endif
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    m = size( harvest, 1 )
    n = size( harvest, 2 )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( m<=0_i4b .or. n<=0_i4b )  return
!
#ifdef _OPENMP
    i           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    ompnested   = omp_get_nested() .and. ompparallel
    test_par    = ( .not.ompparallel .or. ompnested ) .and.      &
                  m*n>=omp_limit                      .and.      &
                  i>1_i4b
#endif
!
!$OMP PARALLEL IF(test_par),PRIVATE(i,j,r,theta,harvest1)
!
#ifdef _ALLOC
!
!   ALLOCATE WORK VARIABLE.
!
    allocate( harvest1(m+1_i4b), stat = i )
!
    if ( i/=0_i4b ) then
!
        call merror( name_proc//allocate_error )
!
    end if
#endif
!
!$OMP DO SCHEDULE(STATIC)
!
    do i = 1_i4b, n
!
        call random_r1( harvest1(:m+1_i4b) )
!
        do j = 1_i4b, m-1_i4b, 2_i4b
!
#ifdef _RANDOM_WITH0
            r = sqrt( -two*log( one-harvest1(j) ) )
#else
            r = sqrt( -two*log( harvest1(j) ) )
#endif
!
            theta              = twopi*harvest1(j+1_i4b)
            harvest(j,i)       = r*sin( theta )
            harvest(j+1_i4b,i) = r*cos( theta )
!
        end do
!
        if ( m/=2_i4b*(m/2_i4b) ) then
!
#ifdef _RANDOM_WITH0
            r = sqrt( -two*log( one-harvest1(m) ) )
#else
            r = sqrt( -two*log( harvest1(m) ) )
#endif
!
            theta        = twopi*harvest1(m+1_i4b)
            harvest(m,i) = r*sin( theta )
!
        end if
!
    end do
!
!$OMP END DO
!
#ifdef _ALLOC
!
!   DEALLOCATE WORK VARIABLE.
!
    deallocate( harvest1 )
#endif
!
!$OMP END PARALLEL
!
!
! END OF SUBROUTINE normal_random3_r2
! ___________________________________
!
    end subroutine normal_random3_r2
!
!
! =========================================================================================
!               PSEUDO-RANDOM ORTHOGONAL MATRICES FROM THE HAAR DISTRIBUTION
!               PSEUDO-RANDOM SYMMETRIC MATRICES WITH A PRESCRIBED SPECTRUM
!               PSEUDO-RANDOM MATRICES WITH A PRESCRIBED SINGULAR VALUE DISTRIBUTION
! =========================================================================================
!
!
    subroutine random_qr_cmp( mat, diagr, beta, fillr, initseed )
!
! Purpose
! _______
!
!   RANDOM_QR_CMP generates the first k columns of a pseudo-random QR factorization of a
!   hypothetical real n-by-n matrix MAT, whose elements follow independently a Laplace_Gauss(0;1)
!   distribution (e.g. the standard normal distribution):
!
!                              MAT = Q * R
!
!   where Q is a pseudo-random orthogonal matrix following the Haar distribution from the group
!   of orthogonal matrices and R is upper triangular. The upper-diagonal elements of R follow a
!   Laplace_Gauss(0;1) distribution (e.g. the standard normal distribution) and the squares
!   of the diagonal elements of R, (R(i,i))**(2) follow a chi-squared distribution with n-i+1
!   degrees of freedom.
!
!
! Arguments
! _________
!
!   MAT      (OUTPUT) real(stnd), dimension(:,:)
!            On exit, the elements above the diagonal of the array MAT contain  
!            the corresponding elements of R if FILLR is present and set to true;
!            otherwise the upper-diagonal elements of MAT are not modified.
!            The elements on and below the diagonal, with the arrays BETA
!            and DIAGR, represent the first k columns (with k=size(MAT,2) and k<=n )
!            of a pseudo-random orthogonal matrix Q following the Haar distribution
!            as a product of elementary reflectors and a diagonal matrix,
!            see Further Details.
!
!            The shape of MAT must verify :
!                         size( MAT, 2 ) <= size( MAT, 1 ).
!
!   DIAGR    (OUTPUT) real(stnd), dimension(:) 
!            On exit, the diagonal elements of the matrix R,
!            see Further Details.
!
!            The size of DIAGR must verify :
!                         size( DIAGR ) = size( MAT, 2 ) <= size( MAT, 1 ).
!
!   BETA     (OUTPUT) real(stnd), dimension(:)
!            On exit, the scalars factors of the elementary reflectors,
!            see Further Details.
!
!            The size of BETA must verify :
!                         size( BETA ) = size( DIAGR ) = size( MAT, 2 ) <= size( MAT, 1 ).
!
!   FILLR    (INPUT, OPTIONAL) logical(lgl)
!            on entry, if FILLR is set to true, the super-diagonal elements of R
!            are generated in MAT on exit.
!            If FILLR is set to false, the super-diagonal elements of MAT are not
!            modified or used.
!
!            The default is FILLR = false.
!
!   INITSEED (INPUT, OPTIONAL) logical(lgl)
!            On entry, if INITSEED=true, a call to RANDOM_SEED_() without arguments
!            is done in the subroutine, in order to initiates a non-repeatable reset
!            of the seed used by the STATPACK random generator.
!
!            The default is INITSEED = false.
!
!
! Further Details
! _______________
!
!   RANDOM_QR_CMP uses the method described in the reference (1), based on Householder
!   transformations, for generating the first k columns of a n-by-n pseudo-random orthogonal
!   matrices Q distributed according to the Haar measure over the orthogonal group of order n,
!   in a factored form.
!
!   The pseudo-random orthogonal matrix Q is represented as a product of n elementary reflectors
!   and of a diagonal matrix
!
!            Q = H(1) * H(2) * ... * H(n) * diag(sign(DIAGR))
!
!   Each H(i) has the form
!
!            H(i) = I + beta * ( v * v' ) ,
!                      
!   where beta is a real scalar and v is a real n-element vector with v(1:i-1) = 0.
!   v(i:n) is stored on exit in MAT(i:n,i) and beta in BETA(i).
!
!   diag(sign(DIAGR)) is the n-by-n diagonal matrix with diagonal elements equal to
!   sign( one, DIAGR) (e.g. its ith diagonal element equals to one if DIAGR(i) is positive
!   and -one otherwise).
!
!   It is possible to compute only the first k columns of Q and R by restricting the number
!   of columns of MAT and the sizes of DIAGR and BETA on entry of the subroutine.
!
!   Finally, note that the computations are parallelized if OPENMP is used.
!
!   Q can be generated with the help of subroutine ORTHO_GEN_RANDOM_QR or can be applied
!   to a vector or a matrix with the help of subroutine APPLY_Q_QR.
!
!   For further details on the QR factorization and its use or pseudo-random orthogonal
!   matrices distributed according to the Haar measure over the orthogonal group, see:
!
!   (1) Stewart, G.W., 1980:
!          The efficient generation of random orthogonal matrices with
!          an application to condition estimators.
!          SIAM J. Numer. Anal., 17, 403-409
!
!   (2) Lawson, C.L., and Hanson, R.J., 1974:
!          Solving least square problems.
!          Prentice-Hall.
!
!   (3) Golub, G.H., and Van Loan, C.F., 1996:
!          Matrix Computations.
!          3rd ed. The Johns Hopkins University Press, Baltimore.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,          only : assert, assert_eq
    use Hous_Procedures,    only : hous1
    use Prob_Procedures,    only : pinvn
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(out), dimension(:,:) :: mat
    real(stnd),   intent(out), dimension(:)   :: diagr, beta
!
    logical(lgl), intent(in),  optional :: fillr, initseed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), dimension( size(mat,1) ) :: tmp
!
    integer(i4b) :: n, k, j
!
    logical(lgl) :: fill_r, init_seed
!
#ifndef _USE_PGI
#ifdef _OPENMP
    logical :: test_par, ompparallel, ompnested
#endif
#endif
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='random_qr_cmp'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
!
    k = assert_eq( int(size(beta),i4b),   &
                   int(size(diagr),i4b),  &
                   int(size(mat,2),i4b),  &
                   name_proc              )
!
    call assert( logical( k<=n, lgl ), name_proc )
!
    init_seed = false
!
    if ( present(initseed) ) then
!
        init_seed = initseed
!
    end if
!
    fill_r = false
!
    if ( present(fillr) ) then
!
        fill_r = fillr
!
    end if
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b .or. k<=0_i4b ) return
!
#ifndef _USE_PGI
#ifdef _OPENMP
    j           = omp_get_max_threads()
    ompparallel = omp_in_parallel()
    ompnested   = omp_get_nested() .and. ompparallel
    test_par    = ( .not.ompparallel .or. ompnested ) .and.      &
                  n*k>=omp_limit                      .and.      &
                  j>1_i4b
#endif
#endif
!
!   INITIALIZE THE RANDOM GENERATOR IF REQUIRED.
!
    if ( init_seed ) then
!
#ifndef _USE_PGI
!$OMP SINGLE
#endif
        call random_seed_()
#ifndef _USE_PGI
!$OMP END SINGLE
#endif
!
    end if
!
    if ( fill_r ) then
!
!       GENERATE A RANDOM ORTHOGONAL MATRIX AND AN UPPER
!       TRIANGULAR MATRIX IN FACTORED FORM.
!
#ifndef _USE_PGI
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC)     &
!$OMP            ,PRIVATE(j,tmp)
#endif
!
        do j = 1_i4b, k
!
            call random_r1( tmp(1_i4b:n) )
!
#ifdef _RANDOM_WITH0
            tmp(1_i4b:n) = max( tmp(1_i4b:n), nearest( zero, one ) )
#endif
!    
!           GENERATE UPPER DIAGONAL ELEMENTS OF TRIANGULAR MATRIX.
!
            mat(i4b:n,j) = pinvn( tmp(i4b:n) )
!    
!           GENERATE ELEMENTARY REFLECTOR H(j) TO ANNIHILATE mat(j+1:n,j) .
!    
            call hous1( mat(j:n,j), beta(j), diagr(j) )
!
        end do
!
#ifndef _USE_PGI
!$OMP END PARALLEL DO
#endif
!
    else
!
!       GENERATE A RANDOM ORTHOGONAL MATRIX IN FACTORED FORM.
!
#ifndef _USE_PGI
!$OMP PARALLEL DO IF(test_par), SCHEDULE(STATIC)     &
!$OMP            ,PRIVATE(j,tmp)
#endif
!
        do j = 1_i4b, k
!
            call random_r1( tmp(j:n) )
!
#ifdef _RANDOM_WITH0
            tmp(j:n) =  max( tmp(j:n), nearest( zero, one ) )
#endif
!
            mat(j:n,j) = pinvn( tmp(j:n) )
!    
!           GENERATE ELEMENTARY REFLECTOR H(j) TO ANNIHILATE mat(j+1:n,j) .
!    
            call hous1( mat(j:n,j), beta(j), diagr(j) )
!
        end do
!
#ifndef _USE_PGI
!$OMP END PARALLEL DO
#endif
!
    end if
!
!
! END OF SUBROUTINE random_qr_cmp
! _______________________________
!
    end subroutine random_qr_cmp
!
! =========================================================================================
!
    subroutine ortho_gen_random_qr( mat, diagr, beta )
!
! Purpose
! _______
!
!   ORTHO_GEN_RANDOM_QR generates a n-by-n real pseudo-random orthogonal matrix following
!   the Haar distribution, which is defined as the product of n elementary reflectors of
!   order n and of a n-by-n diagonal matrix with diagonal elements equal to sign( one, DIAGR):
!
!            Q = H(1) * H(2) * ... * H(n) * diag(sign(DIAGR))
!
!   as returned by RANDOM_QR_CMP.
!
!   Optionnally, it is possible to generate only the first k columns of Q by restricting
!   arguments MAT, BETA and DIAGR to the first k columns or elements of the corresponding
!   arguments as returned by RANDOM_QR_CMP.
!
!
! Arguments
! _________
!
!   MAT      (INPUT/OUTPUT) real(stnd), dimension(:,:)
!            On entry, the i-th column must contain the vector which
!            defines the elementary reflector H(i), for i = 1,2,...,k,
!            ( with k<=n and k=size(MAT,2) ) as returned by RANDOM_QR_CMP
!            in its array argument MAT.
!            On exit, the first k columns of the pseudo-random n-by-n orthogonal
!            matrix Q.
!
!            The shape of MAT must verify : size( MAT, 2 ) <= size( MAT, 1 ).
!
!   DIAGR    (INPUT) real(stnd), dimension(:) 
!            On entry, the diagonal elements of the matrix R, as returned
!            by RANDOM_QR_CMP in its argument DIAGR.
!
!            The size of DIAGR must verify :
!                         size( DIAGR ) = size( MAT, 2 ).
!
!   BETA     (INPUT) real(stnd), dimension(:)
!            On entry, BETA(i) must contain the scalar factor of the elementary
!            reflector H(i), as returned by RANDOM_QR_CMP in its argument BETA.
!
!            The size of BETA must verify :
!                         size( BETA ) = size( DIAGR ) = size( MAT, 2 ) .
!
!
! Further Details
! _______________
!
!   This subroutine used a blocked algorithm for agregating the Householder
!   transformations (e.g. the elementary reflectors) stored in the lower triangle
!   of MAT and generating the pseudo-random orthogonal matrix Q of the QR factorization
!   returned by subroutine RANDOM_QR_CMP.
!
!   Furthermore, the computations are parallelized if OPENMP is used.
!
!   For further details on the QR factorization and its use or the blocked
!   algorithm used here, see:
!
!   (1) Lawson, C.L., and Hanson, R.J., 1974:
!           Solving least square problems.
!           Prentice-Hall.
!
!   (2) Golub, G.H., and Van Loan, C.F., 1996:
!           Matrix Computations.
!           3rd ed. The Johns Hopkins University Press, Baltimore.
!
!   (3) Dongarra, J.J., Sorensen, D.C., and Hammarling, S.J., 1989:
!           Block reduction of matrices to condensed form for eigenvalue computations. 
!           J. of Computational and Applied Mathematics, Vol. 27, pp. 215-227.
!
!   (4) Walker, H.F., 1988:
!           Implementation of the GMRES method using Householder transformations.
!           Siam J. Sci. Stat. Comput., Vol. 9, No 1, pp. 152-163.
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Utilities,     only : assert, assert_eq
    use QR_Procedures, only : ortho_gen_qr
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(inout), dimension(:,:) :: mat
    real(stnd),   intent(in),    dimension(:)   :: diagr, beta
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: n, k, j
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='ortho_gen_random_qr'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    n = size( mat, 1 )
!
    k = assert_eq( int(size(beta),i4b),   &
                   int(size(diagr),i4b),  &
                   int(size(mat,2),i4b),  &
                   name_proc              )
!
    call assert( logical( k<=n, lgl ), name_proc )
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b .or. k<=0_i4b ) return
!
!   GENERATE THE RANDOM ORTHOGONAL MATRIX.
!
    call ortho_gen_qr( mat(:n,:k), beta(:k) )
!
    do j = 1_i4b, k
!
        if ( diagr(j) > zero ) cycle
!
        mat(:n,j) = -mat(:n,j)
!
    end do
!
!
! END OF SUBROUTINE ortho_gen_random_qr
! _____________________________________
!
    end subroutine ortho_gen_random_qr
!
! =========================================================================================
!
    subroutine gen_random_sym_mat( eigval, mat, eigvec, initseed )
!
! Purpose
! _______
!
!   GEN_RANDOM_SYM_MAT generates a pseudo-random n-by-n real symmetric matrix
!   with prescribed eigenvalues.
!
!   Optionally, the corresponding eigenvectors of the generated
!   pseudo-random n-by-n real symmetric matrix can be output if required.
!
!
! Arguments
! _________
!
!   EIGVAL   (INPUT) real(stnd), dimension(:)
!            On entry, the prescribed eigenvalues of the pseudo-random n-by-n
!            real symmetric matrix.
!            IF size(EIGVAL)<n, the other eigenvalues are assumed to be zero.
!
!            The size of EIGVAL must verify :
!                            size( EIGVAL ) <= size( MAT, 1 ) = size( MAT, 2 ) . 
!
!   MAT      (OUTPUT) real(stnd), dimension(:,:)
!            On exit, the pseudo-random n-by-n real symmetric matrix.
!
!            The shape of MAT must verify size( MAT, 2 ) = size( MAT, 1 ).
!
!   EIGVEC   (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!            On exit, the pseudo-random eigenvectors corresponding
!            to the eigenvalues prescribed in EIGVAL. The eigenvectors
!            are returned columnwise.
!
!            The shape of EIGVEC must verify:
!
!            - size( EIGVEC, 1 ) = size( MAT, 1 ) = size( MAT, 2 )
!            - size( EIGVEC, 2 ) = size( EIGVAL ) .
!
!   INITSEED (INPUT, OPTIONAL) logical(lgl)
!            On entry, if INITSEED=true, a call to RANDOM_SEED_() without arguments
!            is done in the subroutine, in order to initiates a non-repeatable reset
!            of the seed used by the STATPACK random generator.
!
!            The default is INITSEED = false.
!
!
! Further Details
! _______________
!
!   Pseudo-random eigenvectors are generated as a pseudo-random orthogonal matrix
!   following the Haar distribution from the group of orthogonal matrices and are
!   computed with the help of subroutines RANDOM_QR_CMP and ORTHO_GEN_RANDOM_QR.
!
!   These computed pseudo-random eigenvectors and the corresponding prescribed
!   eigenvalues are then used to generate a pseudo-random n-by-n symmetric matrix
!   whose eigenvalues are the prescribed eigenvalues.
!
!   These computations are parallelized if OPENMP is used.
!
!   For further details, see:
!
!   (1) Stewart, G.W., 1980:
!          The efficient generation of random orthogonal matrices with
!          an application to condition estimators.
!          SIAM J. Numer. Anal., 17, 403-409
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Char_Constants,  only : allocate_error
    use Utilities,       only : assert, assert_eq
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(in),            dimension(:)   :: eigval
    real(stnd),   intent(out),           dimension(:,:) :: mat
    real(stnd),   intent(out), optional, dimension(:,:) :: eigvec
!
    logical(lgl), intent(in),  optional :: initseed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd),  dimension( size(eigval) )           :: beta, diagr
    real(stnd),  dimension(size(mat,1),size(eigval)) :: eigvec00
    real(stnd),          allocatable, dimension(:,:) :: eigvec0
!
    integer      :: iok
    integer(i4b) :: n, k
!
    logical(lgl)  :: test_eigvec
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='gen_random_sym_mat'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    k = size( eigval )
!
    n = assert_eq( int(size(mat,1),i4b),  &
                   int(size(mat,2),i4b),  &
                   name_proc              )
!
    call assert( logical( k<=n, lgl ), name_proc )
!
    test_eigvec = false
!
    if ( present(eigvec) ) then
!
        call assert( logical( int(size(eigvec,1),i4b)==n, lgl ),  &
                     logical( int(size(eigvec,2),i4b)==k, lgl ),  &
                     name_proc )
!
        test_eigvec = true
!
    end if
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( n<=0_i4b .or. k<=0_i4b ) return
!
!   FIRST GENERATE THE PSEUDO-RANDOM EIGENVECTORS USING mat(:,:) AS WORKSPACE.
!
    call random_qr_cmp( eigvec00(:n,:k), diagr(:k), beta(:k), initseed=initseed )
!
    call ortho_gen_random_qr( eigvec00(:n,:k), diagr(:k), beta(:k) )
!
    if ( test_eigvec ) then
!
        eigvec(:n,:k) = eigvec00(:n,:k)*spread(eigval(:k),1,n)
!
!       FORM THE SYMMETRIC MATRIX FROM THE SPECIFIED EIGENVALUE DISTRIBUTION AND RANDOM EIGENVECTORS.
!
#ifdef _BLAS
        call gemm( 'N', 'T', n, n, k, one, eigvec(1_i4b:n,1_i4b:k), n,          &
                   eigvec00(1_i4b:n,1_i4b:k), n, zero, mat(1_i4b:n,1_i4b:n), n  )
#else
        mat(:n,:n) = matmul( eigvec(:n,:k), transpose(eigvec00(:n,:k)) )
#endif
!
    else
!
!       ALLOCATE WORK ARRAY.
!
        allocate( eigvec0(n,k), stat=iok )
!
        if ( iok/=0 ) then
!
            call merror( name_proc//allocate_error )
!
        end if
!
        eigvec0(:n,:k) = eigvec00(:n,:k)*spread(eigval(:k),1,n)
!
!       FORM THE SYMMETRIC MATRIX FROM THE SPECIFIED EIGENVALUE DISTRIBUTION AND RANDOM EIGENVECTORS.
!
#ifdef _BLAS
        call gemm( 'N', 'T', n, n, k, one, eigvec0(1_i4b:n,1_i4b:k), n,         &
                   eigvec00(1_i4b:n,1_i4b:k), n, zero, mat(1_i4b:n,1_i4b:n), n  )
#else
        mat(:n,:n) = matmul( eigvec0(:n,:k), transpose(eigvec00(:n,:k)) )
#endif
!
!       DEALLOCATE WORK ARRAY.
!
        deallocate( eigvec0 )
!
    end if
!
!
! END OF SUBROUTINE gen_random_sym_mat
! ____________________________________
!
    end subroutine gen_random_sym_mat
!
! =========================================================================================
!
    subroutine gen_random_mat( singval, mat, leftvec, rightvec, initseed )
!
! Purpose
! _______
!
!   GEN_RANDOM_MAT generates a pseudo-random m-by-n real matrix
!   with prescribed singular values.
!
!   Optionnally, the corresponding singular vectors of the generated
!   pseudo-random m-by-n real matrix can be output if required.
!
!
! Arguments
! _________
!
!   SINGVAL  (INPUT) real(stnd), dimension(:)
!            On entry , the prescribed singular values of the pseudo-random m-by-n
!            real matrix.
!            IF size(SINGVAL)<min(m,n), the other eigenvalues are assumed to be zero.
!
!            The size of SINGVAL must verify :
!                            size( SINGVAL ) <= min( size( MAT, 1 ), size( MAT, 2 ) ) . 
!
!   MAT      (OUTPUT) real(stnd), dimension(:,:)
!            On exit, the pseudo-random m-by-n real matrix.
!
!   LEFTVEC  (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!            On exit, the pseudo-random left singular vectors corresponding
!            to the singular values prescribed in SINGVAL. The left singular vectors
!            are returned columnwise.
!
!            The shape of LEFTVEC must verify:
!
!            - size( LEFTVEC, 1 ) = size( MAT, 1 )
!            - size( LEFTVEC, 2 ) = size( SINGVAL ) .
!
!   RIGHTVEC (OUTPUT, OPTIONAL) real(stnd), dimension(:,:)
!            On exit, the pseudo-random right singular vectors corresponding
!            to the singular values prescribed in SINGVAL. The right singular vectors
!            are returned columnwise.
!
!            The shape of RIGHTVEC must verify:
!
!            - size( RIGHTVEC, 1 ) = size( MAT, 2 )
!            - size( RIGHTVEC, 2 ) = size( SINGVAL ) .
!
!   INITSEED (INPUT, OPTIONAL) logical(lgl)
!            On entry, if INITSEED=true, a call to RANDOM_SEED_() without arguments
!            is done in the subroutine, in order to initiates a non-repeatable reset
!            of the seed used by the STATPACK random generator.
!
!            The default is INITSEED = false.
!
!
! Further Details
! _______________
!
!   Pseudo-random singular vectors are generated as pseudo-random orthogonal matrices
!   following the Haar distribution from the group of orthogonal matrices and are
!   computed with the help of subroutines RANDOM_QR_CMP and ORTHO_GEN_RANDOM_QR.
!
!   These computed pseudo-random singular vectors and the corresponding prescribed
!   singular values are then used to generate a pseudo-random m-by-n real matrix
!   whose singular values are the prescribed singular values.
!
!   These computations are parallelized if OPENMP is used.
!
!   For further details, see:
!
!   (1) Stewart, G.W., 1980:
!          The efficient generation of random orthogonal matrices with
!          an application to condition estimators.
!          SIAM J. Numer. Anal., 17, 403-409
!
!
! _________________________________________________________________________________________
!
!
! USED MODULES
! ____________
!
    use Char_Constants,  only : allocate_error
    use Utilities,       only : assert
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd),   intent(in),            dimension(:)   :: singval
    real(stnd),   intent(out),           dimension(:,:) :: mat
    real(stnd),   intent(out), optional, dimension(:,:) :: leftvec, rightvec
!
    logical(lgl), intent(in),  optional :: initseed
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd), allocatable, dimension(:)   :: beta, diagr
    real(stnd), allocatable, dimension(:,:) :: leftvec0, rightvec0
!
    integer      :: iok
    integer(i4b) :: m, n, mn, k
!
    logical(lgl)  :: test_leftvec, test_rightvec
!
!
! PARAMETERS
! __________
!
    character(len=*),  parameter :: name_proc='gen_random_mat'
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
!   TEST THE ARGUMENTS.
!
    m  = size( mat, 1 )
    n  = size( mat, 2 )
    mn = min( m, n )
!
    k = size( singval )
!
    call assert( logical( k<=mn, lgl ), name_proc )
!
    test_leftvec = false
!
    if ( present(leftvec) ) then
!
        call assert( logical( int(size(leftvec,1),i4b)==m, lgl ),  &
                     logical( int(size(leftvec,2),i4b)==k, lgl ),  &
                     name_proc )
!
        test_leftvec = true
!
    end if
!
    test_rightvec = false
!
    if ( present(rightvec) ) then
!
        call assert( logical( int(size(rightvec,1),i4b)==n, lgl ),  &
                     logical( int(size(rightvec,2),i4b)==k, lgl ),  &
                     name_proc )
!
        test_rightvec = true
!
    end if
!
!   QUICK RETURN IF POSSIBLE.
!
    if ( mn<=0_i4b .or. k<=0_i4b ) return
!
!   ALLOCATE WORK ARRAYS.
!
    allocate( leftvec0(m,k), rightvec0(n,k), beta(k), diagr(k), stat=iok )
!
    if ( iok/=0 ) then
!
        call merror( name_proc//allocate_error )
!
    end if
!
!   GENERATE THE PSEUDO-RANDOM LEFT SINGULAR VECTORS.
!
    call random_qr_cmp( leftvec0(:m,:k), diagr(:k), beta(:k), initseed=initseed )
!
    call ortho_gen_random_qr( leftvec0(:m,:k), diagr(:k), beta(:k) )
!
!   GENERATE THE PSEUDO-RANDOM RIGHT SINGULAR VECTORS.
!
    call random_qr_cmp( rightvec0(:n,:k), diagr(:k), beta(:k), initseed=initseed )
!
    call ortho_gen_random_qr( rightvec0(:n,:k), diagr(:k), beta(:k) )
!
    if ( test_leftvec ) then
!
        leftvec(:m,:k) = leftvec0(:m,:k)
!
    end if
!
    if ( test_rightvec ) then
!
        rightvec(:n,:k) = rightvec0(:n,:k)
!
    end if
!
!   FORM THE m-BY-n REAL MATRIX FROM THE SPECIFIED SINGULAR VALUE DISTRIBUTION AND RANDOM SINGULAR VECTORS.
!
#ifdef _BLAS
    rightvec0(:n,:k) = rightvec0(:n,:k)*spread(singval(:k),1,n)
!
    call gemm( 'N', 'T', m, n, k, one, leftvec0(1_i4b:m,1_i4b:k), m,          &
                rightvec0(1_i4b:n,1_i4b:k), n, zero, mat(1_i4b:m,1_i4b:n), m  )
#else
    mat(:m,:n) = matmul( leftvec0(:m,:k)*spread(singval(:k),1,m), transpose(rightvec0(:n,:k)) )
#endif
!
!   DEALLOCATE WORK ARRAYS.
!
    deallocate( leftvec0, rightvec0, beta, diagr  )
!
!
! END OF SUBROUTINE gen_random_mat
! ________________________________
!
    end subroutine gen_random_mat
!
!
! =========================================================================================
!                                   BOOTSTRAP AND PERMUTATION PROCEDURES
! =========================================================================================
!
!
    subroutine simple_shuffle_rv( vec )
!
! purpose
! _______
!
!   This subroutine shuffles all the elements of the real vector VEC.
!
!
! Arguments
! _________
!
!   VEC      (INPUT/OUTPUT) real(stnd), dimension(:)
!            On entry, the real vector to be shuffled.
!
!            On exit, the permuted real vector.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Noreen, E.W., 1989:
!           Computer-intensive methods for testing hypotheses: an introduction.
!           Wiley and Sons, New York, USA, ISBN:978-0-471-61136-3
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    real(stnd), intent(inout), dimension(:) :: vec
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    real(stnd)                       :: temp
    real(stnd), dimension(size(vec)) :: tempvec
!
    integer(i4b) :: i, k, nvec
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nvec = size( vec )
!
    if ( nvec<=0_i4b ) return
!
    call random_r1( tempvec )
!
    do i = 1_i4b, nvec-1_i4b
!
        k = i + int( tempvec(i)*real( nvec-i+1_i4b, stnd ) )
        temp   = vec(k)
        vec(k) = vec(i)
        vec(i) = temp
!
    end do
!
!
! END OF SUBROUTINE simple_shuffle_rv
! ___________________________________
!
    end subroutine simple_shuffle_rv
!
! =========================================================================================
!
    subroutine simple_shuffle_cv( vec )
!
! purpose
! _______
!
!   This subroutine shuffles all the elements of the complex vector VEC.
!
!
! Arguments
! _________
!
!   VEC      (INPUT/OUTPUT) complex(stnd), dimension(:)
!            On entry, the complex vector to be shuffled.
!
!            On exit, the permuted complex vector.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Noreen, E.W., 1989:
!           Computer-intensive methods for testing hypotheses: an introduction.
!           Wiley and Sons, New York, USA, ISBN:978-0-471-61136-3
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    complex(stnd), intent(inout), dimension(:) :: vec
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    complex(stnd)                    :: temp
    real(stnd), dimension(size(vec)) :: tempvec
!
    integer(i4b) :: i, k, nvec
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nvec = size( vec )
!
    if ( nvec<=0_i4b ) return
!
    call random_r1( tempvec )
!
    do i = 1_i4b, nvec-1_i4b
!
        k = i + int( tempvec(i)*real( nvec-i+1_i4b, stnd ) )
        temp   = vec(k)
        vec(k) = vec(i)
        vec(i) = temp
!
    end do
!
!
! END OF SUBROUTINE simple_shuffle_cv
! ___________________________________
!
    end subroutine simple_shuffle_cv
!
! =========================================================================================
!
    subroutine simple_shuffle_iv( vec )
!
! purpose
! _______
!
!   This subroutine shuffles all the elements of the integer vector VEC.
!
!
! Arguments
! _________
!
!   VEC      (INPUT/OUTPUT) intger(i4b), dimension(:)
!            On entry, the integer vector to be shuffled.
!
!            On exit, the permuted integer vector.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Noreen, E.W., 1989:
!           Computer-intensive methods for testing hypotheses: an introduction.
!           Wiley and Sons, New York, USA, ISBN:978-0-471-61136-3
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(inout), dimension(:) :: vec
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: temp
!
    real(stnd), dimension(size(vec)) :: tempvec
!
    integer(i4b)  :: i, k, nvec
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nvec = size( vec )
!
    if ( nvec<=0_i4b ) return
!
    call random_r1( tempvec )
!
    do i = 1_i4b, nvec-1_i4b
!
        k = i + int( tempvec(i)*real( nvec-i+1_i4b, stnd ) )
        temp   = vec(k)
        vec(k) = vec(i)
        vec(i) = temp
!
    end do
!
!
! END OF SUBROUTINE simple_shuffle_iv
! ___________________________________
!
    end subroutine simple_shuffle_iv
!
! =========================================================================================
!
    subroutine drawsample( nsample, pop )
!
! purpose
! _______
!
!   This subroutine may be used to draw a sample, without replacement of size NSAMPLE
!   from a population of size SIZE(POP). On output, the integer vector POP(1:NSAMPLE)
!   indicates which observations are included in the sample.
!   
!   The integer vector POP must be dimensioned at least as large as NSAMPLE
!   in the calling program.
!
!
! Arguments
! _________
!
!   NSAMPLE  (INPUT) intger(i4b)
!            On entry, the size of the sample.
!
!   POP      (OUTPUT) integer(i4b), dimension(:)
!            On exit, the indices of the observations belonging to the sample
!            are in POP(1:NSAMPLE) and the indices of the observations, which
!            are not in the sample are in POP(NSAMPLE+1:).
!
!            The size of POP must greater or equal to NSAMPLE. If this condition is not meet
!            POP(:) is set to -1.
!
!
! Further Details
! _______________
!
!   For more details and algorithm, see:
!
!   (1) Noreen, E.W., 1989:
!           Computer-intensive methods for testing hypotheses: an introduction.
!           Wiley and Sons, New York, USA, ISBN:978-0-471-61136-3
!
!
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in)                :: nsample
    integer(i4b), intent(out), dimension(:) :: pop
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: temp, npop
!
    real(stnd), dimension(nsample) :: tempvec
!
    integer(i4b) :: i, k 
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    npop = size( pop )
!
    if ( npop<=0_i4b .or. nsample<=0_i4b ) return
!
    if ( nsample<=npop ) then
!
        pop(1_i4b:npop) = arth( 1_i4b, 1_i4b, npop )
!
        call random_r1( tempvec )
!
        do i = 1_i4b, nsample
!
            k = i + int( tempvec(i)*real( npop-i+1_i4b, stnd ) )
            temp   = pop(k)
            pop(k) = pop(i)
            pop(i) = temp
!
        end do
!
    else
!
        pop(1_i4b:npop) = -1_i4b
!
    end if
!
!
! END OF SUBROUTINE drawsample
! ____________________________
!
    end subroutine drawsample
!
! =========================================================================================
!
    subroutine drawbootsample( npop, sample )
!
! purpose
! _______
!
!   This subroutine may be used to draw a bootstrap random sample of size SIZE(SAMPLE)
!   from a population of size NPOP. On output, the integer vector SAMPLE indicates which
!   observations are included in the bootstrap sample.
!
!
! Arguments
! _________
!
!   NPOP     (INPUT) intger(i4b)
!            On entry, the size of the population.
!
!   SAMPLE   (OUTPUT) integer(i4b), dimension(:)
!            On exit, the indices of the observations belonging to the sample.
!
!
! Further Details
! _______________
!
!   The sampling is done with replacement, meaning that the sample may contain
!   duplicate observations.
!
!   For more details and algorithm, see:
!
!   (1) Noreen, E.W., 1989:
!           Computer-intensive methods for testing hypotheses: an introduction.
!           Wiley and Sons, New York, USA, ISBN:978-0-471-61136-3
!
!   
! __________________________________________________________________________________________
!
!
! SPECIFICATIONS FOR ARGUMENTS
! ____________________________
!
    integer(i4b), intent(in)                :: npop
    integer(i4b), intent(out), dimension(:) :: sample
!
!
! SPECIFICATIONS FOR LOCAL VARIABLES
! __________________________________
!
    integer(i4b) :: nsample
!
    real(stnd), dimension(size(sample)) :: tempvec
!
!
! EXECUTABLE STATEMENTS
! _____________________
!
    nsample = size( sample )
!
    if ( npop<=0_i4b .or. nsample<=0_i4b ) return
!
    call random_r1( tempvec )
!
    sample(:) = int( tempvec(:)*real( npop, stnd ), i4b ) + 1_i4b
!
!
! END OF SUBROUTINE drawbootsample
! ________________________________
!
    end subroutine drawbootsample
!
! =========================================================================================
!
! **********************
! END OF MODULE Random *
! **********************
!    
end module Random
