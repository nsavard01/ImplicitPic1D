module nitparam_mod
    real(kind=8), parameter :: one = 1.0d0, two = 2.0d0, rtfiv = 2.23606797749978981d0
    real(kind=8), parameter :: tenth = 0.10d0, half = 0.50d0, fournines = one - 1.0d-4
    real(kind=8), parameter :: DFLT_CHOICE1_EXP = (one + rtfiv) * half
    real(kind=8), parameter :: DFLT_CHOICE2_EXP = two
    real(kind=8), parameter :: DFLT_CHOICE2_COEF = one
    real(kind=8), parameter :: DFLT_ETA_CUTOFF = tenth
    real(kind=8), parameter :: DFLT_ETA_MAX = fournines
    real(kind=8), parameter :: DFLT_THMIN = tenth
    real(kind=8), parameter :: DFLT_THMAX = half
    real(kind=8), parameter :: DFLT_ETA_FIXED = tenth
    integer, parameter :: DFLT_PRLVL = 0
    integer, parameter :: STDOUT = 6

    real(kind=8) :: choice1_exp, choice2_exp, choice2_coef
    real(kind=8) :: eta_cutoff, etamax
    real(kind=8) :: thmin, thmax, etafixed

    integer :: iplvl, ipunit
    common /nitparam/ choice1_exp, choice2_exp, choice2_coef, eta_cutoff, etamax, thmin, thmax, etafixed
    common /nitprint/ iplvl, ipunit

    integer :: instep, newstep, krystat
    real(kind=8) :: avrate, fcurnrm
    common /nitinfo/ avrate, fcurnrm, instep, newstep, krystat
end module nitparam_mod