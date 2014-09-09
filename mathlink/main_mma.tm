
:Begin:
:Function:      time_evo
:Pattern:       trgEvolution[eta_List, ps_List, Oeta_List, Ok_List, OmegaBulk_List, opts_List]
:Arguments:     {eta, ps, Oeta, Ok, OmegaBulk, opts}
:ArgumentTypes: {RealList, RealList, RealList, RealList, RealList, RealList}
:ReturnType:    Manual
:End:

:Begin:
:Function:      getA
:Pattern:       A[k_Real]
:Arguments:     {k}
:ArgumentTypes: {Real64}
:ReturnType:    Manual
:End:

:Begin:
:Function:      init_A
:Pattern:       InitA[ps_List, options_List]
:Arguments:     {ps, options}
:ArgumentTypes: {RealList, IntegerList}
:ReturnType:    Manual
:End:

:Begin:
:Function:      clean_up_A
:Pattern:       CleanUpA[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Manual
:End:

:Begin:
:Function:      f_ode
:Pattern:       fOde[eta_Real, X_List]
:Arguments:     {eta, X}
:ArgumentTypes: {Real64, RealList}
:ReturnType:    Manual
:End:

:Begin:
:Function:      init_ode
:Pattern:       InitOde[Oeta_List, Ok_List, Omega_List, k_List, options_List]
:Arguments:     {O_eta, O_k, Omega, k, options}
:ArgumentTypes: {RealList, RealList, RealList, RealList, RealList}
:ReturnType:    Manual
:End:

:Begin:
:Function:      clean_up_ode
:Pattern:       CleanUpOde[]
:Arguments:     {}
:ArgumentTypes: {}
:ReturnType:    Manual
:End:

