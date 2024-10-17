Two methods developed for calculation of the upper limit on the number of signal events:
1. Asymptotic method based on likelihood scan (Sen Jia)
2. CLs method implemented with a counting model parametrized as FCs method.

    CLs method flow:

    1. `SetXstatesParams` function returns properties (mass, width, LHCb_error) of a given X state. It includes the library containing the properties of X stated reported by LHCb.
    
    *Note: Only `LHCb_error` parameters are actually calculated before being included in the library. They are calculated by enabling the `LHCb_syst_scan` flag to `true` (while `syst_scan` is set to `false`). The code then performs the sequence of fits with mass and width parameters being deviated according to reported by LHCb errors to provide a final error, which should be manually written in the `SetXstatesParams` library.
    
    2. The `perform_fit` function then performs a Breit-Wigner+Bernstein-Polynomial fit with the parameters of the Breit-Wigner set to those returned by the `SetXstatesParams` function. Fit results are then used to calculate the number of background (b) and total (nobs) numbers of events in the 3 or 5 sigma region (defined by using CL90 or CL95 scale factor).

    3. Counting model (`CountUL`) then calculated the CLs upper limit for the determined `b` and `nobs` parameters. THhe UL is then printed out.

    4. If the `syst_scan` is set to `false`, is varied by the total systematic uncertainty value. CLs calculation is then performed for each  case, and the largest UL value is taken as a final result.

    Note: The total sysematic uncertainty is calculated as a square root of the quadrature sums of syst_error (base systematics uncertainty defined in the `run_both_methods` function), LHCb_error and stat_error (retrieved from the fit). 
