### PhysicalConstants.jl
#
# Authors: Kangkang Zhao
# 
###
module PhysicalConstants
export PhysicalConstants

const Tau60Cobalt = 1.665e8 / (24 * 3600)  # in days
const ESumPeak60Cobalt = 2505.692  # 60Co peak energy
const Qvalue0nbb130Te = 2527.518   # Q-value of 130Te 0νββ
const E208Tlpeak = 2614.511 # keV - this is the gamma line from 208 Tl decay

# define global constants
const N_A = 6.022E23             # Avogadro constant 1/mol
const mass_TeO2 = 159.6       # g/mol
const mass_130Te = 129.9062244E-3 # kg/mol
const Abundance_130Te = 0.341668 # natural isotopic abundance of 130Te

# to find a rough estimate of the signal rate at the 130 Te Qvalue for 0nbb
const Emin0nbbPeak20keV = 2517.5; # keV - if the event has energy above this value it is counted as signal
const Emax0nbbPeak20keV = 2537.5; # keV - if the event has energy below this value it is counted as signal
# proxy for 4 sigma region large 40 keV, around the Qvalue for 0nbb of 130 Te -- FIXME -- this is asymmetric as a region                                                                           
#to find a rough estimate of the rate outside the central region of  130 Te Qvalue for 0nbb - this is an exclusion region 
const Emin0nbbPeak40keV = 2497.5; # keV - if the event has energy below this value it is counted as bkg
const Emax0nbbPeak40keV = 2537.5; # keV - if the event has energy above this value it is counted as bkg
#correct Cobalt mean peak -- decide what to do 




const me_keV = 510.99895000 # PDG
const m_76 = 75.92E-3 # kg/mol




# //fixed quantities - they will not be changed in the 0nbb analysis 
# constexpr double AvogadroConstant = 6.022e23; // 1/mol - Avogadro constant
const NaturalIsotopicAbuTe130::Float64   = 0.341668; # CUORE 130Te Isotopic abundance (per number)
const NaturalIsotopicAbuTe130Err::Float64 = 0.000016;
# // from PhysicalConstants.hh   
# //source : M. A. Fehr, M. Rehkamper, A. N. Halliday, Int. J. Mass Spectrom. 232, 83-94 (2004)
# constexpr double gTeO2Mass = 159.6; // g/mol - Average molecular mass of a TeO2 molecule
# // from PhysicalConstants.hh     
# constexpr double Tau60Cobalt = 2777.59; // 60Cobalt LIFETIME (not half life!) in days
# //constexpr double ESumPeak60Cobalt = 2505.72;// From literature
# constexpr double ESumPeak60Cobalt = 2506.84;// From TAUP19 -> PRL19 (bkg only + 60Co s,f,LS all datasets combined - high stat)
# //constexpr double ESumPeak60Cobalt = 2504.73983704; // Best fit from TAUP official analysis. Alternative is CUORE-0 shift: 2505.692+1.9 keV
# //constexpr double ESumPeak60CobaltErr = 1.3; // Uncertainty on 60Co mean value from 2017 PRL analysis
# constexpr double ESumPeak60CobaltErr = 0.56; // Uncertainty on 60Co mean value from TAUP19 analysis -> PRL analysis (same fit as for ESumPeak60Cobalt)
# constexpr double SumPeak60CobaltRateErr = 0.35; // Relative uncertainty (in % on 60Co rate from 2017 PRL analysis
# constexpr double E208Tlpeak = 2614.511; //keV - this is the gamma line from 208 Tl decay
const Qββ130Te::Float64    = 2527.518; # keV - Qvalue of the neutrinoless double beta decay of 130Te
const Qββ130TeErr::Float64 = 0.013; # keV
# constexpr double E60Co_1173 = 1173.228; // keV
# constexpr double E60Co_1332 = 1332.492; // keV
# constexpr double E40K       = 1460.820; // keV
# constexpr double E214Bi     = 1764.491; // keV
# constexpr double E54Mn      = 834.848;  // keV

# // parameters or quantities we might change in the 0nbb fit - decide how to treat them ?
# // proxy for 2 sigma region large 20 keV, centered at the Qvalue for 0nbb of 130 Te
# //to find a rough estimate of the signal rate at the 130 Te Qvalue for 0nbb
# constexpr double Emin0nbbPeak20keV = 2517.5; // keV - if the event has energy above this value it is counted as signal
# constexpr double Emax0nbbPeak20keV = 2537.5; // keV - if the event has energy below this value it is conted as signal
# // proxy for 4 sigma region large 40 keV, around the Qvalue for 0nbb of 130 Te -- FIXME -- this is asymmetric as a region                                                                           
# //to find a rough estimate of the rate outside the central region of  130 Te Qvalue for 0nbb - this is an exclusion region 
# constexpr double Emin0nbbPeak40keV = 2497.5; // keV - if the event has energy below this value it is counted as bkg
# constexpr double Emax0nbbPeak40keV = 2537.5; //keV - if the event has energy above this value it is counted as bkg
# //correct Cobalt mean peak -- decide what to do 

# constexpr double Emean130TeXrays = 28.0271; //keV - weighted average of ka1..ka3, kb1..kb5 lines, ref: http://nucleardata.nuclear.lu.se/toi/xray.asp?act=list&el=Te

const MCEfficiency::Float64 = 0.88345; # 0nbb events containment efficiency from Monte Carlo simulations -- updated for the PRL19 analysis
const MCEfficiencyErr::Float64 = 0.00085; # 0nbb events containment efficiency error from Monte Carlo simulations -- updated for the PRL19 analysis 





end

