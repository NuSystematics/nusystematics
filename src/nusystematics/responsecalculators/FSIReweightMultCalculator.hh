#ifndef nusystematics_RESPONSE_CALCULATORS_FSIReweightMultCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_FSIReweightMultCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"
#include "TMath.h"

#include <map>
#include <variant>
#include <string>

NEW_SYSTTOOLS_EXCEPT(invalid_FSIMult_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_FSIMult_FILEPATH);
using namespace std;

// ============================================================
// Structs (from FSIReweight.h)
// ============================================================

struct LinearParams {
    double slope, intercept;
};

struct GammaParams {
    double slope, intercept, floor, exp;
};

struct ModelParamsLinear {
    LinearParams hA, hN, INCL, G4;
};

struct ModelParamsGamma {
    GammaParams hA, hN, INCL, G4;
};

enum class ProbeType { Pion, Nucleon };

struct NucleonTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsGamma  SumGamma;
};

struct PionTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsLinear SumMean;
    ModelParamsLinear SumStd;
};

struct ParticleTargetParams {
    ProbeType probeType;
    std::variant<PionTargetParams, NucleonTargetParams> params;
};

// Keyed by [PDG code][target PDG code]
inline std::map<int, std::map<int, ParticleTargetParams>> fitParams;

// ============================================================
// Fit parameter table (from FSIReweight.h)
// ============================================================

inline void makeAllParams() {

    // PDG 2112 (neutron)
    fitParams[2112][1000060120] = {  // C12
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0458, -0.4058},  // hA
                {0.0098, -0.4984},  // hN
                {-0.0605, -0.2798},  // INCL
                {0.0709, -0.4615}  // G4
            },
            // DifStd
            {
                {0.0131, 2.0646},  // hA
                {-0.0168, 1.5610},  // hN
                {-0.0610, 0.9932},  // INCL
                {-0.0460, 1.1354}  // G4
            },
            // SumGamma
            {
                {-0.0004, 0.0176, 0.1985, 0.1000},  // hA
                {-1.3001, 1.3636, 0.0000, 0.1001},  // hN
                {-1.7230, 1.1991, 0.0000, 0.1373},  // INCL
                {-5.9100, 0.3582, 0.2882, 3.0541}  // G4
            }
        }
    };
    fitParams[2112][1000080160] = {  // O16
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0138, -0.4352},  // hA
                {0.0292, -0.5224},  // hN
                {-0.0759, -0.1609},  // INCL
                {0.0528, -0.3374}  // G4
            },
            // DifStd
            {
                {-0.0045, 2.1463},  // hA
                {-0.0204, 1.7025},  // hN
                {-0.0538, 1.1200},  // INCL
                {-0.0425, 1.2290}  // G4
            },
            // SumGamma
            {
                {-0.0131, 0.0066, 0.1803, 0.1002},  // hA
                {-1.5255, 1.4372, 0.0000, 0.1053},  // hN
                {-7.6522, 9.8971, 0.2289, 0.6249},  // INCL
                {-6.3157, 0.4936, 0.2881, 2.7992}  // G4
            }
        }
    };
    fitParams[2112][1000180400] = {  // Ar40
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0498, -1.4574},  // hA
                {-0.1262, -1.0222},  // hN
                {-0.2771, -0.9959},  // INCL
                {-0.3411, -1.1834}  // G4
            },
            // DifStd
            {
                {-0.0045, 2.3063},  // hA
                {0.1143, 2.1339},  // hN
                {0.1603, 1.4664},  // INCL
                {0.1468, 1.5241}  // G4
            },
            // SumGamma
            {
                {-1.7379, 0.6495, 0.0000, 0.1003},  // hA
                {-6.2676, 6.9147, 0.1467, 0.4568},  // hN
                {-6.6095, 6.2220, 0.1407, 0.6269},  // INCL
                {-3.3357, 0.8267, 0.0913, 1.3471}  // G4
            }
        }
    };
    fitParams[2112][1000260560] = {  // Fe56
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0384, -1.1646},  // hA
                {-0.1299, -0.9215},  // hN
                {-0.3301, -0.1915},  // INCL
                {-0.3707, -0.1054}  // G4
            },
            // DifStd
            {
                {0.0084, 2.4672},  // hA
                {0.2016, 2.2740},  // hN
                {0.2428, 1.4377},  // INCL
                {0.1860, 1.5673}  // G4
            },
            // SumGamma
            {
                {-1.7661, 0.5261, 0.0000, 0.1838},  // hA
                {-6.0021, 0.1411, 0.1204, 3.3621},  // hN
                {-4.8448, 0.3332, 0.0985, 2.1657},  // INCL
                {-3.2130, 0.6480, 0.0569, 1.5080}  // G4
            }
        }
    };

    // PDG 2212 (proton)
    fitParams[2212][1000060120] = {  // C12
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {0.0181, 1.6195},  // hA
                {-0.0048, 1.4827},  // hN
                {0.0520, 1.4161},  // INCL
                {-0.0790, 1.6617}  // G4
            },
            // DifStd
            {
                {0.0265, 2.0748},  // hA
                {-0.0193, 1.5583},  // hN
                {-0.0633, 0.9864},  // INCL
                {-0.0515, 1.1469}  // G4
            },
            // SumGamma
            {
                {-0.0021, 0.1255, 0.0917, 0.1041},  // hA
                {-1.5328, 1.6431, 0.0000, 0.1021},  // hN
                {-1.7220, 1.2051, 0.0000, 0.1248},  // INCL
                {-6.1557, 7.3755, 0.2892, 0.6325}  // G4
            }
        }
    };
    fitParams[2212][1000080160] = {  // O16
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0034, 1.6298},  // hA
                {0.0010, 1.4747},  // hN
                {0.0615, 1.4178},  // INCL
                {-0.1070, 1.7680}  // G4
            },
            // DifStd
            {
                {-0.0059, 2.1764},  // hA
                {-0.0112, 1.6833},  // hN
                {-0.0580, 1.1150},  // INCL
                {-0.0348, 1.2180}  // G4
            },
            // SumGamma
            {
                {-0.0001, 0.0100, 0.1771, 0.1000},  // hA
                {-1.6654, 1.5777, 0.0000, 0.1215},  // hN
                {-6.8995, 4.0652, 0.2352, 0.7567},  // INCL
                {-5.6316, 0.4995, 0.2827, 2.7168}  // G4
            }
        }
    };
    fitParams[2212][1000180400] = {  // Ar40
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0069, 0.4224},  // hA
                {-0.1789, 0.9682},  // hN
                {-0.1022, 0.3722},  // INCL
                {-0.4763, 0.7756}  // G4
            },
            // DifStd
            {
                {-0.0101, 2.3873},  // hA
                {0.1621, 2.1070},  // hN
                {0.1980, 1.4204},  // INCL
                {0.0903, 1.5971}  // G4
            },
            // SumGamma
            {
                {-1.8203, 0.7022, 0.0000, 0.1014},  // hA
                {-6.1361, 5.6958, 0.1339, 0.5170},  // hN
                {-6.2255, 4.0756, 0.1401, 0.6827},  // INCL
                {-3.4799, 0.9218, 0.0915, 1.2652}  // G4
            }
        }
    };
    fitParams[2212][1000260560] = {  // Fe56
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0012, 0.7096},  // hA
                {-0.2181, 1.0861},  // hN
                {-0.1576, 1.1421},  // INCL
                {-0.5327, 1.8451}  // G4
            },
            // DifStd
            {
                {0.0088, 2.5039},  // hA
                {0.2508, 2.2552},  // hN
                {0.2552, 1.4138},  // INCL
                {0.1849, 1.5511}  // G4
            },
            // SumGamma
            {
                {-1.7895, 0.5399, 0.0000, 0.1887},  // hA
                {-6.4229, 6.6016, 0.1044, 0.4978},  // hN
                {-4.8679, 0.3490, 0.0993, 2.0361},  // INCL
                {-3.2558, 0.7821, 0.0525, 1.3369}  // G4
            }
        }
    };

    // PDG -211 (piminus)
    fitParams[-211][1000060120] = {  // C12
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.7609, -1.8365},  // hA
                {0.1833, -1.3554},  // hN
                {-0.1156, -1.2158},  // INCL
                {0.0579, -1.0079}  // G4
            },
            // DifStd
            {
                {0.2195, 3.4732},  // hA
                {-0.0569, 1.3899},  // hN
                {-0.1274, 0.9037},  // INCL
                {0.1882, 0.9306}  // G4
            },
            // SumMean
            {
                {0.6658, 4.2362},  // hA
                {2.3593, 4.8110},  // hN
                {1.6460, 3.8559},  // INCL
                {0.4414, 2.0446}  // G4
            },
            // SumStd
            {
                {0.6186, 0.8173},  // hA
                {0.3440, 2.7930},  // hN
                {13.1251, -0.9103},  // INCL
                {2.7875, 2.9880}  // G4
            }
        }
    };
    fitParams[-211][1000080160] = {  // O16
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.7501, -1.9629},  // hA
                {0.0008, -1.1685},  // hN
                {-0.2474, -0.9797},  // INCL
                {0.0322, -0.9268}  // G4
            },
            // DifStd
            {
                {0.2652, 3.4930},  // hA
                {-0.0121, 1.5923},  // hN
                {-0.1428, 1.0599},  // INCL
                {0.2114, 1.0149}  // G4
            },
            // SumMean
            {
                {1.0830, 4.3215},  // hA
                {5.3797, 3.9193},  // hN
                {-0.1347, 2.7288},  // INCL
                {3.9052, 0.3812}  // G4
            },
            // SumStd
            {
                {0.9601, 0.9015},  // hA
                {0.7569, 2.9448},  // hN
                {9.5732, 3.5767},  // INCL
                {3.6217, 2.6746}  // G4
            }
        }
    };
    fitParams[-211][1000180400] = {  // Ar40
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.6508, -1.8136},  // hA
                {-0.2982, -1.8298},  // hN
                {-0.8457, -1.6877},  // INCL
                {-1.3908, -1.3807}  // G4
            },
            // DifStd
            {
                {0.4772, 3.1882},  // hA
                {0.4976, 2.1689},  // hN
                {0.5290, 1.3862},  // INCL
                {0.2569, 1.6061}  // G4
            },
            // SumMean
            {
                {3.4979, 5.0356},  // hA
                {7.3865, 6.2422},  // hN
                {8.5453, 2.4680},  // INCL
                {11.7655, 0.3049}  // G4
            },
            // SumStd
            {
                {2.3536, 1.5691},  // hA
                {1.2768, 4.3182},  // hN
                {1.8897, 4.6942},  // INCL
                {2.3445, 3.6421}  // G4
            }
        }
    };
    fitParams[-211][1000260560] = {  // Fe56
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.1089, -1.9992},  // hA
                {-0.2659, -1.6327},  // hN
                {-0.9714, -0.7995},  // INCL
                {-1.2066, -0.3517}  // G4
            },
            // DifStd
            {
                {0.5256, 3.2859},  // hA
                {0.6522, 2.3551},  // hN
                {0.6318, 1.3908},  // INCL
                {0.4715, 1.4936}  // G4
            },
            // SumMean
            {
                {4.8382, 5.8804},  // hA
                {8.8250, 6.4708},  // hN
                {10.2677, 3.3301},  // INCL
                {14.4090, -0.7734}  // G4
            },
            // SumStd
            {
                {2.5376, 2.2535},  // hA
                {1.2986, 4.4985},  // hN
                {2.4414, 4.1590},  // INCL
                {2.9453, 3.4111}  // G4
            }
        }
    };

    // PDG 111 (pizero)
    fitParams[111][1000060120] = {  // C12
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.9052, -0.1317},  // hA
                {0.0309, 0.4817},  // hN
                {0.0060, 0.5442},  // INCL
                {-0.1184, 0.6118}  // G4
            },
            // DifStd
            {
                {-0.1003, 3.7673},  // hA
                {0.0066, 1.3559},  // hN
                {-0.1776, 0.9585},  // INCL
                {-0.1778, 1.3106}  // G4
            },
            // SumMean
            {
                {0.8463, 4.1349},  // hA
                {1.7051, 5.4571},  // hN
                {0.5703, 5.0099},  // INCL
                {0.0003, 2.0002}  // G4
            },
            // SumStd
            {
                {0.6222, 0.8041},  // hA
                {-0.0207, 3.1001},  // hN
                {4.0335, 3.7402},  // INCL
                {3.8872, 2.0507}  // G4
            }
        }
    };
    fitParams[111][1000080160] = {  // O16
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {3.1404, -0.3688},  // hA
                {0.0321, 0.4539},  // hN
                {-0.0558, 0.6370},  // INCL
                {-0.1467, 0.7228}  // G4
            },
            // DifStd
            {
                {0.1547, 3.6580},  // hA
                {-0.0448, 1.6472},  // hN
                {-0.1356, 1.1183},  // INCL
                {-0.2467, 1.4226}  // G4
            },
            // SumMean
            {
                {1.3092, 4.1830},  // hA
                {3.1107, 5.4940},  // hN
                {2.4582, 1.3847},  // INCL
                {0.9516, 2.1067}  // G4
            },
            // SumStd
            {
                {0.9129, 0.9172},  // hA
                {0.6364, 2.9601},  // hN
                {21.0099, -0.2930},  // INCL
                {2.9113, 2.6140}  // G4
            }
        }
    };
    fitParams[111][1000180400] = {  // Ar40
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.0891, -0.5501},  // hA
                {-0.6140, -0.0622},  // hN
                {-0.6430, -0.3426},  // INCL
                {-1.2156, 0.0380}  // G4
            },
            // DifStd
            {
                {0.4758, 3.2005},  // hA
                {0.4925, 2.1847},  // hN
                {0.4928, 1.4221},  // INCL
                {0.0108, 1.7389}  // G4
            },
            // SumMean
            {
                {3.4947, 4.9809},  // hA
                {6.8813, 6.9028},  // hN
                {8.4174, 2.7382},  // INCL
                {15.0807, -3.5560}  // G4
            },
            // SumStd
            {
                {2.2506, 1.6529},  // hA
                {1.6484, 4.0574},  // hN
                {1.7047, 4.7834},  // INCL
                {2.4725, 4.1813}  // G4
            }
        }
    };
    fitParams[111][1000260560] = {  // Fe56
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.4231, -0.5136},  // hA
                {-0.5220, -0.0200},  // hN
                {-0.7077, 0.4933},  // INCL
                {-0.9195, 0.9232}  // G4
            },
            // DifStd
            {
                {0.3810, 3.4108},  // hA
                {0.5532, 2.4634},  // hN
                {0.6658, 1.3784},  // INCL
                {0.2116, 1.6657}  // G4
            },
            // SumMean
            {
                {4.9143, 5.6408},  // hA
                {8.0224, 7.3538},  // hN
                {10.7279, 3.1948},  // INCL
                {15.1094, -2.1593}  // G4
            },
            // SumStd
            {
                {2.6421, 2.2420},  // hA
                {1.7590, 4.3464},  // hN
                {2.9992, 3.8938},  // INCL
                {3.3541, 3.4803}  // G4
            }
        }
    };

    // PDG 211 (piplus)
    fitParams[211][1000060120] = {  // C12
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.9002, 1.7511},  // hA
                {-0.1257, 2.3318},  // hN
                {0.0878, 2.3290},  // INCL
                {-0.2527, 2.2158}  // G4
            },
            // DifStd
            {
                {-0.1002, 3.7254},  // hA
                {-0.0056, 1.3204},  // hN
                {-0.1019, 0.8634},  // INCL
                {0.1635, 0.9133}  // G4
            },
            // SumMean
            {
                {1.0699, 4.1030},  // hA
                {2.0520, 4.8404},  // hN
                {-1.9430, 6.1858},  // INCL
                {0.2475, 1.9571}  // G4
            },
            // SumStd
            {
                {0.6212, 0.8323},  // hA
                {0.3243, 3.2116},  // hN
                {6.9417, 2.1470},  // INCL
                {2.1653, 3.2065}  // G4
            }
        }
    };
    fitParams[211][1000080160] = {  // O16
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.8132, 1.9100},  // hA
                {-0.0463, 2.1889},  // hN
                {0.1690, 2.2358},  // INCL
                {-0.1429, 2.2621}  // G4
            },
            // DifStd
            {
                {-0.1951, 3.9639},  // hA
                {0.0106, 1.5698},  // hN
                {-0.1041, 1.0261},  // INCL
                {0.1065, 1.0518}  // G4
            },
            // SumMean
            {
                {1.4981, 4.2037},  // hA
                {3.2590, 5.0669},  // hN
                {1.7280, 1.4570},  // INCL
                {1.6754, 1.6791}  // G4
            },
            // SumStd
            {
                {0.9073, 0.9469},  // hA
                {0.7763, 3.4013},  // hN
                {13.2244, 2.3207},  // INCL
                {1.5898, 3.6865}  // G4
            }
        }
    };
    fitParams[211][1000180400] = {  // Ar40
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.2903, 0.8905},  // hA
                {-0.6633, 1.4800},  // hN
                {-0.4911, 1.0835},  // INCL
                {-1.0216, 1.3981}  // G4
            },
            // DifStd
            {
                {0.2987, 3.2604},  // hA
                {0.3182, 2.2811},  // hN
                {0.6505, 1.2921},  // INCL
                {0.2593, 1.3767}  // G4
            },
            // SumMean
            {
                {3.7611, 4.9256},  // hA
                {6.6766, 7.1275},  // hN
                {8.2516, 2.7582},  // INCL
                {12.4768, -0.6278}  // G4
            },
            // SumStd
            {
                {2.1374, 1.6616},  // hA
                {1.2266, 4.5166},  // hN
                {2.3016, 4.5927},  // INCL
                {2.1749, 3.9374}  // G4
            }
        }
    };
    fitParams[211][1000260560] = {  // Fe56
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.8171, 0.9453},  // hA
                {-0.7655, 1.6732},  // hN
                {-0.6432, 1.9197},  // INCL
                {-0.6985, 2.2513}  // G4
            },
            // DifStd
            {
                {0.3591, 3.3994},  // hA
                {0.6809, 2.3669},  // hN
                {0.6811, 1.3115},  // INCL
                {0.3484, 1.4438}  // G4
            },
            // SumMean
            {
                {5.0295, 5.7037},  // hA
                {9.1261, 6.7006},  // hN
                {11.1969, 2.9899},  // INCL
                {13.6643, -0.3371}  // G4
            },
            // SumStd
            {
                {2.5278, 2.2597},  // hA
                {1.7786, 4.6312},  // hN
                {3.0413, 3.9273},  // INCL
                {3.2988, 3.0931}  // G4
            }
        }
    };
}

// ============================================================
// Accessors (from FSIReweight.h)
// ============================================================

inline const auto& getPion(int pdg, int target) {
    return std::get<PionTargetParams>(fitParams.at(pdg).at(target).params);
}

inline const auto& getNucleon(int pdg, int target) {
    return std::get<NucleonTargetParams>(fitParams.at(pdg).at(target).params);
}

// ============================================================
// Visible-energy reweight from template histograms (from FSIReweight.h)
// ============================================================

inline double GetVisEReweight(TH2D* hist_nom, TH2D* hist_alt, double KEini, double Ebias)
{
    if (!hist_nom || !hist_alt) return 1.;

    int idx_KEini = hist_nom->GetXaxis()->FindBin(KEini);
    int idx_Ebias = hist_nom->GetYaxis()->FindBin(Ebias);
    int nY        = hist_nom->GetNbinsY();

    auto findFirstY = [&](TH2D* h, int ix) {
        for (int iy = 1; iy <= nY; iy++)
            if (h->GetBinContent(ix, iy) > 0.00) return iy;
        return 1;
    };
    auto findLastY = [&](TH2D* h, int ix) {
        for (int iy = nY; iy >= 1; iy--)
            if (h->GetBinContent(ix, iy) > 0.00) return iy;
        return nY;
    };

    int yMin = findFirstY(hist_nom, idx_KEini);
    int yMax = findLastY (hist_nom, idx_KEini);

    double norm_nom = hist_nom->Integral(idx_KEini, idx_KEini, yMin, yMax + 1);
    double norm_alt = hist_alt->Integral(idx_KEini, idx_KEini, yMin, yMax + 1);

    double weight_nom = hist_nom->GetBinContent(idx_KEini, idx_Ebias) / norm_nom;
    double weight_alt = hist_alt->GetBinContent(idx_KEini, idx_Ebias) / norm_alt;

    return weight_alt / weight_nom;
}

// ============================================================
// Multiplicity / charge-asymmetry reweight (from FSIReweight.h)
// ============================================================

inline void computeMultDiffWeights(int pdg, int target, int mult, int diff,
                                   double KEini, int max,
                                   const std::string& reweightModel,
                                   bool isPion,
                                   double& hAEst, double& AltEst)
{

   hAEst   = 1.0;
    AltEst  = 1.0;
    
    int A = (target/10) % 1000;
    if (target!=1000060120 && target!=1000180400 && target!=1000080160) target=1000180400;
    max=A;
    if (isPion) {
        const auto& tp = getPion(pdg, target);

const LinearParams& altSumMean = (reweightModel == "hN2018") ? tp.SumMean.hN   :
                                 (reweightModel == "hA2018") ? tp.SumMean.hA   :
                                 (reweightModel == "Geant4")  ? tp.SumMean.G4   :
                                                                tp.SumMean.INCL;
const LinearParams& altSumStd  = (reweightModel == "hN2018") ? tp.SumStd.hN    :
                                 (reweightModel == "hA2018") ? tp.SumStd.hA    :
                                 (reweightModel == "Geant4")  ? tp.SumStd.G4    :
                                                                tp.SumStd.INCL;
const LinearParams& altDifMean = (reweightModel == "hN2018") ? tp.DifMean.hN   :
                                 (reweightModel == "hA2018") ? tp.DifMean.hA   :
                                 (reweightModel == "Geant4")  ? tp.DifMean.G4   :
                                                                tp.DifMean.INCL;
const LinearParams& altDifStd  = (reweightModel == "hN2018") ? tp.DifStd.hN    :
                                 (reweightModel == "hA2018") ? tp.DifStd.hA    :
                                 (reweightModel == "Geant4")  ? tp.DifStd.G4    :
                                                                tp.DifStd.INCL;

        double meanAltSum = altSumMean.slope * KEini + altSumMean.intercept;
        double stdAltSum  = altSumStd.slope  * KEini + altSumStd.intercept;
        double meanhASum  = tp.SumMean.hA.slope * KEini + tp.SumMean.hA.intercept;
        double stdhASum   = tp.SumStd.hA.slope  * KEini + tp.SumStd.hA.intercept;
        //if (reweightModel=="Geant4") std::cout<<meanAltSum<<","<<stdAltSum<<std::endl;
        double meanAltDif = altDifMean.slope * KEini + altDifMean.intercept;
        double stdAltDif  = altDifStd.slope  * KEini + altDifStd.intercept;
        double meanhADif  = tp.DifMean.hA.slope * KEini + tp.DifMean.hA.intercept;
        double stdhADif   = tp.DifStd.hA.slope  * KEini + tp.DifStd.hA.intercept;

        double P_spike = 1.14 * (0.903 - 0.001989 * A) * (1.35 - 4.767 * KEini);
      
        
        if (P_spike < 0) P_spike = 0;
        if (P_spike > 1) P_spike = 1;
        double normNom = TMath::Freq((max + 1.0 - meanhASum) / stdhASum)
                       - TMath::Freq((2.0        - meanhASum) / stdhASum);
        double normAlt = TMath::Freq((max + 1.0 - meanAltSum) / stdAltSum)
                       - TMath::Freq((2.0        - meanAltSum) / stdAltSum);


        double hAEstSum  = (TMath::Freq((mult + 1.0 - meanhASum) / stdhASum)
             - TMath::Freq((mult       - meanhASum) / stdhASum))/normNom;

        double AltEstSum = (TMath::Freq((mult + 1.0 - meanAltSum) / stdAltSum)
              - TMath::Freq((mult - meanAltSum) / stdAltSum))/ normAlt;

        double normNomDiff = TMath::Freq((max/2 + 1.0 - meanhADif) / stdhADif)
                       - TMath::Freq((-max/2        - meanhADif) / stdhADif);
        double normAltDiff = TMath::Freq((max + 1.0 - meanAltDif) / stdAltDif)
                       - TMath::Freq((-max/2        - meanAltDif) / stdAltDif);
        double hAEstDif  = (TMath::Freq((diff + 1.0 - meanhADif) / stdhADif)
             - TMath::Freq((diff       - meanhADif) / stdhADif))/normNomDiff;
        double AltEstDif = (TMath::Freq((diff + 1.0 - meanAltDif) / stdAltDif)
              - TMath::Freq((diff - meanAltDif) / stdAltDif))/ normAltDiff;

        double F_2 = (TMath::Freq((3.0 - meanhASum) / stdhASum)
                    - TMath::Freq((2.0 - meanhASum) / stdhASum)) / normNom;
        double f_2 = P_spike + (1.0 - P_spike) * F_2;
        double qd_rw;
        if (mult == 2) {
            qd_rw = F_2/f_2;

        } else {
            qd_rw = (1.0 - F_2) / (1.0 - f_2);
        }
        if (KEini<0.35){ AltEstDif=1; hAEstDif=1;}
        hAEst  = hAEstSum*hAEstDif;
        AltEst = AltEstSum * qd_rw*AltEstDif;
        if (AltEst/hAEst<0.01){ AltEst=1; hAEst=100;}
        if (AltEst/hAEst>100){ AltEst=100; hAEst=1;}

    } else {
        const auto& tp = getNucleon(pdg, target);

const GammaParams&  altGamma   = (reweightModel == "hN2018") ? tp.SumGamma.hN :
                                 (reweightModel == "hA2018") ? tp.SumGamma.hA :
                                 (reweightModel == "Geant4")  ? tp.SumGamma.G4 :
                                                                tp.SumGamma.INCL;
const LinearParams& altDifMean = (reweightModel == "hN2018") ? tp.DifMean.hN  :
                                 (reweightModel == "hA2018") ? tp.DifMean.hA  :
                                 (reweightModel == "Geant4")  ? tp.DifMean.G4  :
                                                                tp.DifMean.INCL;
const LinearParams& altDifStd  = (reweightModel == "hN2018") ? tp.DifStd.hN   :
                                 (reweightModel == "hA2018") ? tp.DifStd.hA   :
                                 (reweightModel == "Geant4")  ? tp.DifStd.G4   :
                                                                tp.DifStd.INCL;

        double meanAltDif = altDifMean.slope * KEini + altDifMean.intercept;
        double stdAltDif  = altDifStd.slope  * KEini + altDifStd.intercept;
        double meanhADif  = tp.DifMean.hA.slope * KEini + tp.DifMean.hA.intercept;
        double stdhADif   = tp.DifStd.hA.slope  * KEini + tp.DifStd.hA.intercept;

        double normNomDiff = TMath::Freq((max/2 + 1.0 - meanhADif) / stdhADif)
                       - TMath::Freq((-max/2        - meanhADif) / stdhADif);
        double normAltDiff = TMath::Freq((max + 1.0 - meanAltDif) / stdAltDif)
                       - TMath::Freq((-max/2        - meanAltDif) / stdAltDif);
        double hAEstDif  = (TMath::Freq((diff + 1.0 - meanhADif) / stdhADif)
             - TMath::Freq((diff       - meanhADif) / stdhADif))/normNomDiff;
        double AltEstDif = (TMath::Freq((diff + 1.0 - meanAltDif) / stdAltDif)
              - TMath::Freq((diff - meanAltDif) / stdAltDif))/ normAltDiff;


        double gammaAlt = altGamma.intercept          * TMath::Exp(TMath::Power(KEini, altGamma.exp)           * altGamma.slope)           + altGamma.floor;
        double gammaHA  = tp.SumGamma.hA.intercept   * TMath::Exp(TMath::Power(KEini, tp.SumGamma.hA.exp)    * tp.SumGamma.hA.slope)    + tp.SumGamma.hA.floor;

    double normNomSum = 1.0 - TMath::Exp(-gammaHA  * (max - 2));
    double normAltSum = 1.0 - TMath::Exp(-gammaAlt * (max - 2));
    
    double hAEstSum  = (TMath::Exp(-gammaHA  * (mult - 3)) - TMath::Exp(-gammaHA  * (mult + 1 - 3))) / normNomSum;
    double AltEstSum = (TMath::Exp(-gammaAlt * (mult - 3)) - TMath::Exp(-gammaAlt * (mult + 1 - 3))) / normAltSum;

        if (KEini<0.35){ AltEstDif=1; hAEstDif=1;}
        hAEst  = hAEstSum * hAEstDif;
        AltEst = AltEstSum * AltEstDif;
        if (AltEst/hAEst<0.1){ AltEst=1; hAEst=10;}
        if (AltEst/hAEst>10){ AltEst=10; hAEst=1;}
    }
}

namespace nusyst {

  class FSIReweightMultCalculator {

    enum ENuRange {
      LowE = 0,
      HighE = 20,
    };

  protected:

    TH2D *hist_nom_protonPlus;
    TH2D *hist_alt_protonPlus;
    TH2D *hist_nom_neutron;
    TH2D *hist_alt_neutron;
    TH2D *hist_nom_piPlus;
    TH2D *hist_alt_piPlus;
    TH2D *hist_nom_pi0;
    TH2D *hist_alt_pi0;
    TH2D *hist_nom_piMinus;
    TH2D *hist_alt_piMinus;
    TH3D *hist_nom_2p;
    TH3D *hist_alt_2p;

  public:

    FSIReweightMultCalculator(fhicl::ParameterSet const &InputManifest) {
      makeAllParams();
      LoadInputHistograms(InputManifest);
    }
    ~FSIReweightMultCalculator() {}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetFSIMultReweight(double KEini, double Ebias, int mult, int diff, double parameter_value, int parpdg, int target, int max);
    double GetFSIMultReweight_2par(double KEini_0, double KEini_1, double Ebias, double parameter_value, int parpdg);

    std::string GetCalculatorName() const { return "FSIReweightMultCalculator"; }

  };

  inline double FSIReweightMultCalculator::GetFSIMultReweight(double KEini, double Ebias, int mult, int diff, double parameter_value, int parpdg, int target, int max) {
    // --- Visible-energy reweight from template histograms ---
    TH2D *hist_nom, *hist_alt;
    if (parpdg == 2212) {
      hist_nom = hist_nom_protonPlus;
      hist_alt = hist_alt_protonPlus;
    }
    else if (parpdg == 2112) {
      hist_nom = hist_nom_neutron;
      hist_alt = hist_alt_neutron;
    }
    else if (parpdg == 211) {
      hist_nom = hist_nom_piPlus;
      hist_alt = hist_alt_piPlus;
    }
    else if (parpdg == 111) {
      hist_nom = hist_nom_pi0;
      hist_alt = hist_alt_pi0;
    }
    else if (parpdg == -211) {
      hist_nom = hist_nom_piMinus;
      hist_alt = hist_alt_piMinus;
    }
    else {
      return 1.;
    }

    double visERW = GetVisEReweight(hist_nom, hist_alt, KEini, Ebias);

    // --- Multiplicity / charge-asymmetry reweight from analytic fits ---
    bool isPion = (parpdg == 211 || parpdg == -211 || parpdg == 111);

    // Determine the reweight model name from parameter_value for the mult/diff
    // weight: parameter_value in [0,1] interpolates between hA2018 (nom) and alt.
    // The model string carried by hist_alt determines which alt is loaded;
    // we recover it by reading the histogram title set at load time.
    // Rather than re-parsing, we apply the same linear interpolation
    // directly on hAEst (nom) and AltEst (alt) from computeMultDiffWeights,
    // using the reweight_model stored at construction (passed via fhicl).
    double weight_nom = 1.;
    double weight_alt = 1.;

    // reweight_model is embedded in the histogram names; retrieve from title
std::string reweightModel = "Geant4";
      std::string nameHist=std::string(hist_alt->GetName()).substr(0, std::string(hist_alt->GetName()).find('_'));
        if (nameHist=="HG4BertCasc") reweightModel="Geant4";
      else if (nameHist=="INCL" || nameHist=="HINCL" || nameHist=="INCL++") reweightModel="INCL++";
      else if (nameHist=="hA" || nameHist=="hA2018") reweightModel="hA2018";
      else if (nameHist=="hN" || nameHist=="hN2018") reweightModel="hN2018";
      computeMultDiffWeights(parpdg, target, mult, diff, KEini, max,
                           reweightModel, isPion,
                           weight_nom, weight_alt);

    if (weight_nom == 0.) {
      weight_nom=0.01;
    }
    if (weight_alt == 0.) {
      weight_alt = 1.;
    }

    double multDiffWeight = (weight_nom * (1. - parameter_value) + weight_alt * parameter_value) / weight_nom;

    return multDiffWeight;
  }


  inline void FSIReweightMultCalculator::LoadInputHistograms(fhicl::ParameterSet const &ps) {

    std::string const &default_root_file = ps.get<std::string>("input_file", "");

    for (fhicl::ParameterSet const &val_config :
         ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
      std::string hName      = val_config.get<std::string>("name");
      std::string input_hist = val_config.get<std::string>("input_hist");
      std::string input_file = val_config.get<std::string>("input_file", default_root_file);

      // if it does not start with "/", find it under ${nusystematics_ROOT}/data/
      if (input_file.find("/") != 0) {
        std::string tmp_NUSYSTEMATICS_ROOT = std::getenv("nusystematics_ROOT");
        if (tmp_NUSYSTEMATICS_ROOT == "") {
          throw invalid_FSIMult_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set but put relative path:" << input_file;
        }
        input_file = tmp_NUSYSTEMATICS_ROOT + "/data/" + input_file;
      }

      if (hName == "hist_nom_protonPlus") {
        hist_nom_protonPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_alt_protonPlus") {
        hist_alt_protonPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_nom_neutron") {
        hist_nom_neutron = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_alt_neutron") {
        hist_alt_neutron = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_nom_piPlus") {
        hist_nom_piPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_alt_piPlus") {
        hist_alt_piPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_nom_pi0") {
        hist_nom_pi0 = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_alt_pi0") {
        hist_alt_pi0 = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_nom_piMinus") {
        hist_nom_piMinus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if (hName == "hist_alt_piMinus") {
        hist_alt_piMinus = GetHistogram<TH2D>(input_file, input_hist);
      }
    }
  }

} // namespace nusyst

#endif
