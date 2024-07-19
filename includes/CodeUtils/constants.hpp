#ifndef INCLUDES_CODEUTILS_CONSTANTS_HPP
#define INCLUDES_CODEUTILS_CONSTANTS_HPP

#include <string>

namespace constants
{
    const double maxAtomDistanceFromResidueCenter = 5.4; // based on observations of test data as of june 2024
    const double residueDistanceOverlapCutoff     = 2.0 * maxAtomDistanceFromResidueCenter;
    const double maxCutOff                        = 2.0; // roughly the radius of the common largest atoms
    const double overlapTolerance                 = 0.6; // chimera default
    const double clashWeightBase                  = 0.2;
    const double DEFAULT_ANGLE =
        116.8; // Glycam06 J Comput Chem. 2008 Mar; 29(4): 622–655. Used for C-O-C in glycosidic linkages.
    const double CARBON_SURFACE_AREA  = 36.31681103;
    const double MAX_RESIDUE_DIAMETER = 9.0;
    const double PI_RADIAN            = 3.141592653589793;
    const double PI_DEGREE            = 180.0;
    const double dNotSet              = 123456789.0;
    const int iNotSet                 = -123456;
    const std::string sNotSet         = "?";
    const double dSulfurCutoff        = 3.0;
    const double defaultBondLength    = 1.4;

    inline double degree2Radian(double d)
    {
        return d / PI_DEGREE * PI_RADIAN;
    }

    const double DIST_EPSILON =
        0.000001; // \todo determine where/if this is used, and if the following would be better.
    // const double Float_Machine_Epsilon = std::numeric_limits<float>::epsilon( );
    // const double Double_Machine_Epsilon = std::numeric_limits<double>::epsilon( );
} // namespace constants
#endif
