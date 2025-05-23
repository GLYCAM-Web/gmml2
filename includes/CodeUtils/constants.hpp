#ifndef INCLUDES_CODEUTILS_CONSTANTS_HPP
#define INCLUDES_CODEUTILS_CONSTANTS_HPP

#include <string>

namespace constants
{
    const double maxAtomDistanceFromResidueCenter = 5.4; // based on observations of test data as of june 2024
    const double residueDistanceOverlapCutoff     = 2.0 * maxAtomDistanceFromResidueCenter;
    const double overlapTolerance                 = 0.26; // enough to capture 0 vdw force for H-H
    const double DEFAULT_ANGLE =
        116.8; // Glycam06 J Comput Chem. 2008 Mar; 29(4): 622–655. Used for C-O-C in glycosidic linkages.
    const double PI_RADIAN     = 3.141592653589793;
    const double PI_DEGREE     = 180.0;
    const double dNotSet       = 123456789.0;
    const int iNotSet          = -123456;
    const std::string sNotSet  = "?";
    const double dSulfurCutoff = 3.0;

    inline double toRadians(double d)
    {
        return (d / PI_DEGREE) * PI_RADIAN;
    }

    inline double toDegrees(double d)
    {
        return (d / PI_RADIAN) * PI_DEGREE;
    }
} // namespace constants
#endif
