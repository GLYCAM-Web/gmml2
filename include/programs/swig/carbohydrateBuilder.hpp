#ifndef INCLUDE_PROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP
#define INCLUDE_PROGRAMS_CARBOHYDRATEBUILDER_CARBOHYDRATEBUILDER_HPP

// This is becoming an interface to the carbohydrate class in gmml for Gems.
#include "include/carbohydrate/carbohydrate.hpp"
#include "include/carbohydrate/parameterManager.hpp"
#include "include/metadata/dihedralangledata.hpp"

#include <string>
#include <vector>

namespace gmml
{ // For specifying a specific shape to be built with GenerateSpecific3DStructure

    struct SingleRotamerInfo
    {
        std::string linkageIndex; // What Dan is calling linkageLabel. Internal index determined at C++ level and given
                                  // to frontend to track.
        std::string linkageName;  // Can be whatever the user wants it to be, default to same as index.
        std::string dihedralName; // omg / phi / psi / chi1 / chi2
        std::string selectedRotamer; // gg / tg / g- etc
        std::string numericValue;    // user entered 64 degrees. Could be a v2 feature.
    };

    typedef std::vector<SingleRotamerInfo> SingleRotamerInfoVector;

    struct DihedralOptions
    { // CONSTRUCTOR

        DihedralOptions() {}

        DihedralOptions(const std::string name, const std::vector<std::string> rotamers)
            : dihedralName_(name), rotamers_(rotamers)
        {}

        // DATA
        std::string dihedralName_;          // omg / phi / psi / chi1 / chi2
        std::vector<std::string> rotamers_; // gg / tg / g- etc
    };

    typedef std::vector<DihedralOptions> DihedralOptionsVector;

    struct LinkageOptions
    { // CONSTRUCTOR

        LinkageOptions() {}

        LinkageOptions(
            std::string name,
            std::string index,
            std::string res1,
            std::string res2,
            DihedralOptionsVector likely,
            DihedralOptionsVector possible)
            : linkageName_(name), indexOrderedLabel_(index), firstResidueNumber_(res1), secondResidueNumber_(res2),
              likelyRotamers_(likely), possibleRotamers_(possible)
        {}

        // DATA
        std::string linkageName_;
        std::string indexOrderedLabel_;
        std::string firstResidueNumber_;
        std::string secondResidueNumber_;
        DihedralOptionsVector likelyRotamers_;
        DihedralOptionsVector possibleRotamers_;
    };

    typedef std::vector<LinkageOptions> LinkageOptionsVector;

    class carbohydrateBuilder
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTORS                   //
        //////////////////////////////////////////////////////////
        carbohydrateBuilder(const std::string& condensedSequence, const ParameterManager& param);
        carbohydrateBuilder(std::string condensedSequence);
        std::string GetNumberOfShapes(bool likelyShapesOnly = false) const;
        LinkageOptionsVector GenerateUserOptionsDataStruct();
        void GenerateSpecific3DStructure(SingleRotamerInfoVector conformerInfo, std::string fileOutputDirectory);

        ParameterManager parameters;
        Molecule carbohydrate;
        std::vector<ResidueLinkage> glycosidicLinkages;
    };
} // namespace gmml
#endif
