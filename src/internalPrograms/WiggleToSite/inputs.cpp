#include "include/internalPrograms/WiggleToSite/inputs.hpp"

#include "include/util/files.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <fstream>
#include <iostream>

namespace gmml
{
    WiggleToSiteInputs::WiggleToSiteInputs(std::string inputFileName)
    {
        //    std::cout << "About to read " << inputFileName << std::endl << std::flush;
        std::ifstream infile(inputFileName);
        if (!infile)
        {
            std::string message = "Uh oh, input file: " + inputFileName + ", could not be opened for reading!\n";
            util::log(__LINE__, __FILE__, util::ERR, message);
            throw std::runtime_error(message);
        }
        while (infile) // While there's still stuff left to read
        {
            std::string strInput;
            getline(infile, strInput);
            if (util::startsWith(strInput, "Substrate:"))
            {
                substrateFile_ = util::split(strInput, ':').at(1);
            }
            if (util::startsWith(strInput, "Carbohydrate:"))
            {
                carbohydrateSequence_ = util::split(strInput, ':').at(1);
            }
            if (util::startsWith(strInput, "SuperimpositionTargetResidue:"))
            {
                std::string inputPortion = util::split(strInput, ':').at(1);
                superimpositionTargetResidue_ = pdb::ResidueId(util::split(inputPortion, '_'));
                //            std::cout << "superimpositionTargetResidue_ in input file is: " <<
                //            superimpositionTargetResidue_ << std::endl;
            }
            if (util::startsWith(strInput, "WiggleTargetResidue:"))
            {
                std::string inputPortion = util::split(strInput, ':').at(1);
                wigglingTargetResidue_ = pdb::ResidueId(util::split(inputPortion, '_'));
                //            std::cout << "wigglingTargetResidue_ in input file is: " << wigglingTargetResidue_ <<
                //            std::endl;
            }
            if (util::startsWith(strInput, "TargetModelNumber:"))
            {
                substrateModelNumber_ = std::stoi(util::split(strInput, ':').at(1));
            }
            if (util::startsWith(strInput, "CarbohydrateSuperimpositionResidue:"))
            {
                carbohydrateSuperimpositionResidue_ = std::stoi(util::split(strInput, ':').at(1));
            }
            if (util::startsWith(strInput, "CarbohydrateWigglingResidue:"))
            {
                carbohydrateWigglingResidue_ = std::stoi(util::split(strInput, ':').at(1));
            }
            if (util::startsWith(strInput, "persistCycles:"))
            {
                persistCycles_ = std::stoi(util::split(strInput, ':').at(1));
            }
            if (util::startsWith(strInput, "isDeterministic:"))
            {
                if (util::split(strInput, ':').at(1) == "true")
                {
                    isDeterministic_ = true;
                }
            }
        }
        util::log(__LINE__, __FILE__, util::INF, "Reading input file complete.");
    }

    std::string WiggleToSiteInputs::Print()
    {
        std::stringstream ss;
        ss << "substrateFile_ : " << substrateFile_ << "\n";
        ss << "carbohydrateFile_ : " << carbohydrateSequence_ << "\n";
        ss << "superimpositionTargetResidue_ : " << superimpositionTargetResidue_ << "\n";
        ss << "carbohydrateSuperimpositionResidue_ : " << carbohydrateSuperimpositionResidue_ << "\n";
        ss << "wigglingTargetResidue_ : " << wigglingTargetResidue_ << "\n";
        ss << "carbohydrateWigglingResidue_ : " << carbohydrateWigglingResidue_ << "\n";
        ss << "persistCycles_ : " << persistCycles_ << "\n";
        ss << "isDeterministic_ : " << isDeterministic_ << "\n";
        return ss.str();
    }
} // namespace gmml
