#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Editors/glycamResidueCombinator.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// For loading select residues from glycam06.prep files
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " inputPrepFile <inputPrepFile2> <inputPrepFile3>\n";
        std::cout << "Exmpl: " << argv[0] << " ../dat/prep/GLYCAM_06j-1_GAGS_KDN_ABE.prep\n";
        std::cout << "Exmpl: " << argv[0] << " ../dat/prep/GLYCAM_06j-1_GAGS_KDN_ABE.prep tests/inputs/027.0LU.prep\n";
        std::exit(1);
    }
    // std::vector<std::string> residuesToLoadFromPrep = {{"0LD"}, {"0LU"}};
    std::string glycamPrepInputFile                 = argv[1];
    std::vector<std::string> residuesToLoadFromPrep = {
        {"0AA"},  {"0AB"}, {"0AD"}, {"0AU"}, {"0BA"}, {"0BB"}, {"0BC"}, {"0BD"}, {"0BU"}, {"0CA"}, {"0CB"},  {"0CD"},
        {"0CU"},  {"0DA"}, {"0DB"}, {"0DD"}, {"0DU"}, {"0EA"}, {"0EB"}, {"0FA"}, {"0FB"}, {"0GA"}, {"0GB"},  {"0GL"},
        {"0HA"},  {"0HB"}, {"0JA"}, {"0JB"}, {"0JD"}, {"0JU"}, {"0KA"}, {"0KB"}, {"0LA"}, {"0LB"}, {"0MA"},  {"0MB"},
        {"0NA"},  {"0NB"}, {"0OA"}, {"0OB"}, {"0PA"}, {"0PB"}, {"0PD"}, {"0PU"}, {"0QA"}, {"0QB"}, {"0RA"},  {"0RB"},
        {"0RD"},  {"0RU"}, {"0SA"}, {"0SB"}, {"0TA"}, {"0TB"}, {"0TV"}, {"0Tv"}, {"0UA"}, {"0UB"}, {"0VA"},  {"0VB"},
        {"0WA"},  {"0WB"}, {"0XA"}, {"0XB"}, {"0XD"}, {"0XU"}, {"0YA"}, {"0YB"}, {"0ZA"}, {"0ZB"}, {"0aA"},  {"0aB"},
        {"0aD"},  {"0aU"}, {"0bA"}, {"0bB"}, {"0bC"}, {"0bD"}, {"0bU"}, {"0cA"}, {"0cB"}, {"0cD"}, {"0cU"},  {"0dA"},
        {"0dB"},  {"0dD"}, {"0dU"}, {"0eA"}, {"0eB"}, {"0fA"}, {"0fB"}, {"0gA"}, {"0gB"}, {"0gL"}, {"0hA"},  {"0hB"},
        {"0jA"},  {"0jB"}, {"0jD"}, {"0jU"}, {"0kA"}, {"0kB"}, {"0lA"}, {"0lB"}, {"0mA"}, {"0mB"}, {"0nA"},  {"0nB"},
        {"0oA"},  {"0oB"}, {"0pA"}, {"0pB"}, {"0pD"}, {"0pU"}, {"0qA"}, {"0qB"}, {"0rA"}, {"0rB"}, {"0rD"},  {"0rU"},
        {"0sA"},  {"0sB"}, {"0tA"}, {"0tB"}, {"0tV"}, {"0tv"}, {"0uA"}, {"0uB"}, {"0vA"}, {"0vB"}, {"0wA"},  {"0wB"},
        {"0xA"},  {"0xB"}, {"0xD"}, {"0xU"}, {"0yA"}, {"0yB"}, {"0zA"}, {"0zB"}, {"0dR"}, {"045"}, {"0Yn"},  {"0YN"},
        {"0YS"},  {"0Ys"}, {"0yS"}, {"0ys"}, {"0Kn"}, {"0Ko"}, {"0KN"}, {"0KO"}, {"0AE"}, {"0Ae"}, {"0YnP"}, {"0YNP"},
        {"0ZBP"}, {"0LU"}, {"0LD"}, {"0an"}, {"0DH"}, {"0Dh"}, {"0eC"}, {"0ec"}, {"0gf"}, {"0KX"}, {"0LD"},  {"0LG"},
        {"0Lg"},  {"0LH"}, {"0Lh"}, {"0LU"}, {"0mP"}, {"0mp"}, {"0MR"}};

    std::vector<cds::Residue*> allGeneratedResidues;

    for (++argv; *argv; ++argv) // Don't use argv after this.
    {
        char* inputFileName = *argv;
        std::cout << "Loading " << inputFileName << std::endl;
        prep::PrepFile glycamPrepFile(inputFileName, residuesToLoadFromPrep);
        // allTheResidues.reserve(residuesToLoadFromPrep.size() * 20);
        for (auto& residue : glycamPrepFile.getResidues())
        {
            //            std::ofstream outFileStream;
            //            std::string fileName = residue->getName() + "_original.pdb";
            //            outFileStream.open(fileName.c_str());
            //            std::vector<cds::Residue*> vec = {residue};
            //            cds::writeMoleculeToPdb(outFileStream, {0}, {false}, cds::toPdbWriterData({vec}));
            //            outFileStream.close();
            std::cout << "Generating those combos from " << residue->getName() << std::endl;
            residueCombinator::generateResidueCombinations(allGeneratedResidues, residue);
        }
    }
    //    for (auto& combiRes : allGeneratedResidues)
    //    {
    //        std::ofstream outFileStream;
    //        std::string fileName = combiRes->getName() + ".pdb";
    //        outFileStream.open(fileName.c_str());
    //        std::vector<cds::Residue*> vec = {combiRes};
    //        cds::writeMoleculeToPdb(outFileStream, {0}, {false}, cds::toPdbWriterData({vec}));
    //        outFileStream.close();
    //    }

    // Note "CA2" was in the prep file, but Rob said delete. "Perhaps we had it before Amber did"
    // "4YP" is the same as 4YnP, which I'll generate from 0YnP anyway
    // The .Kn and .Ko entries that are missing were mistakes in naming in the GLYCAM_06j-1_GAGS_KDN.prep file.

    // special cases: I don't want the combos of these because they aren't glycans, or I don't have the 0 version of
    // it.
    std::vector<std::string> theSpecialCases = {{"ROH"},
                                                {"TBT"},
                                                {"OME"},
                                                {"ACX"},
                                                {"MEX"},
                                                {"SO3"},
                                                {"NLN"},
                                                {"OLS"},
                                                {"OLT"},
                                                {"ZOLT"},
                                                {"ZOLS"},
                                                // These are residues I don't have the unsubstituted "0" version of:
                                                {"YGa"},
                                                {"YuAP"},
                                                {"0uA1"},
                                                {"4uA1"},
                                                {"4uA2"},
                                                {"4uA3"}};
    prep::PrepFile glycamPrepFile2(glycamPrepInputFile, theSpecialCases);
    for (auto& specialResidue : glycamPrepFile2.getResidues())
    {
        allGeneratedResidues.push_back(specialResidue);
    }

    std::ofstream outFileStream;
    std::string fileName = "GLYCAM_06k.lib";
    outFileStream.open(fileName.c_str());
    cds::serializeResiduesIndividually(allGeneratedResidues);
    cds::WriteResiduesIndividuallyToOffFile(outFileStream, cds::toOffWriterData(allGeneratedResidues));
    outFileStream.close();
}
