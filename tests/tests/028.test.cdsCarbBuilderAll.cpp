#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h> // stat
#include <stdexcept>

int main(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cerr << "Usage:  " << argv[0]
                  << " <List IDs and sequences> <Char delimiter used in list> "
                     "<Folder to put outputs>\n";
        std::cerr << "Example: " << argv[0] << " exampleLibrary.txt _ outputs/ ";
        std::cerr << "Don't use a delimiter that appears in glycan sequences or ids. Like - or , or [] etc\n";
        std::exit(EXIT_FAILURE);
    }
    // Convert command line inputs to legible variables
    std::ifstream infile(argv[1]);
    char delimiter = argv[2][0]; // The second [0] gets me the first element of the argv which is type char**
    std::string outputFolderName = argv[3];
    struct stat info;
    if (stat(argv[3], &info) != 0)
    {
        std::cerr << "Folder " << outputFolderName << "/ does not exist and it isn't my job to make it.\n";
        std::exit(EXIT_FAILURE);
    }
    std::string line;
    while (std::getline(infile, line))
    {
        if (line.empty())
        {
            continue;
        }
        std::vector<std::string> splitLine = codeUtils::split(line, delimiter);
        if (splitLine.size() != 2)
        {
            std::cerr << "Encountered problem when splitting this line >>>" << line << "<<< from file >>>" << argv[1]
                      << "<<< into an ID and carb string separated by your specified delimiter: >>>" << argv[2]
                      << "<<<\n";
            std::exit(EXIT_FAILURE);
        }
        std::string inputGlycanID = splitLine.at(0);
        std::string inputSequence = splitLine.at(1);
        std::cout << "\n*********************\nBuilding id " << inputGlycanID << ":" << inputSequence
                  << "\n*********************\n";
        try
        { // I want to first build and reorder the sequence like gems does, and pass that to the carbohydrateBuilder
            CondensedSequence::Sequence sequence(inputSequence);
            std::cout << "Input:\n" << sequence.getInputSequence() << "\n";
            std::cout << "Interpreted:\n" << sequence.getInterpretedSequence() << "\n";
            std::cout << "IndexOrdered:\n" << sequence.getIndexOrdered() << "\n";
            std::cout << "IndexOrderedLabeled:\n" << sequence.getIndexLabeled() << "\n";
            std::cout << "Parsed and labelled " << inputGlycanID << " with no exceptions thrown.\n\n";
            cdsCondensedSequence::carbohydrateBuilder carbBuilder(sequence.getIndexOrdered());
            bool likelyShapesOnly = true; // Not sure I like this. Two functions probably more readable.
            std::cout << "Number of residues for this sequence is " << carbBuilder.GetCarbohydrate().GetResidueCount()
                      << "\n";
            std::cout << "Number of likely shapes is " << carbBuilder.GetNumberOfShapes(likelyShapesOnly) << "\n";
            std::cout << "Number of possible shapes is " << carbBuilder.GetNumberOfShapes() << "\n";
            for (auto& residue : carbBuilder.GetCarbohydrate().getResidues())
            {
                std::cout << residue->getStringId() << "\n";
            }
            for (auto& linkageInfo : carbBuilder.GenerateUserOptionsDataStruct())
            {
                std::cout << "Name: " << linkageInfo.linkageName_
                          << ", LinkageIndex: " << linkageInfo.indexOrderedLabel_
                          << ", Res1: " << linkageInfo.firstResidueNumber_
                          << ", Res2: " << linkageInfo.secondResidueNumber_ << "\n";
                for (auto& dihedralInfo : linkageInfo.likelyRotamers_)
                {
                    std::cout << "    LikelyRotamers: " << dihedralInfo.dihedralName_;
                    for (auto& rotamer : dihedralInfo.rotamers_)
                    {
                        std::cout << ", " << rotamer;
                    }
                    std::cout << "\n";
                }
            }
            std::vector<std::string> headerLines {"Produced by GMML carb builder all test"};
            carbBuilder.GetCarbohydrate().Generate3DStructureFiles(outputFolderName, inputGlycanID, headerLines);
        }
        catch (const std::runtime_error& error)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
            std::cerr << "Error thrown by the carbohydrateBuilder in gmml during construction was: " << error.what()
                      << std::endl;
        }
        catch (...)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR,
                      "carbohydrateBuilder class caught a throw that was not anticipated. Curious. Death cometh?");
            std::cerr << "ERROR carbohydrateBuilder caught a throw type that was not anticipated. Pretty please report "
                         "how you got to this to glycam@gmail.com.";
        }
    }
}
