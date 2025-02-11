#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/strings.hpp"
#include <filesystem>
#include <iostream>
#include <string>
#include <fstream>

int main(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cerr << "Usage:   ./buildOligosaccharideLibrary <List IDs and sequences> <Char delimiter used in list> "
                     "<Folder to put outputs>\n";
        std::cerr << "Example: ./buildOligosaccharideLibrary exampleLibrary.txt _ outputs/ ";
        std::cerr << "Don't use a delimiter that appears in glycan sequences or ids. Like - or , or [] etc\n";
        std::exit(EXIT_FAILURE);
    }
    // Convert command line inputs to legible variables
    std::ifstream infile(argv[1]);
    char delimiter = argv[2][0]; // The second [0] gets me the first element of the argv which is type char**
    std::string outputFolderName = argv[3];
    codeUtils::createDirectories(outputFolderName);
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
        std::string inputSequence = splitLine.at(1);
        std::cout << "\n*********************\nBuilding " << inputSequence << "\n*********************\n";
        try
        {
            cdsCondensedSequence::carbohydrateBuilder carbBuilder(inputSequence);
            std::string inputGlycanID = splitLine.at(0);
            carbBuilder.GenerateSingle3DStructureDefaultFiles(outputFolderName, inputGlycanID);
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
