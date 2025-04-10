#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryFile.hpp"
#include "includes/CentralDataStructure/FileFormats/offFileWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/offWriter.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/molecule.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

#include <iostream>
#include <ostream>
#include <stdexcept>

int main()
{
    std::string libFilePath = "../dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib";
    // std::string condensed_sequence =
    // "LIdopAa1-4DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-8DNeup5Aca2-8DNeup5Aca2-6DGalpb1-4DGlcpNAc[6S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    //    prep::PrepFile glycamPrepFile(prepFilePath);
    //    for ( auto &prepResidue : glycamPrepFile.getResidues() )
    //    {
    //    	prepResidue->SetConnectivities();
    //    }
    //    std::cout << "*\n*\n*\n*\n*\n*\n*\n*\n*\n";
    // std::vector<std::string> residuesToLoadFromPrep = {"0GA"};
    cds::Molecule molecule;
    lib::parseMolecule(&molecule, libFilePath);
    std::cout << "Finished loading libfile" << std::endl;
    // Need a central place for this:
    std::string fileName = "./libAsPdbFile.pdb";
    try
    {
        cds::GraphIndexData indices = cds::toIndexData({&molecule});
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WritePdb(stream, indices, {});
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing pdbFile class to file:\n" + fileName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
    }
    // OFF molecule
    molecule.setName("MOLECULE");
    fileName = "./libAsOffFile.off";
    try
    {
        cds::GraphIndexData indices = cds::toIndexData({&molecule});
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WriteOff(stream, molecule.getName(), indices);
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // OFF separate residues
    molecule.setName("LIBRARY");
    fileName = "./libAsLibFile.lib";
    try
    {
        cds::GraphIndexData indices     = cds::toIndexData(molecule.getResidues());
        std::vector<bool> includedAtoms = cds::atomVisibility(indices.atoms);
        assembly::Graph graph           = cds::createAssemblyGraph(indices, includedAtoms);
        cds::serializeResiduesIndividually(indices.residues);
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WriteResiduesIndividuallyToOffFile(stream, graph,
                                                                           cds::toOffFileData(indices.residues));
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
}
