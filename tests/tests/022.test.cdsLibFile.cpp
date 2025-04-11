#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Readers/Lib/LibraryFile.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
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
    lib::LibraryData data   = lib::loadLibraryData(libFilePath);
    cdsParameters::ParameterManager params {data};
    std::cout << "Finished loading libfile" << std::endl;
    cds::Molecule molecule = cds::Molecule();
    for (size_t n = 0; n < data.residues.size(); n++)
    {
        const std::string& name = data.residueNames[n];
        cds::Residue* residue   = molecule.addResidue(std::make_unique<cds::Residue>());
        residue->setName(name);
        residue->determineType(name);
        cdsParameters::createAtomsForResidue(params, residue, name);
    }
    cds::GraphIndexData moleculeIndices = cds::toIndexData({&molecule});
    // Need a central place for this:
    std::string fileName                = "./libAsPdbFile.pdb";
    try
    {
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WritePdb(stream, moleculeIndices, {});
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing pdbFile class to file:\n" + fileName);
        throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
    }
    // OFF molecule
    fileName = "./libAsOffFile.off";
    try
    {
        codeUtils::writeToFile(fileName,
                               [&](std::ostream& stream)
                               {
                                   cds::WriteOff(stream, "MOLECULE", moleculeIndices);
                               });
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Error when writing to file:\n" + fileName);
        throw std::runtime_error("Error when writing to file:\n" + fileName);
    }
    // OFF separate residues
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
