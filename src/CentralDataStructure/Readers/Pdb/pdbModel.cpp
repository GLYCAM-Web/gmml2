#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/print.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Editors/amberMdPrep.hpp" //all preprocessing should move to here.
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/cdsFunctions/bondByDistance.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

namespace
{
    std::vector<cds::Residue*> assemblyResidues(pdb::PdbData& data, size_t assemblyId)
    {
        std::vector<size_t> residueAssembly =
            codeUtils::indicesToValues(data.indices.moleculeAssembly, data.indices.residueMolecule);
        std::vector<size_t> residueIds = codeUtils::indicesOfElement(residueAssembly, assemblyId);
        return codeUtils::indicesToValues(data.objects.residues, residueIds);
    }
} // namespace

void pdb::readAssembly(PdbData& data, size_t assemblyId, cds::Assembly& assembly, std::stringstream& stream_block)
{
    int currentModelNumber        = 1;
    std::string previousResidueId = "InitialValue";
    std::string line;
    while (getline(stream_block, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if (recordName == "MODEL")
        {
            try
            {
                currentModelNumber = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(10, 4)));
            }
            catch (...)
            {
                gmml::log(__LINE__, __FILE__, gmml::WAR,
                          "Model issue: this ain't an int: " + codeUtils::RemoveWhiteSpace(line.substr(10, 4)));
                currentModelNumber = 1; // Seems like a reasonable default.
            }
            assembly.setNumber(currentModelNumber);
        }
        else if ((recordName == "ATOM") || (recordName == "HETATM"))
        { // Gimme everything with the same chain, can be everything with no chain.
            // Function that will read from stringstream until chain ID changes or TER or just not ATOM/HETATM
            std::stringstream singleChainSection =
                extractSingleChainFromRecordSection(stream_block, line, extractChainId(line));
            size_t moleculeId = data.objects.molecules.size();
            data.indices.moleculeAssembly.push_back(assemblyId);
            data.moleculeResidueOrder.push_back({});
            std::unique_ptr<cds::Molecule> molecule = std::make_unique<cds::Molecule>(extractChainId(line));
            cds::Molecule* mol                      = assembly.addMolecule(std::move(molecule));
            data.objects.molecules.push_back(mol);
            readChain(data, moleculeId, singleChainSection);
        }
        else if (recordName == "ENDMDL")
        { // Only happens when reading "modelsAsCoordinates", now read the rest of the entry as extra coords for
          // MODEL 1.
            gmml::log(__LINE__, __FILE__, gmml::INF, "PdbFile being read in as trajectory");
            extractCoordinatesFromModel(assembly, stream_block, line);
        }
    }
    return;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
std::string pdb::extractChainId(const std::string& line)
{ // serialNumber can overrun into position 12 in input.
    int shift = checkShiftFromSerialNumberOverrun(line);
    return codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
}

std::stringstream pdb::extractSingleChainFromRecordSection(std::stringstream& stream_block, std::string line,
                                                           const std::string& initialChainID)
{
    std::streampos previousLinePosition = stream_block.tellg(); // Save current line position
    std::stringstream singleChainSection;
    std::string chainID    = initialChainID;
    std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
    while ((chainID == initialChainID) && (recordName != "TER"))
    {
        singleChainSection << line << std::endl;
        previousLinePosition = stream_block.tellg(); // Save current line position.
        if (!std::getline(stream_block, line))
        {
            break; // // If we hit the end, time to leave.
        }
        chainID    = extractChainId(line);
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
    }
    stream_block.seekg(
        previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
    return singleChainSection;
}

void pdb::extractCoordinatesFromModel(cds::Assembly& assembly, std::stringstream& stream_block, std::string line)
{
    const int iPdbLineLength = 80; // repeat for now, fix later
    gmml::log(__LINE__, __FILE__, gmml::INF, "Section to extract coordinates from is\n" + stream_block.str());
    std::vector<Atom*> myAtoms = assembly.getAtoms();
    if (myAtoms.empty())
    {
        std::string message = "No atoms available when extracting coords from multiple models";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw message;
    }
    std::vector<Atom*>::iterator it = myAtoms.begin();
    while ((std::getline(stream_block, line)))
    {
        expandLine(line, iPdbLineLength);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if (recordName == "ATOM" || recordName == "HETATM")
        {
            Atom* atomPtr = *it;
            atomPtr->addCoordinate(checkShiftsAndExtractCoordinate(line));
            it++;
        }
        if (recordName == "ENDMDL")
        { // reset to read next set of coords
            it = myAtoms.begin();
        }
    }
    return;
}

void pdb::ChangeResidueName(PdbData& data, size_t assemblyId, const std::string& selector, const std::string& newName)
{
    for (auto residue : assemblyResidues(data, assemblyId))
    {
        size_t residueId  = codeUtils::indexOf(data.objects.residues, residue);
        std::size_t found = pdbResidueId(data, residueId).print().find(selector);
        if (found != std::string::npos)
        {
            residue->setName(newName);
            return;
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Could not find residue to rename with this selector " + selector);
    return;
}

//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////
void pdb::preProcessCysResidues(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Start CYS preprocessing for this Model\n");
    std::vector<size_t> residueAssembly =
        codeUtils::indicesToValues(data.indices.moleculeAssembly, data.indices.residueMolecule);
    std::vector<size_t> assemblyResidues   = codeUtils::indicesOfElement(residueAssembly, assemblyId);
    std::vector<cds::Residue*> cysResidues = codeUtils::getElementsWithNames(
        codeUtils::indicesToValues(data.objects.residues, assemblyResidues), std::vector<std::string> {"CYS", "CYX"});
    if (cysResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "No CYS or CYX residues detected in this structure\n");
    }
    for (std::vector<cds::Residue*>::iterator it1 = cysResidues.begin(); it1 != cysResidues.end(); ++it1)
    { // I want to go through the list and compare from current item to end. Thus it2 = std::next it1
        cds::Residue* cysRes1 = *it1;
        size_t cysRes1Id      = codeUtils::indexOf(data.objects.residues, cysRes1);
        size_t sgAtom1Id      = findResidueAtom(data, cysRes1Id, "SG");
        cds::Atom* sgAtom1    = data.objects.atoms[sgAtom1Id];
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1, 1); it2 != cysResidues.end(); ++it2)
        {
            cds::Residue* cysRes2 = *it2;
            size_t cysRes2Id      = codeUtils::indexOf(data.objects.residues, cysRes2);
            size_t sgAtom2Id      = findResidueAtom(data, cysRes2Id, "SG");
            cds::Atom* sgAtom2    = data.objects.atoms[sgAtom2Id];
            size_t atomCount      = data.objects.atoms.size();
            if ((sgAtom1Id < atomCount) && (sgAtom2Id < atomCount))
            {
                double distance = cds::distance(sgAtom1->coordinate(), sgAtom2->coordinate());
                if (distance < constants::dSulfurCutoff && distance > 0.001)
                {
                    cysRes1->setName("CYX");
                    cysRes2->setName("CYX");
                    ResidueId Res1Id = pdbResidueId(data, cysRes1Id);
                    ResidueId Res2Id = pdbResidueId(data, cysRes2Id);
                    addBond(data, sgAtom1Id, sgAtom2Id); // I think I want this here. Not 100%.
                    ppInfo.cysBondResidues_.emplace_back(Res1Id, Res2Id, distance);
                    std::stringstream message;
                    message << "Bonding " << Res1Id.print() << " and " << Res2Id.print() << " with distance "
                            << distance;
                    gmml::log(__LINE__, __FILE__, gmml::INF, message.str());
                }
            }
        }
    }
    return;
}

void pdb::preProcessHisResidues(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo,
                                const PreprocessorOptions& inputOptions)
{
    // HIS protonation, user specified:
    gmml::log(__LINE__, __FILE__, gmml::INF, "User His protonation");
    for (auto& userSelectionPair : inputOptions.hisSelections_)
    {
        ChangeResidueName(data, assemblyId, userSelectionPair.first, userSelectionPair.second);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Auto His protonation");
    // HIS protonation, automatic handling.
    for (auto residue : assemblyResidues(data, assemblyId))
    {
        size_t residueId = codeUtils::indexOf(data.objects.residues, residue);
        if (residue->getName() == "HIE" || residue->getName() == "HID" || residue->getName() == "HIP")
        {
            ppInfo.hisResidues_.emplace_back(pdbResidueId(data, residueId));
        }
        else if (residue->getName() == "HIS")
        {
            size_t atomCount = data.objects.atoms.size();
            size_t atomHE2   = findResidueAtom(data, residueId, "HE2");
            size_t atomHD1   = findResidueAtom(data, residueId, "HD1");
            if ((atomHE2 == atomCount) && (atomHD1 < atomCount))
            {
                residue->setName("HID");
            }
            else if ((atomHE2 < atomCount) && (atomHD1 < atomCount))
            {
                residue->setName("HIP");
            }
            else // HIE is default
            {
                residue->setName("HIE");
            }
            gmml::log(__LINE__, __FILE__, gmml::INF, "About to emplaceBack Id");
            gmml::log(__LINE__, __FILE__, gmml::INF, pdbResidueId(data, residueId).print());
            ppInfo.hisResidues_.emplace_back(pdbResidueId(data, residueId));
        }
    }
    return;
}

void pdb::preProcessChainTerminals(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo,
                                   const PreprocessorOptions& inputOptions)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain terminations");
    std::vector<size_t> moleculeIds = codeUtils::indicesOfElement(data.indices.moleculeAssembly, assemblyId);
    for (size_t moleculeId : moleculeIds)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Chain termination processing started for this chain");
        // Do the thing
        size_t nTerResidue = getNTerminal(data, moleculeId);
        if (nTerResidue >= data.objects.residues.size())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Could not modify terminals of this chain.");
        }
        else
        {
            ModifyTerminal(data, nTerResidue, inputOptions.chainNTermination_);
            size_t cTerResidue = getCTerminal(data, moleculeId);
            ModifyTerminal(data, cTerResidue, inputOptions.chainCTermination_);
            gmml::log(__LINE__, __FILE__, gmml::INF, "N term : " + pdbResidueId(data, nTerResidue).print());
            gmml::log(__LINE__, __FILE__, gmml::INF, "C term : " + pdbResidueId(data, cTerResidue).print());
            // Report the thing
            ppInfo.chainTerminals_.emplace_back(data.residues.chainIds[nTerResidue],
                                                getNumberAndInsertionCode(data, nTerResidue),
                                                getNumberAndInsertionCode(data, cTerResidue),
                                                inputOptions.chainNTermination_, inputOptions.chainCTermination_);
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocessing complete for this chain");
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain termination processing complete");
}

void pdb::preProcessGapsUsingDistance(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo,
                                      const PreprocessorOptions& inputOptions)
{
    // Missing Residues (gaps); If two sequential protein residues in the same molecule aren't close enough to bond:
    // this is a gap regardless of residue number/insertion code. User will want caps(ACE/NME) or zwitterionic, we can't
    // know ourselves without knowledge of the system, but most of the time caps.
    gmml::log(__LINE__, __FILE__, gmml::INF, "Gaps");
    std::vector<size_t> moleculeIds = codeUtils::indicesOfElement(data.indices.moleculeAssembly, assemblyId);
    for (size_t moleculeId : moleculeIds)
    {
        cds::Molecule* cdsMolecule = data.objects.molecules[moleculeId];
        gmml::log(__LINE__, __FILE__, gmml::INF, "Gap detection started for chain " + cdsMolecule->GetChainId());
        std::vector<cds::Residue*> proteinResidues =
            cdsSelections::selectResiduesByType(cdsMolecule->getResidues(), cds::ResidueType::Protein);
        if (proteinResidues.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "No protein residues found in chain with id: " + cdsMolecule->GetChainId());
            break;
        }
        std::vector<cds::Residue*>::iterator it2;
        for (std::vector<cds::Residue*>::iterator it1 = proteinResidues.begin(); it1 != proteinResidues.end() - 1;
             ++it1)
        {
            it2              = std::next(it1);
            size_t res1      = codeUtils::indexOf(data.objects.residues, *it1);
            size_t res2      = codeUtils::indexOf(data.objects.residues, *it2);
            size_t atomCount = data.objects.atoms.size();
            size_t res1AtomC = findResidueAtom(data, res1, "C");
            size_t res2AtomN = findResidueAtom(data, res2, "N");
            if ((res1AtomC < atomCount) && (res2AtomN < atomCount) &&
                (!isWithinBondingDistance(data.objects.atoms[res1AtomC], data.objects.atoms[res2AtomN])))
            { // GAP detected
                // Look for non-natural protein residues within bonding distance, they fall under ResidueType
                // Undefined, this indicates it's not gap.
                if (!amberMdPrep::checkForNonNaturalProteinResidues(
                        data,
                        cdsSelections::selectResiduesByType(cdsMolecule->getResidues(), cds::ResidueType::Undefined),
                        data.objects.atoms[res1AtomC], ppInfo))
                {
                    // Log it
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              inputOptions.gapCTermination_ + " cap for : " + pdbResidueId(data, res1).print());
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              inputOptions.gapNTermination_ + " cap for : " + pdbResidueId(data, res2).print());
                    // Do it
                    InsertCap(data, moleculeId, res1, inputOptions.gapCTermination_);
                    InsertCap(data, moleculeId, res2, inputOptions.gapNTermination_);
                    // Record it
                    ppInfo.missingResidues_.emplace_back(data.residues.chainIds[res1],
                                                         getNumberAndInsertionCode(data, res1),
                                                         getNumberAndInsertionCode(data, res2),
                                                         inputOptions.gapCTermination_, inputOptions.gapNTermination_);
                }
            }
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Gap detection completed for chain " + cdsMolecule->GetChainId());
    }
    return;
}

void pdb::preProcessMissingUnrecognized(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo,
                                        const cdsParameters::ParameterManager& parmManager)
{
    for (auto residue : assemblyResidues(data, assemblyId))
    {
        ResidueId residueId = pdbResidueId(data, codeUtils::indexOf(data.objects.residues, residue));
        size_t index        = codeUtils::indexOf(parmManager.lib.residueNames, residue->GetParmName());
        // Unrecognized residue->
        if (index == parmManager.lib.residueNames.size())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "ParmManager did not recognize residue: " + residue->GetParmName());
            ppInfo.unrecognizedResidues_.emplace_back(residueId);
        }
        else // Recognized residue->
        {
            const lib::ResidueData& parmResidue = parmManager.lib.residues[index];
            std::vector<size_t> parmHeavyAtoms  = cdsSelections::FindHeavyAtoms(parmResidue.atoms.elements);
            std::vector<std::string> parmHeavyAtomNames =
                codeUtils::indicesToValues(parmResidue.atoms.names, parmHeavyAtoms);
            std::vector<std::string> pdbAtomNames = residue->getAtomNames();
            for (auto& parmHeavyAtomName : parmHeavyAtomNames) // What heavy atoms are missing from the pdb residue?
            {
                if (!codeUtils::contains(pdbAtomNames, parmHeavyAtomName))
                { // Residue missing a heavy atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Atom named " + parmHeavyAtomName + " missing from " + residueId.print());
                    ppInfo.missingHeavyAtoms_.emplace_back(parmHeavyAtomName, residueId);
                }
            }
            for (auto& pdbAtomName : pdbAtomNames) // What atoms in the pdb residue are unrecognized?
            {
                if (!codeUtils::contains(parmResidue.atoms.names, pdbAtomName))
                {
                    // Residue contains unrecognized atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Unrecognized atom named " + pdbAtomName + " in " + residueId.print());
                    ppInfo.unrecognizedAtoms_.emplace_back(pdbAtomName, residueId);
                }
            }
        }
    }
    return;
}

void pdb::Print(const cds::Assembly& assembly, std::ostream& out)
{
    for (auto& residue : assembly.getResidues())
    {
        cds::print(out, *residue);
    }
}

void pdb::Write(const PdbData& data, const std::vector<std::vector<size_t>>& moleculeResidues, std::ostream& stream)
{
    for (auto& residueIds : moleculeResidues)
    {
        for (size_t residueId : residueIds)
        {
            std::vector<size_t> residueAtoms   = codeUtils::indicesOfElement(data.indices.atomResidue, residueId);
            cds::Residue* residue              = data.objects.residues[residueId];
            cds::GraphIndexData residueIndices = cds::toIndexData({residue});
            assembly::Graph graph              = cds::createCompleteAssemblyGraph(residueIndices);
            std::vector<cds::Atom*> atoms      = residue->getAtoms();
            cds::PdbFileResidueData residueData {{residue->getNumber()},
                                                 {cds::truncatedResidueName(residue)},
                                                 {data.residues.chainIds[residueId]},
                                                 {data.residues.insertionCodes[residueId]}};
            cds::PdbFileAtomData atomData =
                cds::toPdbFileAtomData(atoms, codeUtils::indicesToValues(data.atoms.recordNames, residueAtoms),
                                       codeUtils::indicesToValues(data.atoms.occupancies, residueAtoms),
                                       codeUtils::indicesToValues(data.atoms.temperatureFactors, residueAtoms));
            cds::PdbFileFormat format;
            cds::PdbFileData writerData {format, {}, residueData, atomData};
            cds::writeAssemblyToPdb(stream, graph, {{0}}, {{data.residues.hasTerCard[residueId]}}, {}, writerData);
        }
        if (!residueIds.empty())
        { // Sometimes you get empty chains after things have been deleted I guess.
            stream << "TER\n";
        }
    }
}
