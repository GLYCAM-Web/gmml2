#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"
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
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

namespace
{
    std::vector<size_t> assemblyResidues(pdb::PdbData& data, size_t assemblyId)
    {
        std::vector<size_t> residueAssembly =
            codeUtils::indicesToValues(data.indices.moleculeAssembly, data.indices.residueMolecule);
        return codeUtils::indicesOfElement(residueAssembly, assemblyId);
    }

    std::string residueParmName(const pdb::PdbData& data, size_t residueId)
    {
        const std::string& name = data.residues.names[residueId];
        if (data.residues.isNTerminal[residueId])
        {
            return "N" + name;
        }
        else if (data.residues.isCTerminal[residueId])
        {
            return "C" + name;
        }
        return name;
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
            data.assemblies.numbers[assemblyId] = currentModelNumber;
            assembly.setNumber(currentModelNumber);
        }
        else if ((recordName == "ATOM") || (recordName == "HETATM"))
        { // Gimme everything with the same chain, can be everything with no chain.
            // Function that will read from stringstream until chain ID changes or TER or just not ATOM/HETATM
            std::string chainId                  = extractChainId(line);
            std::stringstream singleChainSection = extractSingleChainFromRecordSection(stream_block, line, chainId);
            size_t moleculeId                    = data.indices.moleculeCount;
            data.indices.moleculeAssembly.push_back(assemblyId);
            data.indices.moleculeCount++;
            data.molecules.chainIds.push_back(chainId);
            data.molecules.residueOrder.push_back({});
            std::unique_ptr<cds::Molecule> molecule = std::make_unique<cds::Molecule>();
            cds::Molecule* mol                      = assembly.addMolecule(std::move(molecule));
            data.objects.molecules.push_back(mol);
            readChain(data, moleculeId, singleChainSection);
        }
        else if (recordName == "ENDMDL")
        { // Only happens when reading "modelsAsCoordinates", now read the rest of the entry as extra coords for
          // MODEL 1.
            gmml::log(__LINE__, __FILE__, gmml::INF, "PdbFile being read in as trajectory");
            extractCoordinatesFromModel(data, stream_block, line);
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

void pdb::extractCoordinatesFromModel(PdbData& data, std::stringstream& stream_block, std::string line)
{
    const int iPdbLineLength = 80; // repeat for now, fix later
    gmml::log(__LINE__, __FILE__, gmml::INF, "Section to extract coordinates from is\n" + stream_block.str());
    size_t atomCount = data.indices.atomCount;
    if (atomCount == 0)
    {
        std::string message = "No atoms available when extracting coords from multiple models";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw message;
    }
    const std::vector<cds::Coordinate>& initial = data.atoms.coordinates;
    data.trajectory.coordinates                 = {initial};
    std::vector<cds::Coordinate> current        = initial;
    size_t atomId                               = 0;
    while ((std::getline(stream_block, line)))
    {
        expandLine(line, iPdbLineLength);
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
        if (recordName == "ATOM" || recordName == "HETATM")
        {
            if (atomId >= atomCount)
            {
                std::string message = "Trajectory model contains more atoms than initial model";
                gmml::log(__LINE__, __FILE__, gmml::ERR, message);
                throw message;
            }
            current[atomId] = checkShiftsAndExtractCoordinate(line);
            atomId++;
        }
        if (recordName == "ENDMDL")
        { // reset to read next set of coords
            atomId = 0;
            data.trajectory.coordinates.push_back(current);
            current = initial;
        }
    }
    return;
}

void pdb::ChangeResidueName(PdbData& data, size_t assemblyId, const std::string& selector, const std::string& newName)
{
    for (size_t residueId : assemblyResidues(data, assemblyId))
    {
        std::size_t found = pdbResidueId(data, residueId).print().find(selector);
        if (found != std::string::npos)
        {
            data.residues.names[residueId] = newName;
            data.objects.residues[residueId]->setName(newName);
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
    std::function<bool(const size_t&)> isCYS = [&](size_t n)
    {
        return codeUtils::contains({"CYS", "CYX"}, data.residues.names[n]);
    };
    gmml::log(__LINE__, __FILE__, gmml::INF, "Start CYS preprocessing for this Model\n");
    std::vector<size_t> residueAssembly =
        codeUtils::indicesToValues(data.indices.moleculeAssembly, data.indices.residueMolecule);
    std::vector<size_t> assemblyResidues = codeUtils::indicesOfElement(residueAssembly, assemblyId);
    std::vector<size_t> cysResidues      = codeUtils::vectorFilter(isCYS, assemblyResidues);
    if (cysResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "No CYS or CYX residues detected in this structure\n");
    }
    for (size_t n = 0; n < cysResidues.size(); n++)
    { // I want to go through the list and compare from current item to end. Thus it2 = std::next it1
        size_t cysRes1Id   = cysResidues[n];
        size_t sgAtom1Id   = findResidueAtom(data, cysRes1Id, "SG");
        cds::Atom* sgAtom1 = data.objects.atoms[sgAtom1Id];
        for (size_t k = n + 1; k < cysResidues.size(); k++)
        {
            size_t cysRes2Id   = cysResidues[k];
            size_t sgAtom2Id   = findResidueAtom(data, cysRes2Id, "SG");
            cds::Atom* sgAtom2 = data.objects.atoms[sgAtom2Id];
            size_t atomCount   = data.indices.atomCount;
            if ((sgAtom1Id < atomCount) && (sgAtom2Id < atomCount))
            {
                double distance = cds::distance(sgAtom1->coordinate(), sgAtom2->coordinate());
                if (distance < constants::dSulfurCutoff && distance > 0.001)
                {
                    data.residues.names[cysRes1Id] = "CYX";
                    data.residues.names[cysRes2Id] = "CYX";
                    data.objects.residues[cysRes1Id]->setName("CYX");
                    data.objects.residues[cysRes2Id]->setName("CYX");
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
    for (size_t residueId : assemblyResidues(data, assemblyId))
    {
        const std::string& name = data.residues.names[residueId];
        if (codeUtils::contains({"HIE", "HID", "HIP"}, name))
        {
            ppInfo.hisResidues_.emplace_back(pdbResidueId(data, residueId));
        }
        else if (name == "HIS")
        {
            size_t atomCount    = data.indices.atomCount;
            size_t atomHE2      = findResidueAtom(data, residueId, "HE2");
            size_t atomHD1      = findResidueAtom(data, residueId, "HD1");
            auto newResidueName = [&]()
            {
                if ((atomHE2 == atomCount) && (atomHD1 < atomCount))
                {
                    return "HID";
                }
                else if ((atomHE2 < atomCount) && (atomHD1 < atomCount))
                {
                    return "HIP";
                }
                else // HIE is default
                {
                    return "HIE";
                }
            };
            std::string newName            = newResidueName();
            data.residues.names[residueId] = newName;
            data.objects.residues[residueId]->setName(newName);
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
        if (nTerResidue >= data.indices.residueCount)
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
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Gap detection started for chain " + data.molecules.chainIds[moleculeId]);
        std::vector<size_t> proteinResidues = codeUtils::boolsToIndices(
            codeUtils::vectorAnd(codeUtils::vectorEquals(data.residues.types, cds::ResidueType::Protein),
                                 codeUtils::vectorEquals(data.indices.residueMolecule, moleculeId)));
        if (proteinResidues.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "No protein residues found in chain with id: " + data.molecules.chainIds[moleculeId]);
            break;
        }
        for (size_t n = 0; n < proteinResidues.size() - 1; n++)
        {
            size_t res1      = proteinResidues[n];
            size_t res2      = proteinResidues[n + 1];
            size_t atomCount = data.indices.atomCount;
            size_t res1AtomC = findResidueAtom(data, res1, "C");
            size_t res2AtomN = findResidueAtom(data, res2, "N");
            if ((res1AtomC < atomCount) && (res2AtomN < atomCount) &&
                (!isWithinBondingDistance(data, res1AtomC, res2AtomN)))
            { // GAP detected
                // Look for non-natural protein residues within bonding distance, they fall under ResidueType
                // Undefined, this indicates it's not gap.
                if (!amberMdPrep::checkForNonNaturalProteinResidues(
                        data, codeUtils::indicesOfElement(data.residues.types, cds::ResidueType::Undefined), res1AtomC,
                        ppInfo))
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
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Gap detection completed for chain " + data.molecules.chainIds[moleculeId]);
    }
    return;
}

void pdb::preProcessMissingUnrecognized(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo,
                                        const cdsParameters::ParameterManager& parmManager)
{
    for (size_t residueId : assemblyResidues(data, assemblyId))
    {
        std::string parmName   = residueParmName(data, residueId);
        ResidueId residueIdObj = pdbResidueId(data, residueId);
        size_t index           = codeUtils::indexOf(parmManager.lib.residueNames, parmName);
        // Unrecognized residue->
        if (index == parmManager.lib.residueNames.size())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "ParmManager did not recognize residue: " + parmName);
            ppInfo.unrecognizedResidues_.emplace_back(residueIdObj);
        }
        else // Recognized residue->
        {
            const lib::ResidueData& parmResidue = parmManager.lib.residues[index];
            std::vector<size_t> parmHeavyAtoms  = cdsSelections::FindHeavyAtoms(parmResidue.atoms.elements);
            std::vector<std::string> parmHeavyAtomNames =
                codeUtils::indicesToValues(parmResidue.atoms.names, parmHeavyAtoms);
            std::vector<size_t> atomIds           = residueAtoms(data, residueId);
            std::vector<std::string> pdbAtomNames = codeUtils::indicesToValues(data.atoms.names, atomIds);
            for (auto& parmHeavyAtomName : parmHeavyAtomNames) // What heavy atoms are missing from the pdb residue?
            {
                if (!codeUtils::contains(pdbAtomNames, parmHeavyAtomName))
                { // Residue missing a heavy atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Atom named " + parmHeavyAtomName + " missing from " + residueIdObj.print());
                    ppInfo.missingHeavyAtoms_.emplace_back(parmHeavyAtomName, residueIdObj);
                }
            }
            for (auto& pdbAtomName : pdbAtomNames) // What atoms in the pdb residue are unrecognized?
            {
                if (!codeUtils::contains(parmResidue.atoms.names, pdbAtomName))
                {
                    // Residue contains unrecognized atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Unrecognized atom named " + pdbAtomName + " in " + residueIdObj.print());
                    ppInfo.unrecognizedAtoms_.emplace_back(pdbAtomName, residueIdObj);
                }
            }
        }
    }
    return;
}

void pdb::Write(const PdbData& data, const std::vector<std::vector<size_t>>& moleculeResidues, std::ostream& stream)
{
    std::function<bool(const size_t&)> atomAlive = [&](size_t n)
    {
        return data.atomGraph.nodeAlive[n];
    };
    std::function<std::string(const size_t&)> elementString = [&](size_t n)
    {
        return data.atoms.names[n].empty() ? "" : data.atoms.names[n].substr(0, 1);
    };
    std::function<std::string(const size_t&)> truncatedResidueNames = [&](size_t n)
    {
        return codeUtils::truncate(3, data.residues.names[n]);
    };

    assembly::Graph graph = cds::createAssemblyGraph(data.indices, data.atomGraph);
    cds::PdbFileResidueData residueData {
        data.residues.numbers,
        codeUtils::vectorMap(truncatedResidueNames, codeUtils::indexVector(data.indices.residueCount)),
        data.residues.chainIds, data.residues.insertionCodes};
    cds::PdbFileAtomData atomData {data.atoms.coordinates,
                                   data.atoms.numbers,
                                   data.atoms.names,
                                   codeUtils::vectorMap(elementString, codeUtils::indexVector(data.indices.atomCount)),
                                   data.atoms.recordNames,
                                   data.atoms.occupancies,
                                   data.atoms.temperatureFactors};
    cds::PdbFileFormat format;
    cds::PdbFileData writerData {format, {}, residueData, atomData};
    for (auto& residueIds : moleculeResidues)
    {
        for (size_t residueId : residueIds)
        {
            cds::writeAssemblyToPdb(stream, graph, {{residueId}}, {{data.residues.hasTerCard[residueId]}}, {},
                                    writerData);
        }
        if (!residueIds.empty())
        { // Sometimes you get empty chains after things have been deleted I guess.
            stream << "TER\n";
        }
    }
}
