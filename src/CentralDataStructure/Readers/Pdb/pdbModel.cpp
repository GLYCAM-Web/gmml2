#include "includes/CentralDataStructure/Readers/Pdb/pdbModel.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileData.hpp"
#include "includes/CentralDataStructure/FileFormats/pdbFileWriter.hpp"
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
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
#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"

using pdb::PdbModel;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModel::PdbModel()
{}

PdbModel::PdbModel(std::stringstream& stream_block)
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
            this->setNumber(currentModelNumber);
        }
        else if ((recordName == "ATOM") || (recordName == "HETATM"))
        { // Gimme everything with the same chain, can be everything with no chain.
            // Function that will read from stringstream until chain ID changes or TER or just not ATOM/HETATM
            std::stringstream singleChainSection =
                this->extractSingleChainFromRecordSection(stream_block, line, this->extractChainId(line));
            this->addMolecule(std::make_unique<PdbChain>(singleChainSection, this->extractChainId(line)));
        }
        else if (recordName == "ENDMDL")
        { // Only happens when reading "modelsAsCoordinates", now read the rest of the entry as extra coords for
          // MODEL 1.
            gmml::log(__LINE__, __FILE__, gmml::INF, "PdbFile being read in as trajectory");
            this->extractCoordinatesFromModel(stream_block, line);
        }
    }
    // gmml::log(__LINE__, __FILE__, gmml::INF, "PdbModel Constructor Complete Captain");
    return;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
std::string PdbModel::extractChainId(const std::string& line)
{ // serialNumber can overrun into position 12 in input.
    int shift = pdb::checkShiftFromSerialNumberOverrun(line);
    return codeUtils::RemoveWhiteSpace(line.substr(21 + shift, 1));
}

std::stringstream PdbModel::extractSingleChainFromRecordSection(std::stringstream& stream_block, std::string line,
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
        chainID    = this->extractChainId(line);
        recordName = codeUtils::RemoveWhiteSpace(line.substr(0, 6));
    }
    stream_block.seekg(
        previousLinePosition); // Go back to previous line position. E.g. was reading HEADER and found TITLE.
                               //    gmml::log(__LINE__, __FILE__, gmml::INF,
    //              "Single chain section is:\n" + singleChainSection.str() + "\nEnd of single chain section.");
    return singleChainSection;
}

void PdbModel::extractCoordinatesFromModel(std::stringstream& stream_block, std::string line)
{
    const int iPdbLineLength = 80; // repeat for now, fix later
    gmml::log(__LINE__, __FILE__, gmml::INF, "Section to extract coordinates from is\n" + stream_block.str());
    std::vector<Atom*> myAtoms = this->getAtoms();
    if (myAtoms.empty())
    {
        std::string message = "No atoms available when extracting coords from multiple models";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw message;
    }
    std::vector<Atom*>::iterator it = myAtoms.begin();
    while ((std::getline(stream_block, line)))
    {
        pdb::expandLine(line, iPdbLineLength);
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

void PdbModel::addConectRecord(const cds::Atom* atom1, const cds::Atom* atom2)
{
    conectRecords_.emplace_back(std::vector<const cds::Atom*> {atom1, atom2});
    return;
}

void PdbModel::ChangeResidueName(const std::string& selector, const std::string& newName)
{
    for (auto& residue : this->getResidues())
    {
        std::size_t found = codeUtils::erratic_cast<PdbResidue*>(residue)->printId().find(selector);
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
void PdbModel::preProcessCysResidues(pdb::PreprocessorInformation& ppInfo)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Start CYS preprocessing for this Model\n");
    std::vector<cds::Residue*> cysResidues =
        codeUtils::getElementsWithNames(this->getResidues(), std::vector<std::string> {"CYS", "CYX"});
    if (cysResidues.empty())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "No CYS or CYX residues detected in this structure\n");
    }
    for (std::vector<cds::Residue*>::iterator it1 = cysResidues.begin(); it1 != cysResidues.end(); ++it1)
    { // I want to go through the list and compare from current item to end. Thus it2 = std::next it1
        PdbResidue* cysRes1 = codeUtils::erratic_cast<PdbResidue*>(*it1);
        cds::Atom* sgAtom1  = cysRes1->FindAtom("SG");
        for (std::vector<cds::Residue*>::iterator it2 = std::next(it1, 1); it2 != cysResidues.end(); ++it2)
        {
            PdbResidue* cysRes2 = codeUtils::erratic_cast<PdbResidue*>(*it2);
            cds::Atom* sgAtom2  = cysRes2->FindAtom("SG");
            if ((sgAtom1 != nullptr) && (sgAtom2 != nullptr))
            {
                // gmml::log(__LINE__, __FILE__, gmml::INF, "Found SG ATOMS");
                double distance = cds::distance(sgAtom1->coordinate(), sgAtom2->coordinate());
                if (distance < constants::dSulfurCutoff && distance > 0.001)
                {
                    // gmml::log(__LINE__, __FILE__, gmml::INF, "Distance less than cutoff");
                    cysRes1->setName("CYX");
                    cysRes2->setName("CYX");
                    // gmml::log(__LINE__, __FILE__, gmml::INF, "Names set");
                    addBond(sgAtom1, sgAtom2); // I think I want this here. Not 100%.
                    this->addConectRecord(sgAtom1, sgAtom2);
                    ppInfo.cysBondResidues_.emplace_back(cysRes1->getId(), cysRes2->getId(), distance);
                    // gmml::log(__LINE__, __FILE__, gmml::INF, "ThisNoHappen?");
                    std::stringstream message;
                    message << "Bonding " << cysRes1->printId() << " and " << cysRes2->printId() << " with distance "
                            << distance;
                    gmml::log(__LINE__, __FILE__, gmml::INF, message.str());
                }
            }
        }
    }
    return;
}

void PdbModel::preProcessHisResidues(pdb::PreprocessorInformation& ppInfo, const pdb::PreprocessorOptions& inputOptions)
{
    // HIS protonation, user specified:
    gmml::log(__LINE__, __FILE__, gmml::INF, "User His protonation");
    for (auto& userSelectionPair : inputOptions.hisSelections_)
    {
        this->ChangeResidueName(userSelectionPair.first, userSelectionPair.second);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Auto His protonation");
    // HIS protonation, automatic handling.
    for (auto& cdsresidue : this->getResidues())
    {
        PdbResidue* residue = codeUtils::erratic_cast<PdbResidue*>(cdsresidue);
        if (residue->getName() == "HIE" || residue->getName() == "HID" || residue->getName() == "HIP")
        {
            ppInfo.hisResidues_.emplace_back(residue->getId());
        }
        else if (residue->getName() == "HIS")
        {
            if ((residue->FindAtom("HE2") == nullptr) && (residue->FindAtom("HD1") != nullptr))
            {
                residue->setName("HID");
            }
            else if ((residue->FindAtom("HE2") != nullptr) && (residue->FindAtom("HD1") != nullptr))
            {
                residue->setName("HIP");
            }
            else // HIE is default
            {
                residue->setName("HIE");
            }
            gmml::log(__LINE__, __FILE__, gmml::INF, "About to emplaceBack Id");
            gmml::log(__LINE__, __FILE__, gmml::INF, residue->printId());
            ppInfo.hisResidues_.emplace_back(residue->getId());
        }
    }
    return;
}

void PdbModel::preProcessChainTerminals(pdb::PreprocessorInformation& ppInfo,
                                        const pdb::PreprocessorOptions& inputOptions)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain terminations");
    for (auto& cdsMolecule : this->getMolecules())
    {
        PdbChain* chain = codeUtils::erratic_cast<PdbChain*>(cdsMolecule);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Chain termination processing started for this chain");
        // Do the thing
        PdbResidue* nTerResidue = chain->getNTerminal();
        if (nTerResidue == nullptr)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Could not modify terminals of this chain.");
        }
        else
        {
            chain->ModifyTerminal(inputOptions.chainNTermination_, nTerResidue);
            PdbResidue* cTerResidue = chain->getCTerminal();
            chain->ModifyTerminal(inputOptions.chainCTermination_, cTerResidue);
            gmml::log(__LINE__, __FILE__, gmml::INF, "N term : " + nTerResidue->printId());
            gmml::log(__LINE__, __FILE__, gmml::INF, "C term : " + cTerResidue->printId());
            // Report the thing
            ppInfo.chainTerminals_.emplace_back(nTerResidue->getChainId(), nTerResidue->getNumberAndInsertionCode(),
                                                cTerResidue->getNumberAndInsertionCode(),
                                                inputOptions.chainNTermination_, inputOptions.chainCTermination_);
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Preprocessing complete for this chain");
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain termination processing complete");
    return;
}

void PdbModel::preProcessGapsUsingDistance(pdb::PreprocessorInformation& ppInfo,
                                           const pdb::PreprocessorOptions& inputOptions)
{
    // Missing Residues (gaps); If two sequential protein residues in the same molecule aren't close enough to bond:
    // this is a gap regardless of residue number/insertion code. User will want caps(ACE/NME) or zwitterionic, we can't
    // know ourselves without knowledge of the system, but most of the time caps.
    gmml::log(__LINE__, __FILE__, gmml::INF, "Gaps");
    for (auto& cdsMolecule : this->getMolecules())
    {
        PdbChain* chain = codeUtils::erratic_cast<PdbChain*>(cdsMolecule);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Gap detection started for chain " + chain->GetChainId());
        std::vector<cds::Residue*> proteinResidues =
            cdsSelections::selectResiduesByType(chain->getResidues(), cds::ResidueType::Protein);
        if (proteinResidues.empty())
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "No protein residues found in chain with id: " + chain->GetChainId());
            break;
        }
        std::vector<cds::Residue*>::iterator it2;
        for (std::vector<cds::Residue*>::iterator it1 = proteinResidues.begin(); it1 != proteinResidues.end() - 1;
             ++it1)
        {
            it2                        = std::next(it1);
            PdbResidue* res1           = codeUtils::erratic_cast<PdbResidue*>(*it1);
            PdbResidue* res2           = codeUtils::erratic_cast<PdbResidue*>(*it2);
            //                            std::cout << "res1 is " + res1->getNumberAndInsertionCode() + "_" +
            //                            res1->getChainId()
            //                            << std::endl; std::cout << "res2 is " + res2->getNumberAndInsertionCode() +
            //                            "_" + res2->getChainId() << std::endl;
            const cds::Atom* res1AtomC = res1->FindAtom("C");
            const cds::Atom* res2AtomN = res2->FindAtom("N");
            if ((res1AtomC != nullptr) && (res2AtomN != nullptr) && (!isWithinBondingDistance(res1AtomC, res2AtomN)))
            { // GAP detected
                // Look for non-natural protein residues within bonding distance, they fall under ResidueType
                // Undefined, this indicates it's not gap.
                if (!amberMdPrep::checkForNonNaturalProteinResidues(
                        cdsSelections::selectResiduesByType(chain->getResidues(), cds::ResidueType::Undefined),
                        res1AtomC, ppInfo))
                {
                    // Log it
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              inputOptions.gapCTermination_ + " cap for : " + res1->printId());
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              inputOptions.gapNTermination_ + " cap for : " + res2->printId());
                    // Do it
                    chain->InsertCap(*res1, inputOptions.gapCTermination_);
                    chain->InsertCap(*res2, inputOptions.gapNTermination_);
                    // Record it
                    ppInfo.missingResidues_.emplace_back(res1->getChainId(), res1->getNumberAndInsertionCode(),
                                                         res2->getNumberAndInsertionCode(),
                                                         inputOptions.gapCTermination_, inputOptions.gapNTermination_);
                }
            }
        }
        gmml::log(__LINE__, __FILE__, gmml::INF, "Gap detection completed for chain " + chain->GetChainId());
    }
    return;
}

void PdbModel::preProcessMissingUnrecognized(pdb::PreprocessorInformation& ppInfo,
                                             const cdsParameters::ParameterManager& parmManager)
{
    for (auto& cdsResidue : this->getResidues())
    {
        PdbResidue* residue            = codeUtils::erratic_cast<PdbResidue*>(cdsResidue);
        cds::Residue* parameterResidue = parmManager.findParameterResidue(residue->GetParmName());
        // Unrecognized residue->
        if (parameterResidue == nullptr)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "ParmManager did not recognize residue: " + residue->GetParmName());
            ppInfo.unrecognizedResidues_.emplace_back(residue->getId());
        }
        else // Recognized residue->
        {
            std::vector<std::string> parmHeavyAtomNames =
                cdsSelections::FindNamesOfAtoms(cdsSelections::FindHeavyAtoms(parameterResidue->getAtoms()));
            std::vector<std::string> parmAtomNames = parameterResidue->getAtomNames();
            std::vector<std::string> pdbAtomNames  = residue->getAtomNames();
            for (auto& parmHeavyAtomName : parmHeavyAtomNames) // What heavy atoms are missing from the pdb residue?
            {
                if (!codeUtils::contains(pdbAtomNames, parmHeavyAtomName))
                { // Residue missing a heavy atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Atom named " + parmHeavyAtomName + " missing from " + residue->printId());
                    ppInfo.missingHeavyAtoms_.emplace_back(parmHeavyAtomName, residue->getId());
                }
            }
            for (auto& pdbAtomName : pdbAtomNames) // What atoms in the pdb residue are unrecognized?
            {
                if (!codeUtils::contains(parmAtomNames, pdbAtomName))
                {
                    // Residue contains unrecognized atom.
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Unrecognized atom named " + pdbAtomName + " in " + residue->printId());
                    ppInfo.unrecognizedAtoms_.emplace_back(pdbAtomName, residue->getId());
                }
            }
        }
    }
    return;
}

// void PdbModel::bondAtomsByDistance()
//{
//     cds::bondAtomsByDistance(this->getAtoms());
// }

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbModel::Print(std::ostream& out) const
{
    for (auto& residue : this->getResidues())
    {
        residue->Print(out);
    }
}

void PdbModel::Write(std::ostream& stream) const
{
    for (auto& cdsMolecule : this->getMolecules())
    {
        PdbChain* pdbChain                  = codeUtils::erratic_cast<PdbChain*>(cdsMolecule);
        std::vector<cds::Residue*> residues = pdbChain->getResidues();
        for (auto& residue : residues)
        {
            PdbResidue* pdbResidue        = codeUtils::erratic_cast<PdbResidue*>(residue);
            std::vector<cds::Atom*> atoms = pdbResidue->getAtoms();
            std::vector<std::string> recordName;
            std::vector<double> occupancy;
            std::vector<double> temperatureFactor;
            for (auto& atom : atoms)
            {
                PdbAtom* pdbAtom = codeUtils::erratic_cast<PdbAtom*>(atom);
                recordName.push_back(pdbAtom->GetRecordName());
                occupancy.push_back(pdbAtom->GetOccupancy());
                temperatureFactor.push_back(pdbAtom->GetTemperatureFactor());
            }
            cds::PdbFileResidueData residueData {{codeUtils::indexVector(atoms)},
                                                 {pdbResidue->getNumber()},
                                                 {cds::truncatedResidueName(pdbResidue)},
                                                 {pdbResidue->getChainId()},
                                                 {pdbResidue->getInsertionCode()}};
            cds::PdbFileAtomData atomData = cds::toPdbFileAtomData(atoms, recordName, occupancy, temperatureFactor);
            cds::PdbFileFormat format;
            cds::PdbFileData writerData {format, {}, residueData, atomData};
            cds::writeAssemblyToPdb(stream, {{0}}, {{pdbResidue->HasTerCard()}}, {}, writerData);
        }
        if (!residues.empty())
        { // Sometimes you get empty chains after things have been deleted I guess.
            stream << "TER\n";
        }
    }
}
