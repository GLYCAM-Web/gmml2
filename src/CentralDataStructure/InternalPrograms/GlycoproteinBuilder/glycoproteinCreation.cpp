#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinCreation.hpp"
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkageCreation.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/containerTypes.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <string>
#include <functional>
#include <sstream>
#include <fstream>
#include <stdexcept>

using cds::Atom;
using cds::Residue;
using cdsCondensedSequence::Carbohydrate;

namespace
{
    struct ChainAndResidue
    {
        std::string chain;
        std::string residue;
    };

    struct Attachment
    {
        Carbohydrate* glycan;
        Residue* aglycone;
        Residue* reducing;
    };

    ChainAndResidue selectionChainAndResidue(const std::string userSelection)
    { // Chain_residueNumber_insertionCode* *optional.
        std::vector<std::string> split = codeUtils::split(userSelection, '_');
        size_t splitCount              = split.size();
        if (splitCount == 2)
        {
            return {split[0], split[1]};
        }
        else if (splitCount == 3)
        {
            return {split[0], split[1] + split[2]};
        }
        else
        {
            throw std::runtime_error(
                "userSelection (" + userSelection +
                ") for residue to glycosylate is incorrect format.\nMust be "
                "chain_residueNumber_insertionCode.\nInsertionCode is optional. Chain can be ? if no "
                "chain numbers are in input.\nExamples: ?_24_? or ?_24 will use the first residue it "
                "encounters numbered 24. A_24_B is A chain, residue 24, insertion code B");
        }
    }

    size_t selectResidueFromInput(const std::vector<std::string>& residueIds, const std::vector<std::string>& chainIds,
                                  const std::vector<std::string>& numberAndInsertionCode,
                                  const ChainAndResidue& selection)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "We working with " + selection.chain + "_" + selection.residue);
        for (size_t n = 0; n < residueIds.size(); n++)
        {
            if ((chainIds[n] == selection.chain) && (numberAndInsertionCode[n] == selection.residue))
            {
                gmml::log(__LINE__, __FILE__, gmml::INF, "Id of selected glycosite: " + residueIds[n]);
                return n;
            }
        }
        return residueIds.size();
    }

    // Step 1. If the C1 atom has a neighbor that isn't in the queryResidue, return C1.
    // Step 2. If the C2 atom has a neighbor that isn't in the queryResidue, return C2.
    // Step 3. Panic.
    Atom* guessAnomericAtomByForeignNeighbor(cds::Residue* queryResidue)
    {
        std::vector<Atom*> atoms               = queryResidue->getAtoms();
        std::vector<std::string> atomNames     = cds::atomNames(atoms);
        std::vector<std::string> usualSuspects = {"C1", "C2"};
        for (auto& suspectName : usualSuspects)
        {
            size_t potentialAnomer = codeUtils::indexOf(atomNames, suspectName);
            if (potentialAnomer < atomNames.size() &&
                cdsSelections::selectNeighborNotInAtomVector(atoms[potentialAnomer], atoms) != nullptr)
            { // If atom has a foreign neighbor.
                return atoms[potentialAnomer];
            }
        }
        std::string message =
            "Did not find a C1 or C2 with a foreign neighbor in residue: " + cds::residueStringId(queryResidue) +
            ", thus no anomeric atom was found.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
        return nullptr;
    }

    Attachment toAttachment(Carbohydrate* glycan)
    { // Tagging during construction would be better. Ano-ano linkages won't work,
        std::vector<cds::Residue*> residues = glycan->getResidues();
        std::vector<cds::ResidueType> types = cds::residueTypes(residues);
        size_t foundAglycone                = codeUtils::indexOf(types, cds::ResidueType::Aglycone);
        size_t aglycone                     = foundAglycone == residues.size() ? 0 : foundAglycone;
        size_t foundSugar                   = codeUtils::indexOf(types, cds::ResidueType::Sugar);
        size_t reducing                     = foundSugar == residues.size() ? 1 : foundSugar;
        if (residues.size() >= 2)
        {
            return {glycan, residues[aglycone], residues[reducing]};
        }
        else
        {
            std::string message = "Reducing residue and aglycone requested for Carbohydrate with name " +
                                  glycan->getName() + ", but it doesn't have more than 1 residue";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
    }

    size_t glycosylationTableIndex(const glycoproteinMetadata::GlycosylationTable& table, const std::string& aminoAcid,
                                   const std::string& proteinResidueId)
    {
        size_t index = codeUtils::indexOf(table.residueNames, aminoAcid);
        if (index == table.residueNames.size())
        {
            std::string message = "Problem creating glycosylation site. The amino acid requested: " + proteinResidueId +
                                  " has name (" + aminoAcid +
                                  ") that isn't supported. Currently you can glycosylate ASN, THR, SER or TYR. "
                                  "Email us to request " +
                                  "others. Ideally include examples of 3D structures we can use as a template.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        return index;
    }

    // This function prepares the glycan molecule in the glycan_ assembly for superimpostion onto an amino acid in the
    // protein It does this by "growing" the atoms of the amino acid side chain (e.g. Asn, Thr or Ser) out from the
    // glycan reducing terminal Another function will use these additional atoms to superimpose the glycan onto residue
    void prepareGlycansForSuperimpositionToParticularResidue(const glycoproteinMetadata::GlycosylationTable& table,
                                                             size_t index, Attachment& attachment)
    {
        // Want: residue.FindAtomByTag("anomeric-carbon"); The below is risky as it uses atoms names, i.e. would break
        // for Sialic acid. ToDo Ok so I reckon the below is just assuming alpha or beta depending on the concext. Need
        // to fix a lot, but need to reproduce functionality after refactor first.
        //    Atom* anomericAtom = cdsSelections::guessAnomericAtom(reducing_Residue);
        //   This won't work as sometimes want alpha, sometimes beta. i.e. a coordinateOppositeToNeighborAverage
        //   function
        // This needs to be abstracted so it works for C2 reducing residues:
        // Delete aglycon atoms from glycan.
        Residue* aglycone = attachment.aglycone;
        for (auto& atom : aglycone->getAtoms())
        {
            aglycone->deleteAtom(atom);
        }
        aglycone->setName("SUP");

        const std::vector<std::string>& names                                  = table.atomNames[index];
        const std::vector<glycoproteinMetadata::SuperimpositionValues>& values = table.values[index];
        Residue* reducing                                                      = attachment.reducing;
        Atom* anomericAtom                                                     = reducing->FindAtom("C1");
        Coordinate coordC5                                                     = reducing->FindAtom("C5")->coordinate();
        Coordinate coordO5                                                     = reducing->FindAtom("O5")->coordinate();
        Coordinate coordC1                                                     = anomericAtom->coordinate();
        std::vector<Coordinate> coords {coordC5, coordO5, coordC1};
        for (size_t n = 0; n < values.size(); n++)
        {
            coords.push_back(cds::calculateCoordinateFromInternalCoords(
                coords[n], coords[n + 1], coords[n + 2], values[n].angle, values[n].dihedral, values[n].distance));
        }
        std::vector<cds::Atom*> atoms;
        atoms.reserve(names.size());
        for (size_t n = 0; n < names.size(); n++)
        {
            atoms.push_back(aglycone->addAtom(std::make_unique<Atom>(names[n], coords[n + 3])));
        }
        cds::addBond(anomericAtom, atoms[0]);
    }

    void superimposeGlycanToGlycosite(const pdb::PdbData& pdbData, Residue* glycosite, Attachment& attachment)
    {
        // superimposition_atoms_ points to three atoms that were added to the glycan. Based on their names e.g. CG,
        // ND2, we will superimpose them onto the correspoinding "target" atoms in the protein residue
        // (glycosite_residue).
        size_t glycositeId                         = codeUtils::indexOf(pdbData.indices.residues, glycosite);
        std::string residueId                      = pdb::residueStringId(pdbData, glycositeId);
        std::vector<cds::Atom*> proteinAtoms       = glycosite->getAtoms();
        std::vector<std::string> proteinAtomNames  = cds::atomNames(proteinAtoms);
        std::vector<cds::Atom*> aglyconeAtoms      = attachment.aglycone->mutableAtoms();
        std::vector<std::string> aglyconeAtomNames = cds::atomNames(aglyconeAtoms);
        // Sanity check
        if (aglyconeAtoms.size() < 3)
        {
            throw std::runtime_error("The aglycone does not contain enough atoms to perform the requested "
                                     "superimposition to " +
                                     residueId + ".\nCheck your input structure!\n");
        }
        // Get the 3 target atoms from protein residue.
        std::vector<cds::Coordinate> targetCoords;
        for (auto name : aglyconeAtomNames)
        {
            size_t index = codeUtils::indexOf(proteinAtomNames, name);
            if (index == proteinAtoms.size())
            {
                throw std::runtime_error("Error: atom '" + name + "' not found in residue " + residueId +
                                         " with atoms: " + codeUtils::join(", ", proteinAtomNames));
            }
            targetCoords.push_back(proteinAtoms[index]->coordinate());
        }
        std::vector<cds::Atom*> glycanAtoms          = attachment.glycan->mutableAtoms();
        std::vector<cds::Coordinate> aglyconeCoords  = cds::atomCoordinates(aglyconeAtoms);
        std::vector<cds::Coordinate> glycanCoords    = cds::atomCoordinates(glycanAtoms);
        cds::AffineTransform transform               = cds::affineTransform(targetCoords, aglyconeCoords);
        std::vector<cds::Coordinate> updatedAglycone = cds::matrixCoordinates(transform.affine * transform.moving);
        std::vector<cds::Coordinate> updatedGlycan =
            cds::matrixCoordinates(transform.affine * cds::generateMatrix(glycanCoords));
        for (size_t n = 0; n < aglyconeAtoms.size(); n++)
        {
            aglyconeAtoms[n]->setCoordinate(updatedAglycone[n]);
        }
        for (size_t n = 0; n < glycanCoords.size(); n++)
        {
            glycanAtoms[n]->setCoordinate(updatedGlycan[n]);
        }
    }

    std::vector<size_t> linkagesContainingResidue(const std::vector<cds::ResidueLinkage>& glycosidicLinkages,
                                                  cds::Residue* residue)
    {
        std::function<bool(const size_t&)> contains = [&](const size_t& n)
        {
            auto& linkageResidues = glycosidicLinkages[n].link.residues;
            return (linkageResidues.first == residue || linkageResidues.second == residue);
        };
        return codeUtils::vectorFilter(contains, codeUtils::indexVector(glycosidicLinkages));
    }

    cds::ResidueLinkage replaceAglycone(const GlycamMetadata::DihedralAngleDataTable& metadataTable,
                                        cds::Residue* reducingResidue, cds::Residue* newAglycone)
    {
        // Old aglycone atoms are connected to reducing residue, are found during creation of rotatable dihedrals.
        newAglycone->addNeighbor(newAglycone->getName() + "-" + reducingResidue->getName(), reducingResidue);
        cds::ResidueLink link = cds::findResidueLink({reducingResidue, newAglycone});
        return cds::createResidueLinkage(metadataTable, link);
    }
} // namespace

namespace glycoproteinBuilder
{

    std::vector<GlycosylationSite> createGlycosites(const pdb::PdbData& pdbData, cds::Assembly* glycoprotein,
                                                    const std::vector<GlycositeInput>& glycositesInputVector)
    {
        std::vector<GlycosylationSite> result;
        std::vector<cds::Residue*> residues = glycoprotein->getResidues();
        std::vector<std::string> numberAndInsertionCodes;
        std::vector<std::string> chainIds;
        std::vector<std::string> residueIds;
        for (auto residue : residues)
        {
            size_t residueId = codeUtils::indexOf(pdbData.indices.residues, residue);
            numberAndInsertionCodes.push_back(pdb::getNumberAndInsertionCode(pdbData, residueId));
            chainIds.push_back(pdbData.residues.chainIds[residueId]);
            residueIds.push_back(pdb::residueStringId(pdbData, residueId));
        }
        for (auto& glycositeInput : glycositesInputVector)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                          glycositeInput.glycanInput);
            ChainAndResidue selection = selectionChainAndResidue(glycositeInput.proteinResidueId);
            size_t residueIndex = selectResidueFromInput(residueIds, chainIds, numberAndInsertionCodes, selection);
            if (residueIndex == residues.size())
            {
                throw std::runtime_error("Error: Did not find a residue with id matching " +
                                         glycositeInput.proteinResidueId + "\n");
            }
            cds::Residue* glycositeResidue = residues[residueIndex];
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "About to push to glycosites with: " + glycositeInput.proteinResidueId + " and glycan " +
                          glycositeInput.glycanInput);
            result.push_back({glycositeResidue, glycositeInput});
        }
        return result;
    }

    std::vector<cdsCondensedSequence::Carbohydrate*>
    addGlycansToProtein(const cdsParameters::ParameterManager& parameterManager,
                        const codeUtils::SparseVector<double>& elementRadii,
                        const GlycamMetadata::DihedralAngleDataTable& metadataTable, const pdb::PdbData& pdbData,
                        cds::Assembly* glycoprotein, const std::vector<GlycosylationSite>& glycosites)
    {
        const glycoproteinMetadata::GlycosylationTable glycosylationTable =
            glycoproteinMetadata::defaultGlycosylationTable();
        std::vector<Carbohydrate*> result;
        for (auto& glycosite : glycosites)
        {
            Carbohydrate* glycan = codeUtils::erratic_cast<Carbohydrate*>(glycoprotein->addMolecule(
                std::make_unique<Carbohydrate>(parameterManager, elementRadii, glycosite.input.glycanInput)));
            result.push_back(glycan);
            Attachment attachment                = toAttachment(glycan);
            cds::Residue* aglycone               = attachment.aglycone;
            cds::Residue* reducingResidue        = attachment.reducing;
            cds::Residue* glycositeResidue       = glycosite.residue;
            std::string glycositeResidueName     = glycositeResidue->getName();
            const GlycositeInput& glycositeInput = glycosite.input;
            size_t tableIndex =
                glycosylationTableIndex(glycosylationTable, glycositeResidueName, glycositeInput.proteinResidueId);
            prepareGlycansForSuperimpositionToParticularResidue(glycosylationTable, tableIndex, attachment);
            superimposeGlycanToGlycosite(pdbData, glycositeResidue, attachment);
            cds::addBond(glycositeResidue->FindAtom(glycosylationTable.connectingAtomNames[tableIndex]),
                         guessAnomericAtomByForeignNeighbor(attachment.reducing));
            glycositeResidue->setName(glycosylationTable.renamedResidues[tableIndex]);
            std::vector<cds::ResidueLinkage>& linkages = glycan->GetGlycosidicLinkages();
            std::vector<size_t> aglyconeLinkages       = linkagesContainingResidue(linkages, aglycone);
            if (aglyconeLinkages.size() != 1)
            {
                throw std::runtime_error("Error: found more than 1 linkage to aglycone in " +
                                         glycositeInput.glycanInput);
            }
            glycan->cds::Molecule::deleteResidue(aglycone);
            cds::ResidueLinkage newLinkage = replaceAglycone(metadataTable, reducingResidue, glycositeResidue);
            size_t aglyconeLinkage         = aglyconeLinkages[0];
            linkages[aglyconeLinkage]      = newLinkage;
            cds::setShapeToPreference(newLinkage, cds::defaultShapePreference(metadataTable, newLinkage.rotamerType,
                                                                              newLinkage.dihedralMetadata));
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "Completed creating glycosite on residue " + glycositeInput.proteinResidueId + " with glycan " +
                          glycositeInput.glycanInput);
        }
        return result;
    }
} // namespace glycoproteinBuilder