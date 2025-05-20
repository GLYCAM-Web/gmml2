#include "includes/CentralDataStructure/Readers/Pdb/pdbPreProcess.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbChain.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbData.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/bondByDistance.hpp"
#include "includes/CentralDataStructure/Geometry/geometryFunctions.hpp"
#include "includes/CodeUtils/casting.hpp"
#include "includes/CodeUtils/constants.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp"
#include "includes/CentralDataStructure/Selections/residueSelections.hpp"
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Editors/amberMdPrep.hpp" //all preprocessing should move to here.
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/cdsFunctions/graphInterface.hpp"
#include "includes/Assembly/assemblyGraph.hpp"

namespace
{
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

    void addResidueAtoms(pdb::PdbData& data, size_t residueId,
                         const std::vector<std::pair<std::string, cds::Coordinate>>& atoms)
    {
        for (auto& atom : atoms)
        {
            const std::string& name      = atom.first;
            const cds::Coordinate& coord = atom.second;
            pdb::addAtom(data, residueId, name, coord);
        }
    }

    void addResidueBonds(pdb::PdbData& data, size_t residueId,
                         const std::vector<std::array<std::variant<size_t, std::string>, 2>>& bonds)
    {
        auto atomIndex = [&](const std::variant<size_t, std::string>& key)
        {
            if (std::holds_alternative<size_t>(key))
            {
                return std::get<size_t>(key);
            }
            else
            {
                return findResidueAtom(data, residueId, std::get<std::string>(key));
            }
        };
        for (auto& bond : bonds)
        {
            pdb::addBond(data, atomIndex(bond[0]), atomIndex(bond[1]));
        }
    }
} // namespace

void pdb::changeResidueName(PdbData& data, size_t assemblyId, const std::string& selector, const std::string& newName)
{
    for (size_t residueId : assemblyResidues(data.indices, assemblyId))
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
    std::vector<size_t> assemblyResidues = assembly::assemblyResidues(data.indices, assemblyId);
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
        changeResidueName(data, assemblyId, userSelectionPair.first, userSelectionPair.second);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Auto His protonation");
    // HIS protonation, automatic handling.
    for (size_t residueId : assemblyResidues(data.indices, assemblyId))
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

void pdb::modifyTerminal(PdbData& data, size_t residueId, const std::string& type)
{
    if (type == "NH3+") // For now, leaving it to tleap to add the correct H's
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying N Terminal of : " + pdbResidueId(data, residueId).print());
        size_t atomId = findResidueAtom(data, residueId, "H");
        deleteAtom(data, residueId, atomId);
    }
    else if (type == "CO2-")
    {
        size_t atomCount = data.indices.atomCount;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Modifying C Terminal of : " + pdbResidueId(data, residueId).print());
        size_t atomOXT = findResidueAtom(data, residueId, "OXT");
        if (atomOXT < atomCount)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF,
                      "OXT atom already exists: " + data.atoms.names[atomOXT] + "_" +
                          std::to_string(data.atoms.numbers[atomOXT]));
            return;
        }
        // I don't like this, but at least it's somewhat contained:
        size_t atomCA = findResidueAtom(data, residueId, "CA");
        size_t atomC  = findResidueAtom(data, residueId, "C");
        size_t atomO  = findResidueAtom(data, residueId, "O");
        if (atomCA == atomCount || atomC == atomCount || atomO == atomCount)
        {
            gmml::log(
                __LINE__, __FILE__, gmml::WAR,
                "Cterminal residue missing an atoms named CA, C or O, cannot create an OXT atom for this residue.");
            return;
        }
        auto coord = [&](size_t n)
        {
            return data.atoms.coordinates[n];
        };
        cds::Coordinate oxtCoord =
            cds::calculateCoordinateFromInternalCoords(coord(atomCA), coord(atomC), coord(atomO), 120.0, 180.0, 1.25);
        size_t oxtAtom = addAtom(data, residueId, "OXT", oxtCoord);
        addBond(data, oxtAtom, atomC);
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Created new atom named OXT after " + data.atoms.names[atomO] + "_" +
                      std::to_string(data.atoms.numbers[atomO]));
        return;
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Cannot handle this type of terminal option: " + type);
    return;
}

void pdb::preProcessChainTerminals(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo,
                                   const PreprocessorOptions& inputOptions)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Chain terminations");
    for (size_t moleculeId : assemblyMolecules(data.indices, assemblyId))
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
            modifyTerminal(data, nTerResidue, inputOptions.chainNTermination_);
            size_t cTerResidue = getCTerminal(data, moleculeId);
            modifyTerminal(data, cTerResidue, inputOptions.chainCTermination_);
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

void pdb::insertCap(PdbData& data, size_t moleculeId, size_t refResidueId, const std::string& type)
{
    // This approach is bad. When parameter manager is good we can use that to remove the get_carestian stuff
    using cds::Coordinate;
    auto coord = [&](size_t n)
    {
        return data.atoms.coordinates[n];
    };
    if (type == "NHCH3") // NME
    {
        //        NME resid numbers. Otherwise good.
        Coordinate cCoordProtein  = coord(findResidueAtom(data, refResidueId, "C"));
        Coordinate caCoordProtein = coord(findResidueAtom(data, refResidueId, "CA"));
        Coordinate oCoordProtein  = coord(findResidueAtom(data, refResidueId, "O"));
        Coordinate nCoordNME =
            cds::calculateCoordinateFromInternalCoords(oCoordProtein, caCoordProtein, cCoordProtein, 120.0, 180.0, 1.4);
        Coordinate hCoordNME =
            cds::calculateCoordinateFromInternalCoords(oCoordProtein, caCoordProtein, nCoordNME, 109.0, 180.0, 1.0);
        Coordinate ch3CoordNME =
            cds::calculateCoordinateFromInternalCoords(caCoordProtein, cCoordProtein, nCoordNME, 125.0, 180.0, 1.48);
        Coordinate hh31CoordNME =
            cds::calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 180.0, 1.09);
        Coordinate hh32CoordNME =
            cds::calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 60.0, 1.09);
        Coordinate hh33CoordNME =
            cds::calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, -60.0, 1.09);
        size_t residueId =
            readResidue(data, moleculeId, "NME", cds::ResidueType::ProteinCappingGroup, true, refResidueId);
        addResidueAtoms(data, residueId,
                        {
                            {   "N",    nCoordNME},
                            {   "H",    hCoordNME},
                            { "CH3",  ch3CoordNME},
                            {"HH31", hh31CoordNME},
                            {"HH32", hh32CoordNME},
                            {"HH33", hh33CoordNME}
        });
        addResidueBonds(data, residueId,
                        {
                            {"N", findResidueAtom(data, refResidueId, "C")},
                            {"N", "H"},
                            {"N", "CH3"},
                            {"CH3", "HH31"},
                            {"CH3", "HH32"},
                            {"CH3", "HH33"}
        });
    }
    else if (type == "COCH3") // ACE
    {
        //        NME resid numbers. Otherwise good.
        // These are the atoms in residue that I use to build the ACE out from.
        Coordinate cCoordProtein  = coord(findResidueAtom(data, refResidueId, "C"));
        Coordinate caCoordProtein = coord(findResidueAtom(data, refResidueId, "CA"));
        Coordinate nCoordProtein  = coord(findResidueAtom(data, refResidueId, "N"));
        // This is bad, should use templates loaded from lib/prep file instead.
        Coordinate cCoordACE = cds::calculateCoordinateFromInternalCoords(cCoordProtein, caCoordProtein, nCoordProtein,
                                                                          120.0, -130.0, 1.4);
        Coordinate oCoordACE =
            cds::calculateCoordinateFromInternalCoords(caCoordProtein, nCoordProtein, cCoordACE, 120.0, 0.0, 1.23);
        Coordinate ch3CoordACE =
            cds::calculateCoordinateFromInternalCoords(caCoordProtein, nCoordProtein, cCoordACE, 125.0, 180.0, 1.48);
        Coordinate hh31CoordACE =
            cds::calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
        Coordinate hh32CoordACE =
            cds::calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
        Coordinate hh33CoordACE =
            cds::calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
        // Ok this next bit is convoluted, but I look up the position of the first atom in the protein residue and
        // insert the new Atom before it, and get passed back the position of the newly created atom, so I can use that
        // when creating the next one and so on. With ACE we want to insert before the residue, so I'm finding the
        // residue before here:
        const std::vector<size_t>& order = data.molecules.residueOrder[moleculeId];
        size_t position                  = codeUtils::indexOf(order, refResidueId) - 1;
        size_t residueId =
            readResidue(data, moleculeId, "ACE", cds::ResidueType::ProteinCappingGroup, false, order[position]);
        addResidueAtoms(data, residueId,
                        {
                            {   "C",    cCoordACE},
                            {   "O",    oCoordACE},
                            { "CH3",  ch3CoordACE},
                            {"HH31", hh31CoordACE},
                            {"HH32", hh32CoordACE},
                            {"HH33", hh33CoordACE}
        });
        addResidueBonds(data, residueId,
                        {
                            {"C", findResidueAtom(data, refResidueId, "N")},
                            {"C", "O"},
                            {"C", "CH3"},
                            {"CH3", "HH31"},
                            {"CH3", "HH32"},
                            {"CH3", "HH33"}
        });
        gmml::log(__LINE__, __FILE__, gmml::INF, "Created ACE residue: " + pdbResidueId(data, residueId).print());
    }
}

void pdb::preProcessGapsUsingDistance(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo,
                                      const PreprocessorOptions& inputOptions)
{
    // Missing Residues (gaps); If two sequential protein residues in the same molecule aren't close enough to bond:
    // this is a gap regardless of residue number/insertion code. User will want caps(ACE/NME) or zwitterionic, we can't
    // know ourselves without knowledge of the system, but most of the time caps.
    gmml::log(__LINE__, __FILE__, gmml::INF, "Gaps");
    for (size_t moleculeId : assemblyMolecules(data.indices, assemblyId))
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Gap detection started for chain " + data.molecules.chainIds[moleculeId]);
        std::vector<size_t> proteinResidues = codeUtils::boolsToIndices(
            codeUtils::vectorAnd(codeUtils::vectorEquals(data.residues.types, cds::ResidueType::Protein),
                                 isMoleculeResidue(data.indices, moleculeId)));
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
                    insertCap(data, moleculeId, res1, inputOptions.gapCTermination_);
                    insertCap(data, moleculeId, res2, inputOptions.gapNTermination_);
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
    for (size_t residueId : assemblyResidues(data.indices, assemblyId))
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
            std::vector<size_t> atomIds           = residueAtoms(data.indices, residueId);
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