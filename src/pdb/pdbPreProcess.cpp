#include "include/pdb/pdbPreProcess.hpp"

#include "include/CentralDataStructure/Selections/atomSelections.hpp"
#include "include/CentralDataStructure/Selections/residueSelections.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/graphInterface.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/assembly/assemblyGraph.hpp"
#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/geometry/measurements.hpp"
#include "include/pdb/amberMdPrep.hpp" //all preprocessing should move to here.
#include "include/pdb/bondByDistance.hpp"
#include "include/pdb/pdbChain.hpp"
#include "include/pdb/pdbData.hpp"
#include "include/pdb/pdbFunctions.hpp"
#include "include/pdb/pdbResidue.hpp"
#include "include/pdb/pdbResidueId.hpp"
#include "include/util/constants.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        namespace
        {
            std::string residueParmName(const PdbData& data, size_t residueId)
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

            void addResidueAtoms(
                PdbData& data, size_t residueId, const std::vector<std::pair<std::string, Coordinate>>& atoms)
            {
                for (auto& atom : atoms)
                {
                    const std::string& name = atom.first;
                    const Coordinate& coord = atom.second;
                    addAtom(data, residueId, name, coord);
                }
            }

            void addResidueBonds(
                PdbData& data,
                size_t residueId,
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
                    addBond(data, atomIndex(bond[0]), atomIndex(bond[1]));
                }
            }
        } // namespace

        PreprocessorInformation preProcess(
            PdbFile& file, const ParameterManager& parameterManager, PreprocessorOptions inputOptions)
        {
            util::log(__LINE__, __FILE__, util::INF, "Preprocesssing has begun");
            PreprocessorInformation ppInfo;
            PdbData& data = file.data;
            for (size_t assemblyId = 0; assemblyId < file.assemblies.size();
                 assemblyId++) // Now we do all, but maybe user can select at some point.
            {
                Assembly* assembly = &file.assemblies[assemblyId];
                preProcessCysResidues(data, assemblyId, ppInfo);
                preProcessHisResidues(data, assemblyId, ppInfo, inputOptions);
                preProcessChainTerminals(data, assemblyId, ppInfo, inputOptions);
                preProcessGapsUsingDistance(data, assemblyId, ppInfo, inputOptions);
                preProcessMissingUnrecognized(data, assemblyId, ppInfo, parameterManager);
                setAtomChargesForResidues(parameterManager, assembly->getResidues());
            }
            util::log(__LINE__, __FILE__, util::INF, "Preprocessing completed");
            return ppInfo;
        }

        void changeResidueName(
            PdbData& data, size_t assemblyId, const std::string& selector, const std::string& newName)
        {
            for (size_t residueId : assemblyResidues(data.indices, assemblyId))
            {
                std::size_t found = toString(pdbResidueId(data, residueId)).find(selector);
                if (found != std::string::npos)
                {
                    data.residues.names[residueId] = newName;
                    data.objects.residues[residueId]->setName(newName);
                    return;
                }
            }
            util::log(__LINE__, __FILE__, util::WAR, "Could not find residue to rename with this selector " + selector);
            return;
        }

        //////////////////////////////////////////////////////////
        //                      FUNCTIONS                       //
        //////////////////////////////////////////////////////////
        void preProcessCysResidues(PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo)
        {
            std::function<bool(const size_t&)> isCYS = [&](size_t n) {
                return util::contains({"CYS", "CYX"}, data.residues.names[n]);
            };
            util::log(__LINE__, __FILE__, util::INF, "Start CYS preprocessing for this Model\n");
            std::vector<size_t> assemblyResidues = assembly::assemblyResidues(data.indices, assemblyId);
            std::vector<size_t> cysResidues = util::vectorFilter(isCYS, assemblyResidues);
            if (cysResidues.empty())
            {
                util::log(__LINE__, __FILE__, util::INF, "No CYS or CYX residues detected in this structure\n");
            }
            for (size_t n = 0; n < cysResidues.size(); n++)
            { // I want to go through the list and compare from current item to end. Thus it2 = std::next it1
                size_t cysRes1Id = cysResidues[n];
                size_t sgAtom1Id = findResidueAtom(data, cysRes1Id, "SG");
                Atom* sgAtom1 = data.objects.atoms[sgAtom1Id];
                for (size_t k = n + 1; k < cysResidues.size(); k++)
                {
                    size_t cysRes2Id = cysResidues[k];
                    size_t sgAtom2Id = findResidueAtom(data, cysRes2Id, "SG");
                    Atom* sgAtom2 = data.objects.atoms[sgAtom2Id];
                    size_t atomCount = data.indices.atomCount;
                    if ((sgAtom1Id < atomCount) && (sgAtom2Id < atomCount))
                    {
                        double dist = distance(sgAtom1->coordinate(), sgAtom2->coordinate());
                        if (dist < constants::dSulfurCutoff && dist > 0.001)
                        {
                            data.residues.names[cysRes1Id] = "CYX";
                            data.residues.names[cysRes2Id] = "CYX";
                            data.objects.residues[cysRes1Id]->setName("CYX");
                            data.objects.residues[cysRes2Id]->setName("CYX");
                            ResidueId Res1Id = pdbResidueId(data, cysRes1Id);
                            ResidueId Res2Id = pdbResidueId(data, cysRes2Id);
                            addBond(data, sgAtom1Id, sgAtom2Id); // I think I want this here. Not 100%.
                            ppInfo.cysBondResidues.push_back({Res1Id, Res2Id, dist});
                            std::stringstream message;
                            message << "Bonding " << toString(Res1Id) << " and " << toString(Res2Id)
                                    << " with distance " << dist;
                            util::log(__LINE__, __FILE__, util::INF, message.str());
                        }
                    }
                }
            }
            return;
        }

        void preProcessHisResidues(
            PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const PreprocessorOptions& inputOptions)
        {
            // HIS protonation, user specified:
            util::log(__LINE__, __FILE__, util::INF, "User His protonation");
            for (auto& userSelectionPair : inputOptions.hisSelections)
            {
                changeResidueName(data, assemblyId, userSelectionPair.first, userSelectionPair.second);
            }
            util::log(__LINE__, __FILE__, util::INF, "Auto His protonation");
            // HIS protonation, automatic handling.
            for (size_t residueId : assemblyResidues(data.indices, assemblyId))
            {
                const std::string& name = data.residues.names[residueId];
                if (util::contains({"HIE", "HID", "HIP"}, name))
                {
                    ppInfo.hisResidues.emplace_back(pdbResidueId(data, residueId));
                }
                else if (name == "HIS")
                {
                    size_t atomCount = data.indices.atomCount;
                    size_t atomHE2 = findResidueAtom(data, residueId, "HE2");
                    size_t atomHD1 = findResidueAtom(data, residueId, "HD1");
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
                    std::string newName = newResidueName();
                    data.residues.names[residueId] = newName;
                    data.objects.residues[residueId]->setName(newName);
                    util::log(__LINE__, __FILE__, util::INF, "About to emplaceBack Id");
                    util::log(__LINE__, __FILE__, util::INF, toString(pdbResidueId(data, residueId)));
                    ppInfo.hisResidues.emplace_back(pdbResidueId(data, residueId));
                }
            }
            return;
        }

        void modifyTerminal(PdbData& data, size_t residueId, const std::string& type)
        {
            if (type == "NH3+") // For now, leaving it to tleap to add the correct H's
            {
                util::log(
                    __LINE__,
                    __FILE__,
                    util::INF,
                    "Modifying N Terminal of : " + toString(pdbResidueId(data, residueId)));
                size_t atomId = findResidueAtom(data, residueId, "H");
                deleteAtom(data, residueId, atomId);
            }
            else if (type == "CO2-")
            {
                size_t atomCount = data.indices.atomCount;
                util::log(
                    __LINE__,
                    __FILE__,
                    util::INF,
                    "Modifying C Terminal of : " + toString(pdbResidueId(data, residueId)));
                size_t atomOXT = findResidueAtom(data, residueId, "OXT");
                if (atomOXT < atomCount)
                {
                    util::log(
                        __LINE__,
                        __FILE__,
                        util::INF,
                        "OXT atom already exists: " + data.atoms.names[atomOXT] + "_" +
                            std::to_string(data.atoms.numbers[atomOXT]));
                    return;
                }
                // I don't like this, but at least it's somewhat contained:
                size_t atomCA = findResidueAtom(data, residueId, "CA");
                size_t atomC = findResidueAtom(data, residueId, "C");
                size_t atomO = findResidueAtom(data, residueId, "O");
                if (atomCA == atomCount || atomC == atomCount || atomO == atomCount)
                {
                    util::log(
                        __LINE__,
                        __FILE__,
                        util::WAR,
                        "Cterminal residue missing an atoms named CA, C or O, cannot create an OXT atom for this "
                        "residue.");
                    return;
                }
                auto coord = [&](size_t n) { return data.atoms.coordinates[n]; };
                Coordinate oxtCoord = calculateCoordinateFromInternalCoords(
                    coord(atomCA), coord(atomC), coord(atomO), 120.0, 180.0, 1.25);
                size_t oxtAtom = addAtom(data, residueId, "OXT", oxtCoord);
                addBond(data, oxtAtom, atomC);
                util::log(
                    __LINE__,
                    __FILE__,
                    util::INF,
                    "Created new atom named OXT after " + data.atoms.names[atomO] + "_" +
                        std::to_string(data.atoms.numbers[atomO]));
                return;
            }
            util::log(__LINE__, __FILE__, util::WAR, "Cannot handle this type of terminal option: " + type);
            return;
        }

        void preProcessChainTerminals(
            PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const PreprocessorOptions& inputOptions)
        {
            util::log(__LINE__, __FILE__, util::INF, "Chain terminations");
            for (size_t moleculeId : assemblyMolecules(data.indices, assemblyId))
            {
                util::log(__LINE__, __FILE__, util::INF, "Chain termination processing started for this chain");
                // Do the thing
                size_t nTerResidue = getNTerminal(data, moleculeId);
                if (nTerResidue >= data.indices.residueCount)
                {
                    util::log(__LINE__, __FILE__, util::INF, "Could not modify terminals of this chain.");
                }
                else
                {
                    modifyTerminal(data, nTerResidue, inputOptions.chainNTermination);
                    size_t cTerResidue = getCTerminal(data, moleculeId);
                    modifyTerminal(data, cTerResidue, inputOptions.chainCTermination);
                    util::log(__LINE__, __FILE__, util::INF, "N term : " + toString(pdbResidueId(data, nTerResidue)));
                    util::log(__LINE__, __FILE__, util::INF, "C term : " + toString(pdbResidueId(data, cTerResidue)));
                    // Report the thing
                    ppInfo.chainTerminals.push_back(
                        {data.residues.chainIds[nTerResidue],
                         getNumberAndInsertionCode(data, nTerResidue),
                         getNumberAndInsertionCode(data, cTerResidue),
                         inputOptions.chainNTermination,
                         inputOptions.chainCTermination});
                }
                util::log(__LINE__, __FILE__, util::INF, "Preprocessing complete for this chain");
            }
            util::log(__LINE__, __FILE__, util::INF, "Chain termination processing complete");
        }

        void insertCap(PdbData& data, size_t moleculeId, size_t refResidueId, const std::string& type)
        {
            // This approach is bad. When parameter manager is good we can use that to remove the get_carestian stuff
            auto coord = [&](size_t n) { return data.atoms.coordinates[n]; };
            if (type == "NHCH3") // NME
            {
                //        NME resid numbers. Otherwise good.
                Coordinate cCoordProtein = coord(findResidueAtom(data, refResidueId, "C"));
                Coordinate caCoordProtein = coord(findResidueAtom(data, refResidueId, "CA"));
                Coordinate oCoordProtein = coord(findResidueAtom(data, refResidueId, "O"));
                Coordinate nCoordNME = calculateCoordinateFromInternalCoords(
                    oCoordProtein, caCoordProtein, cCoordProtein, 120.0, 180.0, 1.4);
                Coordinate hCoordNME =
                    calculateCoordinateFromInternalCoords(oCoordProtein, caCoordProtein, nCoordNME, 109.0, 180.0, 1.0);
                Coordinate ch3CoordNME =
                    calculateCoordinateFromInternalCoords(caCoordProtein, cCoordProtein, nCoordNME, 125.0, 180.0, 1.48);
                Coordinate hh31CoordNME =
                    calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 180.0, 1.09);
                Coordinate hh32CoordNME =
                    calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 60.0, 1.09);
                Coordinate hh33CoordNME =
                    calculateCoordinateFromInternalCoords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, -60.0, 1.09);
                size_t residueId =
                    readResidue(data, moleculeId, "NME", ResidueType::ProteinCappingGroup, true, refResidueId);
                addResidueAtoms(
                    data,
                    residueId,
                    {
                        {   "N",    nCoordNME},
                        {   "H",    hCoordNME},
                        { "CH3",  ch3CoordNME},
                        {"HH31", hh31CoordNME},
                        {"HH32", hh32CoordNME},
                        {"HH33", hh33CoordNME}
                });
                addResidueBonds(
                    data,
                    residueId,
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
                Coordinate cCoordProtein = coord(findResidueAtom(data, refResidueId, "C"));
                Coordinate caCoordProtein = coord(findResidueAtom(data, refResidueId, "CA"));
                Coordinate nCoordProtein = coord(findResidueAtom(data, refResidueId, "N"));
                // This is bad, should use templates loaded from lib/prep file instead.
                Coordinate cCoordACE = calculateCoordinateFromInternalCoords(
                    cCoordProtein, caCoordProtein, nCoordProtein, 120.0, -130.0, 1.4);
                Coordinate oCoordACE =
                    calculateCoordinateFromInternalCoords(caCoordProtein, nCoordProtein, cCoordACE, 120.0, 0.0, 1.23);
                Coordinate ch3CoordACE =
                    calculateCoordinateFromInternalCoords(caCoordProtein, nCoordProtein, cCoordACE, 125.0, 180.0, 1.48);
                Coordinate hh31CoordACE =
                    calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
                Coordinate hh32CoordACE =
                    calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
                Coordinate hh33CoordACE =
                    calculateCoordinateFromInternalCoords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
                // Ok this next bit is convoluted, but I look up the position of the first atom in the protein residue
                // and insert the new Atom before it, and get passed back the position of the newly created atom, so I
                // can use that when creating the next one and so on. With ACE we want to insert before the residue, so
                // I'm finding the residue before here:
                const std::vector<size_t>& order = data.molecules.residueOrder[moleculeId];
                size_t position = util::indexOf(order, refResidueId) - 1;
                size_t residueId =
                    readResidue(data, moleculeId, "ACE", ResidueType::ProteinCappingGroup, false, order[position]);
                addResidueAtoms(
                    data,
                    residueId,
                    {
                        {   "C",    cCoordACE},
                        {   "O",    oCoordACE},
                        { "CH3",  ch3CoordACE},
                        {"HH31", hh31CoordACE},
                        {"HH32", hh32CoordACE},
                        {"HH33", hh33CoordACE}
                });
                addResidueBonds(
                    data,
                    residueId,
                    {
                        {"C", findResidueAtom(data, refResidueId, "N")},
                        {"C", "O"},
                        {"C", "CH3"},
                        {"CH3", "HH31"},
                        {"CH3", "HH32"},
                        {"CH3", "HH33"}
                });
                util::log(
                    __LINE__, __FILE__, util::INF, "Created ACE residue: " + toString(pdbResidueId(data, residueId)));
            }
        }

        void preProcessGapsUsingDistance(
            PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const PreprocessorOptions& inputOptions)
        {
            // Missing Residues (gaps); If two sequential protein residues in the same molecule aren't close enough to
            // bond: this is a gap regardless of residue number/insertion code. User will want caps(ACE/NME) or
            // zwitterionic, we can't know ourselves without knowledge of the system, but most of the time caps.
            util::log(__LINE__, __FILE__, util::INF, "Gaps");
            for (size_t moleculeId : assemblyMolecules(data.indices, assemblyId))
            {
                util::log(
                    __LINE__,
                    __FILE__,
                    util::INF,
                    "Gap detection started for chain " + data.molecules.chainIds[moleculeId]);
                std::vector<size_t> proteinResidues = util::boolsToIndices(util::vectorAnd(
                    util::vectorEquals(data.residues.types, ResidueType::Protein),
                    isMoleculeResidue(data.indices, moleculeId)));
                if (proteinResidues.empty())
                {
                    util::log(
                        __LINE__,
                        __FILE__,
                        util::INF,
                        "No protein residues found in chain with id: " + data.molecules.chainIds[moleculeId]);
                    break;
                }
                for (size_t n = 0; n < proteinResidues.size() - 1; n++)
                {
                    size_t res1 = proteinResidues[n];
                    size_t res2 = proteinResidues[n + 1];
                    size_t atomCount = data.indices.atomCount;
                    size_t res1AtomC = findResidueAtom(data, res1, "C");
                    size_t res2AtomN = findResidueAtom(data, res2, "N");
                    if ((res1AtomC < atomCount) && (res2AtomN < atomCount) &&
                        (!isWithinBondingDistance(data, res1AtomC, res2AtomN)))
                    { // GAP detected
                        // Look for non-natural protein residues within bonding distance, they fall under ResidueType
                        // Undefined, this indicates it's not gap.
                        if (!checkForNonNaturalProteinResidues(
                                data,
                                util::indicesOfElement(data.residues.types, ResidueType::Undefined),
                                res1AtomC,
                                ppInfo))
                        {
                            // Log it
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::INF,
                                inputOptions.gapCTermination + " cap for : " + toString(pdbResidueId(data, res1)));
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::INF,
                                inputOptions.gapNTermination + " cap for : " + toString(pdbResidueId(data, res2)));
                            // Do it
                            insertCap(data, moleculeId, res1, inputOptions.gapCTermination);
                            insertCap(data, moleculeId, res2, inputOptions.gapNTermination);
                            // Record it
                            ppInfo.missingResidues.push_back(
                                {data.residues.chainIds[res1],
                                 getNumberAndInsertionCode(data, res1),
                                 getNumberAndInsertionCode(data, res2),
                                 inputOptions.gapCTermination,
                                 inputOptions.gapNTermination});
                        }
                    }
                }
                util::log(
                    __LINE__,
                    __FILE__,
                    util::INF,
                    "Gap detection completed for chain " + data.molecules.chainIds[moleculeId]);
            }
            return;
        }

        void preProcessMissingUnrecognized(
            PdbData& data, size_t assemblyId, PreprocessorInformation& ppInfo, const ParameterManager& parmManager)
        {
            for (size_t residueId : assemblyResidues(data.indices, assemblyId))
            {
                std::string parmName = residueParmName(data, residueId);
                ResidueId residueIdObj = pdbResidueId(data, residueId);
                size_t index = util::indexOf(parmManager.lib.residueNames, parmName);
                // Unrecognized residue->
                if (index == parmManager.lib.residueNames.size())
                {
                    util::log(__LINE__, __FILE__, util::INF, "ParmManager did not recognize residue: " + parmName);
                    ppInfo.unrecognizedResidues.emplace_back(residueIdObj);
                }
                else // Recognized residue->
                {
                    const lib::ResidueData& parmResidue = parmManager.lib.residues[index];
                    std::vector<size_t> parmHeavyAtoms = FindHeavyAtoms(parmResidue.atoms.elements);
                    std::vector<std::string> parmHeavyAtomNames =
                        util::indicesToValues(parmResidue.atoms.names, parmHeavyAtoms);
                    std::vector<size_t> atomIds = residueAtoms(data.indices, residueId);
                    std::vector<std::string> pdbAtomNames = util::indicesToValues(data.atoms.names, atomIds);
                    for (auto& parmHeavyAtomName :
                         parmHeavyAtomNames) // What heavy atoms are missing from the pdb residue?
                    {
                        if (!util::contains(pdbAtomNames, parmHeavyAtomName))
                        { // Residue missing a heavy atom.
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::INF,
                                "Atom named " + parmHeavyAtomName + " missing from " + toString(residueIdObj));
                            ppInfo.missingHeavyAtoms.push_back({parmHeavyAtomName, residueIdObj});
                        }
                    }
                    for (auto& pdbAtomName : pdbAtomNames) // What atoms in the pdb residue are unrecognized?
                    {
                        if (!util::contains(parmResidue.atoms.names, pdbAtomName))
                        {
                            // Residue contains unrecognized atom.
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::INF,
                                "Unrecognized atom named " + pdbAtomName + " in " + toString(residueIdObj));
                            ppInfo.unrecognizedAtoms.push_back({pdbAtomName, residueIdObj});
                        }
                    }
                }
            }
            return;
        }
    } // namespace pdb
} // namespace gmml
