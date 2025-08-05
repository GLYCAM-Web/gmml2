#include "include/fileType/pdb/atomicConnectivity.hpp"

#include "include/assembly/assemblyIndices.hpp"
#include "include/assembly/assemblyTypes.hpp"
#include "include/fileType/pdb/bondByDistance.hpp"
#include "include/fileType/pdb/pdbData.hpp"
#include "include/fileType/pdb/pdbFunctions.hpp"
#include "include/fileType/pdb/pdbResidue.hpp"
#include "include/geometry/geometryFunctions.hpp"
#include "include/graph/graphManipulation.hpp"
#include "include/metadata/aminoAcids.hpp"
#include "include/metadata/atomicBonds.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <stdexcept>
#include <vector>

namespace gmml
{
    namespace pdb
    {
        namespace
        {
            void addBonds(
                PdbData& data, size_t proteinRes, const std::vector<std::pair<std::string, std::string>>& bonds)
            {
                size_t atomCount = data.indices.atomCount;
                for (auto& bondPair : bonds)
                {
                    size_t firstAtom = findResidueAtom(data, proteinRes, bondPair.first);
                    size_t secondAtom = findResidueAtom(data, proteinRes, bondPair.second);
                    if (firstAtom < atomCount && secondAtom < atomCount)
                    {
                        addBond(data, firstAtom, secondAtom);
                    }
                }
            }

            void bondHydrogenAtomsToClosestNonHydrogen(PdbData& data, size_t residueId)
            {
                std::vector<bool> isResidueAtom = assembly::isResidueAtom(data.indices, residueId);
                std::vector<bool> isHydrogen = util::vectorEquals(data.atoms.elements, Element::H);
                std::vector<size_t> hydrogenAtoms =
                    indicesOfLivingAtoms(data.indices, util::vectorAnd(isResidueAtom, isHydrogen));
                std::vector<size_t> nonHydrogenAtoms =
                    indicesOfLivingAtoms(data.indices, util::vectorAnd(isResidueAtom, util::vectorNot(isHydrogen)));
                if (nonHydrogenAtoms.empty())
                {
                    return;
                }
                double hydrogenCutoff = hydrogenCovalentBondMaxLength();
                for (size_t hydrogen : hydrogenAtoms)
                {
                    Coordinate coord = data.atoms.coordinates[hydrogen];
                    size_t closest = 0;
                    double minDistance = distance(coord, data.atoms.coordinates[nonHydrogenAtoms[0]]);
                    for (size_t n = 1; n < nonHydrogenAtoms.size(); n++)
                    {
                        double dist = distance(coord, data.atoms.coordinates[nonHydrogenAtoms[n]]);
                        if (dist < minDistance)
                        {
                            if (minDistance < hydrogenCutoff)
                            {
                                util::log(
                                    __LINE__,
                                    __FILE__,
                                    util::WAR,
                                    "Hydrogen atom within " + std::to_string(minDistance) +
                                        " Å of multiple non-hydrogen atoms in residue");
                            }
                            minDistance = dist;
                            closest = n;
                        }
                    }
                    if (minDistance < hydrogenCutoff)
                    {
                        addBond(data, nonHydrogenAtoms[closest], hydrogen);
                    }
                    else
                    {
                        util::log(
                            __LINE__,
                            __FILE__,
                            util::WAR,
                            "Hydrogen atom not within " + std::to_string(hydrogenCutoff) +
                                " Å of any non-hydrogen atom in residue");
                    }
                }
            }

            void setBondingForAminoAcid(const AminoAcidTable& table, PdbData& data, size_t residueId)
            {
                size_t aminoAcid = aminoAcidIndex(table, data.residues.names[residueId]);
                addBonds(data, residueId, table.bonds[aminoAcid]);
                addBonds(data, residueId, carboxylBonds());
                bondHydrogenAtomsToClosestNonHydrogen(data, residueId);
            }

            bool autoConnectSuccessiveResidues(PdbData& data, size_t cTermRes, size_t nTermRes)
            {
                size_t atomCount = data.indices.atomCount;
                size_t cAtom = findResidueAtom(data, cTermRes, "C");
                size_t oxtAtom = findResidueAtom(data, cTermRes, "OXT");
                size_t nAtom = findResidueAtom(data, nTermRes, "N");
                if ((cAtom < atomCount) && (nAtom < atomCount) && (oxtAtom == atomCount) &&
                    isWithinBondingDistance(data, cAtom, nAtom))
                {
                    addBond(data, cAtom, nAtom);
                    return true;
                }
                return false;
            }

            void setProteinInterConnectivity(
                const AminoAcidTable& table, PdbData& data, const std::vector<size_t>& residueIds)
            {
                if (residueIds.empty())
                {
                    return;
                }
                setBondingForAminoAcid(table, data, residueIds[0]); // does the first one, and handles when size is 1.
                for (size_t n = 0; n < residueIds.size() - 1; n++)
                {
                    size_t previousRes = residueIds[n];
                    if (!autoConnectSuccessiveResidues(data, previousRes, residueIds[n + 1]))
                    { // Automatically bond the N and C atoms of successive residues
                        util::log(
                            __LINE__,
                            __FILE__,
                            util::WAR,
                            "Gap detected between " + residueStringId(data, previousRes) + " and " +
                                residueStringId(data, residueIds[n + 1]));
                    }
                }
            }
        } // namespace

        void setIntraConnectivity(const AminoAcidTable& table, PdbData& data)
        {
            std::vector<bool> isProtein = util::vectorEquals(data.residues.types, ResidueType::Protein);
            for (size_t residueId : util::boolsToIndices(isProtein))
            {
                setBondingForAminoAcid(table, data, residueId);
            }
            distanceBondIntra(data, util::boolsToIndices(util::vectorNot(isProtein)));
        }

        void setInterConnectivity(const AminoAcidTable& table, PdbData& data, const assembly::Bounds& bounds)
        {
            std::vector<bool> isProtein = util::vectorEquals(data.residues.types, ResidueType::Protein);
            setProteinInterConnectivity(table, data, util::boolsToIndices(isProtein));
            distanceBondInter(data, bounds, util::boolsToIndices(util::vectorNot(isProtein)));
        }
    } // namespace pdb
} // namespace gmml
