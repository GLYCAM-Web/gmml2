#include "include/fileType/prep/prepAtom.hpp"

#include "include/fileType/prep/prepFunctions.hpp"
#include "include/geometry/measurements.hpp" //get_cartesian_point_from_internal_coords()
#include "include/graph/graphManipulation.hpp"
#include "include/util/casting.hpp"
#include "include/util/containers.hpp"
#include "include/util/logging.hpp"
#include "include/util/strings.hpp"

#include <iomanip>
#include <ios>
#include <sstream>
#include <stdexcept>

namespace gmml
{
    namespace prep
    {
        namespace
        {
            TopologicalType extractAtomTopologicalType(std::istream& ss)
            {
                std::string s;
                ss >> s;
                size_t index = util::indexOf(topologicalTypeNames, s);
                if (index == topologicalTypeNames.size())
                {
                    throw std::runtime_error("Unknown topological type");
                }
                return TopologicalType(index);
            }

            std::vector<size_t> findDihedralAtoms(const PrepData& data, size_t index)
            {
                const std::vector<bool>& nodeAlive = data.atomGraph.nodeAlive;
                std::vector<size_t> result = {index};
                size_t lastId = index;
                auto isParent = [&](const std::array<size_t, 2>& nodes)
                { return (nodes[1] == lastId) && nodeAlive[nodes[0]] && nodeAlive[nodes[1]]; };
                for (int currentDepth = 0; currentDepth < 3; currentDepth++)
                {
                    auto it = std::find_if(data.atomGraph.edgeNodes.begin(), data.atomGraph.edgeNodes.end(), isParent);
                    size_t edgeId =
                        it - data.atomGraph.edgeNodes.begin(); // Go up the first parent only. Loops may create
                                                               // another parent, but they should be ignored.
                    size_t parent = data.atomGraph.edgeNodes[edgeId][0];
                    result.push_back(parent);
                    lastId = parent;
                }
                return result;
            }
        } // namespace

        void initializePrepAtom(PrepData& data, size_t residueId, const std::string& line)
        {
            AtomData& atoms = data.atoms;
            size_t atomId = addAtom(data, residueId);
            std::stringstream ss(line);

            atoms.number[atomId] = util::extractFromStream(ss, uint());
            atoms.name[atomId] = util::extractFromStream(ss, std::string());
            atoms.type[atomId] = util::extractFromStream(ss, std::string());
            atoms.topologicalType[atomId] = extractAtomTopologicalType(ss);
            atoms.bondIndex[atomId] = util::extractFromStream(ss, uint());
            atoms.angleIndex[atomId] = util::extractFromStream(ss, uint());
            atoms.dihedralIndex[atomId] = util::extractFromStream(ss, uint());
            atoms.bondLength[atomId] = util::extractFromStream(ss, double());
            atoms.angle[atomId] = util::extractFromStream(ss, double());
            atoms.dihedral[atomId] = util::extractFromStream(ss, double());
            atoms.charge[atomId] = util::extractFromStream(ss, double());
            atoms.coordinate[atomId] = {0, 0, 0};
        }

        Coordinate determine3dCoordinate(PrepData& data, size_t index)
        {
            std::vector<size_t> foundAtoms = findDihedralAtoms(data, index);
            return calculateCoordinateFromInternalCoords(
                data.atoms.coordinate[foundAtoms[3]],
                data.atoms.coordinate[foundAtoms[2]],
                data.atoms.coordinate[foundAtoms[1]],
                data.atoms.angle[index],
                data.atoms.dihedral[index],
                data.atoms.bondLength[index]);
        }

        void printAtom(const PrepData& data, size_t index, std::ostream& out)
        {
            const AtomData& atoms = data.atoms;
            out << std::setw(3) << atoms.number[index] << std::setw(6) << atoms.name[index] << std::setw(6)
                << atoms.type[index] << std::setw(3) << topologicalTypeNames[atoms.topologicalType[index]]
                << std::setw(4) << atoms.bondIndex[index] << std::setw(4) << atoms.angleIndex[index] << std::setw(4)
                << atoms.dihedralIndex[index] << std::setw(10) << atoms.bondLength[index] << std::setw(10)
                << atoms.angle[index] << std::setw(10) << atoms.dihedral[index] << std::setw(10) << atoms.charge[index];
        }

        void writeAtom(const PrepData& data, size_t index, std::ostream& stream)
        {
            const AtomData& atoms = data.atoms;
            stream << std::right << std::setw(2) << atoms.number[index] << " " << std::left << std::setw(4)
                   << atoms.name[index] << " " << std::left << std::setw(3) << atoms.type[index] << " " << std::setw(1)
                   << topologicalTypeNames[atoms.topologicalType[index]] << " " << std::right << std::setw(2)
                   << atoms.bondIndex[index] << " " << std::right << std::setw(2) << atoms.angleIndex[index] << " "
                   << std::right << std::setw(2) << atoms.dihedralIndex[index] << " " << std::right << std::setw(8)
                   << std::fixed << std::setprecision(3) << atoms.bondLength[index] << " " << std::right << std::setw(8)
                   << std::fixed << std::setprecision(3) << atoms.angle[index] << " " << std::right << std::setw(8)
                   << std::fixed << std::setprecision(3) << atoms.dihedral[index] << "    " << std::right
                   << std::setw(8) << std::fixed << std::setprecision(4) << atoms.charge[index] << std::endl;
        }
    } // namespace prep
} // namespace gmml
