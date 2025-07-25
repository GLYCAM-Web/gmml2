#ifndef INCLUDE_METADATA_DIHEDRALANGLEDATA_HPP
#define INCLUDE_METADATA_DIHEDRALANGLEDATA_HPP

#include "include/CentralDataStructure/residueTypes.hpp"

#include <string>
#include <variant>
#include <vector>

namespace gmml
{
    enum class RotamerType
    {
        conformer,
        permutation
    };

    struct AngleLimit
    {
        double lowerDeviationLimit;
        double upperDeviationLimit;
    };

    struct AngleStd
    {
        double lowerDeviationStd;
        double upperDeviationStd;
    };

    typedef std::variant<AngleLimit, AngleStd> AngleDeviation;

    struct DihedralAngleData
    {
        std::string linking_atom1_;
        std::string linking_atom2_;
        std::string dihedral_angle_name_;
        double default_angle;
        AngleDeviation angle_deviation;
        double weight_;
        RotamerType rotamer_type_; // permutation or conformer
        std::string rotamer_name_;
        unsigned int number_of_bonds_from_anomeric_carbon_;
        int index_; // Used to indicate whether multiple entries are meant to overwrite each other or generate
                    // an additional angle
        std::vector<std::string> residue1_conditions_;
        std::vector<std::string> residue2_conditions_;
        std::string atom1_;
        std::string atom2_;
        std::string atom3_;
        std::string atom4_;
    };

    struct DihedralAngleDataTable
    {
        std::vector<DihedralAngleData> entries;
        std::vector<double> weights;
    };

    const DihedralAngleDataTable& dihedralAngleDataTable();

    std::vector<std::vector<size_t>> getDihedralAngleDataEntriesForLinkage(
        const std::string& atom1Name,
        const ResidueAttributes& residue1Attributes,
        const std::string& atom2Name,
        const ResidueAttributes& residue2Attributes);

    std::vector<size_t> likelyMetadata(const DihedralAngleDataTable& table, const std::vector<size_t>& entries);
} // namespace gmml

#endif
