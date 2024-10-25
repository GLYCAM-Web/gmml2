#ifndef INCLUDES_MOLECULARMETADATA_GLYCAM_DIHEDRALANGLEDATA_HPP
#define INCLUDES_MOLECULARMETADATA_GLYCAM_DIHEDRALANGLEDATA_HPP

#include <string>
#include <functional>
#include <variant>
#include <vector>
#include <stdexcept>

namespace GlycamMetadata
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

        inline std::string print()
        {
            return atom1_ + "_" + atom2_ + "_" + atom3_ + "_" + atom4_ + " : " + rotamer_name_;
        }

        //////////////////////////////////////////////////////////
        //                       OPERATORS                      //
        //////////////////////////////////////////////////////////
        inline bool operator==(const DihedralAngleData& other)
        {
            return (this->index_ == other.index_ &&
                    this->number_of_bonds_from_anomeric_carbon_ == other.number_of_bonds_from_anomeric_carbon_);
        }

        inline bool operator!=(const DihedralAngleData& other)
        {
            return !(operator==(other));
        }
    };

    typedef std::vector<DihedralAngleData> DihedralAngleDataVector;

    // Pass in the two atoms on either side the residue-residue linkage
    std::vector<DihedralAngleDataVector> getDihedralAngleDataEntriesForLinkage(const std::string atom1Name,
                                                                               const std::string residue1Name,
                                                                               const std::string atom2Name,
                                                                               const std::string residue2Name);
    std::vector<double> dihedralAngleDataWeights(const DihedralAngleDataVector& metadataVector);
    DihedralAngleDataVector likelyMetadata(const DihedralAngleDataVector& entries);

    template<typename T>
    T onAngleDeviation(std::function<T(const AngleLimit&)>& onLimit, std::function<T(const AngleStd&)>& onStd,
                       const AngleDeviation& deviation)
    {
        if (std::holds_alternative<AngleLimit>(deviation))
        {
            return onLimit(std::get<AngleLimit>(deviation));
        }
        else if (std::holds_alternative<AngleStd>(deviation))
        {
            return onStd(std::get<AngleStd>(deviation));
        }
        else
        {
            throw std::runtime_error("unhandled angle deviation in GlycamMetadata::onAngleDeviation");
        }
    };
} // namespace GlycamMetadata
#endif
