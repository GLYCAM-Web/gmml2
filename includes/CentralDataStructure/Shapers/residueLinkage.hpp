#ifndef GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_RESIDUE_LINKAGE_HPP
#define GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_RESIDUE_LINKAGE_HPP
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a RotatableDihedral object.
 */
#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CodeUtils/logging.hpp"

#include <string>
#include <array>
#include <vector>
#include <utility>

namespace cds
{
    struct ResidueLink
    {
        std::pair<cds::Residue*, cds::Residue*> residues;
        std::pair<cds::Atom*, cds::Atom*> atoms;
    };

    struct ResidueLinkNames
    {
        std::pair<std::string, std::string> residues;
        std::pair<std::string, std::string> atoms;
    };

    struct DihedralAtoms
    {
        bool isBranching;
        std::array<cds::Atom*, 4> atoms;
    };

    typedef std::vector<gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector> DihedralAngleMetadata;

    ResidueLink findResidueLink(std::pair<cds::Residue*, cds::Residue*> residues);

    class ResidueLinkage
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ResidueLinkage(ResidueLink link);
        //~ResidueLinkage() {std::cout << "Linkage dtor for " << this->GetFromThisResidue1()->getId() << " -Link- "  <<
        // this->GetToThisResidue2()->getId() << "\n";}
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::vector<RotatableDihedral>& GetRotatableDihedralsRef();
        std::vector<RotatableDihedral> GetRotatableDihedrals() const;
        std::vector<RotatableDihedral> GetRotatableDihedralsWithMultipleRotamers() const;
        int GetNumberOfShapes(const bool likelyShapesOnly = false) const;

        inline cds::Residue* GetFromThisResidue1() const
        {
            return link_.residues.first;
        }

        inline cds::Residue* GetToThisResidue2() const
        {
            return link_.residues.second;
        }

        bool CheckIfConformer() const;

        inline unsigned long long GetIndex() const
        {
            return index_;
        }

        std::string GetName() const;

        void AddNonReducingOverlapResidues(std::vector<cds::Residue*> extraResidues);
        std::vector<cds::Residue*>& GetNonReducingOverlapResidues();
        std::vector<cds::Residue*>& GetReducingOverlapResidues();

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void SetRotatableDihedrals(std::vector<RotatableDihedral> rotatableDihedrals)
        {
            rotatableDihedrals_ = rotatableDihedrals;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void SetDefaultShapeUsingMetadata();
        void SetRandomShapeUsingMetadata();
        void SetSpecificShapeUsingMetadata(int shapeNumber);
        void SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
        void SetShapeToPrevious();
        void DetermineAtomsThatMove();
        void SimpleWiggleCurrentRotamers(std::vector<cds::Atom*>& overlapAtomSet1,
                                         std::vector<cds::Atom*>& overlapAtomSet2, const int angleIncrement = 5);
        void SimpleWiggleCurrentRotamers(const std::array<std::vector<cds::Residue*>, 2>& residues,
                                         const int angleIncrement = 5);

        inline void SetIndex(unsigned long long index)
        {
            index_ = index;
        }

        // void AddResiduesForOverlapCheck(std::vector<cds::Residue*> extraResidues, bool reducingSide = true);
        void DetermineResiduesForOverlapCheck();
        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        std::string Print() const;

        //////////////////////////////////////////////////////////
        //                  OPERATOR OVERLOADING                //
        //////////////////////////////////////////////////////////
        bool operator==(const ResidueLinkage& rhs) const
        {
            return (this->GetIndex() == rhs.GetIndex());
        }

        bool operator!=(const ResidueLinkage& rhs) const
        {
            return (this->GetIndex() != rhs.GetIndex());
        }

      private:
        //////////////////////////////////////////////////////////
        //                    PRIVATE FUNCTIONS                 //
        //////////////////////////////////////////////////////////

        void InitializeClass();
        void CreateHydrogenForPsiAngles(std::vector<DihedralAtoms>& dihedralAtoms,
                                        const DihedralAngleMetadata& metadata);
        unsigned long long GenerateIndex();

        inline void SetName(std::string name)
        {
            name_ = name;
        }

        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        ResidueLink link_;
        std::vector<RotatableDihedral> rotatableDihedrals_;
        unsigned long long index_ = 0;
        std::string name_         = ""; // e.g. "DGalpb1-6DGlcpNAc". It being empty works with GetName();
        std::vector<cds::Residue*> nonReducingOverlapResidues_; // overlap speedups
        std::vector<cds::Residue*> reducingOverlapResidues_;    // overlap speedups
    };
} // namespace cds
#endif // GMML_INCLUDES_GEOMETRYTOPOLOGY_RESIDUELINKAGES_RESIDUE_LINKAGE_HPP
