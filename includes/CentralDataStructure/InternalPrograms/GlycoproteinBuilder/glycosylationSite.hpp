#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP

#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/CentralDataStructure/Overlaps/overlaps.hpp"

enum MoleculeType
{
    PROTEIN,
    GLYCAN,
    ALL
};

using cds::Atom;
using cds::Residue;
using cds::ResidueLinkage;
using cds::RotatableDihedral;
using cdsCondensedSequence::Carbohydrate;

class GlycosylationSite
{
  public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    GlycosylationSite(Residue* residue, Carbohydrate* carbohydrate, std::vector<Residue*> otherProteinResidues,
                      unsigned int glycanStartResidueNumber);

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline std::string GetResidueId()
    {
        return this->GetResidue()->getStringId();
    }

    inline Carbohydrate* GetGlycan()
    {
        return glycan_;
    }

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    inline void SetOtherGlycosites(std::vector<GlycosylationSite*> glycosites)
    {
        other_glycosites_ = glycosites;
    }

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    cds::Overlap CountOverlaps(MoleculeType moleculeType);
    cds::Overlap CountOverlapsFast();
    //    void StashCoordinates();
    //    void SetStashedCoordinates();
    void Wiggle(bool firstLinkageOnly, int interval);
    void SetRandomDihedralAnglesUsingMetadata();
    void AddOtherGlycositesToLinkageOverlapAtoms();
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
    void PrintOverlaps();
    void Print(std::string type = "All");

  private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSOR                    //
    //////////////////////////////////////////////////////////
    inline std::vector<GlycosylationSite*>& GetOtherGlycosites()
    {
        return other_glycosites_;
    }

    inline std::vector<Residue*>& GetOtherProteinResidues()
    {
        return otherProteinResidues_;
    }

    inline Residue* GetResidue() const
    {
        return residue_;
    }

    inline ResidueLinkage& GetProteinGlycanLinkage()
    {
        return this->GetGlycan()->GetGlycosidicLinkages().front();
    }

    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void AttachGlycan(unsigned int glycanResidueStartNumber);
    void RenumberGlycanToMatch(unsigned int startNumber);
    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue* glycosite_residue);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    Atom* GetConnectingProteinAtom(const std::string residue_name) const;
    void WiggleOneLinkage(ResidueLinkage& linkage, int interval);
    cds::Overlap CountOverlaps(const std::vector<Residue*>& residuesA, const std::vector<Residue*>& residuesB);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    Residue* residue_; /*!< A pointer back to the residue for this glycosite >*/
    Carbohydrate* glycan_;
    std::vector<GlycosylationSite*> other_glycosites_;
    std::vector<Residue*> otherProteinResidues_;
};
#endif
