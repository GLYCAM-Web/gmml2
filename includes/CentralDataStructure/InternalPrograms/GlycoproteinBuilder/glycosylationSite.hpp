#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOSYLATIONSITE_HPP

#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/assembly.hpp"

using cds::Atom;
using cds::Residue;
using cdsCondensedSequence::Carbohydrate;

class GlycosylationSite
{
  public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    GlycosylationSite(Residue* residue, Carbohydrate* carbohydrate, unsigned int glycanStartResidueNumber);

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

    inline Residue* GetResidue() const
    {
        return residue_;
    }

  private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    void AttachGlycan(unsigned int glycanResidueStartNumber);
    void RenumberGlycanToMatch(unsigned int startNumber);
    void Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name);
    void Superimpose_Glycan_To_Glycosite(Residue* glycosite_residue);
    void Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    Atom* GetConnectingProteinAtom(const std::string residue_name) const;
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    Residue* residue_; /*!< A pointer back to the residue for this glycosite >*/
    Carbohydrate* glycan_;
};
#endif
