#include <bits/std_abs.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip> // For setting precision and formating in output
#include <stdexcept>

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CentralDataStructure/cdsFunctions/cdsFunctions.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // calculateCoordinateFromInternalCoords
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp" // For the ClearAtomLabels sillyness.
#include "includes/CentralDataStructure/Selections/atomSelections.hpp"
#include "includes/CentralDataStructure/Geometry/geometryTypes.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/references.hpp"
#include "includes/MolecularMetadata/glycoprotein.hpp"

namespace glycoproteinBuilder
{
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    GlycosylationSite::GlycosylationSite(Residue* residue, Carbohydrate* carbohydrate, GlycositeInput input,
                                         unsigned int glycanStartResidueNumber)
        : residue_(residue), glycan_(carbohydrate), input_(input)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching glycan!");
        this->AttachGlycan(glycanStartResidueNumber);
        cdsSelections::ClearAtomLabels(carbohydrate->GetReducingResidue()); // jfc
        cdsSelections::ClearAtomLabels(this->GetResidue());
    }

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void GlycosylationSite::AttachGlycan(unsigned int glycanResidueStartNumber)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Start of AttachGlycan. Residue ID is: " + this->GetResidue()->getStringId());
        this->Prepare_Glycans_For_Superimposition_To_Particular_Residue(this->GetResidue()->getName());
        gmml::log(__LINE__, __FILE__, gmml::INF, "Superimpose prep done");
        this->Superimpose_Glycan_To_Glycosite(this->GetResidue());
        gmml::log(__LINE__, __FILE__, gmml::INF, "SuperimposedGlycanToGlycosite");
        gmml::log(__LINE__, __FILE__, gmml::INF, "Setting internal bond count to check if more form later");
        this->RenumberGlycanToMatch(glycanResidueStartNumber);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attach glycan done");
    }

    void GlycosylationSite::RenumberGlycanToMatch(unsigned int startNumber)
    {
        for (auto& glycanResidue : this->GetGlycan()->getResidues())
        {
            glycanResidue->setNumber(++startNumber);
        }
        return;
    }

    // This function prepares the glycan molecule in the glycan_ assembly for superimpostion onto an amino acid in the
    // protein It does this by "growing" the atoms of the amino acid side chain (e.g. Asn, Thr or Ser) out from the
    // glycan reducing terminal Another function will use these additional atoms to superimpose the glycan onto residue
    void GlycosylationSite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name)
    {
        // Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB ,
        // not CB, CA, N.
        Residue* reducing_Residue = this->GetGlycan()->GetReducingResidue();
        // Want: residue.FindAtomByTag("anomeric-carbon"); The below is risky as it uses atoms names, i.e. would break
        // for Sialic acid. ToDo Ok so I reckon the below is just assuming alpha or beta depending on the concext. Need
        // to fix a lot, but need to reproduce functionality after refactor first.
        //    Atom* anomericAtom = cdsSelections::guessAnomericAtom(reducing_Residue);
        //   This won't work as sometimes want alpha, sometimes beta. i.e. a coordinateOppositeToNeighborAverage
        //   function
        // This needs to be abstracted so it works for C2 reducing residues:
        Coordinate coordC5        = reducing_Residue->FindAtom("C5")->coordinate();
        Coordinate coordO5        = reducing_Residue->FindAtom("O5")->coordinate();
        Coordinate coordC1        = reducing_Residue->FindAtom("C1")->coordinate();
        Atom* anomericAtom        = reducing_Residue->FindAtom("C1"); // For adding bond.
        // Delete aglycon atoms from glycan.
        Residue* aglycon          = this->GetGlycan()->GetAglycone();
        for (auto& atom : aglycon->getAtoms())
        {
            aglycon->deleteAtom(atom);
        }
        // Ok so going to set it so that the new "superimposition residue" is the old aglycon residue
        // This avoids having to delete the algycon residue object from assembly and adding the super residue to
        // assembly.
        Residue* superimposition_residue = aglycon; // "renaming" so the below reads better.
        superimposition_residue->setName("SUP");
        // I put both the regular name and the O/N-linked glycam name here, as I'm not sure when it will be renamed.
        if ((amino_acid_name == "ASN") || (amino_acid_name == "NLN"))
        {
            Atom* atomND2 = superimposition_residue->addAtom(std::make_unique<Atom>(
                "ND2", (cds::calculateCoordinateFromInternalCoords(coordC5, coordO5, coordC1, 109.3, 180, 1.53))));
            Atom* atomCG  = superimposition_residue->addAtom(
                std::make_unique<Atom>("CG", (cds::calculateCoordinateFromInternalCoords(
                                                 coordO5, coordC1, atomND2->coordinate(), 109.3, 261, 1.325))));
            superimposition_residue->addAtom(std::make_unique<Atom>(
                "OD1", (cds::calculateCoordinateFromInternalCoords(coordC1, atomND2->coordinate(), atomCG->coordinate(),
                                                                   126, 0, 1.22))));
            cds::addBond(anomericAtom,
                         atomND2); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
        }
        else if ((amino_acid_name == "THR") || (amino_acid_name == "SER") || (amino_acid_name == "OLT") ||
                 (amino_acid_name == "OLS"))
        {
            Atom* atomOG1 = superimposition_residue->addAtom(std::make_unique<Atom>(
                "OG", (cds::calculateCoordinateFromInternalCoords(coordC5, coordO5, coordC1, 112, 68, 1.46))));
            Atom* atomCB  = superimposition_residue->addAtom(
                std::make_unique<Atom>("CB", (cds::calculateCoordinateFromInternalCoords(
                                                 coordO5, coordC1, atomOG1->coordinate(), 109.3, 75, 1.53))));
            superimposition_residue->addAtom(std::make_unique<Atom>(
                "CA", (cds::calculateCoordinateFromInternalCoords(coordC1, atomOG1->coordinate(), atomCB->coordinate(),
                                                                  109.3, 125, 1.53))));
            if ((amino_acid_name == "THR") || (amino_acid_name == "OLT"))
            {
                atomOG1->setName("OG1"); // It's OG in Ser.
            }
            cds::addBond(anomericAtom,
                         atomOG1); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
        }
        else if ((amino_acid_name == "TYR") || (amino_acid_name == "OLY"))
        {
            Atom* atomOH = superimposition_residue->addAtom(std::make_unique<Atom>(
                "OH", (cds::calculateCoordinateFromInternalCoords(coordC5, coordO5, coordC1, 112, 68, 1.46))));
            Atom* atomCZ = superimposition_residue->addAtom(std::make_unique<Atom>(
                "CZ",
                (cds::calculateCoordinateFromInternalCoords(coordO5, coordC1, atomOH->coordinate(), 117, 60, 1.35))));
            superimposition_residue->addAtom(std::make_unique<Atom>(
                "CE1", (cds::calculateCoordinateFromInternalCoords(coordC1, atomOH->coordinate(), atomCZ->coordinate(),
                                                                   120, 180, 1.37))));
            cds::addBond(anomericAtom,
                         atomOH); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
        }
        else
        {
            std::string message =
                "Problem creating glycosylation site. The amino acid requested: " + GetInput().proteinResidueId +
                " has name (" + amino_acid_name +
                ") that isn't supported. Currently you can glycosylate ASN, THR, SER or TYR. "
                "Email us to request " +
                "others. Ideally include examples of 3D structures we can use as a template.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        return;
    }

    void GlycosylationSite::Superimpose_Glycan_To_Glycosite(Residue* glycosite_residue)
    {
        // Get the 3 target atoms from protein residue.
        std::vector<cds::CoordinateReference> targetCoords;
        // superimposition_atoms_ points to three atoms that were added to the glycan. Based on their names e.g. CG,
        // ND2, we will superimpose them onto the correspoinding "target" atoms in the protein residue
        // (glycosite_residue).
        for (auto& superimposition_atom : this->GetGlycan()->GetAglycone()->getAtoms())
        {
            for (auto& protein_atom : glycosite_residue->getAtoms())
            {
                if (protein_atom->getName() == superimposition_atom->getName())
                {
                    targetCoords.push_back(protein_atom->coordinateReference());
                }
            }
        }
        auto aglyconeAtoms                                   = this->GetGlycan()->GetAglycone()->mutableAtoms();
        auto glycanAtoms                                     = this->GetGlycan()->mutableAtoms();
        std::vector<cds::CoordinateReference> aglyconeCoords = cds::atomCoordinateReferences(aglyconeAtoms);
        std::vector<cds::CoordinateReference> glycanCoords   = cds::atomCoordinateReferences(glycanAtoms);
        // Sanity checks:
        if (aglyconeCoords.size() < 3)
        {
            throw std::runtime_error("The aglycone does not contain enough atoms to perform the requested "
                                     "superimposition to " +
                                     this->GetResidueId() + ".\nCheck your input structure!\n");
        }
        if (targetCoords.size() < 3)
        {
            throw std::runtime_error(
                "Did not find the correctly named atoms in target residue to perform the requested "
                "superimposition to " +
                this->GetResidueId() + ".\nCheck your input structure!\n");
        }
        cds::Superimpose(aglyconeCoords, targetCoords, glycanCoords);
        // Connect the glycan and protein atoms to each other.
        Atom* protein_connection_atom = this->GetConnectingProteinAtom(glycosite_residue->getName());
        cds::addBond(protein_connection_atom, this->GetGlycan()->GetAnomericAtom()); // Atom connectivity
        this->Rename_Protein_Residue_To_GLYCAM_Nomenclature();                       // e.g. ASN to NLN
        this->GetGlycan()->replaceAglycone(this->GetResidue());
        return;
    }

    void GlycosylationSite::Rename_Protein_Residue_To_GLYCAM_Nomenclature()
    {
        this->GetResidue()->setName(
            glycoproteinMetadata::ConvertGlycosylatedResidueName(this->GetResidue()->getName()));
    }

    //////////////////////////////////////////////////////////
    //                   PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    Atom* GlycosylationSite::GetConnectingProteinAtom(const std::string residue_name) const
    {
        std::string connectionAtomName = glycoproteinMetadata::GetGlycositeConnectionAtomName(residue_name);
        if (connectionAtomName == "")
        {
            std::string message = "Problem in GetConnectingProteinAtom. The amino acid requested: " + residue_name +
                                  " has name that isn't supported. Currently you can glycosylate ASN, THR, SER or TYR. "
                                  "Email us to request "
                                  "others. Ideally include examples of 3D structures we can use as a template.";
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        return this->GetResidue()->FindAtom(connectionAtomName);
    }

} // namespace glycoproteinBuilder
