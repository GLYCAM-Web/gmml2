#ifndef ONTOLOGYVOCABULARY_HPP
#define ONTOLOGYVOCABULARY_HPP

#include <string>

namespace Ontology
{
    const std::string ONT_DOMAIN = "<http://gmmo.uga.edu/#";
    const std::string ONT_PREFIX = "gmmo:";
    const std::string ENTITY_COMMENT = "\n ###  http://gmmo.uga.edu#";

    const std::string TYPE = "rdf:type";
    const std::string LABEL = "rdfs:label";

    /* Classes */
    const std::string Atom = "gmmo:Atom";
    const std::string RingAtom = "gmmo:RingAtom";
    const std::string SideAtom = "gmmo:SideAtom";
    const std::string Linkage = "gmmo:Linkage";
    const std::string Monosaccharide = "gmmo:Monosaccharide";
    const std::string Oligosaccharide = "gmmo:Oligosaccharide";
    const std::string PDB = "gmmo:PDB";
    const std::string Residue = "gmmo:Residue";
    const std::string SugarName = "gmmo:SugarName";

    /* Object Properties */
    const std::string hasAtom = "gmmo:hasAtom";
    const std::string hasChild = "gmmo:hasChild";
    const std::string hasChildAtomLinkage = "gmmo:hasChildAtomIndex";
    const std::string hasGlycosidicLinkage = "gmmo:hasGlycosidicLinkage";
    const std::string hasNeighbor = "gmmo:hasNeighbor";
    const std::string hasOligo = "gmmo:hasOligo";
    const std::string hasParent = "gmmo:hasParent";
    const std::string hasParentAtomLinkage = "gmmo:hasParentAtomLinkage";
    const std::string hasResidue = "gmmo:hasResidue";
    const std::string hasRingAtom = "gmmo:hasRingAtom";
    const std::string hasRoot = "gmmo:hasRoot";
    const std::string hasSideAtom = "gmmo:hasSideAtom";
    const std::string hasSugarName = "gmmo:hasSugarName";
    const std::string hasmemberOfResidue = "gmmo:memberOfResidue"; ///remove maybe?

    /* Datatype Properties */
    const std::string anomeric_status = "gmmo:anomericStatus";
    const std::string chemical_code_str = "gmmo:stringChemicalCode";
    const std::string configuration = "gmmo:configuration";
    const std::string cycle_atom_str = "gmmo:cycleAtomsString";
    const std::string derivative = "gmmo:derivative";
    const std::string glycosidic_linkage_str = "gmmo:glycosidicLinkageString";
    const std::string group_name = "gmmo:group_name"; ///remove?
    const std::string id = "gmmo:identifier";
    const std::string isomer = "gmmo:isomer";
    const std::string linkage_str = "gmmo:linkageString";
    const std::string mono_name = "gmmo:monosaccharideName";
    const std::string mono_short_name = "gmmo:monosaccharideShortName";
    const std::string mono_stereo_name = "gmmo:monosaccharideStereochemName";
    const std::string mono_stereo_short_name = "gmmo:monosaccharideStereochemShortName";
    const std::string oligo_name = "gmmo:oligoName";
    const std::string orientation = "gmmo:orientation";
    const std::string path = "gmmo:path";
    const std::string ring_index = "gmmo:ringIndex";
    const std::string ring_type = "gmmo:ringType";
    const std::string side_index = "gmmo:sideIndex";
    const std::string x = "gmmo:xCoordinate";
    const std::string y = "gmmo:yCoordinate";
    const std::string z = "gmmo:zCoordinate";
}

#endif // ONTOLOGYVOCABULARY_HPP