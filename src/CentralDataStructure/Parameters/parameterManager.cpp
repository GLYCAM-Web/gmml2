#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/cdsFunctions/atomicBonding.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/filesystem.hpp"

// How exactly this happens can be improved, but the information should only ever be loaded into gmml in one place.
namespace cdsParameters
{
    //    static const std::vector<std::string> prepFilesToLoad = {"/dat/prep/GLYCAM_06j-1_GAGS_KDN.prep"};
    static const std::vector<std::string> prepFilesToLoad;
    static const std::vector<std::string> libFilesToLoad = {
        "/dat/CurrentParams/GLYCAM_06k.lib",
        "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib",
        "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib",
        "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib",
        "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib",
        "/dat/CurrentParams/other/solvents.lib",
        "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib",
        "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib",
        "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib",
    };
} // namespace cdsParameters

using cdsParameters::ParameterManager;

ParameterManager cdsParameters::loadParameters(const std::string& baseDir)
{ // Library files of 3D structures with parameters for simulations
    ParameterManager result;
    for (auto& libFilePath : cdsParameters::libFilesToLoad)
    {
        lib::LibraryData libData = lib::loadLibraryData(baseDir + "/" + libFilePath);
        result.lib.residueNames  = codeUtils::vectorAppend(result.lib.residueNames, libData.residueNames);
        result.lib.residues      = codeUtils::vectorAppend(result.lib.residues, libData.residues);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished construction of ParameterManager.");
    return result;
}

void cdsParameters::setAtomChargesForResidues(const ParameterManager& parameters,
                                              std::vector<cds::Residue*> queryResidues)
{
    for (auto& residue : queryResidues)
    {
        setAtomChargesForResidue(parameters, residue);
    }
}

void cdsParameters::createAtomsForResidue(const ParameterManager& parameters, cds::Residue* queryResidue,
                                          const std::string glycamNameForResidue)
{
    size_t index = codeUtils::indexOf(parameters.lib.residueNames, glycamNameForResidue);
    if (index == parameters.lib.residueNames.size())
    {
        std::string message = "Did not find a parameter residue for " + queryResidue->getName() +
                              " with this glycam residue code: " + glycamNameForResidue;
        throw std::runtime_error(message);
    }
    const std::string& name = parameters.lib.residueNames[index];
    queryResidue->setName(name);
    queryResidue->determineType(name);
    const lib::ResidueData& residue = parameters.lib.residues[index];
    const lib::AtomData& atoms      = residue.atoms;
    std::vector<cds::Atom*> atomVec;
    atomVec.reserve(atoms.names.size());
    for (size_t n = 0; n < atoms.names.size(); n++)
    {
        cds::Atom* atom = queryResidue->addAtom(std::make_unique<cds::Atom>());
        atomVec.push_back(atom);
        atom->setName(atoms.names[n]);
        atom->setType(atoms.types[n]);
        atom->setCharge(atoms.charges[n]);
        atom->setNumber(atoms.numbers[n]);
        if (residue.hasCoordinates)
        {
            atom->setCoordinate(atoms.coordinates[n]);
        }
    }
    for (auto& bond : residue.bonds)
    {
        cds::addBond(atomVec[bond[0]], atomVec[bond[1]]);
    }
}

bool cdsParameters::setAtomChargesForResidue(const ParameterManager& parameters, cds::Residue* queryResidue)
{
    bool allAtomsPresent = true;
    size_t index         = codeUtils::indexOf(parameters.lib.residueNames, queryResidue->GetParmName());
    if (index == parameters.lib.residueNames.size())
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Did not find parameters and so cannot set charges for residue named: " +
                      queryResidue->GetParmName());
        return false;
    }

    const lib::AtomData& atoms = parameters.lib.residues[index].atoms;
    for (auto& queryAtom : queryResidue->getAtoms())
    {
        const std::string& atomName = queryAtom->getName();
        size_t atomIndex            = codeUtils::indexOf(atoms.names, atomName);
        if (atomIndex < atoms.names.size())
        {
            queryAtom->setCharge(atoms.charges[atomIndex]);
            queryAtom->setType(atoms.types[atomIndex]);
        }
        else
        {
            gmml::log(__LINE__, __FILE__, gmml::WAR, "No charges found for atom named " + atomName);
            allAtomsPresent = false;
        }
    }
    return allAtomsPresent;
}
