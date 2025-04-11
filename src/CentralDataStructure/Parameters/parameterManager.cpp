#include "includes/CentralDataStructure/Parameters/parameterManager.hpp"
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

namespace
{
    void initializeResidueMap(ParameterManager& parameters, std::vector<cds::Residue*> incomingResidues)
    {
        for (auto& residue : incomingResidues)
        {
            parameters.residueNames.push_back(residue->getName());
            parameters.residues.push_back(residue);
            gmml::log(__LINE__, __FILE__, gmml::INF, "Added this to map: " + residue->getName());
        }
    }

    bool setChargeForAtom(cds::Atom* queryAtom, std::vector<cds::Atom*> referenceAtoms)
    {
        for (auto& refAtom : referenceAtoms)
        {
            if (queryAtom->getName() == refAtom->getName())
            {
                queryAtom->setCharge(refAtom->getCharge());
                queryAtom->setType(refAtom->getType());
                return true;
            }
        }
        gmml::log(__LINE__, __FILE__, gmml::WAR, "No charges found for atom named " + queryAtom->getName());
        return false;
    }

    cds::Residue copyParameterResidue(const ParameterManager& parameters, const std::string name)
    {
        cds::Residue* reference = findParameterResidue(parameters, name);
        if (reference == nullptr)
        {
            std::string message = "Did not find and therefore cannot copy a parameter residue with this name: " + name;
            gmml::log(__LINE__, __FILE__, gmml::ERR, message);
            throw std::runtime_error(message);
        }
        return *reference; // copy
    }
} // namespace

ParameterManager cdsParameters::loadParameters(const std::string& baseDir)
{ // Library files of 3D structures with parameters for simulations.
    gmml::log(__LINE__, __FILE__, gmml::INF, "gmmlhome is: " + baseDir);
    ParameterManager result;
    for (auto& prepFilePath : cdsParameters::prepFilesToLoad)
    {
        auto& file = result.prepFiles.emplace_back(baseDir + "/" + prepFilePath);
        initializeResidueMap(result, file.getResidues());
    }
    for (auto& libFilePath : cdsParameters::libFilesToLoad)
    {
        cds::Molecule& molecule = result.libFiles.emplace_back(cds::Molecule());
        lib::parseMolecule(&molecule, baseDir + "/" + libFilePath);
        initializeResidueMap(result, molecule.getResidues());
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished construction of ParameterManager.");
    return result;
}

cds::Residue* cdsParameters::findParameterResidue(const ParameterManager& parameters, const std::string name)
{
    size_t index = codeUtils::indexOf(parameters.residueNames, name);
    if (index < parameters.residueNames.size())
    {
        return parameters.residues[index];
    }
    gmml::log(__LINE__, __FILE__, gmml::WAR, "Did not find parameters for residue named: " + name);
    return nullptr;
}

void cdsParameters::setAtomChargesForResidues(const ParameterManager& parameters,
                                              std::vector<cds::Residue*> queryResidues)
{
    for (auto& residue : queryResidues)
    {
        setAtomChargesForResidue(parameters, residue);
    }
}

void cdsParameters::createAtomsForResidue(const ParameterManager& parameters,
                                          cdsCondensedSequence::ParsedResidue* queryResidue,
                                          const std::string glycamNameForResidue)
{
    try
    {
        cds::Residue parameterResidue = copyParameterResidue(parameters, glycamNameForResidue);
        queryResidue->setName(parameterResidue.getName());
        queryResidue->setAtoms(parameterResidue.extractAtoms());
        return;
    }
    catch (const std::runtime_error& error)
    { // I just want to throw a nicer error as this happens a lot:
        std::string message = "Did not find a parameter residue for " + queryResidue->getName() +
                              " with this glycam residue code: " + glycamNameForResidue;
        throw std::runtime_error(message);
    }
}

bool cdsParameters::setAtomChargesForResidue(const ParameterManager& parameters, cds::Residue* queryResidue)
{
    bool allAtomsPresent           = true;
    cds::Residue* parameterResidue = findParameterResidue(parameters, queryResidue->GetParmName());
    if (parameterResidue == nullptr)
    {
        gmml::log(__LINE__, __FILE__, gmml::WAR,
                  "Did not find parameters and so cannot set charges for residue named: " +
                      queryResidue->GetParmName());
        return false;
    }
    std::vector<cds::Atom*> parameterAtoms = parameterResidue->getAtoms();
    for (auto& queryAtom : queryResidue->getAtoms())
    {
        if (!setChargeForAtom(queryAtom, parameterAtoms))
        {
            allAtomsPresent = false;
        }
    }
    return allAtomsPresent;
}
