#include "include/readers/parameterManager.hpp"

#include "include/CentralDataStructure/atom.hpp"
#include "include/CentralDataStructure/cdsFunctions.hpp"
#include "include/CentralDataStructure/residue.hpp"
#include "include/util/containers.hpp"
#include "include/util/filesystem.hpp"
#include "include/util/logging.hpp"

// How exactly this happens can be improved, but the information should only ever be loaded into gmml in one place.
namespace gmml
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

    ParameterManager loadParameters(const std::string& baseDir)
    { // Library files of 3D structures with parameters for simulations
        ParameterManager result;
        for (auto& libFilePath : libFilesToLoad)
        {
            lib::LibraryData libData = lib::loadLibraryData(baseDir + "/" + libFilePath);
            result.lib.residueNames = util::vectorAppend(result.lib.residueNames, libData.residueNames);
            result.lib.residues = util::vectorAppend(result.lib.residues, libData.residues);
        }
        util::log(__LINE__, __FILE__, util::INF, "Finished construction of ParameterManager.");
        return result;
    }

    void setAtomChargesForResidues(const ParameterManager& parameters, std::vector<Residue*> queryResidues)
    {
        for (auto& residue : queryResidues)
        {
            setAtomChargesForResidue(parameters, residue);
        }
    }

    void createAtomsForResidue(
        const ParameterManager& parameters, Residue* queryResidue, const std::string glycamNameForResidue)
    {
        size_t index = util::indexOf(parameters.lib.residueNames, glycamNameForResidue);
        if (index == parameters.lib.residueNames.size())
        {
            std::string message = "Did not find a parameter residue for " + queryResidue->getName() +
                                  " with this glycam residue code: " + glycamNameForResidue;
            throw std::runtime_error(message);
        }
        const std::string& glycamCode = parameters.lib.residueNames[index];
        queryResidue->setName(glycamCode);
        queryResidue->determineType(glycamCode);
        const lib::ResidueData& residue = parameters.lib.residues[index];
        const lib::AtomData& atoms = residue.atoms;
        std::vector<Atom*> atomVec;
        atomVec.reserve(atoms.names.size());
        for (size_t n = 0; n < atoms.names.size(); n++)
        {
            Atom* atom = queryResidue->addAtom(std::make_unique<Atom>());
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
            addBond(atomVec[bond[0]], atomVec[bond[1]]);
        }
    }

    bool setAtomChargesForResidue(const ParameterManager& parameters, Residue* queryResidue)
    {
        bool allAtomsPresent = true;
        size_t index = util::indexOf(parameters.lib.residueNames, queryResidue->GetParmName());
        if (index == parameters.lib.residueNames.size())
        {
            util::log(
                __LINE__,
                __FILE__,
                util::WAR,
                "Did not find parameters and so cannot set charges for residue named: " + queryResidue->GetParmName());
            return false;
        }

        const lib::AtomData& atoms = parameters.lib.residues[index].atoms;
        for (auto& queryAtom : queryResidue->getAtoms())
        {
            const std::string& atomName = queryAtom->getName();
            size_t atomIndex = util::indexOf(atoms.names, atomName);
            if (atomIndex < atoms.names.size())
            {
                queryAtom->setCharge(atoms.charges[atomIndex]);
                queryAtom->setType(atoms.types[atomIndex]);
            }
            else
            {
                util::log(__LINE__, __FILE__, util::WAR, "No charges found for atom named " + atomName);
                allAtomsPresent = false;
            }
        }
        return allAtomsPresent;
    }
} // namespace gmml
