#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_LINKAGEBRANCHES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_LINKAGEBRANCHES_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"

#include <vector>

namespace cdsSelections
{
    // Helper class
    class Branch // This should all be replaced by better code in the Template Graph Library. Looking at you Preston.
    {
      public:
        Branch(cds::Atom* rootAtom, int depth = 0) : depth_(depth), maxDepth_(depth), rootAtom_(rootAtom)
        {
            branchFound = false;
        };

        inline cds::Atom* GetEnd()
        {
            return endAtom_;
        }

        inline cds::Atom* GetRoot()
        {
            return rootAtom_;
        }

        inline int GetDepth()
        {
            return depth_;
        }

        inline bool IsBranchFound()
        {
            return branchFound;
        }

        // inline void AddToPath(Atom* atom) {path_.push_back(atom);}
        inline void SetRoot(cds::Atom* atom)
        {
            rootAtom_ = atom;
        }

        inline void SetEnd(cds::Atom* atom)
        {
            endAtom_    = atom;
            branchFound = true;
        }

        inline void SetDepth(int depth)
        {
            depth_ = depth;
        }

        inline void ChangeDepth(int delta)
        {
            depth_ += delta;
            if (depth_ > maxDepth_)
            {
                maxDepth_ = depth_;
            }
        }

        inline bool AtMaxDepth()
        {
            return maxDepth_ == depth_;
        }

      private:
        int depth_;
        int maxDepth_;
        cds::Atom* rootAtom_ = nullptr;
        cds::Atom* endAtom_  = nullptr;
        bool branchFound;
        // AtomVector path_;
    };

    void FindEndsOfBranchesFromLinkageAtom(cds::Atom* currentAtom, cds::Atom* previousAtom, cds::Residue* residue,
                                           Branch* branch);
    std::tuple<std::vector<cds::DihedralAtoms>, std::vector<cds::Atom*>>
    findRotatableDihedralsinBranchesConnectingResidues(const cds::ResidueLink& link,
                                                       const std::vector<cds::Atom*>& residueCyclePoints);
} // namespace cdsSelections
#endif
