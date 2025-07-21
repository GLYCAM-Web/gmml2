#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_LINKAGEBRANCHES_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_LINKAGEBRANCHES_HPP

#include "include/CentralDataStructure/Shapers/residueLinkageTypes.hpp"
#include "include/CentralDataStructure/cdsTypes.hpp"

#include <vector>

namespace gmml
{
    // Helper class
    class Branch // This should all be replaced by better code in the Template Graph Library. Looking at you Preston.
    {
      public:
        Branch(Atom* rootAtom, int depth = 0) : depth_(depth), maxDepth_(depth), rootAtom_(rootAtom)
        {
            branchFound = false;
        };

        inline Atom* GetEnd() { return endAtom_; }

        inline Atom* GetRoot() { return rootAtom_; }

        inline int GetDepth() { return depth_; }

        inline bool IsBranchFound() { return branchFound; }

        // inline void AddToPath(Atom* atom) {path_.push_back(atom);}
        inline void SetRoot(Atom* atom) { rootAtom_ = atom; }

        inline void SetEnd(Atom* atom)
        {
            endAtom_ = atom;
            branchFound = true;
        }

        inline void SetDepth(int depth) { depth_ = depth; }

        inline void ChangeDepth(int delta)
        {
            depth_ += delta;
            if (depth_ > maxDepth_)
            {
                maxDepth_ = depth_;
            }
        }

        inline bool AtMaxDepth() { return maxDepth_ == depth_; }

      private:
        int depth_;
        int maxDepth_;
        Atom* rootAtom_ = nullptr;
        Atom* endAtom_ = nullptr;
        bool branchFound;
        // AtomVector path_;
    };

    void FindEndsOfBranchesFromLinkageAtom(Atom* currentAtom, Atom* previousAtom, Residue* residue, Branch* branch);

    std::tuple<std::vector<DihedralAtoms>, std::vector<Atom*>> findRotatableDihedralsinBranchesConnectingResidues(
        const ResidueLink& link, const std::vector<Atom*>& residueCyclePoints);
} // namespace gmml

#endif
