#ifndef INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP

#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/randomDescent.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Geometry/overlap.hpp"
#include <string>

namespace glycoproteinBuilder
{
    using cds::Assembly;

    class GlycoproteinBuilder
    {
      public:
        GlycoproteinBuilder(GlycoproteinBuilderInputs inputStruct,
                            pdb::PreprocessorOptions preprocessingOptions = pdb::PreprocessorOptions());

        inline Assembly* getGlycoprotein()
        {
            return glycoprotein_;
        }

        void ResolveOverlaps(const std::string& outputDir, const std::vector<std::string>& headerLines, int numThreads);

      private:
        pdb::PdbFile pdbFile;
        std::vector<GlycosylationSite> glycosites_; // Info about each glycosylation site. See the class.
        std::vector<Residue*> proteinResidues_;
        Assembly* glycoprotein_; // Generated by this code.
        GlycoproteinBuilderInputs settings;
    };
} // namespace glycoproteinBuilder

#endif
