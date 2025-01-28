#include <vector>
#include <string>
#include <iostream>

#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/graphVizDotConfig.hpp"

int main()
{
    std::string s1 = "DTalp[2S,3Me]a1-6DManpa1-6[DAllpb1-3][DNeup5Aca2-6DGalpb1-4DGlcp[3S]b1-2DAltpa1-4]DManpb1-4DGulp["
                     "6Me]b1-4DGlcpNAcb1-OH";
    std::string s2 = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-"
                     "4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string s3 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH";
    std::string s4 = "dUA[2S]1-4DGlcpNAc[3S,6S]a1-4LIdopA(2SO)[2S]a1-4DGlcpNSa1-4DGlcpA[2S]b1-4DGlcpAb1-3DGalpb1-"
                     "3DGalpb1-4DXylpb1-OH";
    std::string s5 = "dUA[2S]a1-4DGlcpNSa1-4LIdopA[2S]a1-4DGlcpNSa1-4LIdopA(4C1)a1-4DGlcpNS[6S]a1-OH";
    std::string s6 = "DGlcpa1-2DFrufb";
    std::string s7 = "DFrufb2-1DGlcpa";
    std::vector<std::string> sequences {s1, s2, s3, s4, s5, s6, s7};
    int count = 0;
    for (auto& sequence : sequences)
    {
        //    	CondensedSequence::SequenceManipulator manipulator(sequence);
        cdsCondensedSequence::GraphVizDotConfig config;
        config.file_name_ = std::to_string(count) + ".dot";
        CondensedSequence::drawGlycan(config, sequence);
        //    	std::cout << manipulator.PrintGraphViz(config) << std::endl;
        count++;
    }
    return 0;
}
