#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequencePrinter.hpp"
#include "includes/CodeUtils/containers.hpp"
#include "includes/CodeUtils/strings.hpp"

#include <vector>
#include <string>
#include <memory>
#include <iostream>

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
    std::string s8 =
        "DNeup5Ac&Label=residue-9;a2-6&Label=link-7;DGalp&Label=residue-8;b1-4&Label=link-6;DGlcpNAc&Label=residue-6;["
        "3&Label=link-5;S&Label=residue-7;]b1-2&Label=link-4;DManp&Label=residue-5;a1-3&Label=link-3;[DGlcpNAc&Label="
        "residue-10;b1-4&Label=link-8;][DManp&Label=residue-12;[2&Label=link-11;S&Label=residue-13;,3&Label=link-12;Me&"
        "Label=residue-14;]a1-6&Label=link-10;DManp&Label=residue-11;a1-6&Label=link-9;]DManp&Label=residue-4;b1-4&"
        "Label=link-2;DGlcpNAc&Label=residue-3;[6&Label=link-13;Me&Label=residue-15;]b1-4&Label=link-1;DGlcpNAc&Label="
        "residue-2;b1-1&Label=link-0;-OH&Label=residue-1;";
    std::string s9  = "";
    std::string s10 = "There will be cake.";
    std::string s11 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH ";
    std::string s12 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3]DGalpb1-4DXylpb1-OH";
    std::string s13 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalp[Boo]b1-3DGalpb1-4DXylpb1-OH";
    std::string s14 = "dUA[2S]1-4DGlcpNAc[3S,6S]a1-4LIdopA(2SO)[2S]a1-4LIdopA(2SO)a1-4DGlcpNSa1-4DGlcpA[2S]b1-OH";
    std::string s15 = "DManpa1-2DManpa1-6[DManpa1-2DManpa1-3]DManpa1-6[DManpa1-2DManpa1-2DManpa1-3]DManpb1-4DGlcpNAcb1-"
                      "4DGlcpNAc[2Ac,3Me,6S]b1-OH";
    std::vector<std::string> sequences {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15};
    int count = 1;
    for (auto& sequence : sequences)
    {
        std::cout << "s" << count << ":\n" << sequence << "\n";
        try
        {
            cdsCondensedSequence::SequenceData sequenceData = cdsCondensedSequence::parseAndReorder(sequence);
            std::cout << "Index Ordered:\n"
                      << cdsCondensedSequence::printSequence(cdsCondensedSequence::defaultConfig(), sequenceData)
                      << "\n";
            std::cout << "Iupac:\n"
                      << cdsCondensedSequence::printSequence(cdsCondensedSequence::iupacConfig(), sequenceData) << "\n";
            std::cout << "Labeled:\n"
                      << cdsCondensedSequence::printSequence(cdsCondensedSequence::defaultConfigLabelled(),
                                                             sequenceData)
                      << "\n\n";
        }
        catch (std::runtime_error& error)
        {
            std::stringstream ss;
            ss << "Runtime error thrown for sequence s" << count << " : " << sequence << "\n" << error.what() << "\n\n";
            std::cout << ss.str();
        }
        catch (...)
        {
            std::stringstream ss;
            ss << "Unknown error type thrown for sequence s" << count << " : " << sequence << "\n\n";
            std::cout << ss.str();
        }
        count++;
    }
    std::cout << "Finished test\n\n";
    return 0;
}
