#include "include/readers/Pdb/SectionClasses/remarkRecord.hpp"

#include "include/util/constants.hpp" //dNotSet
#include "include/util/logging.hpp"
#include "include/util/strings.hpp" //trim

namespace gmml
{
    namespace pdb
    {
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        RemarkRecord::RemarkRecord()
        {
            this->SetResolution(constants::dNotSet);
            this->SetBFactor(constants::dNotSet);
        }

        RemarkRecord::RemarkRecord(std::stringstream& stream_block)
        {
            this->SetResolution(constants::dNotSet);
            this->SetBFactor(constants::dNotSet);
            std::string line;
            getline(stream_block, line);
            std::string temp = line;
            while (!util::Trim(temp).empty())
            {
                if (line.find("REMARK") != std::string::npos)
                {
                    if (line.find("2 RESOLUTION.") != std::string::npos)
                    {
                        std::string tmp_resolution = line.substr(23, 7);
                        util::Trim(tmp_resolution);
                        try
                        {
                            this->SetResolution(std::stof(tmp_resolution));
                        }
                        catch (const std::invalid_argument& error)
                        {
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::ERR,
                                "RESOLUTION is not a valid float value. Value:\t" + tmp_resolution);
                        }
                    }
                    if (line.find("MEAN B VALUE") != std::string::npos)
                    {
                        int start = line.find(":") + 1;
                        std::string tmp_b_factor = line.substr(start, 80 - start);
                        util::Trim(tmp_b_factor);
                        try
                        {
                            this->SetBFactor(std::stof(tmp_b_factor));
                        }
                        catch (const std::invalid_argument& error)
                        {
                            util::log(
                                __LINE__,
                                __FILE__,
                                util::ERR,
                                "MEAN B VALUE is not a valid float value. Value:\t" + tmp_b_factor);
                        }
                    }
                }
                getline(stream_block, line);
                temp = line;
            }
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void RemarkRecord::SetResolution(const float resolution) { this->resolution_ = resolution; }

        void RemarkRecord::SetBFactor(const float b_factor) { this->b_factor_ = b_factor; }

        //////////////////////////////////////////////////////////
        //                      DISPLAY FUNCTION                //
        //////////////////////////////////////////////////////////
        void RemarkRecord::Print(std::ostream& out) const
        {
            out << "Resolution: " << this->GetResolution() << ". BFactor: " << this->GetBFactor() << "\n";
        }

        void RemarkRecord::Write(std::ostream& stream) const
        {
            stream << "REMARK   2\n";
            stream << "REMARK   2 RESOLUTION.  " << this->GetResolution() << " ANGSTROMS.\n";
            stream << "REMARK   3\n";
            stream << "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : " << this->GetBFactor() << "\n";
        }
    } // namespace pdb
} // namespace gmml
