// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBRESIDUEMODIFICATION_HPP
#define PDBRESIDUEMODIFICATION_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbResidueModification
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbResidueModification();
            /*! \fn
              * Constructor with required parameters
              * @param id_code
              * @param residue_name
              * @param chain_identifier
              * @param sequence_number
              * @param insertion_code
              * @param standard_residue_name
              * @param dscr
              */
            PdbResidueModification(const std::string& id_code, const std::string& residue_name, char chain_identifier, int sequence_number,
                                   char insertion_code, const std::string& standard_residue_name, const std::string& dscr);
            PdbResidueModification(std::stringstream& stream_block);
            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the id code in a residue modification
              * @return id_code_ attribute of the current object of this class
              */
            std::string GetIdCode();
            /*! \fn
              * An accessor function in order to access to the residue name in a residue modification
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the chain identifier in a residue modification
              * @return chain_identifier_ attribute of the current object of this class
              */
            char GetChainId();
            /*! \fn
              * An accessor function in order to access to the sequence number in a residue modification
              * @return sequence_number_ attribute of the current object of this class
              */
            int GetSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the insertion code in a residue modification
              * @return insertion_code_ attribute of the current object of this class
              */
            char GetInsertionCode();
            /*! \fn
              * An accessor function in order to access to the standard residue name in a residue modification
              * @return standard_residue_name_ attribute of the current object of this class
              */
            std::string GetStandardResidueName();
            /*! \fn
              * An accessor function in order to access to the description in a residue modification
              * @return dscr_ attribute of the current object of this class
              */
            std::string GetDscr();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the id code of the current object
              * Set the id_code_ attribute of the current residue modification
              * @param id_code The id code of the current object
              */
            void SetIdCode(const std::string id_code);
            /*! \fn
              * A mutator function in order to set the residue_name of the current object
              * Set the residue_name_ attribute of the current residue modification
              * @param residue_name The residue name of the current object
              */
            void SetResidueName(const std::string residue_name);
            /*! \fn
              * A mutator function in order to set the chain identifier of the current object
              * Set the chain_identifier_ attribute of the current residue modification
              * @param chain_identifier The chain identifier of the current object
              */
            void SetChainId(char chain_identifier);
            /*! \fn
              * A mutator function in order to set the sequence number of the current object
              * Set the sequence_number_ attribute of the current residue modification
              * @param sequence_number The sequence number of the current object
              */
            void SetSequenceNumber(int sequence_number);
            /*! \fn
              * A mutator function in order to set the insertion code of the current object
              * Set the insertion_code_ attribute of the current residue modification
              * @param insertion_code The insertion code of the current object
              */
            void SetInsertionCode(char insertion_code);
            /*! \fn
              * A mutator function in order to set the standard residue name of the current object
              * Set the standard_residue_name_ attribute of the current residue modification
              * @param standard_residue_name The standard residue name of the current object
              */
            void SetStandardResidueName(const std::string standard_residue_name);
            /*! \fn
              * A mutator function in order to set the describtion of the current object
              * Set the dscr_ attribute of the current residue modification
              * @param dscr The dscribtion of the current object
              */
            void SetDscr(const std::string dscr);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string id_code_;
            std::string residue_name_;
            char chain_identifier_;
            int sequence_number_;
            char insertion_code_;
            std::string standard_residue_name_;
            std::string dscr_;
    };
}

#endif // PDBRESIDUEMODIFICATION_HPP