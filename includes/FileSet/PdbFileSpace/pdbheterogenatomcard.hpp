// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENATOMCARD_HPP
#define PDBHETEROGENATOMCARD_HPP

#include <string>
#include <map>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAtom;

    class PdbHeterogenAtomCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef PdbAtom PdbHeterogenAtom;
            typedef std::map<int, PdbHeterogenAtom*> PdbHeterogenAtomMap;            

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeterogenAtomCard();
            PdbHeterogenAtomCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a heterogen atom card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the heterogen atoms in a heterogen atom card
              * @return heterogen_atoms_ attribute of the current object of this class
              */
            PdbHeterogenAtomMap GetHeterogenAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current heterogen atom card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            PdbHeterogenAtomMap heterogen_atoms_;

    };
}

#endif // PDBHETEROGENATOMCARD_HPP