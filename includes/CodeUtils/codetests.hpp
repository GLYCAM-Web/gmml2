#ifndef INCLUDES_CODEUTILS_CODETESTS_HPP
#define INCLUDES_CODEUTILS_CODETESTS_HPP
// ToDo This belongs in gmml/tests, doesn't need to include common and can be a function instead of a class.
#include <vector>
#include <string>

namespace CodeUtils
{
    class CodeTests
    {
      public:
        //////////////////////////////////////////////////////////
        //                       Constructor                    //
        //////////////////////////////////////////////////////////
        /*! \fn
         * Default constructor
         */
        CodeTests();

        //////////////////////////////////////////////////////////
        //                           ACCESSOR                   //
        //////////////////////////////////////////////////////////
        /** \addtogroup Code_Utils
         * @{
         */
        /*! \fn
         * Return a list of available tests
         */
        std::vector<std::string> ListCodeTests();
        /** @}*/
        //////////////////////////////////////////////////////////
        //                           MUTATOR                    //
        //////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////
        //                         FUNCTIONS                    //
        //////////////////////////////////////////////////////////

        void ProduceSegmentationFault();

        //////////////////////////////////////////////////////////
        //                     DISPLAY FUNCTIONS                //
        //////////////////////////////////////////////////////////

      private:
        //////////////////////////////////////////////////////////
        //                         ATTRIBUTES                   //
        //////////////////////////////////////////////////////////
    };
} // namespace CodeUtils

#endif
