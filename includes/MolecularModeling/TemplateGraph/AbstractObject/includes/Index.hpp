#ifndef ABSTRACTOBJECT_INCLUDES_INDEX_HPP
#define ABSTRACTOBJECT_INCLUDES_INDEX_HPP

#include <string>
#include <vector>

namespace abstrab
{
    class Index
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        Index() : index_m(generateIndex())
        {}

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline unsigned int getIndex() const
        {
            return index_m;
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void setIndex(unsigned int index)
        {
            index_m = index;
        }

      private:
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        inline unsigned int generateIndex()
        {
            static unsigned int s_NodeIndex =
                0; // static keyword means it is created only once and persists beyond scope of code block.
            return s_NodeIndex++; // makes copy of index, increments the real index, then returns the value in the copy
        }

        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        unsigned int index_m;
    };

} // namespace abstrab
#endif // INDEX_HPP
