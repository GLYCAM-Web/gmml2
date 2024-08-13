#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSFUNCTIONS_CDSFUNCTIONS_HPP

#include <vector>

namespace cds
{

    template<typename T> void serializeNumbers(std::vector<T*> elements)
    {
        unsigned int i = 0;
        for (auto& element : elements)
        {
            element->setNumber(++i);
        }
        return;
    }
} // namespace cds
#endif
