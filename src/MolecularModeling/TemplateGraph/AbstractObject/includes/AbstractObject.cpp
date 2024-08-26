#include "includes/MolecularModeling/TemplateGraph/AbstractObject/includes/AbstractObject.hpp"
#include "includes/CodeUtils/containers.hpp"

#include <string>

namespace abstrab
{

    bool AbstractObject::containsLabel(const std::string query) const
    {
        return codeUtils::contains(labels_m, query);
    }

} // namespace abstrab
