#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_TEMPLATEDSELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_TEMPLATEDSELECTIONS_HPP

#include "includes/CodeUtils/containers.hpp"

#include <vector>

// ToDo move this to CentralDataStructure/Selections and rename the namespace.
namespace codeUtils
{
    template<class T>
    std::vector<T> findElementsNotInVector(const std::vector<T>& inputVector, const std::vector<T>& excludeElements)
    {
        std::vector<T> elementsInInputVectorButNotInQueryElements;
        for (auto element : inputVector)
        { // if element is not in the exclude list
            if (!codeUtils::contains(excludeElements, element))
            {
                elementsInInputVectorButNotInQueryElements.push_back(element);
            }
        }
        return elementsInInputVectorButNotInQueryElements;
    }

    template<class T>
    void findPathBetweenElementsInGraph(const T current, const T target, std::vector<T>& visited, std::vector<T>& path,
                                        bool& targetFound)
    {
        if (current == target)
        {
            targetFound = true;
            path.push_back(current);
            return;
        }
        visited.push_back(current);
        for (auto& neighbor : current->getNeighbors())
        {
            if (!(targetFound || codeUtils::contains(visited, neighbor)))
            {
                findPathBetweenElementsInGraph(neighbor, target, visited, path, targetFound);
            }
        }
        if (targetFound)
        {
            path.push_back(current);
        }
    }
} // namespace codeUtils
#endif
