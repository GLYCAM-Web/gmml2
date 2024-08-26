#ifndef ABSTRACTOBJECT_INCLUDES_GENERIC_OBJECT_HPP
#define ABSTRACTOBJECT_INCLUDES_GENERIC_OBJECT_HPP

#include <vector>
#include <string>

namespace abstrab
{
    /* TODO: Figure out a better name and what data we want to include. Do
     * 			we want/need index counting? This could be a pain when
     * 			deletions and insertions come into play. All my algos (except
     * 			subgraph matching) ignore these signifiers.
     * 			We could name it "generic identifiers" or something of the
     * 			sort.
     */
    class AbstractObject
    {
      public:
        /************************************************
         *  CONSTRUCTORS/DESTRUCTORS
         ***********************************************/
        AbstractObject(const std::string& name, const std::vector<std::string>& labels)
            : index_m(generateIndex()), name_m(name), labels_m(labels)
        {}

        inline unsigned int getIndex() const
        {
            return index_m;
        }

        inline void setIndex(unsigned int index)
        {
            index_m = index;
        }

        inline std::string getName() const
        {
            return name_m;
        }

        inline std::string getLabel() const
        {
            return labels_m.empty() ? "" : labels_m.back();
        }

        inline std::vector<std::string> getLabels() const
        {
            return labels_m;
        }

        inline void setName(std::string name_t)
        {
            name_m = name_t;
        }

        inline void clearLabels()
        {
            labels_m.clear();
        }

        inline void setLabels(std::vector<std::string> labels_t)
        {
            labels_m = labels_t;
        }

        inline void addLabel(std::string label_t)
        {
            labels_m.push_back(label_t);
        }

        bool containsLabel(const std::string query) const;

      private:
        inline unsigned int generateIndex()
        {
            static unsigned int s_NodeIndex =
                0; // static keyword means it is created only once and persists beyond scope of code block.
            return s_NodeIndex++; // makes copy of index, increments the real index, then returns the value in the copy
        }

        unsigned int index_m;
        std::string name_m;
        std::vector<std::string> labels_m;
    };

} // namespace abstrab
#endif // ABSTRACT_OBJECT_HPP
