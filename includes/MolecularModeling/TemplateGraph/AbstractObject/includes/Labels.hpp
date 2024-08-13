#ifndef ABSTRACTOBJECT_INCLUDES_LABELS_HPP
#define ABSTRACTOBJECT_INCLUDES_LABELS_HPP

#include <algorithm>
#include <string>
#include <vector>

// TODO: Pick a better namespace name for this, figure out what we should use instead of `unsigned int`
namespace abstrab
{
    class Labels
    {
      public:
        /************************************************
         *  CONSTRUCTORS/DESTRUCTORS
         ***********************************************/
        inline Labels(const std::string name_t, const std::vector<std::string> labels_t)
            : name_m(name_t), labels_m(labels_t)
        {}

        /************************************************
         *  GETTER/SETTER
         ***********************************************/
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

        /************************************************
         *  MUTATORS
         ***********************************************/

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

        /************************************************
         *  FUNCTIONS
         ***********************************************/
        inline bool containsLabel(const std::string query) const
        {
            return std::find(labels_m.begin(), labels_m.end(), query) != labels_m.end();
        }

      private:
        /************************************************
         *  ATTRIBUTES
         ***********************************************/
        std::string name_m;
        std::vector<std::string> labels_m;
    };
} // namespace abstrab

#endif // LABELS_HPP
