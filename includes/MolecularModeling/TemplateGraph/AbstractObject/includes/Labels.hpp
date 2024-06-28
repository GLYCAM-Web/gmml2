#ifndef ABSTRACTOBJECT_INCLUDES_LABELS_HPP
#define ABSTRACTOBJECT_INCLUDES_LABELS_HPP

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
        inline Labels() : name_m(""), labels_m({""})
        {}

        inline Labels(std::string name_t) : name_m(name_t), labels_m({name_t})
        {}

        inline Labels(std::vector<std::string> labels_t) : name_m(""), labels_m(labels_t)
        {}

        inline Labels(std::string name_t, std::string label_t) : name_m(name_t), labels_m({label_t})
        {}

        inline Labels(std::string name_t, std::vector<std::string> labels_t) : name_m(name_t), labels_m(labels_t)
        {}

        inline ~Labels()
        {}

        // copy constructor
        inline Labels(const Labels& rhs) : name_m(rhs.name_m), labels_m(rhs.labels_m)
        {}

        // move constructor
        inline Labels(Labels&& rhs) : name_m(rhs.name_m), labels_m(rhs.labels_m)
        {}

        // copy assignment
        inline Labels& operator=(const Labels& rhs)
        {
            // this->name_m = rhs.name_m;
            // this->labels_m = rhs.labels_m;
            return *this = Labels(rhs);
        }

        // move assignment
        Labels& operator=(Labels&& rhs)
        {
            this->name_m   = rhs.name_m;
            this->labels_m = rhs.labels_m;
            return *this;
        }

        /************************************************
         *  GETTER/SETTER
         ***********************************************/
        inline std::string getName() const
        {
            return this->name_m;
        }

        inline std::string getLabel() const
        {
            if (this->labels_m.empty())
            {
                return "";
            }
            return this->labels_m.back();
        }

        inline std::vector<std::string> getLabels() const
        {
            return this->labels_m;
        }

        inline void setName(std::string name_t)
        {
            this->name_m = name_t;
        }

        inline void setLabels(std::string label_t)
        {
            this->labels_m.clear();
            this->addLabel(label_t);
        }

        inline void setLabels(std::vector<std::string> labels_t)
        {
            this->labels_m = labels_t;
        }

        inline void clearLabels()
        {
            this->labels_m.clear();
        }

        /************************************************
         *  MUTATORS
         ***********************************************/
        inline void addLabel(std::string label_t)
        {
            this->labels_m.push_back(label_t);
        }

        inline void addLabels(std::vector<std::string> labels)
        {
            this->labels_m.insert(this->labels_m.end(), labels.begin(), labels.end());
        }

        /************************************************
         *  FUNCTIONS
         ***********************************************/
        inline bool containsLabel(const std::string query) const
        {
            for (std::string currLabel : this->labels_m)
            {
                // TODO: Why was regex used in your version?
                if (currLabel == query)
                {
                    return true;
                }
            }
            return false;
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
