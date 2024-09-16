#ifndef INCLUDES_CODEUTILS_REFERENCES_HPP
#define INCLUDES_CODEUTILS_REFERENCES_HPP

#include <cstddef>
#include <vector>

namespace codeUtils
{
    template<typename T> class reference
    {
      public:
        reference(size_t index_, std::vector<T>* vector_) : index(index_), vector(vector_) {};

        T& get() const
        {
            return (*vector)[index];
        }

        void set(const T& value)
        {
            (*vector)[index] = value;
        }

        bool invalid() const
        {
            return (vector == nullptr) || (index >= vector->size());
        }

        bool operator==(const reference<T>& a) const
        {
            return (vector == a.vector) && (index == a.index);
        }

        bool operator!=(const reference<T>& a) const
        {
            return !operator==(a);
        }

      private:
        size_t index;
        std::vector<T>* vector;
    };

} // namespace codeUtils
#endif
