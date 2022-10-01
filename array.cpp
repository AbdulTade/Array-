#include <Windows.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <variant>
#include <random>
#include <cmath>
#include <time.h>
#include <memory>
#include <typeinfo>

namespace maths_suite
{
    template <typename _T>
    struct apply_t
    {
        /* function to apply to each array element */
        virtual _T apply(_T x);
    };

    /**
     * @brief templatized fixed-size array class. Template type parameter should be an intergral type
     * choices : int, bool e.t.c
     *
     * @tparam _T
     */
    template <typename _T = int>
    class array
    {

    private:
        size_t len = 0;
        _T *arr_ptr = nullptr;
        const std::map<std::string, std::string> type2name = {{"i", "int"}, {"c", "char"}, {"f", "float"}, {"d", "double"}, {"b", "bool"}};
        std::string _type;

        /* Swap two values */
        void swap(_T *_x, _T *_y)
        {
            _T tmp = *_x;
            *_x = *_y;
            *_y = tmp;
        }

        /* partition function for quicksort algorithm */
        int partition(_T *arr, int low, int high)
        {
            _T pivot = arr[high]; // pivot
            int i = (low - 1);

            for (int j = 0; j <= (high - 1); j++)
            {
                if (arr[j] < pivot)
                {
                    i++;
                    swap(&arr[i], &arr[j]);
                }
            }

            swap(&arr[i + 1], &arr[high]);
            return (i + 1);
        }

        void quicksort(_T *arr, int low, int high)
        {
            if (low < high)
            {
                /* p_indx is partitioning index, arr[p] is now at right place */
                int p_indx = partition(arr, low, high);
                quicksort(arr, low, p_indx - 1);
                quicksort(arr, p_indx + 1, high);
            }
        }

        void PopulateMapUnique()
        {
            int pos = 0;
            for (size_t i = 0; i < len; i++)
            {
                if (!dictionary.count(arr_ptr[i]))
                {
                    dictionary[arr_ptr[i]] = pos;
                    pos++;
                }
            }
        }

    public:
        std::map<_T, int> dictionary; /* to categorical dictionary */

        /* constant iterator for array container class */
        struct const_iterator
        {
            using reference_type  = _T&;
            using pointer_type    = _T*;
            using value_type      = _T;
            using difference_type = std::ptrdiff_t;

            const_iterator(_T *ptr) : a_ptr(ptr){};
            _T &operator*() const { return *a_ptr; }
            _T *operator->() const { return a_ptr; }
            const_iterator &operator++()
            {
                a_ptr++;
                return *this;
            }
            const_iterator operator++(int)
            {
                const_iterator tmp = *this;
                ++(*this);
                return tmp;
            }

            friend bool operator==(const const_iterator &a, const const_iterator &b)
            {
                return a.a_ptr == b.a_ptr;
            }

            friend bool operator!=(const const_iterator &a, const const_iterator &b)
            {
                return a.a_ptr != b.a_ptr;
            }

        private:
            _T *a_ptr; /* iterator pointer */
        };

        struct maths
        {
            /* Natural Logarithm */
            array<_T> &ln(const array<_T> &a)
            {
                bool inValid = false;
                 array<_T> r = array<_T>(a.len);
                for (size_t i = 0; i < a.len; i++)
                {
                    if (a.arr_ptr[i] <= 0)
                    {
                        a.arr_ptr[i] = (_T)(-1);
                        inValid = true;
                    }
                    r.arr_ptr[i] = log(a.arr_ptr[i]);
                }
                if (inValid)
                    std::cout << "warning: encountered some inconsistent values for ln(x)" << std::endl;
                return r;
            };

            /* exponential function */
            array<_T> &_exp(const array<_T> &a)
            {
                 array<_T> r = array<_T>(a.len);
                for (size_t i = 0; i < a.len; i++)
                {
                    r.arr_ptr[i] = exp(a.arr_ptr[i]);
                }
                return r;
            };

            /* Logarithm to base 10 */
            array<_T> &_log10(const array<_T> &a)
            {
                bool inValid = false;
                 array<_T> r = array<_T>(a.len);
                for (size_t i = 0; i < a.len; i++)
                {
                    if (a.arr_ptr[i] <= 0)
                    {
                        a.arr_ptr[i] = (_T)(-1);
                        inValid = true;
                    }
                    r.arr_ptr[i] = log10(a.arr_ptr[i]);
                }
                if (inValid)
                    std::cout << "warning: encountered some inconsistent values for log10(x)" << std::endl;
                return r;
            };

            array<_T> &_sqrt(const array<_T> &a)
            {
                bool inValid = false;
                 array<_T> r = array<_T>(a.len);
                for (size_t i = 0; i < a.len; i++)
                {
                    if (a.arr_ptr[i] < 0)
                    {
                        a.arr_ptr[i] = 0;
                        inValid = true;
                    }
                    r.arr_ptr[i] = sqrt(a.arr_ptr[i]);
                }
                if (inValid)
                    std::cout << "warning: encountered some inconsistent values for sqrt(x)" << std::endl;
                return r;
            };

            array<_T> &_abs(const array<_T> &a)
            {
                 array<_T> r = array<_T>(a.len);
                for (size_t i = 0; i < a.len; i++)
                    r.arr_ptr[i] = abs(a.arr_ptr[i]);
                return r;
            }
        };

        /* iterator pointer to first array element */
        const_iterator begin() { return const_iterator(arr_ptr); }
        /* iterator pointer to one element past last array element */
        const_iterator end() { return const_iterator(arr_ptr + len); }

        /* default constructor class */
        array() = default;

        /**
         * @brief Construct a new array object
         *
         * @param size
         */
        array(size_t size) : len{size}
        {
            arr_ptr = new _T[size];
            _type = typeid(arr_ptr[0]).name();
        }

        /* Array pointer should be allocated with new[] operator */
        array(_T *arr, size_t size) : len{size}
        {
            if (arr == nullptr)
                throw std::runtime_error("Pointer argument cannot be null");
            arr_ptr = arr;
            _type = typeid(arr_ptr[0]).name();
        };

        /**
         * @brief Construct a new array object from a std::vector<_T> container
         *
         * @param std::vector<_T>& vec
         */
        array(std::vector<_T> &vec) : len{vec.size()}
        {
            arr_ptr = new _T[len];
            for (size_t i = 0; i < len; i++)
                arr_ptr[i] = vec[i];
            _type = typeid(arr_ptr[0]).name();
        }

        /**
         * @brief Construct a new array object from another array object;
         *
         * @param array<_T> arr;
         */
        array(const array<_T> &a) : len{a.len}
        {
            arr_ptr = new _T[len];
            memcpy(arr_ptr, a.arr_ptr, sizeof(_T) * len);
            _type = typeid(arr_ptr[0]).name();
        }

        _T operator[](size_t indx)
        {
            if (indx > len - 1)
                throw std::runtime_error("Out of bounds access to array");
            return arr_ptr[indx];
        }

        /* Get an array slice. range should be delimeted by ":" */
        array<_T> &operator[](const std::string &_Range)
        {
            std::vector<std::string> str_vector;
            std::string token;
            char delimeter = ':';

            std::stringstream ss{_Range};

            while (std::getline(ss, token, delimeter))
                str_vector.push_back(token);

            long long start = atoll(str_vector[0].c_str()); /* start of */
            long long end = atoll(str_vector[1].c_str());

            if (start < 0 || end > len)
                throw std::runtime_error("Out of bound slice");

             array<_T> arr = array<_T>(end - start);
            for (size_t i = 0; i < (end - start); i++)
                arr.arr_ptr[i] = arr_ptr[start + i];
            return arr;
        }

        /* make sure type passed supports + operation */
        array<_T> &operator+(const array<_T> &a)
        {
            if (a.len != len)
                throw std::logic_error("cannot add arrays of unequal dimensions");
             array<_T> arr = array<_T>(a.len);
            for (size_t j = 0; j < a.len; j++)
                arr.arr_ptr[j] = a.arr_ptr[j] + arr_ptr[j];
            return arr;
        }

        array<_T> &operator+(const _T s)
        {
             array<_T> a = array<_T>(len);
            for (size_t i = 0; i < len; i++)
                a.arr_ptr[i] = a.arr_ptr[i] + s;
            return a;
        }

        friend array<_T> &operator+(const _T s, const array<_T> &a)
        {
             array<_T> sum = array<_T>(a.len);
            for (size_t i = 0; i < a.len; i++)
                sum.arr_ptr[i] = a.arr_ptr[i] + s;
            return sum;
        }

        /* Make sure type passed supports "-" operation */
        array<_T> &operator-(const array<_T> &a)
        {
            if (a.len != len)
                throw std::logic_error("cannot add arrays of unequal dimensions");
             array<_T> arr = array<_T>(a.len);
            for (size_t j = 0; j < a.len; j++)
                arr.arr_ptr[j] = arr_ptr[j] - a.arr_ptr[j];
            return arr;
        }

        /* make sure type passed supports the * operation */
        array<_T> &operator*(_T factor)
        {
             array<_T> arr = array<_T>(len);
            for (size_t j = 0; j < len; j++)
                arr.arr_ptr[j] = arr_ptr[j] * factor;
            return arr;
        }

        /* make sure type supports / operation */
        array<_T> &operator/(_T divisor)
        {
             array<_T> arr = array<_T>(len);
            for (size_t j = 0; j < len; j++)
                arr.arr_ptr[j] = arr_ptr[j] / divisor;
            return arr;
        }

        /* Dot product of two vectors */
        _T dot(const array<_T> &a)
        {
            _T d = 0;
            for (size_t j = 0; j < a.len; j++)
                d += arr_ptr[j] * a.arr_ptr[j];
            return d;
        }

        /* Hadamard Product of two vectors */
        array<_T> &operator*(const array<_T> &a)
        {
            if (a.len != len)
                throw std::logic_error("cannot add arrays of unequal dimensions");
             array<_T> arr = array<_T>(a.len);
            for (size_t j = 0; j < a.len; j++)
                arr.arr_ptr[j] = arr_ptr[j] * a.arr_ptr[j];
            return arr;
        }

        /* How to find the nth norm of a number */
        _T norm(int n)
        {
            _T sum = 0;
            if (n <= 0)
                throw std::invalid_argument("invalid: there cannot be a zero norm");
            if (n == 1)
                return this->sum();

            for (size_t i = 0; i < len; i++)
                sum += pow(arr_ptr[i], n);

            return sqrt(sum);
        }

        friend array<_T> operator*(const _T m, const array<_T> &a)
        {
            array<_T> arr = array<_T>(a.len);
            return arr * m;
        }

        /* Equality operator */
        bool operator==(const array<_T> &a)
        {
            if (a.len != len)
                return false;
            for (size_t k = 0; k < len; k++)
            {
                if (a.arr_ptr[k] != arr_ptr[k])
                    return false;
            }
            return true;
        }

        /* inequality operator */
        bool operator!=(const array<_T> &a)
        {
            return !operator==(a);
        }

        /* type name as a string */
        std::string dtype()
        {
            return type2name.at(_type);
        }

        /* Length of the array */
        size_t length()
        {
            return len;
        }

        /* Assignment operator */
        array<_T> &operator=(const array<_T> &a)
        {
            return *this;
        }

        /* sum of all elements of array container class */
        _T sum()
        {
            _T s = 0;
            for (size_t i = 0; i < len; i++)
                s += arr_ptr[i];
            return s;
        }

        /* arithmetic mean */
        _T mean()
        {
            _T s = this->sum();
            return s / len;
        }

        array<_T> &mode()
        {
            std::vector<_T> frequencies;
            std::vector<_T> modes;

            PopulateMapUnique();
            for (auto &[k, v] : dictionary)
                v = 0;

            for (size_t n = 0; n < len; n++)
                dictionary[arr_ptr[n]]++;

            for (auto &[k, v] : dictionary)
                frequencies.push_back(v);

            auto a = array<_T>(frequencies);
            _T m = a._max();

            for (auto &[k, v] : dictionary)
            {
                if (v == m)
                    modes.push_back(k);
            }

             array<_T> _mode = array<_T>(modes);
            return _mode;
        }

        /* vector power */
        array<_T> &Pow(int _p)
        {
             array<_T> a = array<_T>(len);
            for (size_t j = 0; j < len; j++)
                a.arr_ptr[j] = (_T)(pow((double)arr_ptr[j], (double)_p));
            return a;
        }

        /* Reverse the array */
        array<_T> &reverse()
        {
             array<_T> arr = array<_T>(len);
            for (size_t n = 0; n < len; n++)
                arr.arr_ptr[n] = arr_ptr[len - n - 1];
            return arr;
        }

        /**
         * @brief applies a function to every element and returns an array class with the result;
         *
         */
        array<_T> &apply(apply_t<_T> &_func)
        {
             array<_T> arr = array<_T>(len);
            for (size_t i = 0; i < len; i++)
                arr.arr_ptr[i] = _func.apply(arr_ptr[i]);
            return arr;
        }

        /* Converts an array to its categorical representation using the map given */
        array<int> &toCategorical(std::map<_T, int> &_map)
        {
             array<int> a = array<int>(len);
            for (size_t m = 0; m < len; m++)
                a.insert(m, _map[arr_ptr[m]]);
            return a;
        }

        /* Converts an array to its categorical representation using an internally generated map */
        array<int> &toCategorical()
        {
             array<int> a = array<int>(len);
            PopulateMapUnique();
            for (size_t m = 0; m < len; m++)
                a.insert(m, dictionary[arr_ptr[m]]);
            return a;
        }

        /* perform a quicksort of the array */
        void qsort()
        {
            quicksort(arr_ptr, 0, len - 1);
        }

        /* Zeros the memory */
        void zeros()
        {
            for (size_t n = 0; n < len; n++)
                arr_ptr[n] = 0;
        }

        /* Returns c-style representation of array */
        _T *c_arr()
        {
            _T *new_arr = new _T[len];
            memcpy(new_arr, arr_ptr, sizeof(_T) * len);
            return new_arr;
        }

        /* concatenate two arrays */
        array<_T> &concat(const array<_T> &a)
        {
             array<_T> c = array<_T>(len + a.len);
            size_t pos = 0;
            for (; pos < len; pos++)
                c.arr_ptr[pos] = arr_ptr[pos];

            for (; pos < (a.len + len); pos++)
                c.arr_ptr[pos] = a.arr_ptr[pos - len];

            return c;
        }

        /**
         * @brief insert a value at a given position in the array.
         *        Usage:
         *              array.insert(size_t pos,_T value);
         * @param pos   position index
         * @param value value to insert at position
         */
        void insert(size_t pos, _T value)
        {
            if (pos >= len)
                throw std::runtime_error("Out of bound access to array");
            arr_ptr[pos] = value;
        }

        /* insertion operator */
        friend std::ostream &operator<<(std::ostream &oss, const array<_T> &a)
        {
            std::string _type = typeid(a.arr_ptr[0]).name();
            oss << "array {\n\t"
                << "dtype : " << a.type2name.at(_type) << "\n"
                << "\tsize : " << a.len << "\n}" << std::endl;
            return oss;
        }

        /* print first few elements of array */
        void head(size_t count = 10)
        {
            size_t tcount = (count > 10) ? 10 : count;
            std::cout << "array "
                      << "{"
                      << " ";

            for (size_t c = 0; c < ((len < tcount) ? len : tcount); c++)
                std::cout << arr_ptr[c] << " , ";

            std::cout << (((int)(len - tcount) < 0) ? 0 : (len - tcount)) << " more items ... "
                      << "}" << std::endl;
        }

        /* print last few elements of array */
        void tail(size_t count = 10)
        {
            size_t tcount = (count > 10) ? 10 : count;
            std::cout << "array "
                      << "{"
                      << " ... ";
            std::cout << (((int)(len - tcount) < 0) ? len : (len - tcount)) << " items, ";

            for (size_t c = (((int)(len - tcount) < 0) ? 0 : (len - tcount)); c < len; c++)
                std::cout << arr_ptr[c] << " , ";

            std::cout << "}" << std::endl;
        }

        void *operator new[](size_t size) noexcept
        {
            void *ret_ptr = malloc(size * sizeof(array<_T>));
            if (!ret_ptr)
                return nullptr;
            return ret_ptr;
        }

        void operator delete[](void *ptr)
        {
            if (ptr == nullptr)
                return;
            free(ptr);
        }

        /* minimum value of array */
        _T _min()
        {
            _T minimum = arr_ptr[0];
            for (size_t i = 1; i < len; i++)
            {
                if (arr_ptr[i] < minimum)
                    minimum = arr_ptr[i];
            }

            return minimum;
        }

        /* maximum value of array */
        _T _max()
        {
            _T maximum = arr_ptr[0];
            for (size_t i = 1; i < len; i++)
            {
                if (arr_ptr[i] > maximum)
                    maximum = arr_ptr[i];
            }
            return maximum;
        }

        /* class destructor */
        ~array() noexcept
        {
            delete[] arr_ptr;
        }
    };

    /**
     * @brief Matrix class
     *
     */
    template <typename T>
    class Matrix
    {
    private:
        size_t n_rows = 0;
        size_t n_cols = 0;
        std::vector<array<T>> mat;

    public:
        array<size_t> &shape()
        {
             auto a = array<size_t>(2);
            a.insert(0, n_rows);
            a.insert(1, n_cols);
            return a;
        }

        Matrix(size_t n_rows, size_t n_cols) : n_rows(n_rows), n_cols(n_cols)
        {
            for (size_t i = 0; i < n_rows; i++)
            {
                 auto a = array<T>(n_cols);
                a.zeros();
                mat.push_back(a);
            }
        };

        Matrix<T> &operator+(const Matrix<T> &m)
        {
            if ((m.n_cols != n_cols) || (m.n_rows != n_rows))
                throw std::runtime_error("Cannot add two matrices of unequal dimension");
             auto m_ret = Matrix<T>(n_cols, n_rows);
            for (size_t k = 0; k < n_rows; k++)
            {
                for (size_t j = 0; j < n_cols; j++)
                    m_ret.mat[k].insert(j, m.mat[k][j] + mat[k][j]);
            }
            return m_ret;
        }

        Matrix<T> &operator-(const Matrix<T> &m)
        {
            if ((m.n_cols != n_cols) || (m.n_rows != n_rows))
                throw std::runtime_error("Cannot add two matrices of unequal dimension");
             auto m_ret = Matrix<T>(n_cols, n_rows);
            for (size_t k = 0; k < n_rows; k++)
            {
                for (size_t j = 0; j < n_cols; j++)
                    m_ret.mat[k].insert(j, m.mat[k][j] - mat[k][j]);
            }
            return m_ret;
        }

        Matrix<T> &operator*(const Matrix<T> &m)
        {
            if ((n_cols != m.n_rows))
                throw std::runtime_error("Cannot multiply matrices of unmatching row and column dimensions");
             auto m_ret = Matrix<T>(n_rows, m.n_cols);
            for (size_t i = 0; i < n_rows; i++)
            {
                T tmp = 0;
                for (size_t j = 0; j < m.n_cols; j++)
                {
                    for (size_t k = 0; k < n_cols; k++)
                        tmp += mat[i][k] * m.mat[k][j];
                    m_ret.mat[i].insert(j, tmp);
                }
            }
            return m_ret;
        }

        /* Returns a submatrix in the specified column and row range exclusive*/
        Matrix<T> &submatrix(size_t rbegin_indx, size_t rend_indx, size_t cbegin_indx, size_t cend_indx)
        {
            bool is_subscriptable = ((rend_indx <= n_rows) && (cend_indx <= n_rows)) && ((rbegin_indx < rend_indx) && (cbegin_indx < rend_indx));
            if (!is_subscriptable)
                throw std::runtime_error("cannot get submatrix for index out of range");
            size_t mat_rows = rend_indx - rbegin_indx;
            size_t mat_cols = cend_indx - cbegin_indx;
             auto m = Matrix<T>(mat_rows, mat_cols);
            for (size_t i = 0; i < mat_rows; i++)
            {
                for (size_t j = 0; j < mat_cols; j++)
                    m.mat[i].insert(j, mat[rbegin_indx + i][cbegin_indx + j]);
            }
            return m;
        }

        /**
         * @brief erases an entire row.
         *
         * @param r_pos
         * @return Matrix<T>&
         */
        Matrix<T> &rerase(size_t r_pos)
        {
            if (r_pos >= n_rows)
                throw std::runtime_error("Out of bounds access to matrix");
             Matrix<T> m = Matrix<T>(n_rows-1, n_cols);
            size_t r_m = 0;
            for (size_t i = 0; i < n_rows; i++)
            {
                if (i == r_pos)
                    continue;
                for (size_t j = 0; j < n_cols; j++)
                    m[r_m].insert(j, mat[i][j]);
                r_m++;
            }

            return m;
        }

        /**
         * @brief erase an entire column
         * 
         * @param c_pos 
         * @return Matrix<T>& 
         */
        Matrix<T> &cerase(size_t c_pos)
        {
            if(c_pos >= n_cols)
                throw std::runtime_error("Out of bound access to array");
            size_t c_m = 0;
             Matrix<T> m = Matrix<T>(n_rows,n_cols-1);
            for(size_t i = 0; i < n_rows; i++)
            {
                for(size_t j = 0; j < n_cols; j++)
                {
                    if(j == c_m)
                        continue;
                    m.mat[i].insert(c_m,mat[i][j]);
                    c_m++;
                }
                c_m = 0;
            }
            return m;
        }

        /**
         * @brief Returns a new matrix with the row and column shown removed
         *
         * @param r_pos
         * @param c_pos
         * @return reference to a newly created matrix
         */
        Matrix<T> &erase(size_t r_pos, size_t c_pos)
        {
            auto c_matrix = this->rerase(r_pos);
             auto m = c_matrix.cerase(c_pos);
            return m;
        }

        /** @brief  the determinant of a matrix
         * Use formula
         * det(A) = (-1)^(i+j) * aij * Mij
         * @return T
         */
        T det()
        {
            T d = 0;
            if (n_cols != n_rows)
                throw std::runtime_error("Matrix must be square.");
            return d;
        }

        void _rand()
        {
            time_t t = time(NULL);
            srand(t);
            zeros();
            for (size_t i = 0; i < n_rows; i++)
            {
                for (size_t j = 0; j < n_cols; j++)
                {
                    mat[i].insert(j, (T)rand());
                }
            }
        }

        void zeros()
        {
            for (size_t k = 0; k < n_rows; k++)
            {
                for (size_t j = 0; j < n_cols; j++)
                    mat[k].insert(j, 0);
            }
        }

        array<T> &operator[](size_t indx)
        {
            if (indx >= n_rows)
                throw std::runtime_error("Can not access matrix row out of bounds");
             auto a = array<T>(mat[indx]);
            return a;
        }

        friend std::ostream &operator<<(std::ostream &oss, Matrix<T> &m)
        {
            oss << "Matrix {\n\t"
                << "dtype : " << m.mat[0].dtype() << "\n\tshape : ( " << m.n_rows << " , " << m.n_cols << " )\n}";
            return oss;
        }
    };

}


template <typename T>
struct iterator_traits
{
    typedef T value_type;
    typedef T &reference_type;
    typedef T *pointer_type;
    typedef std::ptrdiff_t difference_type;
};

int main()
{
    auto m = maths_suite::Matrix<double>(10,10);
    m._rand();
    std::vector<double> v = {1,2,4,3,4,3,5,3,5,4,3,5,5,3,3,4};
    auto a = maths_suite::array<double>(v);
    auto it = a.begin();

    while(++it != a.end())
    {
        maths_suite::array<double>::const_iterator::value_type x = *it;
        std::cout << x << std::endl;
    }
}