// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

#pragma once

#ifdef SAFE_PTR_DEBUG
    #define SAFE_PTR_DEBUG_BOOL 1
#else
    #define SAFE_PTR_DEBUG_BOOL 0
#endif

#ifdef SAFE_PTR_TEST
    #define SAFE_PTR_TEST_BOOL 1
#else
    #define SAFE_PTR_TEST_BOOL 0
#endif

#define SAFE_PTR_WARNING(msg) _warning(msg, __FILE__, __LINE__, __func__)

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <iterator>
#include <type_traits>
#if SAFE_PTR_DEBUG_BOOL
    #include <unordered_map>
    #include <mutex>
#endif

namespace fz {

#if SAFE_PTR_TEST_BOOL
    struct _SafePtrWarning {};
#endif

template<typename T>
class SafePtr
{
public:
    // constructor
    SafePtr() {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = true;
        #endif
    }

    // constructor
    SafePtr(const size_t size) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = false;
        #endif
        _begin = new T[size];
        _end = _begin + size;
    }

    // constructor
    SafePtr(const size_t size, const T value) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = false;
        #endif
        _begin = new T[size];
        _end = _begin + size;
        fill(value);
    }

    // constructor
    SafePtr(const std::initializer_list<T>& il) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = false;
        #endif
        _begin = new T[il.size()];
        _end = _begin + il.size();
        std::copy(il.begin(), il.end(), this->_begin);
    }

    // constructor
    template<
        typename InputIt,
        typename std::enable_if<
            !std::is_integral<InputIt>::value, int
        >::type = 0
    >
    SafePtr(InputIt first, InputIt last) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            _memory_id = _get_new_memory_id();
            _ref_count[_memory_id] = 1;
            _is_deleted[_memory_id] = false;
        #endif
        _construct_from_range(first, last, _sp_has_subtraction<InputIt>{});
    }

    // destructor
    ~SafePtr() noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            if (_get_is_view()) {
                return;
            }
            --_get_ref_count();
            if (_get_ref_count() == 0) {
                if(!_get_is_deleted()) {
                    SAFE_PTR_WARNING("Memory was leaked.");
                }
                _ref_count.erase(_memory_id);
                _is_deleted.erase(_memory_id);
            }
        #endif
    }

    // copy constructor
    SafePtr(const SafePtr& other) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_use_after_free();
            this->_memory_id = _get_new_memory_id();
            _ref_count[this->_memory_id] = 1;
            _is_deleted[this->_memory_id] = false;
        #endif
        this->_begin = new T[other.size()];
        this->_end = this->_begin + other.size();
        std::copy(other.begin(), other.end(), this->_begin);
    }
    
    // move constructor
    SafePtr(SafePtr&& other) noexcept(!SAFE_PTR_TEST_BOOL) {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_use_after_free();
            this->_memory_id = other._memory_id;
            if (!_get_is_view()) {
                ++_get_ref_count();
            }
        #endif
        this->_begin = other._begin;
        this->_end = other._end;
    }

    // copy assignment operator
    SafePtr& operator=(const SafePtr& other) {
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            if (this != &other) {
        #endif
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_use_after_free();
            if (!_get_is_view()) {
                --_get_ref_count();
                if (_get_ref_count()==0 && !_get_is_deleted()) {
                    SAFE_PTR_WARNING("Memory was leaked.");
                }
            }
            this->_memory_id = _get_new_memory_id();
            _ref_count[this->_memory_id] = 1;
            _is_deleted[this->_memory_id] = false;
        #endif
        this->_begin = new T[other.size()];
        this->_end = this->_begin + other.size();
        std::copy(other.begin(), other.end(), this->_begin);
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            }
        #endif
        return *this;
    }

    // move assignment operator
    SafePtr& operator=(SafePtr&& other) noexcept(!SAFE_PTR_TEST_BOOL) {
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            if (this != &other) {
        #endif
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            other._check_for_use_after_free();
            if (!_get_is_view()) {
                --_get_ref_count();
                if (_get_ref_count()==0 && !_get_is_deleted()) {
                    SAFE_PTR_WARNING("Memory was leaked.");
                }
            }
            this->_memory_id = other._memory_id;
            if (!_get_is_view()) {
                ++_get_ref_count();
            }
        #endif
        this->_begin = other._begin;
        this->_end = other._end;
        #ifndef SAFE_PTR_DISABLE_SELF_ASSIGNING_CHECKING
            }
        #endif
        return *this;
    }

    void free() const {
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            if (_get_is_deleted()) {
                throw std::logic_error(
                    "it was tried to free the memory of a SafePtr that does "
                    "not own data."
                );
            }
            if (_get_is_view()) {
                throw std::logic_error(
                    "it was tried to free the memory of a view"
                );
            }
            _get_is_deleted() = true;
        #endif
        delete[] _begin;
    }

    static SafePtr<T> make_view(T* const data, const size_t size) {
        SafePtr<T> safe_ptr;
        #if SAFE_PTR_DEBUG_BOOL
            std::lock_guard<std::mutex> lock(_mtx);
            safe_ptr._memory_id = 0;
        #endif
        safe_ptr._begin = data;
        safe_ptr._end = data + size;
        return safe_ptr;
    }

    size_t size() const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return _end - _begin;
    }

    const T* begin() const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return _begin;
    }

    T* begin() {
        return const_cast<T*>(
            const_cast<const SafePtr<T>&>(*this).begin()
        );
    }

    const T* cbegin() const {
        return begin();
    }

    const T* cbegin() {
        return const_cast<const SafePtr<T>&>(*this).cbegin();
    }

    const T* end() const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return _end;
    }

    T* end() {
        return const_cast<T*>(
            const_cast<const SafePtr<T>&>(*this).end()
        );
    }

    const T* cend() const {
        return end();
    }

    const T* cend() {
        return const_cast<const SafePtr<T>&>(*this).cend();
    }

    const T& operator[](const size_t index) const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return *(_begin + index);
    }

    T& operator[](const size_t index) {
        return const_cast<T&>(
            const_cast<const SafePtr<T>&>(*this)[index]
        );
    }

    const T& at(const size_t index) const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        if (index >= this->size()) {
            throw std::out_of_range(
                "tried to access SafePtr element out of range"
            );
        }
        return *(_begin + index);
    }

    T& at(const size_t index) {
        return const_cast<T&>(
            const_cast<const SafePtr<T>&>(*this).at(index)
        );
    }

    bool empty() const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return this->size() == 0;
    }

    const T* data() const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return _begin;
    }

    T* data() {
        return const_cast<T*>(
            const_cast<const SafePtr<T>&>(*this).data()
        );
    }

    const T& front() const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return *(this->_begin);
    }

    T& front() {
        return const_cast<T&> (
            const_cast<const SafePtr<T>&>(*this).front()
        );
    }

    const T& back() const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        return *(this->_end-1);
    }

    T& back() {
        return const_cast<T&>(
            const_cast<const SafePtr<T>&>(*this).back()
        );
    }

    void fill(const T& value) {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        for (auto& p : *this) {
            p = value;
        }
    }

    void
    print_all(const char* const variable_name = "SafePtr::print_all") const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif
        std::cout << variable_name << ": {\n";
        if (!empty()) {
            size_t idx = 0;
            for (const T* it = cbegin(); it != cend()-1; ++it) {
                std::cout << "    " << idx << ": " << *it << ",\n";
                ++idx;
            }
            std::cout << "    " << idx << ": " << back() << "\n";
        }
        std::cout << "}\n";
    }

    void print(const char* const variable_name = "SafePtr::print") const {
        #if SAFE_PTR_DEBUG_BOOL
            _check_for_use_after_free();
        #endif

        constexpr size_t MAX_ELEMENTS_TO_PRINT = 26;
        if (size() <= MAX_ELEMENTS_TO_PRINT) {
            print_all(variable_name);
            return;
        }
        
        constexpr size_t ELEMENTS_BEFORE_DOTS = 13;
        constexpr size_t ELEMENTS_AFTER_DOTS = 12;
        std::cout << variable_name << ": {\n";
        size_t idx = 0;
        for (const T* it=cbegin(); it!=cbegin()+ELEMENTS_BEFORE_DOTS; ++it) {
            std::cout << "    " << idx << ": " << *it << ",\n";
            ++idx;
        }
        idx = size() - ELEMENTS_AFTER_DOTS;
        std::cout << "    ...\n";
        for (const T* it = cend()-ELEMENTS_AFTER_DOTS; it != cend()-1; ++it) {
            std::cout << "    " << idx << ": " << *it << ",\n";
            ++idx;
        }
        std::cout << "    " << idx << ": " << back() << "\n";
        std::cout << "}\n";
    }

private:
    T* _begin; // points to the first element
    T* _end;   // points to the byte after the last byte of the last element

    // C++11-compatible replacement for std::void_t (which is C++17)
    template<typename...>
    struct _sp_void { using type = void; };
    template<typename... Ts>
    using _sp_void_t = typename _sp_void<Ts...>::type;

    // Detects whether subtracting the iterators to get the size is valid.
    template<typename It, typename = void>
    struct _sp_has_subtraction : std::false_type {};

    template<typename It>
    struct _sp_has_subtraction<
        It,
        _sp_void_t<decltype(std::declval<It>() - std::declval<It>())>
    > : std::true_type {};

    // Allocates memory by subtracting the iterators to get the size (faster).
    template<typename InputIt>
    void _construct_from_range(
        InputIt first, InputIt last, std::true_type
    ) {
        const size_t n = static_cast<size_t>(last - first);
        _begin = new T[n];
        _end = _begin + n;
        std::copy(first, last, _begin);
    }

    // Allocates memory by iterating along the range to count size (slower).
    template<typename InputIt>
    void _construct_from_range(
        InputIt first, InputIt last, std::false_type
    ) {
        size_t n = 0;
        for (InputIt it = first; it != last; ++it) {
            ++n;
        }
        _begin = new T[n];
        _end = _begin + n;
        std::copy(first, last, _begin);
    }

    #if SAFE_PTR_DEBUG_BOOL
        size_t _memory_id; // 0 is for if the ptr is a view
        static size_t _next_available_memory_id;
        static std::unordered_map<size_t,size_t> _ref_count;
        static std::unordered_map<size_t,bool> _is_deleted;
        static bool _id_overflow_occurred;
        static std::mutex _mtx;

        static size_t& _get_new_memory_id() {
            ++_next_available_memory_id;
            if (_next_available_memory_id == 0) {
                _id_overflow_occurred = true;
            }
            if (!_id_overflow_occurred) {
                return _next_available_memory_id;
            }
            while (true) { // handle overflow
                if (_ref_count.count(_next_available_memory_id) == 0) {
                    return _next_available_memory_id;
                }
                ++_next_available_memory_id;
            }
        }

        void _check_for_use_after_free() const noexcept(!SAFE_PTR_TEST_BOOL) {
            if (_get_is_deleted() == true) {
                SAFE_PTR_WARNING(
                    "Tried to access data after free() was called."
                );
            }
        }

        bool _get_is_view() const {
            return _memory_id == 0;
        }

        size_t& _get_ref_count() const {
            return _ref_count.at(_memory_id);
        }

        bool& _get_is_deleted() const {
            return _is_deleted.at(_memory_id);
        }

        void _warning(
            const char* const msg,
            const char* const file,
            int line,
            const char* const func
        ) const {
        #if SAFE_PTR_TEST_BOOL
            throw _SafePtrWarning();
        #endif
            std::cerr << "\033[33m"
                    << "SafePtr warning at " << file << ":" << line
                    << " in function " << func << ": "
                    << "\033[0m" << msg << "\n";
        }
    #endif
};

#if SAFE_PTR_DEBUG_BOOL
    template<typename T>
    size_t SafePtr<T>::_next_available_memory_id = 1;

    template<typename T>
    std::unordered_map<size_t,size_t> SafePtr<T>::_ref_count;

    template <typename T>
    std::unordered_map<size_t, bool> SafePtr<T>::_is_deleted = [] {
        std::unordered_map<size_t, bool> m;
        m[0] = false;
        return m;
    }();

    template<typename T>
    bool SafePtr<T>::_id_overflow_occurred = false;
    
    template<typename T>
    std::mutex SafePtr<T>::_mtx;
#endif

} // namespace fz
