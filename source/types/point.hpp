// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef TYPES_POINT_HPP
#define TYPES_POINT_HPP

#pragma pack(push,1)
template <typename type, unsigned int dimensions>
class point final {
public:
    type data[dimensions];

public:
    point() = default;

    const type& operator[](unsigned int index) const {
        return this->data[index];
    }

    type& operator[](unsigned int index) {
        return this->data[index];
    }
};
#pragma pack(pop)

#pragma pack(push,1)
template <typename type>
class point<type, 3> final {
public:
    #if defined(_MSC_VER)
        #pragma warning(push)
        #pragma warning(disable : 4201)
    #endif
    #if defined(__GNUC__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wgnu-anonymous-struct"
        #pragma GCC diagnostic ignored "-Wnested-anon-types"
    #endif
    union {
        type data[3];
        struct {
            type x, y, z;
        };
    };
    #if defined(__GNUC__)
        #pragma GCC diagnostic pop
    #endif
    #if defined(_MSC_VER)
        #pragma warning(pop)
    #endif
public:
    point() = default;

    const type& operator[](unsigned int index) const {
        return this->data[index];
    }

    type& operator[](unsigned int index) {
        return this->data[index];
    }
};
#pragma pack(pop)

#pragma pack(push,1)
template <typename type>
class point<type, 2> final {
public:
    #if defined(_MSC_VER)
        #pragma warning(push)
        #pragma warning(disable : 4201)
    #endif
    #if defined(__GNUC__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wgnu-anonymous-struct"
        #pragma GCC diagnostic ignored "-Wnested-anon-types"
    #endif
    union {
        type data[2];
        struct {
            type x, y;
        };
    };
    #if defined(__GNUC__)
        #pragma GCC diagnostic pop
    #endif
    #if defined(_MSC_VER)
        #pragma warning(pop)
    #endif
public:
    point() = default;

    const type& operator[](unsigned int index) const {
        return this->data[index];
    }

    type& operator[](unsigned int index) {
        return this->data[index];
    }
};
#pragma pack(pop)

template <typename type>
using point2 = point<type, 2>;

using point2i = point<int, 2>;
using point2f = point<float, 2>;
using point2d = point<double, 2>;

template <typename type>
using point3 = point<type, 3>;

using point3i = point<int, 3>;
using point3f = point<float, 3>;
using point3d = point<double, 3>;



#endif // TYPES_POINT_HPP
