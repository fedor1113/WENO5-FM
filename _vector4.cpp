/*
-----------------------------------------------------------------------------
Adapted from Ogre3D
-----------------------------------------------------------------------------
*/

#include <array>

#include "_vector4.h"

template class Vector4<double>;

template<>
const Vector4<int> Vector4<int>::ZERO( 0, 0, 0, 0 );

template<>
const Vector4<double> Vector4<double>::ZERO( 0., 0., 0., 0. );

template<>
const Vector4<long double> Vector4<long double>::ZERO(std::valarray<long double>(0., 4));
