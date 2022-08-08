/*
-----------------------------------------------------------------------------
Adapted from Ogre3D
-----------------------------------------------------------------------------
*/

#ifndef __Vector4_H__
#define __Vector4_H__

#include <assert.h>
// #include <concepts>
// #include <cstddef>
#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

/** 4-dimensional homogeneous vector.
*/



//template<class N>
//concept SignedNumber = std::is_signed_v<N>
//	&& std::is_arithmetic<N>
//	&& std::is_floating_point<N>;


template <typename N>
//requires SignedNumber<N>
class Vector4
{
public:

	N x, y, z, w;

	static const Vector4<N> ZERO;

	Vector4()
		: x( ZERO.x ), y( ZERO.y ), z( ZERO.z ), w( ZERO.w )
	{
	}

	inline Vector4( const N fX, const N fY, const N fZ, const N fW )
        : x( fX ), y( fY ), z( fZ ), w( fW)
    {
    }

	inline explicit Vector4( const N scalar )
        : x( scalar )
        , y( scalar )
        , z( scalar )
        , w( scalar )
    {
    }

	inline Vector4(std::valarray<N> V) : x(V[0]), y(V[1]), z(V[2]), w(V[3]) {}
	inline Vector4(std::vector<N> V) : x(V[0]), y(V[1]), z(V[2]), w(V[3]) {}

	inline N square() const
	{
		return x*x + y*y + z*z + w*w;
	}

	inline N operator [] ( const std::size_t i ) const
    {
        assert( i < 4 );

        return *(&x+i);
    }

	inline N& operator [] ( const std::size_t i )
    {
        assert( i < 4 );

        return *(&x+i);
    }

    inline Vector4& operator = ( const Vector4& rkVector )
    {
        x = rkVector.x;
        y = rkVector.y;
        z = rkVector.z;
        w = rkVector.w;

        return *this;
    }

	inline Vector4& operator = ( const N fScalar)
	{
		x = fScalar;
		y = fScalar;
		z = fScalar;
		w = fScalar;
		return *this;
	}

    inline bool operator == ( const Vector4& rkVector ) const
    {
        return ( x == rkVector.x &&
            y == rkVector.y &&
            z == rkVector.z &&
            w == rkVector.w );
    }

    inline bool operator != ( const Vector4& rkVector ) const
    {
        return ( x != rkVector.x ||
            y != rkVector.y ||
            z != rkVector.z ||
            w != rkVector.w );
    }

	inline bool operator < ( const Vector4& rkVector ) const
	{
		return (square() < rkVector.square());
	}

    // arithmetic operations
    inline Vector4 operator + ( const Vector4& rkVector ) const
    {
        return Vector4(
            x + rkVector.x,
            y + rkVector.y,
            z + rkVector.z,
            w + rkVector.w);
    }

    inline Vector4 operator - ( const Vector4& rkVector ) const
    {
        return Vector4(
            x - rkVector.x,
            y - rkVector.y,
            z - rkVector.z,
            w - rkVector.w);
    }

	inline Vector4 operator * ( const N fScalar ) const
    {
        return Vector4(
            x * fScalar,
            y * fScalar,
            z * fScalar,
            w * fScalar);
    }

    inline Vector4 operator * ( const Vector4& rhs) const
    {
        return Vector4(
            rhs.x * x,
            rhs.y * y,
            rhs.z * z,
            rhs.w * w);
    }

	inline Vector4 operator / ( const N fScalar ) const
    {
        assert( fScalar != 0.0 );

		N fInv = 1.0 / fScalar;

        return Vector4(
            x * fInv,
            y * fInv,
            z * fInv,
            w * fInv);
    }

    inline Vector4 operator / ( const Vector4& rhs) const
    {
        return Vector4(
            x / rhs.x,
            y / rhs.y,
            z / rhs.z,
            w / rhs.w);
    }

    inline const Vector4& operator + () const
    {
        return *this;
    }

    inline Vector4 operator - () const
    {
        return Vector4(-x, -y, -z, -w);
    }

	inline friend Vector4 operator * ( const N fScalar, const Vector4& rkVector )
    {
        return Vector4(
            fScalar * rkVector.x,
            fScalar * rkVector.y,
            fScalar * rkVector.z,
            fScalar * rkVector.w);
    }

	inline friend Vector4 operator / ( const N fScalar, const Vector4& rkVector )
    {
        return Vector4(
            fScalar / rkVector.x,
            fScalar / rkVector.y,
            fScalar / rkVector.z,
            fScalar / rkVector.w);
    }

	inline friend Vector4 operator + (const Vector4& lhs, const N rhs)
    {
        return Vector4(
            lhs.x + rhs,
            lhs.y + rhs,
            lhs.z + rhs,
            lhs.w + rhs);
    }

	inline friend Vector4 operator + (const N lhs, const Vector4& rhs)
    {
        return Vector4(
            lhs + rhs.x,
            lhs + rhs.y,
            lhs + rhs.z,
            lhs + rhs.w);
    }

	inline friend Vector4 operator - (const Vector4& lhs, N rhs)
    {
        return Vector4(
            lhs.x - rhs,
            lhs.y - rhs,
            lhs.z - rhs,
            lhs.w - rhs);
    }

	inline friend Vector4 operator - (const N lhs, const Vector4& rhs)
    {
        return Vector4(
            lhs - rhs.x,
            lhs - rhs.y,
            lhs - rhs.z,
            lhs - rhs.w);
    }

    // arithmetic updates
    inline Vector4& operator += ( const Vector4& rkVector )
    {
        x += rkVector.x;
        y += rkVector.y;
        z += rkVector.z;
        w += rkVector.w;

        return *this;
    }

    inline Vector4& operator -= ( const Vector4& rkVector )
    {
        x -= rkVector.x;
        y -= rkVector.y;
        z -= rkVector.z;
        w -= rkVector.w;

        return *this;
    }

	inline Vector4& operator *= ( const N fScalar )
    {
        x *= fScalar;
        y *= fScalar;
        z *= fScalar;
        w *= fScalar;
        return *this;
    }

	inline Vector4& operator += ( const N fScalar )
    {
        x += fScalar;
        y += fScalar;
        z += fScalar;
        w += fScalar;
        return *this;
    }

	inline Vector4& operator -= ( const N fScalar )
    {
        x -= fScalar;
        y -= fScalar;
        z -= fScalar;
        w -= fScalar;
        return *this;
    }

    inline Vector4& operator *= ( const Vector4& rkVector )
    {
        x *= rkVector.x;
        y *= rkVector.y;
        z *= rkVector.z;
        w *= rkVector.w;

        return *this;
    }

	inline Vector4& operator /= ( const N fScalar )
    {
        assert( fScalar != 0.0 );

		N fInv = 1.0 / fScalar;

        x *= fInv;
        y *= fInv;
        z *= fInv;
        w *= fInv;

        return *this;
    }

    inline Vector4& operator /= ( const Vector4& rkVector )
    {
        x /= rkVector.x;
        y /= rkVector.y;
        z /= rkVector.z;
        w /= rkVector.w;

        return *this;
    }

    /** Calculates the dot (scalar) product of this vector with another.
    */
	inline N dotProduct(const Vector4& vec) const
    {
        return x * vec.x + y * vec.y + z * vec.z + w * vec.w;
	}

	/** Vectorizing common std functions.
	*/
	inline friend Vector4 abs(const Vector4& rkVector)
	{
		return Vector4(
			std::abs(rkVector.x),
			std::abs(rkVector.y),
			std::abs(rkVector.z),
			std::abs(rkVector.w)
		);
	}


	inline friend Vector4 exp(const Vector4& rkVector)
	{
		return Vector4(
			std::exp(rkVector.x),
			std::exp(rkVector.y),
			std::exp(rkVector.z),
			std::exp(rkVector.w)
		);
	}

	inline friend Vector4 log(const Vector4& rkVector)
	{
		return Vector4(
			std::log(rkVector.x),
			std::log(rkVector.y),
			std::log(rkVector.z),
			std::log(rkVector.w)
		);
	}

	inline friend Vector4 log10(const Vector4& rkVector)
	{
		return Vector4(
			std::log10(rkVector.x),
			std::log10(rkVector.y),
			std::log10(rkVector.z),
			std::log10(rkVector.w)
		);
	}

	inline friend Vector4 pow(const Vector4& base, const Vector4& exp)
	{
		return Vector4(
			std::pow(base.x, exp.x),
			std::pow(base.y, exp.y),
			std::pow(base.z, exp.z),
			std::pow(base.w, exp.w)
		);
	}

	inline friend Vector4 pow(const Vector4& base, const N exp)
	{
		return Vector4(
			std::pow(base.x, exp),
			std::pow(base.y, exp),
			std::pow(base.z, exp),
			std::pow(base.w, exp)
		);
	}

	inline friend Vector4 pow(const N base, const Vector4& exp)
	{
		return Vector4(
			std::pow(base, exp.x),
			std::pow(base, exp.y),
			std::pow(base, exp.z),
			std::pow(base, exp.w)
		);
	}

	inline friend Vector4 sqrt(const Vector4& rkVector)
	{
		return Vector4(
			std::sqrt(rkVector.x),
			std::sqrt(rkVector.y),
			std::sqrt(rkVector.z),
			std::sqrt(rkVector.w)
		);
	}

	inline friend Vector4 sin(const Vector4& rkVector)
	{
		return Vector4(
			std::sin(rkVector.x),
			std::sin(rkVector.y),
			std::sin(rkVector.z),
			std::sin(rkVector.w)
		);
	}

	inline friend Vector4 cos(const Vector4& rkVector)
	{
		return Vector4(
			std::cos(rkVector.x),
			std::cos(rkVector.y),
			std::cos(rkVector.z),
			std::cos(rkVector.w)
		);
	}

	inline friend Vector4 tan(const Vector4& rkVector)
	{
		return Vector4(
			std::tan(rkVector.x),
			std::tan(rkVector.y),
			std::tan(rkVector.z),
			std::tan(rkVector.w)
		);
	}

	inline friend Vector4 asin(const Vector4& rkVector)
	{
		return Vector4(
			std::asin(rkVector.x),
			std::asin(rkVector.y),
			std::asin(rkVector.z),
			std::asin(rkVector.w)
		);
	}

	inline friend Vector4 acos(const Vector4& rkVector)
	{
		return Vector4(
			std::acos(rkVector.x),
			std::acos(rkVector.y),
			std::acos(rkVector.z),
			std::acos(rkVector.w)
		);
	}

	inline friend Vector4 atan(const Vector4& rkVector)
	{
		return Vector4(
			std::atan(rkVector.x),
			std::atan(rkVector.y),
			std::atan(rkVector.z),
			std::atan(rkVector.w)
		);
	}

//	inline friend Vector4 atan2(const Vector4& rkVector)
//	{
//		return Vector4(
//			std::atan2(rkVector.x),
//			std::atan2(rkVector.y),
//			std::atan2(rkVector.z),
//			std::atan2(rkVector.w)
//		);
//	}

	inline friend Vector4 sinh(const Vector4& rkVector)
	{
		return Vector4(
			std::sinh(rkVector.x),
			std::sinh(rkVector.y),
			std::sinh(rkVector.z),
			std::sinh(rkVector.w)
		);
	}

	inline friend Vector4 cosh(const Vector4& rkVector)
	{
		return Vector4(
			std::cosh(rkVector.x),
			std::cosh(rkVector.y),
			std::cosh(rkVector.z),
			std::cosh(rkVector.w)
		);
	}

	inline friend Vector4 tanh(const Vector4& rkVector)
	{
		return Vector4(
			std::tanh(rkVector.x),
			std::tanh(rkVector.y),
			std::tanh(rkVector.z),
			std::tanh(rkVector.w)
		);
	}

    /** Function for writing to a stream.
    */
	inline friend std::ostream& operator <<
		( std::ostream& o, const Vector4& v )
    {
        o << "Vector4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
        return o;
    }
};



#endif
