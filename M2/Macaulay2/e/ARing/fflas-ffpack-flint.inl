/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// ==========================================================================
// Copyright(c)'2021 by Mahrud Sayrafi
// This file is NOT a part of Givaro, but the majority of it is based
// on modifying givaro/src/kernel/ring/modular-balanced-int64.inl.
// This file is distributed under the GPLv2 license.
// ==========================================================================

#ifndef __fflas_ffpack_flint_INL
#define __fflas_ffpack_flint_INL

#include <cmath> // fmod

#define NORMALISE(x)				\
{						\
    if (x < _mhalfp) x += _p;	\
    else if (x > _halfp) x -= _p;	\
}

#define NORMALISE_HI(x)				\
{						\
    if (x > _halfp) x -= _p;		\
}

// TODO: remove this clause
namespace Givaro {
    int64_t Integer::operator % (unsigned long int l) const
    {
        int64_t res ;
        if (l>0) {
            res = static_cast<int64_t>(this->operator%( static_cast<uint64_t>(l) ) );
        }
        else {
            res = static_cast<int64_t>(this->operator%( static_cast<uint64_t>( -l ) ) );
        }
        return res;
    }
} // namespace Givaro

namespace FFLAS {
    namespace Protected {
        template<>
        inline int WinogradThreshold (const M2::ModularFlint<fmpz> & F) {return __FFLASFFPACK_WINOTHRESHOLD_BAL;}
    } // namespace Protected

    template <typename Element>
    struct ModeTraits<M2::ModularFlint<Element>> {typedef typename ModeCategories::DelayedTag value;};

    template <>
    inline size_t bitsize(const M2::ModularFlint<fmpz>& F, size_t M, size_t N, const typename M2::ModularFlint<fmpz>::ConstElement_ptr A, size_t lda){
        fmpz max = FLINT_MAX(F.maxElement(), FLINT_ABS(F.minElement()));
	return fmpz_bits(&max);
    }

    template <> struct ModeTraits<M2::ModularFlint<fmpz>> {typedef typename ModeCategories::ConvertTo<ElementCategories::MachineFloatTag> value;};
} // namespace FFLAS

namespace ARing
{

    //----- Classic arithmetic

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::mul(Element& r, const Element& a, const Element& b) const
    {
        Element q = static_cast<Element>(double(a) * double(b) * _dinvp);
        r = static_cast<Element>(a * b - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::div(Element& r, const Element& a, const Element& b) const
    {
        Element tmp;
        return mul (r, a, inv(tmp, b));
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::add(Element& r, const Element& a, const Element& b) const
    {
        r = a + b;
        NORMALISE(r);
        return r;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::sub(Element& r, const Element& a, const Element& b) const
    {
        r = a - b;
        NORMALISE(r);
        return r;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::neg(Element& r, const Element& a) const
    {
        return r = -a;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::inv(Element& r, const Element& a) const
    {
	fmpz_invmod(&r, &a, &_p);
        NORMALISE(r);
        return r;
    }

    inline bool ModularFlint<fmpz>::isUnit(const Element& a) const
    {
        fmpz_t d;
	fmpz_init(d);
        fmpz_gcd(d, &a, &_p);
        return fmpz_is_pm1(d);
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::mulin(Element& r, const Element& a) const
    {
        return mul(r, r, a);
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::divin(Element& r, const Element& a) const
    {
        return div(r, r, a);
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::addin(Element& r, const Element& a) const
    {
        return add(r, r, a);
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::subin(Element& r, const Element& a) const
    {
        return sub(r, r, a);
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::negin(Element& r) const
    {
        return neg(r, r);
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::invin(Element& r) const
    {
        return inv(r, r);
    }

    //----- Special arithmetic

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::axpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) + (double) y) * _dinvp);
        r = static_cast<Element>(a * x + y - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::axpyin(Element& r, const Element& a, const Element& x) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) + (double) r) * _dinvp);
        r = static_cast<Element>(a * x + r - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::axmy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) - (double) y) * _dinvp);
        r = static_cast<Element>(a * x - y - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::axmyin(Element& r, const Element& a, const Element& x) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) - (double) r) * _dinvp);
        r = static_cast<Element>(a * x - r - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::maxpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return negin(axmy(r, a, x, y));
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>:: maxpyin(Element& r, const Element& a, const Element& x) const
    {
        return negin(axmyin(r, a, x));
    }

    //----- Init

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::init(Element& x) const
    {
        return x = 0;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::init(Element& x, const float y) const
    {
        x = static_cast<Element>(fmod(y, double(_p)));
        NORMALISE(x);
        return x;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::init(Element& x, const double y) const
    {
        x = static_cast<Element>(fmod(y, double(_p)));
        NORMALISE(x);
        return x;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::init(Element& x, const int64_t y) const
    {
        x = static_cast<Element>(y % _p);
        NORMALISE_HI(x);
        return x;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::init(Element& x, const uint64_t y) const
    {
        x = static_cast<Element>(y % _p);
        NORMALISE_HI(x);
        return x;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::assign(Element& x, const Element& y) const
    {
        return x = y;
    }

    //----- Reduce

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::reduce(Element& x, const Element& y) const
    {
        x = y % _p;
        NORMALISE(x);
        return x;
    }

    inline ModularFlint<fmpz>::Element&
    ModularFlint<fmpz>::reduce(Element& x) const
    {
        x %= _p;
        NORMALISE(x);
        return x;
    }

    //----- IO

    inline std::ostream&
    ModularFlint<fmpz>::write(std::ostream& os) const
    {
        return os << "ModularFlint<fmpz> modulo " << _p;
    }

    inline std::ostream&
    ModularFlint<fmpz>::write(std::ostream& os, const Element& x) const
    {
        return os << x;
    }

    inline std::istream&
    ModularFlint<fmpz>::read(std::istream& is, Element& x) const
    {
        Element tmp;
        is >> tmp;
        init(x, tmp);
        return is;
    }

} // namespace ARing

#undef NORMALISE
#undef NORMALISE_HI

#endif // __fflas_ffpack_flint_INL
