/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// ==========================================================================
// Copyright(c)'2021 by Mahrud Sayrafi
// This file is NOT a part of Givaro, but the majority of it is based
// on modifying givaro/src/kernel/ring/modular-balanced-int64.h.
// This file is distributed under the GPLv2 license.
// ==========================================================================

#ifndef __fflas_ffpack_flint_H
#define __fflas_ffpack_flint_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>

#include <fflas-ffpack/fflas-ffpack.h>

#include <iostream>

namespace ARing
{
    // compare with Givaro::ModularBalanced<int64_t>
    template<class TAG> class ModularFlint;

    template <>
    class ModularFlint<fmpz>
    {
    public:

        // ----- Exported types
        using Self_t = ModularFlint<fmpz>;
        using Element = fmpz;
        using Element_ptr = Element*;
        using ConstElement = const Element;
        using ConstElement_ptr = const Element*;
        using Residu_t = long; // FIXME
        enum { size_rep = sizeof(Element) };

        // ----- Constantes
        const Element zero = 0;
        const Element one  = 1;
        const Element mOne = -1;

        // ----- Flint
        // TODO: use fq_nmod arithmetic more
        fq_nmod_ctx_t _ctx;

        // ----- Constructors
        ModularFlint()
        : _p(zero), _halfp(zero), _mhalfp(zero), _dinvp(0.0)
        {}

        ModularFlint(Element p)
	: _p(p), _dinvp(1. / fmpz_get_d(&p))
        {
	    assert(fmpz_get_ui(&_p) >= minCardinality());
	    assert(fmpz_get_ui(&_p) <= maxCardinality());
	    fmpz_tdiv_q_2exp(&_halfp, &_p, 1);
	    fmpz_sub(&_mhalfp, &_halfp, &_p);
	    fmpz_sub_ui(&_mhalfp, &_mhalfp, 1);
	    fq_nmod_ctx_init(_ctx, &_p, 1, "a");
        }

        ModularFlint(const Self_t& F)
        : _p(F._p), _halfp(F._halfp), _mhalfp(F._mhalfp), _dinvp(F._dinvp)
        {}

        // ----- Accessors
        inline Element minElement() const { return _mhalfp; }
        inline Element maxElement() const { return _halfp; }

        // ----- Access to the modulus
        inline Residu_t residu() const { return _p; }
        inline Residu_t size() const { return _p; }
        inline Residu_t characteristic() const { return _p; }
        inline Residu_t cardinality() const { return _p; }
        template<class T> inline T& characteristic(T& p) const { return p = _p; }
        template<class T> inline T& cardinality(T& p) const { return p = _p; }

        // TODO: maxCardinality of Flint is much higher
        static inline Residu_t maxCardinality() { return 6074000999ull; } // p=floor(2^32.5) s.t. a*b+c fits in int64_t, with abs(a,b,c) <= (p-1)/2
        static inline Residu_t minCardinality() { return 3; }

        // ----- Checkers
        inline bool isZero(const Element& a) const { return a == zero; }
        inline bool isOne (const Element& a) const { return a == one; }
        inline bool isMOne(const Element& a) const { return a == mOne; }
        inline bool isUnit(const Element& a) const;
        inline bool areEqual(const Element& a, const Element& b) const { return a == b; }
        inline size_t length(const Element a) const { return size_rep; }

        // ----- Ring-wise operators
        inline bool operator==(const Self_t& F) const { return _p == F._p; }
        inline bool operator!=(const Self_t& F) const { return _p != F._p; }
        inline Self_t& operator=(const Self_t& F)
        {
            F.assign(const_cast<Element&>(one),  F.one);
            F.assign(const_cast<Element&>(zero), F.zero);
            F.assign(const_cast<Element&>(mOne), F.mOne);
            _p = F._p;
            _halfp = F._halfp;
            _mhalfp = F._mhalfp;
            _dinvp = F._dinvp;
            return *this;
        }

        // ----- Initialisation
        Element& init(Element& a) const;
        Element& init(Element& r, const float a) const;
        Element& init(Element& r, const double a) const;
        Element& init(Element& r, const int64_t a) const;
        Element& init(Element& r, const uint64_t a) const;
        Element& init(Element& r, const fmpz_t a) const;
        // template<typename T> Element& init(Element& r, const T& a) const
        // { r = Caster<Element>(a); return reduce(r); }

        Element& assign(Element& r, const Element& a) const;

        // ----- Convert
        template<typename T> T& convert(T& r, const Element& a) const
        { return r = static_cast<T>(a); }

        Element& reduce(Element& r, const Element& a) const;
        Element& reduce(Element& r) const;

        // ----- Classic arithmetic
        Element& mul(Element& r, const Element& a, const Element& b) const;
        Element& div(Element& r, const Element& a, const Element& b) const;
        Element& add(Element& r, const Element& a, const Element& b) const;
        Element& sub(Element& r, const Element& a, const Element& b) const;
        Element& neg(Element& r, const Element& a) const;
        Element& inv(Element& r, const Element& a) const;

        Element& mulin(Element& r, const Element& a) const;
        Element& divin(Element& r, const Element& a) const;
        Element& addin(Element& r, const Element& a) const;
        Element& subin(Element& r, const Element& a) const;
        Element& negin(Element& r) const;
        Element& invin(Element& r) const;

        // -- axpy:   r <- a * x + y
        // -- axpyin: r <- a * x + r
        Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const;
        Element& axpyin(Element& r, const Element& a, const Element& x) const;

        // -- axmy:   r <- a * x - y
        // -- axmyin: r <- a * x - r
        Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const;
        Element& axmyin(Element& r, const Element& a, const Element& x) const;

        // -- maxpy:   r <- y - a * x
        // -- maxpyin: r <- r - a * x
        Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const;
        Element& maxpyin(Element& r, const Element& a, const Element& x) const;

        // TODO
        // // -- type_string
        // static const std::string type_string () {
        //     return "ModularFlint<" + TypeString<Element>::get() +  ">";
        // }

        // TODO
        // // ----- Random generators
        // typedef ModularRandIter<Self_t> RandIter;
        // typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        // template< class Random > Element& random(Random& g, Element& r) const
        // { return init(r, g()); }
        // template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        // { while (isZero(init(a, g())))
        //     ;
        //     return a; }

        // --- IO methods
        std::ostream& write(std::ostream& s) const;
        std::istream& read (std::istream& s, Element& a) const;
        std::ostream& write(std::ostream& s, const Element& a) const;

    protected:

        Residu_t _p;
        Element _halfp;
        Element _mhalfp;
        double _dinvp;
    };
} // namespace ARing

#include "fflas-ffpack-flint.inl"

#endif // __fflas_ffpack_flint_H
