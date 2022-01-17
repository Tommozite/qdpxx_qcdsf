#ifndef GENERIC_SPIN_RECON_H
#define GENERIC_SPIN_RECON_H

namespace QDP
{

    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 4> Spin4;
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 2> Spin2;

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir0Minus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir0Minus>::Type_t
    spinReconstructDir0Minus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir0Minus>::Type_t d;

        inlineSpinReconDir0Minus(&(s1.elem(0).elem(0).real()),
                                 &(d.elem(0).elem(0).real()),
                                 1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir0Plus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir0Plus>::Type_t
    spinReconstructDir0Plus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir0Plus>::Type_t d;

        inlineSpinReconDir0Plus(&(s1.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir1Minus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir1Minus>::Type_t
    spinReconstructDir1Minus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir1Minus>::Type_t d;

        inlineSpinReconDir1Minus(&(s1.elem(0).elem(0).real()),
                                 &(d.elem(0).elem(0).real()),
                                 1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir1Plus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir1Plus>::Type_t
    spinReconstructDir1Plus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir1Plus>::Type_t d;

        inlineSpinReconDir1Plus(&(s1.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir2Minus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir2Minus>::Type_t
    spinReconstructDir2Minus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir2Minus>::Type_t d;

        inlineSpinReconDir2Minus(&(s1.elem(0).elem(0).real()),
                                 &(d.elem(0).elem(0).real()),
                                 1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir2Plus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir2Plus>::Type_t
    spinReconstructDir2Plus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir2Plus>::Type_t d;

        inlineSpinReconDir2Plus(&(s1.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir3Minus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir3Minus>::Type_t
    spinReconstructDir3Minus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir3Minus>::Type_t d;

        inlineSpinReconDir3Minus(&(s1.elem(0).elem(0).real()),
                                 &(d.elem(0).elem(0).real()),
                                 1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin2, FnSpinReconstructDir3Plus>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin2, FnSpinReconstructDir3Plus>::Type_t
    spinReconstructDir3Plus(const Spin2 &s1)
    {
        UnaryReturn<Spin2, FnSpinReconstructDir3Plus>::Type_t d;

        inlineSpinReconDir3Plus(&(s1.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        return d;
    }
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir0MinusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir0MinusFull>::Type_t
    spinReconstructDir0MinusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir0MinusFull>::Type_t d;

        inlineSpinReconDir0MinusFull(&(s1.elem(0).elem(0).real()),
                                     &(d.elem(0).elem(0).real()),
                                     1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir0PlusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir0PlusFull>::Type_t
    spinReconstructDir0PlusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir0PlusFull>::Type_t d;

        inlineSpinReconDir0PlusFull(&(s1.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir1MinusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir1MinusFull>::Type_t
    spinReconstructDir1MinusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir1MinusFull>::Type_t d;

        inlineSpinReconDir1MinusFull(&(s1.elem(0).elem(0).real()),
                                     &(d.elem(0).elem(0).real()),
                                     1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir1PlusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir1PlusFull>::Type_t
    spinReconstructDir1PlusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir1PlusFull>::Type_t d;

        inlineSpinReconDir1PlusFull(&(s1.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir2MinusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir2MinusFull>::Type_t
    spinReconstructDir2MinusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir2MinusFull>::Type_t d;

        inlineSpinReconDir2MinusFull(&(s1.elem(0).elem(0).real()),
                                     &(d.elem(0).elem(0).real()),
                                     1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir2PlusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir2PlusFull>::Type_t
    spinReconstructDir2PlusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir2PlusFull>::Type_t d;

        inlineSpinReconDir2PlusFull(&(s1.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir3MinusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir3MinusFull>::Type_t
    spinReconstructDir3MinusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir3MinusFull>::Type_t d;

        inlineSpinReconDir3MinusFull(&(s1.elem(0).elem(0).real()),
                                     &(d.elem(0).elem(0).real()),
                                     1);

        return d;
    }

    template <>
    struct UnaryReturn<Spin4, FnSpinReconstructDir3PlusFull>
    {
        typedef Spin4 Type_t;
    };

    template <>
    inline UnaryReturn<Spin4, FnSpinReconstructDir3PlusFull>::Type_t
    spinReconstructDir3PlusFull(const Spin4 &s1)
    {
        UnaryReturn<Spin4, FnSpinReconstructDir3PlusFull>::Type_t d;

        inlineSpinReconDir3PlusFull(&(s1.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        return d;
    }

} // namespace QDP

#endif
