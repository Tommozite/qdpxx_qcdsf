#ifndef GENERIC_FUSED_SPIN_PROJ_H
#define GENERIC_FUSED_SPIN_PROJ_H

namespace QDP
{

    // Convenience Types
    typedef PColorVector<RComplex<REAL>, 3> ColVec;
    typedef PSpinVector<ColVec, 4> Spin4;
    typedef PSpinVector<ColVec, 2> Spin2;
    typedef PColorMatrix<RComplex<REAL>, 3> ColMat;

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir0Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir0Plus(y)
    struct FnAdjMultSprojDir0Plus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir0Plus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Plus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir0Plus" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir0Plus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir0Plus(r)
    // to:    FnAdjMultSprojDir0Plus(l,r)
    //
    // The spinProjectDir0Plus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir0Plus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir0Plus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir0Plus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir0Plus, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir0Plus(l,r) <- adj(l)*FnSpinProjectDir0Plus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir0Plus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir0Plus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Plus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0Plus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Plus>::Type_t
        adjMultSprojDir0Plus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Plus>::Type_t tmp;

        tmp = spinProjectDir0Plus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0Plus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Plus>::Type_t
    adjMultSprojDir0Plus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Plus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Plus>::Type_t d;

        inlineSpinProjDir0Plus(&(b.elem(0).elem(0).real()),
                               &(d.elem(0).elem(0).real()),
                               1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir0Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

    // This is a struct for adj(x)*spinProjectDir0Minus(y)
    struct FnAdjMultSprojDir0Minus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir0Minus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Minus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Call the appropriate match
            return (adjMultSprojDir0Minus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir0Minus(r)
    // to:    FnAdjMultSprojDir0Minus(l,r)
    //
    // The spinProjectDir0Minus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir0Minus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir0Minus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir0Minus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir0Minus, T2>, C2> &r)
    {
        // Print Diagnostics
        //  cout << "FnAdjMultSprojDir0Minus(l,r) <- adj(l)*FnSpinProjectDir0Minus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir0Minus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir0Minus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Minus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0Minus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Minus>::Type_t
        adjMultSprojDir0Minus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Minus>::Type_t tmp;

        tmp = spinProjectDir0Minus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0Minus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Minus>::Type_t
    adjMultSprojDir0Minus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Minus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Minus>::Type_t d;

        inlineSpinProjDir0Minus(&(b.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir1Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

    // This is a struct for adj(x)*spinProjectDir1Plus(y)
    struct FnAdjMultSprojDir1Plus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir1Plus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Plus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Call the appropriate match
            return (adjMultSprojDir1Plus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir1Plus(r)
    // to:    FnAdjMultSprojDir1Plus(l,r)
    //
    // The spinProjectDir1Plus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir1Plus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir1Plus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir1Plus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir1Plus, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir1Plus(l,r) <- adj(l)*FnSpinProjectDir1Plus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir1Plus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir1Plus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Plus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1Plus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Plus>::Type_t
        adjMultSprojDir1Plus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Plus>::Type_t tmp;

        tmp = spinProjectDir1Plus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1Plus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Plus>::Type_t
    adjMultSprojDir1Plus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Plus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Plus>::Type_t d;

        inlineSpinProjDir1Plus(&(b.elem(0).elem(0).real()),
                               &(d.elem(0).elem(0).real()),
                               1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir1Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

    // This is a struct for adj(x)*spinProjectDir1Minus(y)
    struct FnAdjMultSprojDir1Minus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir1Minus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Minus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {

            // Call the appropriate match
            return (adjMultSprojDir1Minus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir1Minus(r)
    // to:    FnAdjMultSprojDir1Minus(l,r)
    //
    // The spinProjectDir1Minus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir1Minus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir1Minus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir1Minus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir1Minus, T2>, C2> &r)
    {
        // Print Diagnostics
        //  cout << "FnAdjMultSprojDir1Minus(l,r) <- adj(l)*FnSpinProjectDir1Minus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir1Minus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir1Minus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Minus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1Minus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Minus>::Type_t
        adjMultSprojDir1Minus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Minus>::Type_t tmp;

        tmp = spinProjectDir1Minus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1Minus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Minus>::Type_t
    adjMultSprojDir1Minus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Minus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Minus>::Type_t d;

        inlineSpinProjDir1Minus(&(b.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir2Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

    // This is a struct for adj(x)*spinProjectDir2Plus(y)
    struct FnAdjMultSprojDir2Plus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir2Plus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Plus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Call the appropriate match
            return (adjMultSprojDir2Plus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir2Plus(r)
    // to:    FnAdjMultSprojDir2Plus(l,r)
    //
    // The spinProjectDir2Plus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir2Plus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir2Plus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir2Plus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir2Plus, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir2Plus(l,r) <- adj(l)*FnSpinProjectDir2Plus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir2Plus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir2Plus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Plus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2Plus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Plus>::Type_t
        adjMultSprojDir2Plus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Plus>::Type_t tmp;

        tmp = spinProjectDir2Plus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2Plus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Plus>::Type_t
    adjMultSprojDir2Plus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Plus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Plus>::Type_t d;

        inlineSpinProjDir2Plus(&(b.elem(0).elem(0).real()),
                               &(d.elem(0).elem(0).real()),
                               1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir2Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir2Minus(y)
    struct FnAdjMultSprojDir2Minus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir2Minus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Minus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {

            // Call the appropriate match
            return (adjMultSprojDir2Minus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir2Minus(r)
    // to:    FnAdjMultSprojDir2Minus(l,r)
    //
    // The spinProjectDir2Minus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir2Minus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir2Minus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir2Minus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir2Minus, T2>, C2> &r)
    {
        // Print Diagnostics
        //  cout << "FnAdjMultSprojDir2Minus(l,r) <- adj(l)*FnSpinProjectDir2Minus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir2Minus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir2Minus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Minus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2Minus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Minus>::Type_t
        adjMultSprojDir2Minus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Minus>::Type_t tmp;

        tmp = spinProjectDir2Minus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2Minus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Minus>::Type_t
    adjMultSprojDir2Minus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Minus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Minus>::Type_t d;

        inlineSpinProjDir2Minus(&(b.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir3Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir3Plus(y)
    struct FnAdjMultSprojDir3Plus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir3Plus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Plus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Call the appropriate match
            return (adjMultSprojDir3Plus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir3Plus(r)
    // to:    FnAdjMultSprojDir3Plus(l,r)
    //
    // The spinProjectDir3Plus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir3Plus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir3Plus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir3Plus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir3Plus, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir3Plus(l,r) <- adj(l)*FnSpinProjectDir3Plus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir3Plus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir3Plus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Plus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3Plus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Plus>::Type_t
        adjMultSprojDir3Plus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Plus>::Type_t tmp;

        tmp = spinProjectDir3Plus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3Plus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Plus>::Type_t
    adjMultSprojDir3Plus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Plus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Plus>::Type_t d;

        inlineSpinProjDir3Plus(&(b.elem(0).elem(0).real()),
                               &(d.elem(0).elem(0).real()),
                               1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir3Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

    // This is a struct for adj(x)*spinProjectDir3Minus(y)
    struct FnAdjMultSprojDir3Minus
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir3Minus)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Minus>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {

            // Call the appropriate match
            return (adjMultSprojDir3Minus(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir3Minus(r)
    // to:    FnAdjMultSprojDir3Minus(l,r)
    //
    // The spinProjectDir3Minus gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir3Minus
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir3Minus,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir3Minus>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir3Minus, T2>, C2> &r)
    {
        // Print Diagnostics
        //  cout << "FnAdjMultSprojDir3Minus(l,r) <- adj(l)*FnSpinProjectDir3Minus(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir3Minus,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir3Minus>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin2
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Minus>
    {
        typedef Spin2 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3Minus struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Minus>::Type_t
        adjMultSprojDir3Minus(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Minus>::Type_t tmp;

        tmp = spinProjectDir3Minus(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3Minus struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Minus>::Type_t
    adjMultSprojDir3Minus(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Minus>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Minus>::Type_t d;

        inlineSpinProjDir3Minus(&(b.elem(0).elem(0).real()),
                                &(d.elem(0).elem(0).real()),
                                1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));

        return ret;
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
    //
    //
    //
    //
    //
    //
    //
    //
    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir0PlusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir0PlusFull(y)
    struct FnAdjMultSprojDir0PlusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir0PlusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir0PlusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir0PlusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir0PlusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir0PlusFull(r)
    // to:    FnAdjMultSprojDir0PlusFull(l,r)
    //
    // The spinProjectDir0PlusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir0PlusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir0PlusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir0PlusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir0PlusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir0PlusFull(l,r) <- adj(l)*FnSpinProjectDir0PlusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir0PlusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir0PlusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0PlusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0PlusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0PlusFull>::Type_t
        adjMultSprojDir0PlusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0PlusFull>::Type_t tmp;

        tmp = spinProjectDir0PlusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0PlusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0PlusFull>::Type_t
    adjMultSprojDir0PlusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0PlusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0PlusFull>::Type_t d;

        inlineSpinProjDir0PlusFull(&(b.elem(0).elem(0).real()),
                                   &(d.elem(0).elem(0).real()),
                                   1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }
    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir0MinusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir0MinusFull(y)
    struct FnAdjMultSprojDir0MinusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir0MinusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir0MinusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir0MinusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir0MinusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir0MinusFull(r)
    // to:    FnAdjMultSprojDir0MinusFull(l,r)
    //
    // The spinProjectDir0MinusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir0MinusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir0MinusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir0MinusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir0MinusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir0MinusFull(l,r) <- adj(l)*FnSpinProjectDir0MinusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir0MinusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir0MinusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0MinusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0MinusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0MinusFull>::Type_t
        adjMultSprojDir0MinusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir0MinusFull>::Type_t tmp;

        tmp = spinProjectDir0MinusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir0MinusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0MinusFull>::Type_t
    adjMultSprojDir0MinusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0MinusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0MinusFull>::Type_t d;

        inlineSpinProjDir0MinusFull(&(b.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir1PlusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir1PlusFull(y)
    struct FnAdjMultSprojDir1PlusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir1PlusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir1PlusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir1PlusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir1PlusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir1PlusFull(r)
    // to:    FnAdjMultSprojDir1PlusFull(l,r)
    //
    // The spinProjectDir1PlusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir1PlusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir1PlusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir1PlusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir1PlusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir1PlusFull(l,r) <- adj(l)*FnSpinProjectDir1PlusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir1PlusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir1PlusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1PlusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1PlusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1PlusFull>::Type_t
        adjMultSprojDir1PlusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1PlusFull>::Type_t tmp;

        tmp = spinProjectDir1PlusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1PlusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1PlusFull>::Type_t
    adjMultSprojDir1PlusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1PlusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1PlusFull>::Type_t d;

        inlineSpinProjDir1PlusFull(&(b.elem(0).elem(0).real()),
                                   &(d.elem(0).elem(0).real()),
                                   1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }
    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir1MinusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir1MinusFull(y)
    struct FnAdjMultSprojDir1MinusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir1MinusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir1MinusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir1MinusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir1MinusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir1MinusFull(r)
    // to:    FnAdjMultSprojDir1MinusFull(l,r)
    //
    // The spinProjectDir1MinusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir1MinusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir1MinusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir1MinusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir1MinusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir1MinusFull(l,r) <- adj(l)*FnSpinProjectDir1MinusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir1MinusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir1MinusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1MinusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1MinusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1MinusFull>::Type_t
        adjMultSprojDir1MinusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir1MinusFull>::Type_t tmp;

        tmp = spinProjectDir1MinusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir1MinusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1MinusFull>::Type_t
    adjMultSprojDir1MinusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1MinusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1MinusFull>::Type_t d;

        inlineSpinProjDir1MinusFull(&(b.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir2PlusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir2PlusFull(y)
    struct FnAdjMultSprojDir2PlusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir2PlusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir2PlusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir2PlusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir2PlusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir2PlusFull(r)
    // to:    FnAdjMultSprojDir2PlusFull(l,r)
    //
    // The spinProjectDir2PlusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir2PlusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir2PlusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir2PlusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir2PlusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir2PlusFull(l,r) <- adj(l)*FnSpinProjectDir2PlusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir2PlusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir2PlusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2PlusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2PlusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2PlusFull>::Type_t
        adjMultSprojDir2PlusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2PlusFull>::Type_t tmp;

        tmp = spinProjectDir2PlusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2PlusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2PlusFull>::Type_t
    adjMultSprojDir2PlusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2PlusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2PlusFull>::Type_t d;

        inlineSpinProjDir2PlusFull(&(b.elem(0).elem(0).real()),
                                   &(d.elem(0).elem(0).real()),
                                   1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }
    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir2MinusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir2MinusFull(y)
    struct FnAdjMultSprojDir2MinusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir2MinusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir2MinusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir2MinusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir2MinusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir2MinusFull(r)
    // to:    FnAdjMultSprojDir2MinusFull(l,r)
    //
    // The spinProjectDir2MinusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir2MinusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir2MinusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir2MinusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir2MinusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir2MinusFull(l,r) <- adj(l)*FnSpinProjectDir2MinusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir2MinusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir2MinusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2MinusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2MinusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2MinusFull>::Type_t
        adjMultSprojDir2MinusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir2MinusFull>::Type_t tmp;

        tmp = spinProjectDir2MinusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir2MinusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2MinusFull>::Type_t
    adjMultSprojDir2MinusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2MinusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2MinusFull>::Type_t d;

        inlineSpinProjDir2MinusFull(&(b.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }

    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir3PlusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir3PlusFull(y)
    struct FnAdjMultSprojDir3PlusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir3PlusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir3PlusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir3PlusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir3PlusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir3PlusFull(r)
    // to:    FnAdjMultSprojDir3PlusFull(l,r)
    //
    // The spinProjectDir3PlusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir3PlusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir3PlusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir3PlusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir3PlusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir3PlusFull(l,r) <- adj(l)*FnSpinProjectDir3PlusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir3PlusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir3PlusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3PlusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3PlusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3PlusFull>::Type_t
        adjMultSprojDir3PlusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3PlusFull>::Type_t tmp;

        tmp = spinProjectDir3PlusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3PlusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3PlusFull>::Type_t
    adjMultSprojDir3PlusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3PlusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3PlusFull>::Type_t d;

        inlineSpinProjDir3PlusFull(&(b.elem(0).elem(0).real()),
                                   &(d.elem(0).elem(0).real()),
                                   1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }
    /* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir3MinusFull(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
    // This is a struct for adj(x)*spinProjectDir3MinusFull(y)
    struct FnAdjMultSprojDir3MinusFull
    {
        // Boilerplate
        PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir3MinusFull)

        // OK This is an operator() so we can create instances of this
        // object and treat them as functions()
        template <class T1, class T2>
        inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir3MinusFull>::Type_t
        operator()(const T1 &a, const T2 &b) const
        {
            // Print Diagnostics
            //    cout << "FnAdjMultSprojDir3MinusFull" << endl << flush;

            // Call the appropriate match
            return (adjMultSprojDir3MinusFull(a, b));
        }
    };

    // This is an operator* that rewrites:
    //
    // from:  adj(l)*spinProjectDir3MinusFull(r)
    // to:    FnAdjMultSprojDir3MinusFull(l,r)
    //
    // The spinProjectDir3MinusFull gets grabbed and turned into an op identity
    //     adj(l)              gets grabbed and turned into an op identity
    //  the operation is encoded in the FnAdjMultSprojDir3MinusFull
    //

    template <class T1, class C1, class T2, class C2>
    inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir3MinusFull,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                                          typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>,
                               typename BinaryReturn<C1, C2, FnAdjMultSprojDir3MinusFull>::Type_t>::Expression_t
    operator*(const QDPExpr<UnaryNode<FnAdjoint, T1>, C1> &l,
              const QDPExpr<UnaryNode<FnSpinProjectDir3MinusFull, T2>, C2> &r)
    {
        // Print Diagnostics
        // cout << "FnAdjMultSprojDir3MinusFull(l,r) <- adj(l)*FnSpinProjectDir3MinusFull(r)" << endl;

        typedef UnaryNode<OpIdentity, T1> NewExpr1_t; // Shorthand for the new type for adj(l)
        typedef UnaryNode<OpIdentity, T2> NewExpr2_t; // Shorthand for the new type for spinProj

        // A name for the new node: Tree_t
        typedef BinaryNode<FnAdjMultSprojDir3MinusFull,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T1>, C1>>::Leaf_t,
                           typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity, T2>, C2>>::Leaf_t>
            Tree_t;

        // Create the result:  Tree_t = Binary node
        //                     Return type is the return type defined in structs below
        //                     The CreateLeaf-s do the rewriting -
        //                     unwrap the previous expression and rewrap it as the new type
        return MakeReturn<Tree_t, typename BinaryReturn<C1, C2, FnAdjMultSprojDir3MinusFull>::Type_t>::make(Tree_t(
            CreateLeaf<QDPExpr<NewExpr1_t, C1>>::make(NewExpr1_t(l.expression().child())),
            CreateLeaf<QDPExpr<NewExpr2_t, C2>>::make(NewExpr2_t(r.expression().child()))));
    }

    // Return types Fused su3*spinProj(Spin4)->Spin4
    template <>
    struct BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3MinusFull>
    {
        typedef Spin4 Type_t;
    };

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3MinusFull struct
    template <typename T1, typename T2>
    inline
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3MinusFull>::Type_t
        adjMultSprojDir3MinusFull(const T1 &a, const T2 &b)
    {
        typename BinaryReturn<T1, T2, FnAdjMultSprojDir3MinusFull>::Type_t tmp;

        tmp = spinProjectDir3MinusFull(b);
        return (adj(a) * tmp);
    }

    // This is what you need to specialise now. It'll get called by
    // the operator() of the FnAdjMultSprojDir3MinusFull struct
    template <>
    inline BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3MinusFull>::Type_t
    adjMultSprojDir3MinusFull(const PScalar<ColMat> &a, const Spin4 &b)
    {
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3MinusFull>::Type_t ret;
        BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3MinusFull>::Type_t d;

        inlineSpinProjDir3MinusFull(&(b.elem(0).elem(0).real()),
                                    &(d.elem(0).elem(0).real()),
                                    1);

        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(0), ret.elem(0));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(1), ret.elem(1));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(2), ret.elem(2));
        _inline_mult_adj_su3_mat_vec(a.elem(), d.elem(3), ret.elem(3));

        return ret;
    }
} // namespace QDP;

#endif
