#ifndef QDP_GENERIC_SPIN_PROJECT_EVALUTATES_H
#define QDP_GENERIC_SPIN_PROJECT_EVALUTATES_H

using namespace QDP;
namespace QDP
{

    // Typedefs
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 2> HVec;
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 4> FVec;

    // Four spinor (Ns * Nc * Ncomplex ) Ncomplex fastest
    typedef REAL SpinColFull[4][3][2];

    // Half spinor (Ns/2 * Nc * Ncomplex ) Ncomplex fastest
    typedef REAL SpinColHalf[2][3][2];
// d = SpinProjectDir0Plus(Vec);

////////////////////////////////
// Threading evaluates
//
// by Xu Guo, EPCC, 28 August, 2008
////////////////////////////////

// the wrappers for the functions to be threaded
#include "qdp_generic_spin_project_evaluates_wrapper.h"

    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir0Plus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {

        //  Get at pointer for 4 vec
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir0Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;

            //inlineSpinProjDir0Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir0Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir0Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir1Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir1Plus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir1Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir1Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir1Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir1Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir2Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir2Plus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir2Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir2Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir2Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir2Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir3Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir3Plus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir3Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir3Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir3Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir3Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir0Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir0Minus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {

        //  Get at pointer for 4 vec
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir0Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir0Minus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir0Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir0Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir1Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir1Minus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir1Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir1Minus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir1Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
      for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir1Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir2Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir2Minus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir2Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir2Minus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir2Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*  
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir2Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir3Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<HVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir3Minus,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<HVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {

            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinProjDir3Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir3Minus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinProjDir3Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*  
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir3Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir0Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {

            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir0Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir0Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir0Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir0Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir1Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir1Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir1Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir1Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir1Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir2Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir2Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir2Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir2Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir2Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir3Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir3Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir3Plus(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir3Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir3Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir0Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir0Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir0Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir0Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir0Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir1Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir1Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir1Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir1Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir1Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir2Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir2Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir2Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir2Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir2Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir3Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineSpinReconDir3Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir3Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineSpinReconDir3Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir3Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir0Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {

            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir0Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir0Plus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir0Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir0Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir1Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir1Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir1Plus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir1Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir1Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir2Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir2Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir2Plus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir2Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir2Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir3Plus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3Plus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir3Plus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir3Plus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir3Plus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir3Plus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir0Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir0Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir0Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir0Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir0Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir1Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir1Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir1Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir1Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir1Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir2Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir2Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir2Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir2Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir2Minus(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir3Minus(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3Minus,
                                       Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg arg = {aptr, bptr, inlineAddSpinReconDir3Minus};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir3Minus(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg<A, B> arg = {a, b, tab, inlineAddSpinReconDir3Minus};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir3Minus(aptr, bptr, 1);    
      }*/
        }
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
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir0PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        //  Get at pointer for 4 vec
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir0PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;

            //inlineSpinProjDir0PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir0PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir0PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir1PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir1PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir1PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir1PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir1PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir1PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir2PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir2PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir2PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir2PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir2PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir2PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir3PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir3PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir3PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir3PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir3PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir3PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir0MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir0MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        //  Get at pointer for 4 vec
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir0MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir0MinusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir0MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir0MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir1MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir1MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir1MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir1MinusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir1MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*
      for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir1MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir2MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir2MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir2MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir2MinusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir2MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*  
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir2MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinProjectDir3MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinProjectDir3MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {

            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinProjDir3MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinProjDir3MinusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinProjDir3MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*  
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinProjDir3MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir0PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {

            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir0PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir0PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir0PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir0PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir1PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir1PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir1PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir1PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir1PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir2PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir2PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir2PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir2PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir2PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir3PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir3PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir3PlusFull(aptr, bptr, n_vec);
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir3PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir3PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir0MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir0MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir0MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir0MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir0MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir1MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir1MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir1MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir1MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir1MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir2MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir2MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir2MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir2MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir2MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d = SpinReconstructDir3MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineSpinReconDir3MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineSpinReconDir3MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineSpinReconDir3MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineSpinReconDir3MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir0PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {

            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir0PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir0PlusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir0PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir0PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir1PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir1PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir1PlusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir1PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir1PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir2PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir2PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir2PlusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir2PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir2PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir3PlusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3PlusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir3PlusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir3PlusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir3PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir3PlusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir0MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir0MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir0MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir0MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir0MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir0MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir1MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir1MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir1MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir1MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir1MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir1MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir2MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir2MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir2MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir2MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir2MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir2MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

    // d += SpinReconstructDir3MinusFull(Vec);
    template <class A, class B>
    inline void evaluate(OLattice<FVec> &b,
                         const OpAddAssign &op,
                         const QDPExpr<
                             UnaryNode<FnSpinReconstructDir3MinusFull,
                                       Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {

        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().child());

        if (s.hasOrderedRep())
        {
            REAL *aptr = (REAL *)&(a.elem(s.start()).elem(0).elem(0).real());
            REAL *bptr = (REAL *)&(b.elem(s.start()).elem(0).elem(0).real());

            int total_n_vec = s.end() - s.start() + 1;

            ordered_spin_project_user_arg_full arg = {aptr, bptr, inlineAddSpinReconDir3MinusFull};

            dispatch_to_threads(total_n_vec, arg, ordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            //unsigned int n_vec=s.end() - s.start()+1;
            //inlineAddSpinReconDir3MinusFull(aptr, bptr, n_vec);
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_spin_project_user_arg_full<A, B> arg = {a, b, tab, inlineAddSpinReconDir3MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_spin_project_evaluate_function_full);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      
      REAL *aptr =(REAL *)&(a.elem(i).elem(0).elem(0).real());
      REAL *bptr =(REAL *)&(b.elem(i).elem(0).elem(0).real());
      inlineAddSpinReconDir3MinusFull(aptr, bptr, 1);    
      }*/
        }
    }

} // namespace QDP;

#endif
