#ifndef QDP_GENERIC_FUSED_SPIN_RECON_EVALUATES_H
#define QDP_GENERIC_FUSED_SPIN_RECON_EVALUATES_H

namespace QDP
{

////////////////////////////////
// Threading evaluates
//
// by Xu Guo, EPCC, 28 August, 2008
////////////////////////////////

// the wrappers for the functions to be threaded
#include "qdp_generic_fused_spin_recon_evaluates_wrapper.h"

    // Vec = SpinReconstructDir0Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir0Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
     }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir0Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);      
    }*/
        }
    }

    // Vec = SpinReconstructDir0Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir0Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
 
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
   
				}*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir0Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

      
				}*/
        }
    }

    // Vec = SpinReconstructDir1Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir1Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
 
			       }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir1Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
      
			       }*/
        }
    }

    // Vec = SpinReconstructDir1Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir1Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
      
				}*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir1Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
      
      
				}*/
        }
    }

    // Vec = SpinReconstructDir2Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir2Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
			       }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir2Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
           
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
      
			       }*/
        }
    }

    // Vec = SpinReconstructDir2Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir2Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
				}*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir2Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

 
				}*/
        }
    }

    // Vec = SpinReconstructDir3Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir3Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir3Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);

			       }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir3Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir3Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
      
			       }*/
        }
    }

    // Vec = SpinReconstructDir3Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir3Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {

      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      inlineSpinReconDir3Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

				}*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir3Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      inlineSpinReconDir3Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

      
				}*/
        }
    }

    // Vec += SpinReconstructDir0Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir0Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
   
   
				  }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir0Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
      
 
				  }*/
        }
    }

    // Vec += SpinReconstructDir0Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir0Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
      
   
   
				   }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir0Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

      
				   }*/
        }
    }

    // Vec += SpinReconstructDir1Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir1Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
  
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      

      inlineAddSpinReconDir1Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
   
				  }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir1Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      

      inlineAddSpinReconDir1Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
      
      
				  }*/
        }
    }

    // Vec += SpinReconstructDir1Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir1Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir1Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

  
				   }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir1Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir1Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

 
				   }*/
        }
    }

    // Vec += SpinReconstructDir2Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir2Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
				  }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir2Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
      
      
				  }*/
        }
    }

    // Vec += SpinReconstructDir2Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir2Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
				   }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir2Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
      
      
				   }*/
        }
    }

    // Vec += SpinReconstructDir3Plus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3PlusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir3Plus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
   
				  } */
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir3Plus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3Plus( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
				  }*/
        }
    }

    // Vec += SpinReconstructDir3Minus( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3MinusProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<HVec, OLattice<HVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<HVec> &a = static_cast<const OLattice<HVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir3Minus};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
    HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
    
				   } */
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir3Minus};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      HVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3Minus( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

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
    //
    //
    //
    // Vec = SpinReconstructDir0PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir0PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
     }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir0PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);      
    }*/
        }
    }

    // Vec = SpinReconstructDir0MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir0MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
 
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
   
				}*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir0MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir0MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

      
				}*/
        }
    }

    // Vec = SpinReconstructDir1PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir1PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
 
			       }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir1PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
      
			       }*/
        }
    }

    // Vec = SpinReconstructDir1MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir1MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
      
				}*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir1MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir1MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
      
      
				}*/
        }
    }

    // Vec = SpinReconstructDir2PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir2PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
			       }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir2PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
           
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
      
			       }*/
        }
    }

    // Vec = SpinReconstructDir2MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir2MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);
				}*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir2MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir2MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

 
				}*/
        }
    }

    // Vec = SpinReconstructDir3PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir3PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir3PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);

			       }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir3PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineSpinReconDir3PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
			       (REAL *)&(d.elem(site).elem(0).elem(0).real()),
			       1);
      
      
			       }*/
        }
    }

    // Vec = SpinReconstructDir3MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineSpinReconDir3MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {

      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      inlineSpinReconDir3MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

				}*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineSpinReconDir3MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      inlineSpinReconDir3MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				(REAL *)&(d.elem(site).elem(0).elem(0).real()),
				1);

      
				}*/
        }
    }

    // Vec += SpinReconstructDir0PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir0PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
   
   
				  }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir0PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
      
 
				  }*/
        }
    }

    // Vec += SpinReconstructDir0MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir0MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir0MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
      
   
   
				   }*/
        }
        else
        {

            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir0MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir0MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

      
				   }*/
        }
    }

    // Vec += SpinReconstructDir1PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir1PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /* 
    for(int site = s.start(); site <= s.end(); site++) {
  
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      

      inlineAddSpinReconDir1PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
   
				  }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir1PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      

      inlineAddSpinReconDir1PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
      
      
				  }*/
        }
    }

    // Vec += SpinReconstructDir1MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir1MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir1MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir1MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

  
				   }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir1MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir1MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

 
				   }*/
        }
    }

    // Vec += SpinReconstructDir2PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir2PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
				  }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir2PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
      
      
				  }*/
        }
    }

    // Vec += SpinReconstructDir2MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir2MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir2MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
				   }*/
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir2MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir2MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
      
      
				   }*/
        }
    }

    // Vec += SpinReconstructDir3PlusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3PlusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {

            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir3PlusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
   
				  } */
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir3PlusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3PlusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				  (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				  1);
				  }*/
        }
    }

    // Vec += SpinReconstructDir3MinusFull( u * psi);
    template <>
    inline void evaluate(OLattice<FVec> &d,
                         const OpAddAssign &op,
                         const QDPExpr<
                             BinaryNode<
                                 FnSReconDir3MinusFullProd,
                                 Reference<QDPType<SU3Mat, OLattice<SU3Mat>>>,
                                 Reference<QDPType<FVec, OLattice<FVec>>>>,
                             OLattice<FVec>> &rhs,
                         const Subset &s)
    {
        const OLattice<SU3Mat> &u = static_cast<const OLattice<SU3Mat> &>(rhs.expression().left());
        const OLattice<FVec> &a = static_cast<const OLattice<FVec> &>(rhs.expression().right());

        if (s.hasOrderedRep())
        {
            int totalSize = s.end() - s.start() + 1;

            ordered_fused_spin_recon_user_arg arg = {u, a, d, s.start(), inlineAddSpinReconDir3MinusFull};

            dispatch_to_threads(totalSize, arg, ordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int site = s.start(); site <= s.end(); site++) {
    FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);
    
				   } */
        }
        else
        {
            const int *tab = s.siteTable().slice();

            int totalSize = s.numSiteTable();

            unordered_fused_spin_recon_user_arg arg = {u, a, d, tab, inlineAddSpinReconDir3MinusFull};

            dispatch_to_threads(totalSize, arg, unordered_fused_spin_recon_evaluate_function);

            ///////////////////
            // Original code
            ////////////////////
            /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int site = tab[j];
      
      FVec tmp;
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(0), tmp.elem(0));
      _inline_mult_su3_mat_vec(u.elem(site).elem(), a.elem(site).elem(1), tmp.elem(1));
      
      
      inlineAddSpinReconDir3MinusFull( (REAL *)&(tmp.elem(0).elem(0).real()),
				   (REAL *)&(d.elem(site).elem(0).elem(0).real()),
				   1);

				   }*/
        }
    }

} // namespace QDP;

#endif
