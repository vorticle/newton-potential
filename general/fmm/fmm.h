/*
 * Copyright (C) 2020 Matthias Kirchhart
 *
 * This is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3, or (at your option) any later
 * version.
 *
 * This software is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this software; see the file COPYING.  If not see http://www.gnu.org/licenses.
 */
#ifndef FMM_FMM_H
#define FMM_FMM_H

#include <fmm/tree.h>

////////////////////////////
// ‘Public’ interface     //
////////////////////////////


namespace fmm
{

template <typename strategy, typename source_iter, typename target_iter>
void fmm( strategy& S, source_iter sbegin, source_iter send,
                       target_iter tbegin, target_iter tend );
}








/////////////////////////////////
// Implementation interface.   //
/////////////////////////////////


namespace fmm
{

template <typename strategy, typename source_box>
void upward_pass( source_box *n, strategy *S );

template <typename strategy, typename source_box, typename target_box>
void dual_tree_traversal( target_box *A, source_box *B, strategy *S );

template <typename strategy, typename source_box, typename target_box>
void interact( target_box *A, source_box *B, strategy *S );

template <typename strategy, typename target_box>
void downward_pass( target_box *n, strategy *S );

}




/////////////////////////
//  Implementation.    //
/////////////////////////


namespace fmm
{

template <typename strategy, typename source_iter, typename target_iter>
void fmm( strategy& S, source_iter sbegin, source_iter send,
                       target_iter tbegin, target_iter tend )
{
    using Mcoeff_t = typename strategy::Mcoeff_t;
    using Lcoeff_t = typename strategy::Mcoeff_t;
    using source_tree = tree<Mcoeff_t,source_iter>;
    using target_tree = tree<Lcoeff_t,target_iter>;

    if ( sbegin == send || tbegin == tend ) return;

    source_tree sources { sbegin, send, S };
    target_tree targets { tbegin, tend, S };

    #pragma omp parallel
    #pragma omp single
    {
        upward_pass        ( sources.root.get(), &S );
        dual_tree_traversal( targets.root.get(), sources.root.get(), &S );
        downward_pass( targets.root.get(), &S );
    }
}

template <typename strategy, typename source_box>
void upward_pass( source_box *n, strategy *S )
{
    if ( n->is_leaf() )
    {
        S->p2m( n->coeffs, n->bounds.centre, *(n->begin) );
    }
    else
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                #pragma omp task untied
                upward_pass( n->children[ i ].get(), S );
            }
        }
        #pragma omp taskwait

        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                S->m2m( n->coeffs,  n->children[ i ]->coeffs,
                        n->bounds.centre - n->children[ i ]->bounds.centre );
            }
        }
    }
}

template <typename strategy, typename source_box, typename target_box>
void dual_tree_traversal( target_box *A, source_box *B, strategy *S )
{
    if ( ( A->bounds.radius > B->bounds.radius && ! A->is_leaf() ) ||
         ( ! A->is_leaf() && B->is_leaf() ) )
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( A->children[ i ] )
            {
                #pragma omp task untied
                interact( A->children[ i ].get(), B, S );
            }
        }
        #pragma omp taskwait
    }
    else if ( ! B->is_leaf() )
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( B->children[ i ] )
            {
                interact( A, B->children[ i ].get(), S );
            }
        }
    }
    else interact( A, B, S );
}

template <typename strategy, typename source_box, typename target_box>
void interact( target_box *A, source_box *B, strategy *S )
{
    if ( S->mac( A->bounds, B->bounds ) )
    {
        S->m2l( A->coeffs, B->coeffs, A->bounds.centre - B->bounds.centre );
    }
    else if ( A->is_leaf() && B->is_leaf() )
    {
        S->p2p( *(A->begin), *(B->begin) );
    }
    else
    {
        dual_tree_traversal( A, B, S );
    }
}

template <typename strategy, typename target_box>
void downward_pass( target_box *n, strategy *S )
{
    if ( n->is_leaf() )
    {
        S->l2p( n->coeffs, n->bounds.centre, *(n->begin) );
    }
    else
    {
        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                S->l2l( n->children[ i ]->coeffs, n->coeffs,
                        n->children[ i ]->bounds.centre - n->bounds.centre );
            }
        }

        for ( uint i = 0; i < 8; ++i )
        {
            if ( n->children[ i ] )
            {
                #pragma omp task untied
                downward_pass( n->children[ i ].get(), S );
            }
        }
        #pragma omp taskwait
    }
}

}

#endif

