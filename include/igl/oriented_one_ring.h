// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_ORIENTED_ONE_RING
#define IGL_ORIENTED_ONE_RING

#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

// Computes the oriented the one-ring neighborhood (faces and indices within faces)
// of a given vertex in a mesh. Care is taken to get it right when we are at a boundary.
//
// Inputs:
//   V                #V by 3 list of the vertex positions
//   F                #F by 3 list of the faces (must be triangles)
//   VF               #V list of lists of incident faces (adjacency list), e.g.
//                    as returned by igl::vertex_triangle_adjacency
//   TT               #F by 3 triangle to triangle adjacent matrix (e.g. computed
//                    via igl:triangle_triangle_adjacency)
//   vi               the selected one ring
//
// Output:
//   fi               #numOneRingFaces by 1 list of the sequentially visited faces in the one ring neighborhood.
//                    The one-ring is traversed in COUNTERCLOCKWISE order with respect to the outward normal.
//   fik              #numOneRingFaces by 1, fik[j] is the index of the
//                    *edge* (fi[j],fi[j+1]) inside face fi[j]. Not that
//                    this is not the same as the index of vi inside fi[j]
//                    - depends on orientation of the edge.
//

namespace igl {
    template <typename DerivedF, typename VFType, typename DerivedTT>
    IGL_INLINE void oriented_one_ring(const Eigen::PlainObjectBase<DerivedF> &F,
                                      const std::vector<std::vector<VFType> >& VF,
                                      const Eigen::PlainObjectBase<DerivedTT>& TT,
                                      const int vi,
                                      Eigen::VectorXi &fi,
                                      Eigen::VectorXi &fik);
}

#ifndef IGL_STATIC_LIBRARY
#include "oriented_one_ring.cpp"
#endif

#endif /* defined(IGL_ORIENTED_ONE_RING) */
