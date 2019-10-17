
#include "oriented_one_ring.h"
#include <iostream>

template <typename DerivedF, typename VFType, typename DerivedTT>
IGL_INLINE void igl::oriented_one_ring(const Eigen::PlainObjectBase<DerivedF> &F,
                                       const std::vector<std::vector<VFType> >& VF,
                                       const Eigen::PlainObjectBase<DerivedTT>& TT,
                                       const int vi,
                                       Eigen::VectorXi &fi,
                                       Eigen::VectorXi &fik)
{
    fi.resize(VF[vi].size()+1,1);
    fik.resize(VF[vi].size(),1);
    
    //find a face to start from
    //take care if it's a boundary vertex
    //first, check if the vertex is on a boundary
    //then there must be two faces that are on the boundary
    //(other cases not supported)
    
    int fstart = -1, fstart_k = -1;
    int ind = 0;
    for (int i =0; i<VF[vi].size(); ++i)
    {
        int fi = VF[vi][i];
        for (int  j=0; j<3; ++j)
            if (F(fi,j)==vi && TT(fi,j) == -1)
            {
                ind ++;
                fstart = fi;
                fstart_k = j;
                //        break;
            }
    }
    if (ind >1 )
    {
        std::cerr<<"igl::polyvector_field_one_ring_matchings -- vertex "<<vi<< " is on an unusual boundary"<<std::endl;
        exit(1);
    }
    if (fstart == -1)
    {
        fstart = VF[vi][0];
        for (int  j=0; j<3; ++j)
            if (F(fstart,j)==vi)
                fstart_k = j;
    }
    
    
    int current_face = fstart;
    int current_face_k = fstart_k;
    int i =0;
    fi[i] = current_face;
    int next_face = -1;
    while (next_face != fstart && current_face!=-1)
    {
        // look for the vertex
        int j=-1;
        for (unsigned z=0; z<3; ++z)
            if (F(current_face,(z+1)%3) == vi)
            {
                j=z;
                break;
            }
        assert(j!=-1);
        
        fik[i] = j;
        next_face = TT(current_face, j);
        ++i;
        
        fi[i] = next_face;
        current_face = next_face;
    }
    
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::oriented_one_ring<Eigen::Matrix<int, -1, -1, 0, -1, -1>, int, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::__1::vector<std::__1::vector<int, std::__1::allocator<int> >, std::__1::allocator<std::__1::vector<int, std::__1::allocator<int> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, int, Eigen::Matrix<int, -1, 1, 0, -1, 1>&, Eigen::Matrix<int, -1, 1, 0, -1, 1>&);
#endif

