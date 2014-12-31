#ifndef __EIGEN_HELPERS__
#define __EIGEN_HELPERS__

#include <Eigen/Dense>

namespace Mechanics{

template<class T,int rows,int cols>
void Flatten_Matrix(const std::vector<Triplet<Matrix<T,rows,cols>>>& block_terms,SparseMatrix<T>& flat_matrix)
{
    std::vector<Triplet<T>> terms;
    for(const auto& block_term : block_terms){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                terms.push_back(Triplet<T>(block_term.row()*rows+i,block_term.col()*cols+j,block_term.value()(i,j)));
            }
        }
    }
    flat_matrix.setFromTriplets(terms);
}

}

#endif

