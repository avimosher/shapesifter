#ifndef __EIGEN_HELPERS__
#define __EIGEN_HELPERS__

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen{

template<typename _Scalar>
class Rotation1D : public RotationBase<Rotation1D<_Scalar>,1>
{
public:
    inline Rotation1D()
    {}

    inline Rotation1D(const Rotation1D& other)
    {}

    //inline Rotation1D inverse() const { return *this;}
    inline Rotation1D operator*(const Rotation1D& other) const
    {return Rotation1D();}
};

namespace internal {
template<typename _Scalar> struct traits<Rotation1D<_Scalar> >
{
  typedef _Scalar Scalar;
};
} // end namespace internal


}

namespace Mechanics{

template<class T,int rows,int cols>
void Flatten_Matrix(const std::vector<Eigen::Triplet<Eigen::Matrix<T,rows,cols>>>& block_terms,Eigen::SparseMatrix<T>& flat_matrix)
{
    std::vector<Eigen::Triplet<T>> terms;
    for(const auto& block_term : block_terms){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                terms.push_back(Eigen::Triplet<T>(block_term.row()*rows+i,block_term.col()*cols+j,block_term.value()(i,j)));
            }
        }
    }
    flat_matrix.setFromTriplets(terms.begin(),terms.end());
}

template<class T>
Eigen::Matrix<T,3,3> Cross_Product_Matrix(const Eigen::Matrix<T,3,1>& v)
{
    Eigen::Matrix<T,3,3> matrix;matrix<<0,v[2],-v[1],-v[2],0,v[0],v[1],-v[0],0;return matrix;
}

template<class T>
Eigen::Matrix<T,1,2> Cross_Product_Matrix(const Eigen::Matrix<T,2,1>& v)
{
    Eigen::Matrix<T,1,2> matrix;matrix<<-v[1],v[0];return matrix;
}

template<class T>
Eigen::Matrix<T,0,1> Cross_Product_Matrix(const Eigen::Matrix<T,1,1>& v)
{
    return Eigen::Matrix<T,0,1>();
}

}
#endif

