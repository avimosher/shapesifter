#ifndef __EIGEN_HELPERS__
#define __EIGEN_HELPERS__

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <json/json.h>

// define Rotation1D to match Rotation2D and Quaternion
namespace Eigen{
template<typename _Scalar>
class Rotation1D : public RotationBase<Rotation1D<_Scalar>,1>
{
public:
    inline Rotation1D() {}
    inline Rotation1D(const Rotation1D& other) {}
    inline Rotation1D operator*(const Rotation1D& other) const
    {return Rotation1D();}
    static inline Rotation1D Identity()
    {return Rotation1D();}
};

namespace internal {
template<typename _Scalar> struct traits<Rotation1D<_Scalar> >
{
  typedef _Scalar Scalar;
};
} // end namespace internal

// define cout overloads for Eigen rotation types
template<class T>
std::ostream& operator<<(std::ostream& os,const Rotation1D<T>&)
{}

template<class T>
std::ostream& operator<<(std::ostream& os,const Rotation2D<T>& rotation)
{
    os<<rotation.angle();
    return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os,const Quaternion<T>& rotation)
{
    os<<rotation.coeffs();
    return os;
}
}


// define utilities for manipulating form of Eigen matrices with non-scalar entries
namespace Mechanics{
template<class T,int rows,int cols>
void Flatten_Matrix(const std::vector<Eigen::Triplet<Eigen::Matrix<T,rows,cols>>>& block_terms,Eigen::SparseMatrix<T>& flat_matrix)
{
    std::vector<Eigen::Triplet<T>> terms;
    for(const auto& block_term : block_terms){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                terms.push_back(Eigen::Triplet<T>(block_term.row()*rows+i,block_term.col()*cols+j,block_term.value()(i,j)));}}}
    flat_matrix.setFromTriplets(terms.begin(),terms.end());
}

template<class T>
void Merge_Block_Matrices(const Eigen::Matrix<Eigen::SparseMatrix<T>,Eigen::Dynamic,Eigen::Dynamic>& block_matrix,Eigen::SparseMatrix<T>& matrix)
{
    std::vector<Eigen::Triplet<T>> triplets;
    int rowbase=0,colbase=0;
    int cols=0;
    for(int i=0;i<block_matrix.rows();i++){
        colbase=0;
        for(int j=0;j<block_matrix.cols();j++){
            auto& block=block_matrix(i,j);
            for(int k=0;k<block.outerSize();k++){
                for(typename Eigen::SparseMatrix<T>::InnerIterator it(block,k);it;++it){
                    triplets.push_back(Eigen::Triplet<T>(rowbase+it.row(),colbase+it.col(),it.value()));}}
            colbase+=block_matrix(i,j).cols();}
        cols=std::max(colbase,cols);
        rowbase+=block_matrix(i,0).rows();}
    matrix.resize(rowbase,cols);
    matrix.setFromTriplets(triplets.begin(),triplets.end());
}

template<class T>
void Merge_Block_Vectors(const Eigen::Matrix<Eigen::Matrix<T,Eigen::Dynamic,1>,Eigen::Dynamic,1>& block_matrix,Eigen::Matrix<T,Eigen::Dynamic,1>& matrix)
{
    int rowbase=0;
    for(int i=0;i<block_matrix.rows();i++){
        auto& block=block_matrix(i,0);
        rowbase+=block.rows();}
    matrix.resize(rowbase,1);
    rowbase=0;
    for(int i=0;i<block_matrix.rows();i++){
        auto& block=block_matrix(i,0);
        matrix.block(rowbase,0,block.rows(),1)=block;
        rowbase+=block.rows();}
}

template<class T>
Eigen::Matrix<T,3,3> Cross_Product_Matrix(const Eigen::Matrix<T,3,1>& v)
{
    Eigen::Matrix<T,3,3> matrix;matrix<<0,-v[2],v[1],v[2],0,-v[0],-v[1],v[0],0;return matrix;
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

template<class T,int d>
void Get_Vector(Json::Value& node,Eigen::Matrix<T,d,1>& vector) {
    for(int i=0;i<d;i++){vector[i]=node[i].asDouble();}
}

}
#endif

