#ifndef __EIGEN_HELPERS__
#define __EIGEN_HELPERS__

#include <Utilities/LOG.h>
#include <Utilities/MATH.h>
#include <Eigen/CXX11/Tensor>
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

template<class T>
Matrix<T,3,1> cwiseMinMag(const Matrix<T,3,1>& v1,const Matrix<T,3,1>& v2)
{
    Matrix<T,3,1> result;result<<minmag(v1[0],v2[0]),minmag(v1[1],v2[1]),minmag(v1[2],v2[2]);
    return result;
}
}


using namespace Eigen;

// define utilities for manipulating form of Eigen matrices with non-scalar entries
namespace Mechanics{

template<class T>
using Quadruplet = std::tuple<int,int,int,T>; // row, col, o-index, value


template<class T,int dim>
void Flatten_Term(int row,int col,const Eigen::DiagonalMatrix<T,dim>& term,std::vector<Eigen::Triplet<T>>& flat_terms)
{
    for(int i=0;i<dim;i++){
        flat_terms.push_back(Eigen::Triplet<T>(row*dim+i,col*dim+i,term.diagonal()(i)));}
}

template<class T,int rows,int cols>
void Flatten_Matrix_Term(int row,int col,const Eigen::Matrix<T,rows,cols>& term,std::vector<Eigen::Triplet<T>>& flat_terms)
{
    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            flat_terms.push_back(Eigen::Triplet<T>(row*rows+i,col*cols+j,term(i,j)));}}
}

template<class T,int rows,int cols,int rows_per_block,int cols_per_block>
void Flatten_Matrix_Term(int row,int col,int block_row,int block_col,const Eigen::Matrix<T,rows_per_block,cols_per_block>& term,std::vector<Eigen::Triplet<T>>& flat_terms)
{
    // rows_per_block: rows per sub-block
    // rows: rows per outer block
    for(int i=0;i<rows_per_block;i++){
        for(int j=0;j<cols_per_block;j++){
            flat_terms.push_back(Eigen::Triplet<T>(rows_per_block*block_row+rows*row+i,cols_per_block*block_col+cols*col+j,term(i,j)));}}
}

template<class T,int d1,int d2,int d3,int d1_per_block,int d2_per_block,int d3_per_block>
void Flatten_Quadruplet_Term(int d1_index,int d2_index,int d3_index,int d1_block,int d2_block,int d3_block,const Eigen::TensorFixedSize<T,Eigen::Sizes<d1_per_block,d2_per_block,d3_per_block>>& term,std::vector<Quadruplet<T>>& flat_terms)
{
    // d_per_block: rows per sub-block
    // d: rows per outer block
    for(int i=0;i<d1_per_block;i++){
        for(int j=0;j<d2_per_block;j++){
            for(int k=0;k<d3_per_block;k++){
                flat_terms.push_back(Quadruplet<T>(d1_per_block*d1_block+d1*d1_index+i,d2_per_block*d2_block+d2*d2_index+j,d3_per_block*d3_block+d3*d3_index+k,term(i,j,k)));}}}
}

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

template<class T,int rows1,int rows2,int cols>
void Flatten_Matrices(const std::vector<Eigen::Triplet<Eigen::Matrix<T,rows1,cols>>>& block_terms1,
    const int row_base,
    const std::vector<Eigen::Triplet<Eigen::Matrix<T,rows2,cols>>>& block_terms2,
    Eigen::SparseMatrix<T>& flat_matrix)
{
    std::vector<Eigen::Triplet<T>> terms;
    for(const auto& block_term : block_terms1){
        for(int i=0;i<rows1;i++){
            for(int j=0;j<cols;j++){
                terms.push_back(Eigen::Triplet<T>(block_term.row()*rows1+i,block_term.col()*cols+j,block_term.value()(i,j)));}}}
    for(const auto& block_term : block_terms2){
        for(int i=0;i<rows2;i++){
            for(int j=0;j<cols;j++){
                terms.push_back(Eigen::Triplet<T>(row_base+block_term.row()*rows2+i,block_term.col()*cols+j,block_term.value()(i,j)));}}}
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
    matrix.resize(rowbase,1);matrix.setZero();
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

template<class T>
Eigen::Matrix<T,0,0> Cross_Product_Matrix(const Eigen::Matrix<T,0,1>& v)
{
    return Eigen::Matrix<T,0,0>();
}

// assumption: the columns of m1 should go in index 0, columns of m2 in index 1, and cross product results in index 2
template<class T,int d,typename DerivedLeft>
static TensorFixedSize<T,Sizes<d,d,d>> Cross_Product(const MatrixBase<DerivedLeft>& m1,const Matrix<T,d,d>& m2){
    TensorFixedSize<T,Sizes<d,d,d>> tensor;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Matrix<T,d,1> cross=m1.col(i).cross(m2.col(j));
            for(int k=0;k<3;k++){
                tensor(i,j,k)=cross(k);}}}
    return tensor;
}

template<class T,std::ptrdiff_t d,typename Derived>
// cross product between cross product matrix and dimension 2 of a tensor
static TensorFixedSize<T,Sizes<d,d,d>> Cross_Product(const MatrixBase<Derived>& m,const TensorFixedSize<T,Sizes<d,d,d>>& t){
    TensorFixedSize<T,Sizes<d,d,d>> tensor;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Matrix<T,d,1> tvec;
            for(int k=0;k<3;k++){tvec[k]=t(i,j,k);}
            Matrix<T,d,1> cross=m*tvec;
            for(int k=0;k<3;k++){
                tensor(i,j,k)=cross(k);}}}
    return tensor;
}

template<class T,int d,typename Derived>
static TensorFixedSize<T,Sizes<d,d,d>> Outer_Product(const MatrixBase<Derived>& m,const Matrix<T,d,1>& v,const std::array<int,3>& indices){
    TensorFixedSize<T,Sizes<d,d,d>> tensor;
    std::array<int,3> index;
    for(index[0]=0;index[0]<3;index[0]++){
        for(index[1]=0;index[1]<3;index[1]++){
            for(index[2]=0;index[2]<3;index[2]++){
                tensor(index[0],index[1],index[2])=m(index[indices[0]],index[indices[1]])*v(index[indices[2]]);}}}
    return tensor;
}

        
template<class T,int d,std::ptrdiff_t dt>
static Matrix<T,d,1> Contract(const TensorFixedSize<T,Sizes<dt,dt,dt>>& t,const Matrix<T,d,1>& v1,const Matrix<T,d,1>& v2){
    Matrix<T,d,1> result;result.setZero();
    std::array<int,3> index;
    for(index[0]=0;index[0]<3;index[0]++){
        for(index[1]=0;index[1]<3;index[1]++){
            for(index[2]=0;index[2]<3;index[2]++){
                result(index[2])+=t(index[0],index[1],index[2])*v1(index[0])*v2(index[1]);}}}
    return result;
}

template<class T,int d>
static T Contract(const Matrix<T,d,d>& m,const Matrix<T,d,1>& v1,const Matrix<T,d,1>& v2){
    return v1.transpose()*m*v2;
}

template<class T,int d,std::ptrdiff_t dt>
static Matrix<T,d,d> Contract(const TensorFixedSize<T,Sizes<dt,dt,dt>>& t,const Matrix<T,d,1>& v){
    Matrix<T,d,d> result;result.setZero();
    std::array<int,3> index{};
    for(index[0]=0;index[0]<3;index[0]++){
        for(index[1]=0;index[1]<3;index[1]++){
            for(index[2]=0;index[2]<3;index[2]++){
                result(index[1],index[2])+=t(index[1],index[2],index[0])*v(index[0]);
            }}}
    return result;
}


template<class T>
T Angle_Between(const Eigen::Matrix<T,3,1>& v1,const Eigen::Matrix<T,3,1>& v2)
{T s=v1.cross(v2).norm(),c=v1.dot(v2);return atan2(s,c);}

template<class T,int d>
void Safe_Normalize(Eigen::Matrix<T,d,1>& matrix)
{
    if(matrix.squaredNorm()){matrix.normalize();}
    else{matrix=Eigen::Matrix<T,d,1>::UnitX();}
}
}
#endif

