#ifndef __MATRIX_BUNDLE__
#define __MATRIX_BUNDLE__

#include <Data/DATA.h>
#include <Force/FORCE.h>
#include <Utilities/EIGEN_HELPERS.h>

namespace Mechanics{

template<class TV>
class MATRIX_BUNDLE
{
    typedef typename TV::Scalar T;
public:
    Matrix<SparseMatrix<T>,Dynamic,Dynamic> jacobian_blocks;
    Matrix<Matrix<T,Dynamic,1>,Dynamic,1> right_hand_side_blocks;
    std::vector<std::vector<Triplet<T>>> matrix_block_terms;
    Matrix<std::vector<Triplet<T>>,Dynamic,Dynamic> hessian_block_terms;

    void Initialize(const DATA<TV>& data,const FORCE<TV>& force){
        int full_size=data.size()+force.size();
        jacobian_blocks.resize(full_size,full_size);
        hessian_block_terms.resize(full_size,full_size);
        right_hand_side_blocks.resize(full_size);
        matrix_block_terms.resize(data.size());

        for(int i=0;i<data.size();i++){
            int data_size=data[i]->Velocity_DOF();
            matrix_block_terms[i].clear();
            jacobian_blocks(i,i).resize(data_size,data_size);
            hessian_block_terms(i,i).clear();
            right_hand_side_blocks[i].resize(data_size);
            right_hand_side_blocks[i].setZero();}
    }

    template<class SUBTYPE> int Index(const DATA<TV>& data,const FORCE<TV>& force){
        int data_index=data.template Index<SUBTYPE>();
        if(data_index<data.size()){return data_index;}
        return data_index+force.template Index<SUBTYPE>();
    }

    template<class SUBTYPE1,class SUBTYPE2>
    SparseMatrix<T>& Matrix_Block(const DATA<TV>& data,const FORCE<TV>& force){
        return jacobian_blocks(Index<SUBTYPE1>(data,force),Index<SUBTYPE2>(data,force));
    }

    template<class SUBTYPE1,class SUBTYPE2>
    SparseMatrix<T>& Matrix_Block(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE1& row_object,const SUBTYPE2& column_object){
        return jacobian_blocks(Index<SUBTYPE1>(data,force),Index<SUBTYPE2>(data,force));
    }

    template<class SUBTYPE>
    Matrix<T,Dynamic,1>& RHS(const DATA<TV>& data,const FORCE<TV>& force){
        return right_hand_side_blocks(Index<SUBTYPE>(data,force));
    }

    template<class SUBTYPE>
    std::vector<Triplet<T>>& Matrix_Block_Terms(const DATA<TV>& data,const FORCE<TV>& force){
        return matrix_block_terms[Index<SUBTYPE>(data,force)];
    }

    void Scale_Blocks(const DATA<TV>& data,const FORCE<TV>& force,std::vector<SparseMatrix<T>>& kinematic_projection_matrices,std::vector<SparseMatrix<T>>& inverse_inertia_matrices){
        for(int i=0;i<data.size();i++){
            jacobian_blocks(i,i).setFromTriplets(matrix_block_terms[i].begin(),matrix_block_terms[i].end());
            for(int j=0;j<data.size()+force.size();j++){
                jacobian_blocks(i,j)=inverse_inertia_matrices[i]*jacobian_blocks(i,j);}
            right_hand_side_blocks[i]=kinematic_projection_matrices[i]*inverse_inertia_matrices[i]*right_hand_side_blocks[i];}
        for(int i=0;i<data.size()+force.size();i++){
            for(int j=0;j<data.size()+force.size();j++){
                if(i<data.size()){
                    jacobian_blocks(i,j)=kinematic_projection_matrices[i]*jacobian_blocks(i,j);}
                if(j<data.size()){
                    jacobian_blocks(i,j)=jacobian_blocks(i,j)*kinematic_projection_matrices[j].transpose();}}}
    }

    template<class SUBTYPE1,class SUBTYPE2,class T_MATRIX>
    void Flatten_Jacobian_Block(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE1& row_object,const SUBTYPE2& column_object,const std::vector<Triplet<T_MATRIX>>& terms){
        SparseMatrix<T>& block=Matrix_Block(data,force,row_object,column_object);
        block.resize(row_object.DOF(),column_object.DOF());
        Flatten_Matrix(terms,block);
    }
};
}

#endif

