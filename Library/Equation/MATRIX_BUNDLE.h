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
    Matrix<SparseMatrix<T>,Dynamic,Dynamic> hessian_blocks;
    Matrix<Matrix<T,Dynamic,1>,Dynamic,1> error_blocks;
    std::vector<std::vector<Triplet<T>>> jacobian_block_terms;
    Matrix<std::vector<Triplet<T>>,Dynamic,Dynamic> hessian_block_terms;

    void Initialize(const DATA<TV>& data,const FORCE<TV>& force){
        int full_size=data.size()+force.size();
        jacobian_blocks.resize(full_size,full_size);
        hessian_blocks.resize(full_size,full_size);
        hessian_block_terms.resize(full_size,full_size);
        error_blocks.resize(full_size);
        jacobian_block_terms.resize(data.size());

        for(int i=0;i<data.size();i++){
            int data_size=data[i]->Velocity_DOF();
            jacobian_block_terms[i].clear();
            jacobian_blocks(i,i).resize(data_size,data_size);
            hessian_block_terms(i,i).clear();
            hessian_blocks(i,i).resize(data_size,data_size);
            error_blocks[i].resize(data_size);
            error_blocks[i].setZero();}
    }

    template<class SUBTYPE> int Index(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE& object){
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
        return jacobian_blocks(Index(data,force,row_object),Index(data,force,column_object));
    }

    template<class SUBTYPE1,class SUBTYPE2>
    SparseMatrix<T>& Hessian_Block(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE1& row_object,const SUBTYPE2& column_object){
        return hessian_blocks(Index(data,force,row_object),Index(data,force,column_object));
    }

    template<class SUBTYPE>
    Matrix<T,Dynamic,1>& RHS(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE& object){
        return error_blocks(Index(data,force,object));
    }

    template<class SUBTYPE>
    std::vector<Triplet<T>>& Matrix_Block_Terms(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE& object){
        return jacobian_block_terms[Index(data,force,object)];
    }

    template<class SUBTYPE1,class SUBTYPE2>
    std::vector<Triplet<T>>& Hessian_Block_Terms(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE1& o1,const SUBTYPE2& o2){
        return hessian_block_terms(Index(data,force,o1),Index(data,force,o2));
    }

    void Scale_Errors(const DATA<TV>& data,const FORCE<TV>& force,std::vector<SparseMatrix<T>>& kinematic_projection_matrices,std::vector<SparseMatrix<T>>& inverse_inertia_matrices){
        for(int i=0;i<data.size();i++){
            error_blocks[i]=kinematic_projection_matrices[i]*inverse_inertia_matrices[i]*error_blocks[i];}
    }

    void Scale_Derivatives(const DATA<TV>& data,const FORCE<TV>& force,std::vector<SparseMatrix<T>>& kinematic_projection_matrices,std::vector<SparseMatrix<T>>& inverse_inertia_matrices){
        for(int i=0;i<data.size();i++){
            jacobian_blocks(i,i).setFromTriplets(jacobian_block_terms[i].begin(),jacobian_block_terms[i].end());
            hessian_blocks(i,i).setFromTriplets(hessian_block_terms(i,i).begin(),hessian_block_terms(i,i).end());
            for(int j=0;j<data.size()+force.size();j++){
              jacobian_blocks(i,j)=inverse_inertia_matrices[i]*jacobian_blocks(i,j);}
        }
        /*for(int i=0;i<data.size()+force.size();i++){
            for(int j=0;j<data.size()+force.size();j++){
                if(i<data.size()){
                    jacobian_blocks(i,j)=kinematic_projection_matrices[i]*jacobian_blocks(i,j);}
                if(j<data.size()){
                jacobian_blocks(i,j)=jacobian_blocks(i,j)*kinematic_projection_matrices[j].transpose();}}}*/
    }

    template<class SUBTYPE1,class SUBTYPE2,class T_MATRIX>
    void Flatten_Jacobian_Block(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE1& row_object,const SUBTYPE2& column_object,const std::vector<Triplet<T_MATRIX>>& terms){
        SparseMatrix<T>& block=Matrix_Block(data,force,row_object,column_object);
        block.resize(row_object.DOF(),column_object.DOF());
        Flatten_Matrix(terms,block);
    }

    template<class SUBTYPE1,class SUBTYPE2>
    void Flatten_Hessian_Block(const DATA<TV>& data,const FORCE<TV>& force,const SUBTYPE1& row_object,const SUBTYPE2& column_object){
        SparseMatrix<T>& block=Hessian_Block(data,force,row_object,column_object);
        block.resize(row_object.DOF(),column_object.DOF());
    }
};
}

#endif

