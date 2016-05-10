#include <Data/DATA.h>
#include <Data/RIGID_STRUCTURE_DATA.h>
#include <Driver/SIMULATION.h>
#include <Parsing/PARSER_REGISTRY.h>
using namespace Mechanics;
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Size()
{
    return structures.size();
}
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Velocity_DOF() const
{
    return TWIST<TV>::STATIC_SIZE*structures.size();
}
///////////////////////////////////////////////////////////////////////
template<class TV> int RIGID_STRUCTURE_DATA<TV>::
Position_DOF() const
{
    return FRAME<TV>::STATIC_SIZE*structures.size();
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Identify_DOF(int index) const
{
    int body_index=index/s;
    LOG::cout<<Name()<<": body "<<structures[body_index]->name<<" dof "<<index-body_index*s<<std::endl;
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Pack_Velocities(Block<Matrix<T,Dynamic,1>>& velocities) const
{
    for(int i=0;i<structures.size();i++){
        velocities.template block<TWIST<TV>::STATIC_SIZE,1>(i*TWIST<TV>::STATIC_SIZE,0)=structures[i]->twist.Pack();}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Unpack_Velocities(const Matrix<T,Dynamic,1>& velocities)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->twist.Unpack(velocities.template block<TWIST<TV>::STATIC_SIZE,1>(i*TWIST<TV>::STATIC_SIZE,0));}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Pack_Positions(Block<Matrix<T,Dynamic,1>>& positions) const
{
    for(int i=0;i<structures.size();i++){
        positions.template block<FRAME<TV>::STATIC_SIZE,1>(i*FRAME<TV>::STATIC_SIZE,0)=structures[i]->frame.Pack();}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Unpack_Positions(const Matrix<T,Dynamic,1>& positions)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->frame.Unpack(positions.template block<FRAME<TV>::STATIC_SIZE,1>(i*FRAME<TV>::STATIC_SIZE,0));}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Store_Errors(const Matrix<T,Dynamic,1>& errors)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->error.Unpack(errors.template block<s,1>(i*s,0));}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Step(const DATA<TV>& data)
{
    for(int i=0;i<structures.size();i++){
        //LOG::cout<<structures[i]->name<<" velocity: "<<structures[i]->twist.linear.transpose()<<std::endl;
        structures[i]->frame=Updated_Frame(data,structures[i]->frame,structures[i]->twist);}
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Eliminate_Rotation(const DATA<TV>& data)
{
    for(int i=0;i<structures.size();i++){
        structures[i]->frame.orientation=ROTATION<TV>();
        structures[i]->twist.linear=TV();
        structures[i]->twist.angular=T_SPIN();
    }
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Inertia(const T dt,std::vector<Triplet<T>>& force_terms,SparseMatrix<T>& inverse_inertia,Matrix<T,Dynamic,1>& rhs)
{
    T one_over_dt=1/dt;
    std::vector<Triplet<T>> inverse_inertia_terms;
    for(int i=0;i<structures.size();i++){
        auto structure=structures[i];
        DiagonalMatrix<T,s> inertia_matrix=one_over_dt*structures[i]->Inertia_Matrix();
        Flatten_Term(i,i,inertia_matrix,force_terms);
        rhs.template block<s,1>(s*i,0)=inertia_matrix*structures[i]->twist.Pack();
        DiagonalMatrix<T,s> inverse_inertia_matrix(inertia_matrix.inverse());
        Flatten_Term(i,i,inverse_inertia_matrix,inverse_inertia_terms);}
    inverse_inertia.resize(structures.size()*s,structures.size()*s);
    inverse_inertia.setFromTriplets(inverse_inertia_terms.begin(),inverse_inertia_terms.end());
}
///////////////////////////////////////////////////////////////////////
template<class TV> void RIGID_STRUCTURE_DATA<TV>::
Kinematic_Projection(SparseMatrix<T>& kinematic_projection)
{
    int index=0;
    std::vector<Triplet<T>> projection_terms;
    for(int i=0;i<structures.size();i++){
        auto structure=structures[i];
        if(!structure->kinematic){
            DiagonalMatrix<T,s> projection_matrix;projection_matrix.setIdentity();
            Flatten_Term(index++,i,projection_matrix,projection_terms);
        }}
    kinematic_projection.resize(s*index,s*structures.size());
    kinematic_projection.setFromTriplets(projection_terms.begin(),projection_terms.end());
}
///////////////////////////////////////////////////////////////////////
GENERIC_CEREAL_REGISTRATION(RIGID_STRUCTURE_DATA)
GENERIC_TYPE_DEFINITION(RIGID_STRUCTURE_DATA)
DEFINE_AND_REGISTER_PARSER(RIGID_STRUCTURE_DATA,void)
{
    simulation.data.template Find_Or_Create<RIGID_STRUCTURE_DATA<TV>>();
    return 0;
}
///////////////////////////////////////////////////////////////////////
