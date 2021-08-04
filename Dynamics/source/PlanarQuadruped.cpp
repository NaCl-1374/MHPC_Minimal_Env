#include "PlanarQuadruped.h"
#include "CasadiGen.h"
#include <vector>
#include "orientation_tools.h"


template <typename T>
void PlanarQuadruped<T>::dynamics(VecM<T, xsize_WB> &x, VecM<T, usize_WB> &u, VecM<T, xsize_WB> &x_next, VecM<T, ysize_WB> &y, int mode)
{    
    vector<T *> arg = {x.data(), u.data()};
    vector<T *> res = {_xdot.data(), y.data()};
    switch (mode)
    {
    case 2:
    case 4: // flight
        casadi_interface(arg, res, x.size(), Dyn_FL, Dyn_FL_sparsity_out, Dyn_FL_work);
        break;
    case 1: // back stance
        casadi_interface(arg, res, x.size(), Dyn_BS, Dyn_BS_sparsity_out, Dyn_BS_work);
        break;
    case 3: // front stance
        casadi_interface(arg, res, x.size(), Dyn_FS, Dyn_FS_sparsity_out, Dyn_FS_work);
        break;
    }
   
    x_next = x + _xdot * this->_dt;
}

template<typename T>
void PlanarQuadruped<T>::dynamics_par(VecM<T, xsize_WB> &x, VecM<T, usize_WB> &u, MatMN<T, xsize_WB, xsize_WB> &A, MatMN<T, xsize_WB, usize_WB> &B,
                                      MatMN<T, ysize_WB, xsize_WB> &C, MatMN<T, ysize_WB, usize_WB> &D, int mode)
{
   vector<T *> arg = {x.data(), u.data()};
   vector<T *> res = {_Ac.data(), _Bc.data(), C.data(), D.data()};

    switch (mode)
    {
    case 2:
    case 4: // flight
        casadi_interface(arg, res, _Ac.size(), Dyn_FL_par, Dyn_FL_par_sparsity_out, Dyn_FL_par_work);
        break;
    case 1: // back stance
        casadi_interface(arg, res, _Ac.size(), Dyn_BS_par, Dyn_BS_par_sparsity_out, Dyn_BS_par_work);
        break;
    case 3: // front stance
        casadi_interface(arg, res, _Ac.size(), Dyn_FS_par, Dyn_FS_par_sparsity_out, Dyn_FS_par_work);
        break;
    } 

    A = MatMN<T,xsize_WB,xsize_WB>::Identity() + _Ac * this->_dt;
    B = _Bc * this->_dt;
   
}

template<typename T>
// void PlanarQuadruped<T>::resetmap(DVec<T> &x, DVec<T> &x_next, DVec<T> &y, int mode)

void PlanarQuadruped<T>::resetmap(VecM<T, xsize_WB> &x, VecM<T, xsize_WB> &x_next, VecM<T, ysize_WB> &y, int mode)
{
    vector<T *> arg = {x.data()};
    vector<T *> res = {x_next.data(), y.data()};

    switch (mode)
    {
    case 1:
    case 3:
        x_next = x;
        break;
    case 2: //front impact
        casadi_interface(arg, res, x.size(), Imp_F, Imp_F_sparsity_out, Imp_F_work);
        break;
    case 4: // back impact
        casadi_interface(arg, res, x.size(), Imp_B, Imp_B_sparsity_out, Imp_B_work);
        break;
    }

    
}

template<typename T>
// void PlanarQuadruped<T>::resetmap_par(DVec<T> &x, DMat<T> &Px, int mode)
void PlanarQuadruped<T>::resetmap_par(VecM<T, xsize_WB> &x, MatMN<T, xsize_WB, xsize_WB> &Px, int mode)
{
    vector<T *> arg = {x.data()};
    vector<T *> res = {Px.data()};    

    switch (mode)
    {
    case 1:
    case 3: // smooth transition
        Px = MatMN<T,xsize_WB,xsize_WB>::Identity();      
        break;
    case 2: // front impact
        casadi_interface(arg, res, Px.size(), Imp_F_par, Imp_F_par_sparsity_out, Imp_F_par_work);
        break;
    case 4: // back impact
        casadi_interface(arg, res, Px.size(), Imp_B_par, Imp_B_par_sparsity_out, Imp_B_par_work);
        break;
    }   
}

template<typename T>
void PlanarQuadruped<T>::FootJacobian(VecM<T, xsize_WB> &x, MatMN<T,2,qsize_WB> &J, MatMN<T,2,qsize_WB> &Jd, int foot)
{
   vector<T *> arg = {x.data()};
   vector<T *> res = {J.data(), Jd.data()};

    if (legID::FLEG == foot) // front foot
    {
       casadi_interface(arg, res, x.size(), Jacob_F, Jacob_F_sparsity_out, Jacob_F_work);
    }
    else if (legID::HLEG == foot) // back foot
    {
       casadi_interface(arg, res, x.size(), Jacob_B, Jacob_B_sparsity_out, Jacob_B_work);
    }

}

template<typename T>
void PlanarQuadruped<T>::build_quadruped()
{
    bodyMass = 3.3;
    bodyWidth = 0.05*2;
    bodyLength = 0.19*2;
    bodyCoM.setZero();

    hipLinkMass = 0.634*2;
    hipLinkWidth = 0.03;
    hipLinkLength = 0.209;
    hipLinkCoM.setZero();
    hipLoc << bodyLength/2, 0;

    kneeLinkMass = 0.064*2;
    kneeLinkWidth = 0.03;
    kneeLinkLength = 0.195;
    kneeLinkCoM.setZero();
    kneeLoc << 0, -hipLinkLength;
}


template<typename T>
MatMN<T,4,4> PlanarQuadruped<T>::get_homoTransformation(VecM<T, qsize_WB> &q, size_t linkidx)
{
    MatMN<T,4,4> T_body;
    T_body = ori::homoTransformation(ori::coordinateRotation(ori::CoordinateAxis::Y, -q[2]), VecM<T,3>(q[0], 0, q[1]));
    if(linkidx == linkID2D::body) {return T_body;}

    MatMN<T,4,4> T_Fhip, T_Hhip;
    T_Fhip = T_body *
             ori::homoTransformation(ori::coordinateRotation(ori::CoordinateAxis::Y, -q[3]), VecM<T,3>(hipLoc[0], 0, hipLoc[1]));
    if(linkidx == linkID2D::F_hip) {return T_Fhip;}

    T_Hhip = T_body *
             ori::homoTransformation(ori::coordinateRotation(ori::CoordinateAxis::Y, -q[5]), VecM<T,3>(-hipLoc[0], 0, hipLoc[1]));
    if(linkidx == linkID2D::H_hip) {return T_Hhip;}                     

    MatMN<T,4,4> T_Fknee, T_Hknee;
    T_Fknee = T_Fhip *
             ori::homoTransformation(ori::coordinateRotation(ori::CoordinateAxis::Y, -q[4]), VecM<T,3>(kneeLoc[0], 0, kneeLoc[1]));
    if(linkidx == linkID2D::F_knee) {return T_Fknee;}                     

    T_Hknee = T_Hhip *
             ori::homoTransformation(ori::coordinateRotation(ori::CoordinateAxis::Y, -q[6]), VecM<T,3>(kneeLoc[0], 0, kneeLoc[1]));
    if(linkidx == linkID2D::H_knee) {return T_Hknee;}   

    MatMN<T,4,4> T_Ffoot, T_Hfoot;
    T_Ffoot = T_Fknee * 
             ori::homoTransformation(MatMN<T,3,3>::Identity(), VecM<T,3>(0, 0, -kneeLinkLength));
    if(linkidx == linkID2D::F_foot) {return T_Ffoot;}

    T_Hfoot = T_Hknee *
             ori::homoTransformation(MatMN<T,3,3>::Identity(), VecM<T,3>(0, 0, -kneeLinkLength));
    return T_Hfoot;
        
}

template<typename T>
VecM<T,2> PlanarQuadruped<T>::get_contact_position(VecM<T,qsize_WB> &q, size_t linkidx, VecM<T,2> &contactLoc)
{
    VecM<T,4> p_world;
    VecM<T,2> p_world_2D;
    MatMN<T,4,4> H = get_homoTransformation(q, linkidx);
    p_world = H * VecM<T,4>(contactLoc[0], 0, contactLoc[1], 1);
    p_world_2D << p_world[0], p_world[2];
    return p_world_2D;
}

/*
  @brief: get leg extension (vector from hip joint to foot)
  @params: 
          q: generalized joints of full planar quadruped 
          legidx: 0->front leg, 1->rear leg
*/
template<typename T>
VecM<T,2> PlanarQuadruped<T>::get_leg_ext_vec(VecM<T,qsize_WB> &q, int legidx)
{
    VecM<T, 2> loc;
    loc.setZero();
    if(legID::FLEG == legidx)
    return (get_contact_position(q, linkID2D::F_foot, loc) - get_contact_position(q, linkID2D::F_hip, loc));

    if(legID::HLEG == legidx)
    return (get_contact_position(q, linkID2D::H_foot, loc) - get_contact_position(q, linkID2D::H_hip, loc));

}

template class PlanarQuadruped<double>; 


