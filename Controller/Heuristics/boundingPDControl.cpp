#include"boundingPDControl.h"
template<typename T>
void bounding_PDcontrol(PlanarQuadruped<T>* robot, ModelState<T,xsize_WB,usize_WB,ysize_WB> *ms, int modeidx, size_t N_TIME_STEPS)
{
    VecM<T,4> qjoint_nom(PI/4, -PI*7/12, PI/4, -PI*7/12);
    T legext_nom = 0.2462;
    VecM<T,2>  leg_ext_vec(0,0);
    MatMN<T,4,4> Kp, Kd;
    Kp = 5*VecM<T,4>(8,1,12,10).asDiagonal();
    Kd.setIdentity();
    T Kspring = 2200;
    VecM<T,2> F;
    MatMN<T,2,7> J, Jd;
    VecM<T, 7> q;
    q.setZero();

    for (size_t k = 0; k < N_TIME_STEPS-1; k++)
    {
        F.setZero();
        J.setZero();
        Jd.setZero();
        q = ms[k].x.head(7);

        switch (modeidx)
        {
        case 1:
            robot->FootJacobian(ms[k].x, J, Jd, legID::HLEG);
            leg_ext_vec = robot->get_leg_ext_vec(q, legID::HLEG);
            F = -leg_ext_vec.normalized()*Kspring*(leg_ext_vec.norm() - legext_nom);
            ms[k].u = J.bottomRightCorner(2,4).transpose()*F*3;
            break;
        case 3:
            robot->FootJacobian(ms[k].x, J, Jd, legID::FLEG);
            leg_ext_vec = robot->get_leg_ext_vec(q, legID::FLEG);
            F = -leg_ext_vec.normalized()*Kspring*(leg_ext_vec.norm() - legext_nom);
            ms[k].u = J.bottomRightCorner(2,4).transpose()*F*2.2;
            break;
        case 2:
        case 4:
            ms[k].u = Kp*(qjoint_nom-ms[k].x.segment(3,4)) - Kd*ms[k].x.tail(4);
            break;
                           
        }        
        robot->dynamics(ms[k].x, ms[k].u, ms[k+1].x, ms[k].y, modeidx);      
    }    
}

template void bounding_PDcontrol(PlanarQuadruped<double>* robot, ModelState<double,xsize_WB,usize_WB,ysize_WB> *ms, int modeidx, size_t N_TIME_STEPS);

