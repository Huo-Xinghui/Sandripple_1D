/*reflect boundary conditions for DPM particles*/
/*Edited by Xinghui Huo*/
#include "udf.h"
DEFINE_DPM_BC(reflect_ang_bc,p,t,f,f_normal,dim)
{
    int i;
    real v_mag;
    real alpha;
    real en; /*normal restitution coefficient*/
    real et; /*tangential restitution coefficient*/
    real vn; /*normal velocity of the impact particle*/
    real vt; /*tangential velocity of the impact particle*/
    /*FILE *fp;*/

    v_mag = NV_MAG(P_VEL(p)); /*velocity magnitude of the impact particle*/
    alpha = M_PI/2.0 - acos(MAX(-1.0,MIN(1.0,NV_DOT(f_normal,P_VEL(p))/MAX(v_mag,DPM_SMALL)))); /*angle between the impact particle and the wall normal*/
    en = 0.993 - 1.76*alpha + 1.56*alpha*alpha - 0.49*alpha*alpha*alpha;
    et = 0.988 - 1.66*alpha + 2.11*alpha*alpha - 0.67*alpha*alpha*alpha;
    vn += NV_DOT(P_VEL(p),f_normal);
    /*tangential velocity of the impact particle*/
    for (i=0;i<dim;i++)
        P_VEL(p)[i] -= vn*f_normal[i];
    /*after rebound*/
    for (i=0;i<dim;i++)
        P_VEL(p)[i] *= et;
    /*total velocity*/
    for (i=0;i<dim;i++)
        P_VEL(p)[i] -= en*vn*f_normal[i];
    /*store new velocity in P_VEL0 of the particle*/
    NV_V(P_VEL0(p), =, P_VEL(p));
    /*out put*/
    /*if (myid ==0)
    {
        fp = fopen("impact_angle_blade.txt","a");
        fprintf(fp,"impact angle = %f\n",alpha*180.0/M_PI);
        fclose(fp);
        fp = fopen("impact_velocity_blade.txt","a");
        fprintf(fp,"impact velocity = %f\n",v_mag);
        fclose(fp);
    }*/
    return PATH_ACTIVE;
}