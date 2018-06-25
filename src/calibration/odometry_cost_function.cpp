#include "calibration/odometry_cost_function.h"

#include "projection/jacobian.h"
#include "geometry/geometry.h"

#include "eigen.h"
#include "io.h"



Vector3d odom_zeta_i (const Vector2d & delta, const double * const odom_intrin)
{

    double r1 = odom_intrin [ 0 ] ;
    double r2 = odom_intrin [ 1 ] ;
    double g  = odom_intrin [ 2 ] ;

    //-Converting wheel speed into platform speed and rotation
    Matrix2d joint_to_cartesian;
    joint_to_cartesian <<  ( r1 / 2 ) ,  ( r2 / 2 ) ,
            -( r1 / g ) ,  ( r2 / g ) ;

    Vector2d v_w = joint_to_cartesian * delta;

    //-Updating States
    Matrixd < 3 , 2 > A;
    A <<    1,   0,
            0,   0,
            0,   1;

    //-Updating zeta_i
    Vector3d zeta_i = A * v_w;

    return zeta_i;
}


Matrix3d zeta_i_jacobian(Vector2d & delta, const double * const odom_intrin)
{
    Matrix3d zeta_i_jac_mat;
    double r1 = odom_intrin [ 0 ] ;
    double r2 = odom_intrin [ 1 ] ;
    double g  = odom_intrin [ 2 ] ;

    double d_ql = delta (0);
    double d_qr = delta (1);

    double j11 = d_ql / 2;
    double j12 = d_qr / 2;
    double j13 = 0;

    double j21 = 0;
    double j22 = 0;
    double j23 = 0;

    double j31 = -d_ql / g;
    double j32 =  d_qr / g;
    double j33 = (r1 * d_ql - r2 * d_qr) / (g * g);

    //-Updating zeta_i_jac_mat
    zeta_i_jac_mat <<   j11,	j12,	j13,
            j21,	j22,	j23,
            j31,	j32,	j33;

    return zeta_i_jac_mat;
}


void OdometryCost::tf0n_jac_calc (const double * const odom_intrin,
                                  vector <Matrix3d> & jac_zeta_vec, 
                                  vector <Transf> & tf0_n_vec ) const
{
    Transformation <double> tf_0_i = Transformation <double> (0, 0, 0, 0, 0, 0);        // Initially 0T0
    for(int i = 0; i < _deltaQ.size(); i++)
    {
        Vector2d delta = _deltaQ[i];
        Vector3d zeta_i = odom_zeta_i (delta, odom_intrin);
        Matrix3d zeta_i_jac = zeta_i_jacobian (delta, odom_intrin);

        double delta_x = zeta_i(0);
        double delta_y = zeta_i(1);
        double delta_theta = zeta_i(2);

        Transformation <double> tf_j_i = Transformation <double>(delta_x, delta_y, 0, 0, 0, delta_theta);       // i-1Ti, here j = i-1
        tf_0_i = tf_0_i.compose(tf_j_i);                    // Updating 0Ti at each iteration

        jac_zeta_vec.push_back(zeta_i_jac);
        tf0_n_vec.push_back(tf_0_i);

    }
}

Matrixd<6, 3> OdometryCost::calc_acc(const vector<Transf> & tf0_n_vec, const vector<Matrix3d> & jac_zeta) const
{
    Matrix3d ACC;
    ACC <<  0,  0,  0,
            0,  0,  0,
            0,  0,  0;


    Transformation <double> tf0_n = tf0_n_vec.back();



    for(int i = 0; i < jac_zeta.size(); i++)
    {
        Transformation <double> tf0_j;
        if(i>0)        {
            tf0_j = tf0_n_vec[i-1];                         // 0Ti-1
        }
        else        {
            tf0_j = Transformation <double>(0, 0, 0, 0, 0, 0);    //0T0
        }
        Matrix3<double> rot_0_j_mat = tf0_j.rotMat();             // Rotation Matrix 0Ri-1

        Transformation <double> tf0_i = tf0_n_vec[i];                             // 0Ti
        Transformation <double> tfi_0 = tf0_i.inverse();                         // iT0
        Transformation <double> tfi_n = tfi_0.compose(tf0_n);                   // iTn

        Vector3<double> trans_i_n = tfi_n.trans();


        Matrix3d Jacobian;
        Jacobian << 1,  0,  -trans_i_n(1),
                0,  1,   trans_i_n(0),
                0,  0,      1;

        ACC = ACC + rot_0_j_mat*Jacobian*jac_zeta[i];
    }


    Matrixd<6, 3> acc_jacobian;
    acc_jacobian << ACC(0, 0), ACC(0, 1), ACC(0, 2),
            ACC(1, 0), ACC(1, 1), ACC(1, 2),
            0,          0,          0,
            0,          0,          0,
            0,          0,          0,
            ACC(2, 0), ACC(2, 1), ACC(2, 2);
    return acc_jacobian;
}


OdometryCost::OdometryCost(const double errV, const double errW, const double lambda,
                           const vector <Vector2d> & deltaQ_vec, const double * const intrinsic_prior) :
    _deltaQ(deltaQ_vec),
    _A(Matrix6d::Zero())
{

    //initial motion guess
    vector<Matrix3d> jac_zeta_vec;
    vector<Transf>  tf0_n_vec;
    tf0n_jac_calc (intrinsic_prior, jac_zeta_vec, tf0_n_vec);
    _zetaPrior = tf0_n_vec.back();
    
    const double MIN_SIGMA_V = 0.01;
    const double MIN_SIGMA_W = 0.01;

    const double MIN_DELTA = 0.01;
    const double MIN_L = 0.01;

    const double delta = max(_zetaPrior.rot().norm(), MIN_DELTA);
    const double l = max(_zetaPrior.trans().norm(), MIN_L);

    const double delta2 = delta / 2.;
    const double l2 = l / 2.;

    const double s = sin(delta2);
    const double c = cos(delta2);


    Matrixd<3, 2> dfdu;
    dfdu << c,      l2 * s,
            -s,      l2 * c,
            0,           1;

    Matrix2d Cu;
    Cu <<   errV * errV * l * l,              0,
            0,    errW * errW * delta * delta;
    Cu(0, 0) = max(Cu(0, 0), MIN_SIGMA_V*MIN_SIGMA_V);
    Cu(1, 1) = max(Cu(1, 1), MIN_SIGMA_W*MIN_SIGMA_W);

    Matrix3d Cx = dfdu * Cu * dfdu.transpose() + (lambda * lambda) * Matrix3d::Identity();
    Matrix3d CxInv = Cx.inverse();
    Eigen::LLT<Matrix3d> lltOfCxInv(CxInv);
    Matrix3d U = lltOfCxInv.matrixU();
    _A.topLeftCorner<2, 2>() = U.topLeftCorner<2, 2>();
    _A.topRightCorner<2, 1>() = U.topRightCorner<2, 1>();
    _A(2, 2) = 1. / lambda;

    Matrix3d B = Matrix3d::Zero();
    B(0, 0) = 1. / lambda;
    B(1, 1) = 1. / lambda;
    B(2, 2) = U(2, 2);

    _A.bottomRightCorner<3, 3>() = B /* interOmegaRot(_zetaPrior.rot())*/; //FIXME

    //    _A = Matrix6d::Identity();
}



bool OdometryCost::Evaluate(double const * const * params,
                            double * residual, double ** jacobian) const
{
    Transf xi1(params[0]);
    Transf xi2(params[1]);
    const double * const  odom_intrin = params[2];

    Transf zeta = xi1.inverseCompose(xi2);

    Map<Vector6d> res(residual);
    Vector6d err;


    vector<Transf>  zeta_odo_vec;
    vector <Matrix3d> jac_zeta_vec;

    tf0n_jac_calc (odom_intrin, jac_zeta_vec, zeta_odo_vec);
    Matrixd<6, 3> jac_intrinsic = calc_acc(zeta_odo_vec, jac_zeta_vec);

    Transf  zeta_odo;
    zeta_odo = zeta_odo_vec.back();
    // Make zeta_odo transformation 0Tn
    Transf delta = zeta_odo.inverseCompose(zeta);
    delta.toArray(err.data());
    res = _A * err;

    if (jacobian != NULL)
    {
        //first transformation
        if (jacobian[0] != NULL)
        {
            Matrix3d R10 = xi1.rotMatInv();
            Matrix3d R10_Mwr1 = R10 * interOmegaRot(xi1.rot());
            Matrix6d J1;
            J1 << R10, Matrix3d::Zero(), Matrix3d::Zero(), R10_Mwr1;
            Matrix6d TT = zeta.screwTransfInv();
            Map<Matrix6drm> jac(jacobian[0]);
            jac = -_A * TT * J1;
        }

        //second transformation
        if (jacobian[1] != NULL)
        {
            Matrix3d R20 = xi2.rotMatInv();
            Matrix3d R20_Mwr2 = R20 * interOmegaRot(xi2.rot());
            Matrix6d J2;
            J2 << R20, Matrix3d::Zero(), Matrix3d::Zero(), R20_Mwr2;
            Map<Matrix6drm> jac(jacobian[1]);
            jac = _A * J2;
        }

        //Odometry transformation
        if (jacobian[2] != NULL)
        {
            Matrix3d R31 = zeta_odo.rotMatInv();
            Matrix3d R31_Mwr3 = R31 * interOmegaRot(zeta_odo.rot());
            Matrix6d J3;
            J3 << R31, Matrix3d::Zero(), Matrix3d::Zero(), R31_Mwr3;
            Matrix6d TT = delta.screwTransfInv();
            Map<Matrix <double, 6, 3, RowMajor>> jac(jacobian[2]);
            jac = -_A * TT * J3 * jac_intrinsic;
        }
    }
    return true;
}

