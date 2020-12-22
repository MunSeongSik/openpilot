
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8238155516381470371) {
   out_8238155516381470371[0] = delta_x[0] + nom_x[0];
   out_8238155516381470371[1] = delta_x[1] + nom_x[1];
   out_8238155516381470371[2] = delta_x[2] + nom_x[2];
   out_8238155516381470371[3] = delta_x[3] + nom_x[3];
   out_8238155516381470371[4] = delta_x[4] + nom_x[4];
   out_8238155516381470371[5] = delta_x[5] + nom_x[5];
   out_8238155516381470371[6] = delta_x[6] + nom_x[6];
   out_8238155516381470371[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4318322549209788406) {
   out_4318322549209788406[0] = -nom_x[0] + true_x[0];
   out_4318322549209788406[1] = -nom_x[1] + true_x[1];
   out_4318322549209788406[2] = -nom_x[2] + true_x[2];
   out_4318322549209788406[3] = -nom_x[3] + true_x[3];
   out_4318322549209788406[4] = -nom_x[4] + true_x[4];
   out_4318322549209788406[5] = -nom_x[5] + true_x[5];
   out_4318322549209788406[6] = -nom_x[6] + true_x[6];
   out_4318322549209788406[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_1866736144955826514) {
   out_1866736144955826514[0] = 1.0;
   out_1866736144955826514[1] = 0.0;
   out_1866736144955826514[2] = 0.0;
   out_1866736144955826514[3] = 0.0;
   out_1866736144955826514[4] = 0.0;
   out_1866736144955826514[5] = 0.0;
   out_1866736144955826514[6] = 0.0;
   out_1866736144955826514[7] = 0.0;
   out_1866736144955826514[8] = 0.0;
   out_1866736144955826514[9] = 1.0;
   out_1866736144955826514[10] = 0.0;
   out_1866736144955826514[11] = 0.0;
   out_1866736144955826514[12] = 0.0;
   out_1866736144955826514[13] = 0.0;
   out_1866736144955826514[14] = 0.0;
   out_1866736144955826514[15] = 0.0;
   out_1866736144955826514[16] = 0.0;
   out_1866736144955826514[17] = 0.0;
   out_1866736144955826514[18] = 1.0;
   out_1866736144955826514[19] = 0.0;
   out_1866736144955826514[20] = 0.0;
   out_1866736144955826514[21] = 0.0;
   out_1866736144955826514[22] = 0.0;
   out_1866736144955826514[23] = 0.0;
   out_1866736144955826514[24] = 0.0;
   out_1866736144955826514[25] = 0.0;
   out_1866736144955826514[26] = 0.0;
   out_1866736144955826514[27] = 1.0;
   out_1866736144955826514[28] = 0.0;
   out_1866736144955826514[29] = 0.0;
   out_1866736144955826514[30] = 0.0;
   out_1866736144955826514[31] = 0.0;
   out_1866736144955826514[32] = 0.0;
   out_1866736144955826514[33] = 0.0;
   out_1866736144955826514[34] = 0.0;
   out_1866736144955826514[35] = 0.0;
   out_1866736144955826514[36] = 1.0;
   out_1866736144955826514[37] = 0.0;
   out_1866736144955826514[38] = 0.0;
   out_1866736144955826514[39] = 0.0;
   out_1866736144955826514[40] = 0.0;
   out_1866736144955826514[41] = 0.0;
   out_1866736144955826514[42] = 0.0;
   out_1866736144955826514[43] = 0.0;
   out_1866736144955826514[44] = 0.0;
   out_1866736144955826514[45] = 1.0;
   out_1866736144955826514[46] = 0.0;
   out_1866736144955826514[47] = 0.0;
   out_1866736144955826514[48] = 0.0;
   out_1866736144955826514[49] = 0.0;
   out_1866736144955826514[50] = 0.0;
   out_1866736144955826514[51] = 0.0;
   out_1866736144955826514[52] = 0.0;
   out_1866736144955826514[53] = 0.0;
   out_1866736144955826514[54] = 1.0;
   out_1866736144955826514[55] = 0.0;
   out_1866736144955826514[56] = 0.0;
   out_1866736144955826514[57] = 0.0;
   out_1866736144955826514[58] = 0.0;
   out_1866736144955826514[59] = 0.0;
   out_1866736144955826514[60] = 0.0;
   out_1866736144955826514[61] = 0.0;
   out_1866736144955826514[62] = 0.0;
   out_1866736144955826514[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6537832588653464867) {
   out_6537832588653464867[0] = state[0];
   out_6537832588653464867[1] = state[1];
   out_6537832588653464867[2] = state[2];
   out_6537832588653464867[3] = state[3];
   out_6537832588653464867[4] = state[4];
   out_6537832588653464867[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6537832588653464867[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6537832588653464867[7] = state[7];
}
void F_fun(double *state, double dt, double *out_9172075340303248199) {
   out_9172075340303248199[0] = 1;
   out_9172075340303248199[1] = 0;
   out_9172075340303248199[2] = 0;
   out_9172075340303248199[3] = 0;
   out_9172075340303248199[4] = 0;
   out_9172075340303248199[5] = 0;
   out_9172075340303248199[6] = 0;
   out_9172075340303248199[7] = 0;
   out_9172075340303248199[8] = 0;
   out_9172075340303248199[9] = 1;
   out_9172075340303248199[10] = 0;
   out_9172075340303248199[11] = 0;
   out_9172075340303248199[12] = 0;
   out_9172075340303248199[13] = 0;
   out_9172075340303248199[14] = 0;
   out_9172075340303248199[15] = 0;
   out_9172075340303248199[16] = 0;
   out_9172075340303248199[17] = 0;
   out_9172075340303248199[18] = 1;
   out_9172075340303248199[19] = 0;
   out_9172075340303248199[20] = 0;
   out_9172075340303248199[21] = 0;
   out_9172075340303248199[22] = 0;
   out_9172075340303248199[23] = 0;
   out_9172075340303248199[24] = 0;
   out_9172075340303248199[25] = 0;
   out_9172075340303248199[26] = 0;
   out_9172075340303248199[27] = 1;
   out_9172075340303248199[28] = 0;
   out_9172075340303248199[29] = 0;
   out_9172075340303248199[30] = 0;
   out_9172075340303248199[31] = 0;
   out_9172075340303248199[32] = 0;
   out_9172075340303248199[33] = 0;
   out_9172075340303248199[34] = 0;
   out_9172075340303248199[35] = 0;
   out_9172075340303248199[36] = 1;
   out_9172075340303248199[37] = 0;
   out_9172075340303248199[38] = 0;
   out_9172075340303248199[39] = 0;
   out_9172075340303248199[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_9172075340303248199[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_9172075340303248199[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_9172075340303248199[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_9172075340303248199[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_9172075340303248199[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_9172075340303248199[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_9172075340303248199[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_9172075340303248199[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_9172075340303248199[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_9172075340303248199[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9172075340303248199[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9172075340303248199[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_9172075340303248199[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_9172075340303248199[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_9172075340303248199[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9172075340303248199[56] = 0;
   out_9172075340303248199[57] = 0;
   out_9172075340303248199[58] = 0;
   out_9172075340303248199[59] = 0;
   out_9172075340303248199[60] = 0;
   out_9172075340303248199[61] = 0;
   out_9172075340303248199[62] = 0;
   out_9172075340303248199[63] = 1;
}
void h_25(double *state, double *unused, double *out_2761919617857195544) {
   out_2761919617857195544[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1633799369625853511) {
   out_1633799369625853511[0] = 0;
   out_1633799369625853511[1] = 0;
   out_1633799369625853511[2] = 0;
   out_1633799369625853511[3] = 0;
   out_1633799369625853511[4] = 0;
   out_1633799369625853511[5] = 0;
   out_1633799369625853511[6] = 1;
   out_1633799369625853511[7] = 0;
}
void h_24(double *state, double *unused, double *out_8975958342254413471) {
   out_8975958342254413471[0] = state[4];
   out_8975958342254413471[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5992468293127456264) {
   out_5992468293127456264[0] = 0;
   out_5992468293127456264[1] = 0;
   out_5992468293127456264[2] = 0;
   out_5992468293127456264[3] = 0;
   out_5992468293127456264[4] = 1;
   out_5992468293127456264[5] = 0;
   out_5992468293127456264[6] = 0;
   out_5992468293127456264[7] = 0;
   out_5992468293127456264[8] = 0;
   out_5992468293127456264[9] = 0;
   out_5992468293127456264[10] = 0;
   out_5992468293127456264[11] = 0;
   out_5992468293127456264[12] = 0;
   out_5992468293127456264[13] = 1;
   out_5992468293127456264[14] = 0;
   out_5992468293127456264[15] = 0;
}
void h_30(double *state, double *unused, double *out_2604732949506066399) {
   out_2604732949506066399[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7980099541323329769) {
   out_7980099541323329769[0] = 0;
   out_7980099541323329769[1] = 0;
   out_7980099541323329769[2] = 0;
   out_7980099541323329769[3] = 0;
   out_7980099541323329769[4] = 1;
   out_7980099541323329769[5] = 0;
   out_7980099541323329769[6] = 0;
   out_7980099541323329769[7] = 0;
}
void h_26(double *state, double *unused, double *out_6860041952196644452) {
   out_6860041952196644452[0] = state[7];
}
void H_26(double *state, double *unused, double *out_292363235930205959) {
   out_292363235930205959[0] = 0;
   out_292363235930205959[1] = 0;
   out_292363235930205959[2] = 0;
   out_292363235930205959[3] = 0;
   out_292363235930205959[4] = 0;
   out_292363235930205959[5] = 0;
   out_292363235930205959[6] = 0;
   out_292363235930205959[7] = 1;
}
void h_27(double *state, double *unused, double *out_7763592885518041189) {
   out_7763592885518041189[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2133033255914739710) {
   out_2133033255914739710[0] = 0;
   out_2133033255914739710[1] = 0;
   out_2133033255914739710[2] = 0;
   out_2133033255914739710[3] = 1;
   out_2133033255914739710[4] = 0;
   out_2133033255914739710[5] = 0;
   out_2133033255914739710[6] = 0;
   out_2133033255914739710[7] = 0;
}
void h_29(double *state, double *unused, double *out_7002685030787695411) {
   out_7002685030787695411[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6642845763257389894) {
   out_6642845763257389894[0] = 0;
   out_6642845763257389894[1] = 1;
   out_6642845763257389894[2] = 0;
   out_6642845763257389894[3] = 0;
   out_6642845763257389894[4] = 0;
   out_6642845763257389894[5] = 0;
   out_6642845763257389894[6] = 0;
   out_6642845763257389894[7] = 0;
}
void h_28(double *state, double *unused, double *out_7976698941224527876) {
   out_7976698941224527876[0] = state[5];
   out_7976698941224527876[1] = state[6];
}
void H_28(double *state, double *unused, double *out_5436443302345600486) {
   out_5436443302345600486[0] = 0;
   out_5436443302345600486[1] = 0;
   out_5436443302345600486[2] = 0;
   out_5436443302345600486[3] = 0;
   out_5436443302345600486[4] = 0;
   out_5436443302345600486[5] = 1;
   out_5436443302345600486[6] = 0;
   out_5436443302345600486[7] = 0;
   out_5436443302345600486[8] = 0;
   out_5436443302345600486[9] = 0;
   out_5436443302345600486[10] = 0;
   out_5436443302345600486[11] = 0;
   out_5436443302345600486[12] = 0;
   out_5436443302345600486[13] = 0;
   out_5436443302345600486[14] = 1;
   out_5436443302345600486[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
