/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8238155516381470371);
void inv_err_fun(double *nom_x, double *true_x, double *out_4318322549209788406);
void H_mod_fun(double *state, double *out_1866736144955826514);
void f_fun(double *state, double dt, double *out_6537832588653464867);
void F_fun(double *state, double dt, double *out_9172075340303248199);
void h_25(double *state, double *unused, double *out_2761919617857195544);
void H_25(double *state, double *unused, double *out_1633799369625853511);
void h_24(double *state, double *unused, double *out_8975958342254413471);
void H_24(double *state, double *unused, double *out_5992468293127456264);
void h_30(double *state, double *unused, double *out_2604732949506066399);
void H_30(double *state, double *unused, double *out_7980099541323329769);
void h_26(double *state, double *unused, double *out_6860041952196644452);
void H_26(double *state, double *unused, double *out_292363235930205959);
void h_27(double *state, double *unused, double *out_7763592885518041189);
void H_27(double *state, double *unused, double *out_2133033255914739710);
void h_29(double *state, double *unused, double *out_7002685030787695411);
void H_29(double *state, double *unused, double *out_6642845763257389894);
void h_28(double *state, double *unused, double *out_7976698941224527876);
void H_28(double *state, double *unused, double *out_5436443302345600486);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
