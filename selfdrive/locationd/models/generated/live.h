/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7972485626059406578);
void inv_err_fun(double *nom_x, double *true_x, double *out_6282590760582799899);
void H_mod_fun(double *state, double *out_8806751877630439200);
void f_fun(double *state, double dt, double *out_1928131735915360732);
void F_fun(double *state, double dt, double *out_1553718604023551924);
void h_3(double *state, double *unused, double *out_3777104109957885046);
void H_3(double *state, double *unused, double *out_3573155468537353559);
void h_4(double *state, double *unused, double *out_8209757366902474422);
void H_4(double *state, double *unused, double *out_2464317652654592224);
void h_9(double *state, double *unused, double *out_3032326798753568955);
void H_9(double *state, double *unused, double *out_6904217179851120571);
void h_10(double *state, double *unused, double *out_8611182023020486556);
void H_10(double *state, double *unused, double *out_8352008416246006182);
void h_12(double *state, double *unused, double *out_1998294514485401513);
void H_12(double *state, double *unused, double *out_6688267147719474170);
void h_31(double *state, double *unused, double *out_1498544089605041803);
void H_31(double *state, double *unused, double *out_3032851916058802476);
void h_32(double *state, double *unused, double *out_7683168867542331438);
void H_32(double *state, double *unused, double *out_8843819672581157880);
void h_13(double *state, double *unused, double *out_5819649201521763974);
void H_13(double *state, double *unused, double *out_1530529908191257688);
void h_14(double *state, double *unused, double *out_3032326798753568955);
void H_14(double *state, double *unused, double *out_6904217179851120571);
void h_19(double *state, double *unused, double *out_67499669903076517);
void H_19(double *state, double *unused, double *out_7109847042755999864);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);