
#include <sstream>
#include <cstdio>
#include <vector>
#include <cmath>
#include "checkit.h"
#include "doubledouble.h"

using namespace doubledouble;


#define assert_isnan(test, d) \
    assert_true(test, std::isnan(d.upper), "upper is NAN"); \
    assert_true(test, std::isnan(d.lower), "lower is NAN"); \


struct unary_case {
    const char* label;
    double xhi, xlo, yhi, ylo;
    double reltol;
};

struct atan2_case {
    const char* label;
    double yhi, ylo, xhi, xlo, zhi, zlo;
    double reltol;
};

void assert_dd_close(CheckIt& test, DoubleDouble actual, double expected_hi,
                     double expected_lo, double reltol,
                     const std::string& label)
{
    if (std::fabs(expected_hi) < 1e-90) {
        assert_true(test, std::fabs(actual.upper - expected_hi) < 1e-90, label + " (upper near zero)");
    }
    else {
        assert_equal_fp(test, actual.upper, expected_hi, label + " (upper)");
    }
    if (std::fabs(expected_lo) < 1e-30) {
        assert_true(test, std::fabs(actual.lower - expected_lo) < 1e-30, label + " (lower near zero)");
    }
    else {
        assert_close_fp(test, actual.lower, expected_lo, reltol, label + " (lower)");
    }
}

static const unary_case sin_cases[] = {
    {"0", 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-20+3e-38", 9.99999999999999945e-21, 5.78467285457904316e-37, 9.99999999999999945e-21, 5.78467285457904316e-37, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133522040777182e-29, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.00000000000000002e-08, -3.75892274967951411e-25, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 9.99999998333333298e-05, 3.66021070168176944e-21, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.00000000000000000e+00, 1.24674733385227693e-01, -2.92594749605785818e-18, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.00000000000000000e+00, 4.79425538604203005e-01, -5.10396986055601290e-18, 5e-14},
    {"1", 1.00000000000000000e+00, 0.00000000000000000e+00, 8.41470984807896505e-01, 1.77684509293553611e-18, 5e-14},
    {"pi/24", 1.30899693899574732e-01, -1.34010220806386373e-17, 1.30526192220051601e-01, -9.60766602390152377e-18, 5e-14},
    {"pi/6", 5.23598775598298927e-01, -5.36040883225545492e-17, 5.00000000000000000e-01, 0.00000000000000000e+00, 5e-14},
    {"pi/4", 7.85398163397448279e-01, 3.06161699786838302e-17, 7.07106781186547573e-01, -4.83364665672645673e-17, 5e-14},
    {"pi/3", 1.04719755119659785e+00, -1.07208176645109098e-16, 8.66025403784438597e-01, 5.01754211090345140e-17, 5e-14},
    {"pi/2", 1.57079632679489656e+00, 6.12323399573676604e-17, 1.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"pi", 3.14159265358979312e+00, 1.22464679914735321e-16, 2.86889280175963737e-102, 6.66716889325310485e-119, 5e-14},
    {"10", 1.00000000000000000e+01, 0.00000000000000000e+00, -5.44021110889369774e-01, -3.89498986682235567e-17, 5e-14},
    {"1e6+1e-10", 1.00000000000000012e+06, -1.64153218269348140e-11, -3.49993502077617757e-01, 1.80771655912666880e-17, 5e-10},
};

static const unary_case cos_cases[] = {
    {"0", 0.00000000000000000e+00, 0.00000000000000000e+00, 1.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 1.00000000000000000e+00, -4.99999999999999981e-81, 5e-14},
    {"1e-20+3e-38", 9.99999999999999945e-21, 5.78467285457904316e-37, 1.00000000000000000e+00, -4.99999999999999965e-41, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 1.00000000000000000e+00, -4.99999999999999962e-25, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.00000000000000000e+00, -4.99999999999999990e-17, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 9.99999995000000030e-01, -2.62206881893329720e-17, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.00000000000000000e+00, 9.92197667229329006e-01, 4.75487057518936408e-17, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.00000000000000000e+00, 8.77582561890372759e-01, -4.26231498642799968e-17, 5e-14},
    {"1", 1.00000000000000000e+00, 0.00000000000000000e+00, 5.40302305868139765e-01, -4.76095461260441722e-17, 5e-14},
    {"pi/24", 1.30899693899574732e-01, -1.34010220806386373e-17, 9.91444861373810382e-01, 2.89898987445453692e-17, 5e-14},
    {"pi/6", 5.23598775598298927e-01, -5.36040883225545492e-17, 8.66025403784438597e-01, 5.01754211090345140e-17, 5e-14},
    {"pi/4", 7.85398163397448279e-01, 3.06161699786838302e-17, 7.07106781186547573e-01, -4.83364665672645673e-17, 5e-14},
    {"pi/3", 1.04719755119659785e+00, -1.07208176645109098e-16, 5.00000000000000000e-01, -3.57183559775710932e-102, 5e-14},
    {"pi/2", 1.57079632679489656e+00, 6.12323399573676604e-17, 1.43444640087981869e-102, 3.33358444662655242e-119, 5e-14},
    {"pi", 3.14159265358979312e+00, 1.22464679914735321e-16, -1.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"10", 1.00000000000000000e+01, 0.00000000000000000e+00, -8.39071529076452438e-01, -1.41471199889534177e-17, 5e-14},
    {"1e6+1e-10", 1.00000000000000012e+06, -1.64153218269348140e-11, 9.36752127568144188e-01, -5.12010725356683857e-17, 5e-10},
};

static const unary_case cosm1_cases[] = {
    {"0", 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, -4.99999999999999981e-81, -1.92871978607688389e-97, 5e-14},
    {"1e-20+3e-38", 9.99999999999999945e-21, 5.78467285457904316e-37, -4.99999999999999965e-41, -3.83535603000599262e-57, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, -4.99999999999999962e-25, -3.81497501824819198e-41, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, -4.99999999999999990e-17, -6.28443997131603143e-34, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, -4.99999999583333337e-09, 3.29639839020549711e-26, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.00000000000000000e+00, -7.80233277067094669e-03, -1.56189837468557130e-19, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.00000000000000000e+00, -1.22417438109627283e-01, -9.89786440836627501e-19, 5e-14},
    {"1", 1.00000000000000000e+00, 0.00000000000000000e+00, -4.59697694131860291e-01, 7.90160510521365173e-18, 5e-14},
    {"pi/6", 5.23598775598298927e-01, -5.36040883225545492e-17, -1.33974596215561348e-01, -5.33573012222331378e-18, 5e-14},
    {"pi/4", 7.85398163397448279e-01, 3.06161699786838302e-17, -2.92893218813452483e-01, 7.17468466399326131e-18, 5e-14},
    {"pi/2", 1.57079632679489656e+00, 6.12323399573676604e-17, -1.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"pi", 3.14159265358979312e+00, 1.22464679914735321e-16, -2.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"10", 1.00000000000000000e+01, 0.00000000000000000e+00, -1.83907152907645255e+00, 9.68751824735622394e-17, 5e-14},
};

static const unary_case tan_cases[] = {
    {"0", 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133527040777193e-29, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.00000000000000002e-08, 1.24107725032048620e-25, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 1.00000000333333341e-04, -6.64441504909537898e-21, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.00000000000000000e+00, 1.25655136575130971e-01, -3.43996068035248496e-18, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.00000000000000000e+00, 5.46302489843790484e-01, 2.90965762168371759e-17, 5e-14},
    {"1", 1.00000000000000000e+00, 0.00000000000000000e+00, 1.55740772465490229e+00, -6.18646417603759204e-17, 5e-14},
    {"pi/24", 1.30899693899574732e-01, -1.34010220806386373e-17, 1.31652497587395861e-01, -7.91769915790245361e-18, 5e-14},
    {"pi/6", 5.23598775598298927e-01, -5.36040883225545492e-17, 5.77350269189625731e-01, 3.34502807393563447e-17, 5e-14},
    {"pi/4", 7.85398163397448279e-01, 3.06161699786838302e-17, 1.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"pi/3", 1.04719755119659785e+00, -1.07208176645109098e-16, 1.73205080756887719e+00, 1.00350842218069028e-16, 5e-14},
    {"10", 1.00000000000000000e+01, 0.00000000000000000e+00, 6.48360827459086630e-01, 4.07615160389350079e-17, 5e-14},
    {"pi/2 - 2^-20", 1.57079537312058015e+00, 6.12323399573676604e-17, 1.04857599999968207e+06, 3.88051072563702319e-11, 5e-14},
};

static const unary_case atan_cases[] = {
    {"0", 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133520374110511e-29, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.00000000000000000e+00, 1.24354994546761438e-01, -3.12532414245393831e-18, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.00000000000000000e+00, 4.63647609000806094e-01, 2.26987774529616871e-17, 5e-14},
    {"1", 1.00000000000000000e+00, 0.00000000000000000e+00, 7.85398163397448279e-01, 3.06161699786838302e-17, 5e-14},
    {"-1", -1.00000000000000000e+00, 0.00000000000000000e+00, -7.85398163397448279e-01, -3.06161699786838302e-17, 5e-14},
    {"10", 1.00000000000000000e+01, 0.00000000000000000e+00, 1.47112767430373470e+00, -1.08492227620614239e-16, 5e-14},
    {"-10", -1.00000000000000000e+01, 0.00000000000000000e+00, -1.47112767430373470e+00, 1.08492227620614239e-16, 5e-14},
    {"1e40", 1.00000000000000003e+40, -3.03786028427003667e+23, 1.57079632679489656e+00, 6.12323399573676604e-17, 5e-14},
};

static const unary_case asin_cases[] = {
    {"0", 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133525374110523e-29, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.00000000000000000e+00, 1.25327831168065396e-01, 1.29060104888106167e-18, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.00000000000000000e+00, 5.23598775598298927e-01, -5.36040883225545492e-17, 5e-14},
    {"sqrt(1/2)", 7.07106781186547573e-01, -4.83364665672645673e-17, 7.85398163397448279e-01, 3.06161699786838302e-17, 5e-14},
    {"1 - 2^-54", 1.00000000000000000e+00, -5.55111512312578270e-17, 1.57079631625818439e+00, 1.02948797653591831e-16, 5e-14},
    {"1", 1.00000000000000000e+00, 0.00000000000000000e+00, 1.57079632679489656e+00, 6.12323399573676604e-17, 5e-14},
    {"-1", -1.00000000000000000e+00, 0.00000000000000000e+00, -1.57079632679489656e+00, -6.12323399573676604e-17, 5e-14},
};

static const unary_case acos_cases[] = {
    {"0", 0.00000000000000000e+00, 0.00000000000000000e+00, 1.57079632679489656e+00, 6.12323399573676604e-17, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 1.57079632679489656e+00, 6.12323399573676604e-17, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 1.57079632679389669e+00, -7.19116826266520367e-17, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.00000000000000000e+00, 1.44546849562683133e+00, -1.06591714785286887e-16, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.00000000000000000e+00, 1.04719755119659785e+00, -1.07208176645109098e-16, 5e-14},
    {"sqrt(1/2)", 7.07106781186547573e-01, -4.83364665672645673e-17, 7.85398163397448279e-01, 3.06161699786838302e-17, 5e-14},
    {"1 - 2^-54", 1.00000000000000000e+00, -5.55111512312578270e-17, 1.05367121277235087e-08, -6.71527394826313338e-25, 5e-14},
    {"1", 1.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"-1", -1.00000000000000000e+00, 0.00000000000000000e+00, 3.14159265358979312e+00, 1.22464679914735321e-16, 5e-14},
};

static const atan2_case atan2_cases[] = {
    {"0,1", 0.00000000000000000e+00, 0.00000000000000000e+00, 1.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 5e-14},
    {"0,-1", 0.00000000000000000e+00, 0.00000000000000000e+00, -1.00000000000000000e+00, 0.00000000000000000e+00, 3.14159265358979312e+00, 1.22464679914735321e-16, 5e-14},
    {"1,0", 1.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 1.57079632679489656e+00, 6.12323399573676604e-17, 5e-14},
    {"-1,0", -1.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00, -1.57079632679489656e+00, -6.12323399573676604e-17, 5e-14},
    {"1,1", 1.00000000000000000e+00, 0.00000000000000000e+00, 1.00000000000000000e+00, 0.00000000000000000e+00, 7.85398163397448279e-01, 3.06161699786838302e-17, 5e-14},
    {"1,-1", 1.00000000000000000e+00, 0.00000000000000000e+00, -1.00000000000000000e+00, 0.00000000000000000e+00, 2.35619449019234484e+00, 9.18485099360514844e-17, 5e-14},
    {"-1,-1", -1.00000000000000000e+00, 0.00000000000000000e+00, -1.00000000000000000e+00, 0.00000000000000000e+00, -2.35619449019234484e+00, -9.18485099360514844e-17, 5e-14},
    {"-1,1", -1.00000000000000000e+00, 0.00000000000000000e+00, 1.00000000000000000e+00, 0.00000000000000000e+00, -7.85398163397448279e-01, -3.06161699786838302e-17, 5e-14},
    {"1e-40,1", 9.99999999999999929e-41, 7.07071206001198558e-57, 1.00000000000000000e+00, 0.00000000000000000e+00, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1,1e-40", 1.00000000000000000e+00, 0.00000000000000000e+00, 9.99999999999999929e-41, 7.07071206001198558e-57, 1.57079632679489656e+00, 6.12323399573676604e-17, 5e-14},
    {"3+1e-30,4-2e-30", 3.00000000000000000e+00, 1.00000000000000008e-30, 4.00000000000000000e+00, -2.00000000000000017e-30, 6.43501108793284371e-01, 1.58347850514446868e-17, 5e-14},
    {"-3,4", -3.00000000000000000e+00, 0.00000000000000000e+00, 4.00000000000000000e+00, 0.00000000000000000e+00, -6.43501108793284371e-01, -1.58347850514442862e-17, 5e-14},
    {"3,-4", 3.00000000000000000e+00, 0.00000000000000000e+00, -4.00000000000000000e+00, 0.00000000000000000e+00, 2.49809154479650886e+00, -4.39240759922462181e-18, 5e-14},
};


//
// Some of the equality tests of the lower parts of the results might be
// too optimistic.  They might fail on different combinations of platforms,
// compilers, and compiler options.
//

void test_constructor(CheckIt& test)
{
    DoubleDouble d;
    assert_equal_fp(test, d.upper, 0.0, "DoubleDouble() (upper)");
    assert_equal_fp(test, d.lower, 0.0, "DoubleDouble() (lower)");

    auto three = DoubleDouble(3.0);
    assert_equal_fp(test, three.upper, 3.0, "DoubleDouble(3.0) (upper)");
    assert_equal_fp(test, three.lower, 0.0, "DoubleDouble(3.0) (lower)");

    auto x = DoubleDouble(10.0, 3e-18);
    assert_equal_fp(test, x.upper, 10.0, "DoubleDouble(10.0, 3e-18) (upper)");
    assert_equal_fp(test, x.lower, 3e-18, "DoubleDouble(10.0, 3e-18)  (lower)");

    // Check canonicalization.
    auto z = DoubleDouble(-1.0, 1.25);
    assert_equal_fp(test, z.upper, 0.25, "DoubleDouble(-1.0, 1.25) (upper)");
    assert_equal_fp(test, z.lower, 0.0, "DoubleDouble(-1.0, 1.25) (lower)");

    auto y = DoubleDouble(INFINITY, 0.0);
    assert_equal_fp(test, y.upper, (double) INFINITY, "DoubleDouble(INF, 0) (upper)");
    assert_equal_fp(test, y.lower, 0.0, "DoubleDouble(INF, 0) (lower)");

    auto bad = DoubleDouble((double) INFINITY, (double) -INFINITY);
    assert_true(test, std::isnan(bad.upper), "DoubleDouble(INF, -INF) (upper)");
    assert_true(test, std::isnan(bad.lower), "DoubleDouble(INF, -INF) (lower)");

    auto dd1 = DoubleDouble(NAN, 0.0);
    assert_isnan(test, dd1);

    auto dd2 = DoubleDouble(1.0, NAN);
    assert_isnan(test, dd2);

    auto dd3 = DoubleDouble(NAN);
    assert_isnan(test, dd3);
}

void test_add(CheckIt& test)
{
    auto z = DoubleDouble(1.0, 1e-18) + DoubleDouble(1.5, 5e-19);
    assert_equal_fp(test, z.upper, 2.5, "(1, 1e-18) + (1.5, 5e-19) (upper)");
    assert_close_fp(test, z.lower, 1.5e-18, 5e-16, "(1, 1e-18) + (1.5, 5e-19) (lower)");

    auto s1 = DoubleDouble(4.0, 4e-17) + NAN;
    assert_isnan(test, s1);

    auto s2 = DoubleDouble(4.0, 4e-17) + DoubleDouble(NAN);
    assert_isnan(test, s2);
}

void test_inplace_add(CheckIt& test)
{
    auto z = DoubleDouble(2.5, 1.5e-18);
    z += 3.5;
    assert_equal_fp(test, z.upper, 6.0, "Check z += 3.5 (upper)");
    assert_close_fp(test, z.lower, 1.5e-18, 5e-16, "Check z += 3.5 (lower)");

    z += DoubleDouble(3.0, 1e-18);
    assert_equal_fp(test, z.upper, 9.0, "Check z += (3.5, 1e-18) (upper)");
    assert_close_fp(test, z.lower, 2.5e-18, 5e-16, "Check z += (3.5, 1e-18) (lower)");

    z += NAN;
    assert_isnan(test, z);

    z = DoubleDouble(100.0, 1.3e-15);
    z += DoubleDouble(54.0, NAN);
    assert_isnan(test, z);
}

void test_subtract(CheckIt& test)
{
    DoubleDouble y;
    auto three = DoubleDouble(3.0);
    auto x = DoubleDouble(10.0, 3e-18);

    y = x - three;
    assert_equal_fp(test, y.upper, 7.0, "x - 3 (upper)");
    assert_equal_fp(test, y.lower, 3e-18, "x - 3 (lower)");

    y = three - x;
    assert_equal_fp(test, y.upper, -7.0, "3 - x (upper)");
    assert_equal_fp(test, y.lower, -3e-18, "3 - x (lower)");

    y = x - 10;
    assert_equal_fp(test, y.upper, 3e-18, "x - 10 (upper)");
    assert_equal_fp(test, y.lower, 0.0, "x - 10 (lower)");

    y = x - NAN;
    assert_isnan(test, y);

    y = x - DoubleDouble(5.0, NAN);
    assert_isnan(test, y);
}

void test_inplace_subtract(CheckIt& test)
{
    DoubleDouble y;

    y = DoubleDouble(2.5);
    y -= 4.5e-18;
    assert_equal_fp(test, y.upper, 2.5, "test_inplace_subtract (upper)");
    assert_equal_fp(test, y.lower, -4.5e-18, "test_inplace_subtract (lower)");

    y -= 1.5;
    assert_equal_fp(test, y.upper, 1.0, "test_inplace_subtract (upper)");
    assert_equal_fp(test, y.lower, -4.5e-18, "test_inplace_subtract (lower)");

    y -= 1.0;
    assert_equal_fp(test, y.upper, -4.5e-18, "test_inplace_subtract (upper)");
    assert_equal_fp(test, y.lower, 0.0, "test_inplace_subtract (lower)");

    /*
    // See https://bugs.llvm.org/show_bug.cgi?id=43124
    y = DoubleDouble(17.5, -2.25e-21);
    y -= y;
    assert_equal_fp(test, y.upper, 0.0, "test_inplace_subtract (upper)");
    assert_equal_fp(test, y.lower, 0.0, "test_inplace_subtract (lower)");
    */

    y -= NAN;
    assert_isnan(test, y);

    y = DoubleDouble(100.0, 1.3e-15);
    y -= DoubleDouble(54.0, NAN);
    assert_isnan(test, y);
}

void test_multiply(CheckIt& test)
{
    DoubleDouble y;
    auto x = DoubleDouble(10.0, 3e-18);

    y = x*x;
    assert_equal_fp(test, y.upper, 100.0, "Check x*x (upper)");
    assert_equal_fp(test, y.lower, 6e-17, "Check x*x (lower)");

    y = x*NAN;
    assert_isnan(test, y);

    y = x*DoubleDouble(NAN);
    assert_isnan(test, y);
}

void test_inplace_multiply(CheckIt& test)
{
    auto z = DoubleDouble(9.0, 2.5e-18);
    z *= -3.0;
    assert_equal_fp(test, z.upper, -27.0, "Check z *= -3.0 (upper)");
    assert_equal_fp(test, z.lower, -7.5e-18, "Check z *= -3.0 (lower)");

    z *= DoubleDouble(-2.0, 2e-18);
    assert_equal_fp(test, z.upper, 54.0, "Check z *= (-2.0, 2e-18) (upper)");
    assert_close_fp(test, z.lower, -27*2e-18 + -7.5e-18*-2.0, 5e-16, "Check z *= (-2.0, 2e-18) (lower)");

    z *= NAN;
    assert_isnan(test, z);

    z = DoubleDouble(3.0);
    z *= DoubleDouble(NAN);
    assert_isnan(test, z);
}

void test_divide(CheckIt &test)
{
    auto three = DoubleDouble(3.0);
    auto y = 1.0 / three;
    assert_equal_fp(test, y.upper, 0.3333333333333333, "1/3 (upper)");
    assert_equal_fp(test, y.lower, 1.850371707708594e-17, "1/3 (lower)");

    y = three/NAN;
    assert_isnan(test, y);

    y = three/DoubleDouble(NAN);
    assert_isnan(test, y);
}

void test_inplace_divide(CheckIt& test)
{
    DoubleDouble z;

    z = DoubleDouble(5.0, 4e-21);
    z /= 5.0;
    assert_equal_fp(test, z.upper, 1.0, "Check z /= 5.0 (upper)");
    assert_close_fp(test, z.lower, 8e-22, 5e-16, "Check z /= 5.0 (lower)");

    z = DoubleDouble(5.0, 4e-21);
    z /= DoubleDouble(5.0, 1e-20);
    assert_equal_fp(test, z.upper, 1.0, "Check z /= (5.0, 1e-20) (upper)");
    assert_close_fp(test, z.lower, -1.2e-21, 5e-16, "Check z /= (5.0, 1e-20) (lower)");

    z /= NAN;
    assert_isnan(test, z);

    z = DoubleDouble(3.0);
    z /= DoubleDouble(NAN);
    assert_isnan(test, z);
}

void test_expressions(CheckIt& test)
{
    DoubleDouble y;
    DoubleDouble one{1.0};
    DoubleDouble p6{3.14159, 0.0};

    y = one - DoubleDouble(2.0)/3;
    assert_equal_fp(test, y.upper, 0.3333333333333333, "Check 1 - 2/3 (upper)");
    assert_close_fp(test, y.lower, 1.850371707708594e-17, 5e-16, "Check 1 - 2/3 (lower)");

    y = (dd_pi - p6)*((dd_pi + p6)/2);
    assert_equal_fp(test, y.upper, 8.336494679678182e-06, "Check (pi - p6)*((pi + p6)/2) (upper)");
    assert_close_fp(test, y.lower, -1.4573851200099609e-22, 5e-16, "Check (pi - p6)*((pi + p6)/2) (lower)");
}

void test_comparisons(CheckIt& test)
{
    DoubleDouble w, x, y, z;

    w = DoubleDouble(-3.0);
    x = DoubleDouble(13.5, 1e-20);
    y = DoubleDouble(13.5, -3e-18);
    z = DoubleDouble(99.0, -3e-18);

    assert_true(test, x == x, "x == x");
    assert_true(test, !(x == 2.0), "!(x == 2.0)");
    assert_true(test, !(x == y), "!(x == y)");
    assert_true(test, !(x == z), "!(x == z)");
    assert_true(test, w == -3.0, "w == -3.0");
    assert_true(test, -3.0 == w, "-3.0 == w");

    assert_true(test, x != y, "x != y");
    assert_true(test, x != z, "x != z");
    assert_true(test, x != 1.0, "x != 1.0");
    assert_true(test, 1.0 != x, "1.0 != x");

    assert_true(test, y < x, "y < x");
    assert_true(test, y < 14.0, "y < 14.0");
    assert_true(test, y < 13.5, "y < 13.5");
    assert_true(test, 13.5 < x, "13.5 < x");

    assert_true(test, y <= x, "y <= x");
    assert_true(test, y <= y, "y <= x");
    assert_true(test, y <= 14.0, "y <= 14.0");
    assert_true(test, y <= 13.5, "y <= 13.5");
    assert_true(test, 13.5 <= x, "13.5 <= x");
    assert_true(test, -3.0 <= w, "-3.0 <= w");

    assert_true(test, x > y, "x > y");
    assert_true(test, x > 12.0, "x > 12.0");
    assert_true(test, x > 13.5, "x > 13.5");
    assert_true(test, 100.0 > z, "100.0 > z");

    assert_true(test, x >= y, "x >= y");
    assert_true(test, x >= x, "x >= x");
    assert_true(test, x >= 12.0, "x >= 12.0");
    assert_true(test, x >= 13.5, "x >= 13.5");
    assert_true(test, 100.0 >= z, "100.0 >= z");
    assert_true(test, -3.0 >= w, "-3.0 >= w");

    assert_true(test, DoubleDouble(NAN) != DoubleDouble(NAN), "check NAN != NAN is true");
    assert_false(test, DoubleDouble(NAN) == DoubleDouble(NAN), "NAN == NAN");
    assert_false(test, DoubleDouble(NAN) < DoubleDouble(NAN), "NAN < NAN");
    assert_false(test, DoubleDouble(NAN) > DoubleDouble(NAN), "NAN > NAN");
    assert_false(test, DoubleDouble(NAN) <= DoubleDouble(NAN), "NAN <= NAN");
    assert_false(test, DoubleDouble(NAN) >= DoubleDouble(NAN), "NAN >= NAN");
}


void test_abs(CheckIt& test)
{
    DoubleDouble a;
    auto x = DoubleDouble(-7.0, -3e-18);

    a = x.abs();
    assert_equal_fp(test, a.upper, 7.0, "abs(x) (upper)");
    assert_equal_fp(test, a.lower, 3e-18, "abs(x) (lower)");

    a = DoubleDouble(NAN).abs();
    assert_isnan(test, a);
}

void test_powi(CheckIt& test)
{
    auto x = DoubleDouble(10.0, 3e-18);
    auto y = x.powi(4);
    assert_equal_fp(test, y.upper, 10000.0, "Check x**4 (upper)");
    assert_close_fp(test, y.lower, 1.2e-14, 5e-16, "Check x**4 (lower)");

    y = DoubleDouble(NAN).powi(0);
    assert_true(test, y == 1.0, "NAN**0 is 1");
    y = DoubleDouble(NAN).powi(1);
    assert_isnan(test, y);
    y = DoubleDouble(NAN).powi(2);
    assert_isnan(test, y);
    y = DoubleDouble(NAN).powi(3);
    assert_isnan(test, y);
}

void test_sqrt(CheckIt& test)
{
    DoubleDouble y;
    DoubleDouble x{0.0026, 2e-20};

    y = x.sqrt();
    assert_equal_fp(test, y.upper, 0.05099019513592785, "sqrt((1.0026, 5e-17)) (upper)");
    assert_close_fp(test, y.lower, -8.798476828364278e-20, 5e-16, "sqrt((1.0026, 5e-17)) (lower)");

    y = DoubleDouble(101.0, 7.046e-15).sqrt();
    assert_equal_fp(test, y.upper, 10.04987562112089, "sqrt((101, 7.046e-15)) (upper)");
    assert_close_fp(test, y.lower, 6.749565974288905e-16, 5e-16, "sqrt((101, 7.046e-15)) (lower)");

    y = DoubleDouble(0.0, 0.0).sqrt();
    assert_equal_fp(test, y.upper, 0.0, "sqrt((0, 0)) (upper)");
    assert_equal_fp(test, y.lower, 0.0, "sqrt((0, 0)) (lower)");

    y = DoubleDouble(NAN).sqrt();
    assert_isnan(test, y);
}

void test_log(CheckIt& test)
{
    DoubleDouble y;
    auto x = DoubleDouble(10.0, 3e-18);

    y = x.log();
    assert_equal_fp(test, y.upper, 2.302585092994046, "Check log(x) (upper)");
    assert_equal_fp(test, y.lower, -2.1677562233822494e-16, "Check log(x) (lower)");

    y = dd_e.log();
    assert_equal_fp(test, y.upper, 1.0, "log(e) (upper)");
    assert_equal_fp(test, y.lower, 0.0, "log(e) (lower)");

    // The actual relative error for log(1.000000099005) is approx. 8.7e-10.
    // Can we do better?
    y = DoubleDouble(1.000000099005, 0.0).log();
    assert_equal_fp(test, y.upper, 9.900499507506536e-08, "log(1.000000099005) (upper)");
    assert_close_fp(test, y.lower, 4.563816054961034e-24, 2e-9, "log(1.000000099005) (lower)");

    y = DoubleDouble(1.1, 0.0).log();
    assert_equal_fp(test, y.upper, 0.09531017980432493, "log(1.1) (upper)");
    assert_close_fp(test, y.lower, 5.927240202146761e-18, 8e-16, "log(1.1) (lower)");

    y = DoubleDouble(1e-21, 3.5e-43).log();
    assert_equal_fp(test, y.upper, -48.35428695287496, "log((1e-21, 3.5e-43)) (upper)");
    assert_close_fp(test, y.lower, -1.7511230665702564e-15, 5e-16, "log((1e-21, 3.5e-43)) (lower)");

    y = DoubleDouble(NAN).log();
    assert_isnan(test, y);
}

struct log1p_case {
    // Acceptable relative error in the lower part of the result.
    double reltol;

    double xhi, xlo, yhi, ylo;
};

void test_log1p(CheckIt& test)
{
  struct log1p_case samples[] = {
    {5e-16, 3.1e-50, 7e-68, 3.1e-50, 7e-68},
    {5e-16, 1.4e-30, 5e-48, 1.4e-30, 4.9999999999990196e-48},

    {2e-15, 7.7e-16, 8.6e-34, 7.699999999999997e-16, 2.328394578795871e-34},
    {5e-16, 3e-8, 0.0, 2.9999999550000006e-08, 4.44602137450418e-26},
    {5e-16, 5e-7, 0.0, 4.999998750000417e-07, -1.9801357223218586e-23},
    {5e-16, 2e-6, 0.0, 1.9999980000026667e-06, -1.6536017667488005e-22},
    {4e-16, 2e-5, 0.0, 1.9999800002666627e-05, 1.3694019941860543e-21},
    {5e-16, 7e-4, 0.0, 0.0006997551142733419, 8.045441159182502e-21},
    {5e-16, 2.6e-3, 2e-20, 0.002596625847265978, 1.1029570690454358e-19},
    {5e-16, -2.6e-3, -3e-21, -0.002603385870114881, 1.2070558387143255e-20},

    {5e-16, -1.0, 2e-20, -45.35855467932097, -1.1302366343060957e-15},
    {5e-16, -0.5, -4.7e-18, -0.6931471805599453, -3.2590468138462995e-17},
    {5e-16, 1.5, 5e-19, 0.9162907318741551, -4.121195369011963e-17},
    {5e-16, 3.0, 0.0, 1.3862943611198906, 4.638093627692599e-17},
    };

    for (size_t i = 0; i < sizeof(samples)/sizeof(struct log1p_case); ++i) {
        struct log1p_case sample = samples[i];
        DoubleDouble x{sample.xhi, sample.xlo};
        auto y = x.log1p();
        std::stringstream s1, s2;
        s1 << "log1p case " << i << "  (upper)";
        s2 << "log1p case " << i << "  (lower)";
        assert_equal_fp(test, y.upper, sample.yhi, s1.str());
        assert_close_fp(test, y.lower, sample.ylo, sample.reltol, s2.str());
    }

    auto y = DoubleDouble(NAN).log1p();
    assert_isnan(test, y);
}

void test_exp(CheckIt& test)
{
    DoubleDouble y;
    auto x = DoubleDouble(10.0, 3e-18);
    auto t = DoubleDouble(3.7500950075008e-15, 0.0);

    y = DoubleDouble(0.0).exp();
    assert_equal_fp(test, y.upper, 1.0, "exp(0) (upper)");
    assert_equal_fp(test, y.lower, 0.0, "exp(0) (lower)");

    y = DoubleDouble(1.0).exp();
    assert_equal_fp(test, y.upper, 2.718281828459045, "exp(0) (upper)");
    assert_equal_fp(test, y.lower, 1.4456468917292502e-16, "exp(0) (lower)");

    y = x.exp();
    assert_equal_fp(test, y.upper, 22026.465794806718, "exp(x) (upper)");
    assert_equal_fp(test, y.lower, -1.3119340726673171e-12, "exp(x) (lower)");

    y = dd_pi.exp();
    assert_equal_fp(test, y.upper, 23.14069263277927, "exp(pi) (upper)");
    assert_equal_fp(test, y.lower, -1.3488747091995788e-15, "exp(pi) (lower)");

    y = (x - 10.0).exp();
    assert_equal_fp(test, y.upper, 1.0, "exp(x - 10) (upper)");
    assert_equal_fp(test, y.lower, 3e-18, "exp(x - 10) (lower)");

    y = t.exp();
    assert_equal_fp(test, y.upper, 1.0000000000000038, "exp(t) (upper)");
    assert_close_fp(test, y.lower, -2.4663276224724858e-17, 5e-16, "exp(t) (lower)");

    // The actual relative error for exp(709).lower is approx. 2e-14.
    // Can we do better?
    y = DoubleDouble(709.0).exp();
    assert_equal_fp(test, y.upper, 8.218407461554972e+307, "exp(709) upper");
    assert_close_fp(test, y.lower, -1.955965507696277e+291, 1e-13, "exp(709) lower");

    y = DoubleDouble(710.0).exp();
    assert_true(test, std::isinf(y.upper), "isinf(exp(710).upper)");
    assert_equal_fp(test, y.lower, 0.0, "exp(710).lower == 0");

    y = DoubleDouble(NAN).exp();
    assert_isnan(test, y);

    y = DoubleDouble(-1317.0).exp(); // 1.08e-572
    assert_equal_fp(test, y.upper, 0.0, "exp(-1317) (upper)");
    assert_equal_fp(test, y.lower, 0.0, "exp(-1317) (lower)");

    // >>> em669 = exp(-669).evalf(36)
    // >>> f = float(em669)
    // >>> e = float(em669 - f)
    // >>> f, e, em669
    // (2.864122616676439e-291, -4.0065259918771103e-308, 2.86412261667643910238552476779313774e-291)
    y = DoubleDouble(-669.0).exp();
    assert_equal_fp(test, y.upper, 2.864122616676439e-291, "exp(-669) (upper)");
    assert_close_fp(test, y.lower, -4.0065259918771103e-308, 4e-14, "exp(-669) (lower)");

}

void test_expm1(CheckIt& test)
{
    DoubleDouble y;
    auto t = DoubleDouble(3.7500950075008e-15, 0.0);

    y = DoubleDouble(11.75).expm1();
    assert_equal_fp(test, y.upper, 126752.55900574342, "expm1(11.75) upper");
    assert_close_fp(test, y.lower, 2.550992675486234e-12, 1e-15, "expm1(11.75) lower");

    y = DoubleDouble(1.0).expm1();
    assert_equal_fp(test, y.upper, 1.7182818284590453, "expm1(1) upper");
    assert_equal_fp(test, y.lower, -7.747991575210629e-17, "expm1(1) lower");

    y = DoubleDouble(0.46875).expm1();
    assert_equal_fp(test, y.upper, 0.5979954499506333, "expm1(0.46875) upper");
    assert_equal_fp(test, y.lower, 1.6864630310268093e-17, "expm1(0.46875) lower");

    y = t.expm1();
    assert_equal_fp(test, y.upper, 3.7500950075008074e-15, "expm1(t) upper");
    assert_equal_fp(test, y.lower, -6.814186434788356e-32, "expm1(t) lower");

    y = DoubleDouble(-0.000244140625).expm1();
    assert_equal_fp(test, y.upper, -0.00024411082510278348, "expm1(-2**-12) upper");
    assert_equal_fp(test, y.lower, -7.227720384831839e-21, "expm1(-2**-12) lower");

    y = DoubleDouble(-0.46875).expm1();
    assert_equal_fp(test, y.upper, -0.37421599039540887, "expm1(-0.46875) upper");
    assert_equal_fp(test, y.lower, -7.658883125910196e-18, "expm1(-0.46875) lower");

    y = DoubleDouble(-17.5).expm1();
    assert_equal_fp(test, y.upper, -0.9999999748900085, "expm1(-17.5) upper");
    assert_equal_fp(test, y.lower, 2.6947836187301293e-17, "expm1(-17.5) lower");

    y = DoubleDouble(NAN).expm1();
    assert_isnan(test, y);
}

void test_sin(CheckIt& test)
{
    for (const auto& sample : sin_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, sin(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("sin ") + sample.label);
        assert_dd_close(test, x.sin(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member sin ") + sample.label);
    }

    assert_isnan(test, sin(DoubleDouble(NAN)));
    assert_isnan(test, sin(dd_inf));
}

void test_cos(CheckIt& test)
{
    for (const auto& sample : cos_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, cos(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("cos ") + sample.label);
        assert_dd_close(test, x.cos(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member cos ") + sample.label);
    }

    assert_isnan(test, cos(DoubleDouble(NAN)));
    assert_isnan(test, cos(dd_inf));
}

void test_sincos(CheckIt& test)
{
    for (const auto& sample : sin_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        DoubleDouble s, c;
        sincos(x, s, c);
        assert_dd_close(test, s, sample.yhi, sample.ylo, sample.reltol,
                        std::string("sincos sin ") + sample.label);
    }

    for (const auto& sample : cos_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        DoubleDouble s, c;
        sincos(x, s, c);
        assert_dd_close(test, c, sample.yhi, sample.ylo, sample.reltol,
                        std::string("sincos cos ") + sample.label);
    }

    DoubleDouble s, c;
    sincos(DoubleDouble(NAN), s, c);
    assert_isnan(test, s);
    assert_isnan(test, c);
}

void test_cosm1(CheckIt& test)
{
    for (const auto& sample : cosm1_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, cosm1(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("cosm1 ") + sample.label);
        assert_dd_close(test, x.cosm1(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member cosm1 ") + sample.label);
    }

    assert_isnan(test, cosm1(DoubleDouble(NAN)));
    assert_isnan(test, cosm1(dd_inf));
}

void test_tan(CheckIt& test)
{
    for (const auto& sample : tan_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, tan(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("tan ") + sample.label);
        assert_dd_close(test, x.tan(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member tan ") + sample.label);
    }

    assert_isnan(test, tan(DoubleDouble(NAN)));
    assert_isnan(test, tan(dd_inf));
}

void test_atan(CheckIt& test)
{
    for (const auto& sample : atan_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, atan(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("atan ") + sample.label);
        assert_dd_close(test, x.atan(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member atan ") + sample.label);
    }

    assert_isnan(test, atan(DoubleDouble(NAN)));
}

void test_atan2(CheckIt& test)
{
    for (const auto& sample : atan2_cases) {
        DoubleDouble y{sample.yhi, sample.ylo};
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, atan2(y, x), sample.zhi, sample.zlo, sample.reltol,
                        std::string("atan2 ") + sample.label);
    }

    assert_isnan(test, atan2(DoubleDouble(NAN), dd_one));
    assert_isnan(test, atan2(dd_one, DoubleDouble(NAN)));
    assert_isnan(test, atan2(dd_zero, dd_zero));
}

void test_asin(CheckIt& test)
{
    for (const auto& sample : asin_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, asin(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("asin ") + sample.label);
        assert_dd_close(test, x.asin(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member asin ") + sample.label);
    }

    assert_isnan(test, asin(DoubleDouble(NAN)));
    assert_isnan(test, asin(DoubleDouble(1.0, 1e-30)));
    assert_isnan(test, asin(DoubleDouble(-1.0, -1e-30)));
}

void test_acos(CheckIt& test)
{
    for (const auto& sample : acos_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, acos(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("acos ") + sample.label);
        assert_dd_close(test, x.acos(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member acos ") + sample.label);
    }

    assert_isnan(test, acos(DoubleDouble(NAN)));
    assert_isnan(test, acos(DoubleDouble(1.0, 1e-30)));
    assert_isnan(test, acos(DoubleDouble(-1.0, -1e-30)));
}

void test_hypot(CheckIt& test)
{
    DoubleDouble x, y, h;

    x = DoubleDouble(3.0, 0.0);
    y = DoubleDouble(4.0, 0.0);
    h = hypot(x, y);
    assert_equal_fp(test, h.upper, 5.0, "hypot((3, 0), (4, 0)) (upper)");
    assert_equal_fp(test, h.lower, 0.0, "hypot((3, 0), (4, 0)) (lower)");

    x = DoubleDouble(0.0, 0.0);
    y = DoubleDouble(0.0, 0.0);
    h = hypot(x, y);
    assert_equal_fp(test, h.upper, 0.0, "hypot((0, 0), (0, 0)) (upper)");
    assert_equal_fp(test, h.lower, 0.0, "hypot((0, 0), (0, 0)) (lower)");

    x = DoubleDouble(1.0) / 3.0;
    y = DoubleDouble(7.0) / 5.0;
    h = hypot(x, y);
    assert_equal_fp(test, h.upper, 1.4391355429948602, "hypot(1/3, 7/5) (upper)");
    assert_equal_fp(test, h.lower, -3.962996388323609e-17, "hypot(1/3, 7/5)) (lower)");

    //
    // Reference values for the overflow and underflow tests of hypot
    // were computed with mpmath in Python:
    //
    //    from mpath import mp
    //    mp.dps = 500
    //    x = mp.mpf(5e250) + mp.mpf(5e225)
    //    y = mp.mpf(2e251) + mp.mpf(6e220)
    //    h = mp.sqrt(x**2 + y**2)
    //    h_upper = float(h)
    //    h_lower = float(h - h_upper)
    //
    x = DoubleDouble(5e250, 5e225);
    y = DoubleDouble(2e251, 6e220);
    h = hypot(x, y);
    assert_equal_fp(test, h.upper, 2.0615528128088305e+251, "hypot(BIG1, BIG2) (upper)");
    assert_close_fp(test, h.lower, -9.087579332326068e+234, 5e-16, "hypot(BIG1, BIG2) (lower)");

    x = DoubleDouble(3e-165, 5e-185);
    y = DoubleDouble(7.5e-165, -8e-184);
    h = hypot(x, y);
    assert_equal_fp(test, h.upper, 8.077747210701756e-165, "hypot(TINY1, TINY2) (upper)");
    assert_close_fp(test, h.lower, -5.79126190591365e-182, 1e-15, "hypot(TINY1, TINY2) (lower)");

    x = DoubleDouble(-17.25, -2.25e-19);
    y = DoubleDouble(0.0, 0.0);
    h = hypot(x, y);
    assert_equal_fp(test, h.upper, -x.upper, "hypot(x, (0, 0)) (upper)");
    assert_equal_fp(test, h.lower, -x.lower, "hypot(x, (0, 0)) (lower)");

    x = DoubleDouble(INFINITY);
    y = DoubleDouble(NAN);
    h = hypot(x, y);
    assert_equal_fp(test, h.upper, (double) INFINITY, "hypot(INF, NAN) is INF");
    h = hypot(y, x);
    assert_equal_fp(test, h.upper, (double) INFINITY, "hypot(INF, NAN) is INF");
}

static const unary_case sinh_cases[] = {
    {"0", 0.0, 0.0, 0.0, 0.0, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-20+3e-38", 9.99999999999999945e-21, 5.78467285457904316e-37, 9.99999999999999945e-21, 5.78467285457904316e-37, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133525374110523e-29, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.00000000000000002e-08, -4.25589416346180570e-26, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 1.00000000166666666e-04, 4.74635916281843722e-22, 5e-14},
    {"0.05", 5.00000000000000028e-02, -2.77555756156289151e-18, 5.00208359376550157e-02, -4.31868895548629297e-19, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.0, 1.25325775241115461e-01, -4.31830988622961368e-18, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.0, 5.21095305493747385e-01, -2.33281834764045972e-17, 5e-14},
    {"1", 1.0, 0.0, 1.17520119364380138e+00, 7.84967214228566853e-17, 5e-14},
    {"-1", -1.0, 0.0, -1.17520119364380138e+00, -7.84967214228566853e-17, 5e-14},
    {"10", 10.0, 0.0, 1.10132328747033935e+04, -8.60416647466444226e-14, 5e-14},
    {"-10", -10.0, 0.0, -1.10132328747033935e+04, 8.60416647466444226e-14, 5e-14},
    {"40", 40.0, 0.0, 1.17692633418510000e+17, -7.29605004462548301e+00, 5e-13},
    {"710", 710.0, 0.0, 1.11699738308085546e+308, 5.77253803440148121e+291, 5e-12}
};

static const unary_case cosh_cases[] = {
    {"0", 0.0, 0.0, 1.0, 0.0, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 1.0, 4.99999999999999981e-81, 5e-14},
    {"1e-20+3e-38", 9.99999999999999945e-21, 5.78467285457904316e-37, 1.0, 4.99999999999999965e-41, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 1.0, 4.99999999999999962e-25, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.0, 4.99999999999999990e-17, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 1.00000000499999997e+00, 3.45540215226663072e-17, 5e-14},
    {"0.05", 5.00000000000000028e-02, -2.77555756156289151e-18, 1.00125026043836907e+00, -4.19125468510407076e-17, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.0, 1.00782267782571089e+00, -2.88080034379573318e-17, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.0, 1.12762596520638070e+00, 8.70348011445619159e-17, 5e-14},
    {"1", 1.0, 0.0, 1.54308063481524371e+00, 6.60679677500683305e-17, 5e-14},
    {"-1", -1.0, 0.0, 1.54308063481524371e+00, 6.60679677500683305e-17, 5e-14},
    {"10", 10.0, 0.0, 1.10132329201033226e+04, 5.27017598240763681e-13, 5e-14},
    {"40", 40.0, 0.0, 1.17692633418510000e+17, -7.29605004462548301e+00, 5e-13},
    {"710", 710.0, 0.0, 1.11699738308085546e+308, 5.77253803440148121e+291, 5e-12}
};

static const unary_case tanh_reference_cases[] = {
    {"0", 0.0, 0.0, 0.0, 0.0, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-20+3e-38", 9.99999999999999945e-21, 5.78467285457904316e-37, 9.99999999999999945e-21, 5.78467285457904316e-37, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133520374110511e-29, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.00000000000000002e-08, -5.42558941634618019e-25, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 9.99999996666666683e-05, -2.73265489009813563e-22, 5e-14},
    {"0.05", 5.00000000000000028e-02, -2.77555756156289151e-18, 4.99583749578799696e-02,  2.57390062653096227e-18, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.0, 1.24353001771596208e-01, -2.14518808131414414e-19, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.0, 4.62117157260009737e-01,  2.19166032382609282e-17, 5e-14},
    {"1", 1.0, 0.0, 7.61594155955764851e-01,  3.70902144821649241e-17, 5e-14},
    {"-1", -1.0, 0.0, -7.61594155955764851e-01, -3.70902144821649241e-17, 5e-14},
    {"10", 10.0, 0.0, 9.99999995877692727e-01,  3.69327899661659591e-17, 5e-14},
    {"40", 40.0, 0.0, 1.0, -3.60970277569083060e-35, 5e-14},
    {"710", 710.0, 0.0, 1.0, 0.0, 5e-14}
};

static const unary_case asinh_cases[] = {
    {"0", 0.0, 0.0, 0.0, 0.0, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133522040777182e-29, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.00000000000000002e-08, -3.75892274967951365e-25, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 9.99999998333333298e-05,  4.32687736390399171e-21, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.0, 1.24676746921442747e-01, -2.58168289489360060e-18, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.0, 4.81211825059603471e-01, -2.32578170134627362e-17, 5e-14},
    {"1", 1.0, 0.0, 8.81373587019543048e-01, -2.25054589282586606e-17, 5e-14},
    {"-1", -1.0, 0.0, -8.81373587019543048e-01, 2.25054589282586606e-17, 5e-14},
    {"10", 10.0, 0.0, 2.99822295029796981e+00, -7.44026148092429399e-17, 5e-14},
    {"1e40", 1.00000000000000003e+40, -3.03786028427003667e+23, 9.27965509003217761e+01, -3.44178620965229906e-15, 5e-14}
};

static const unary_case acosh_cases[] = {
    {"1", 1.0, 0.0, 0.0, 0.0, 5e-14},
    {"1+2^-106", 1.0, 1.23259516440783095e-32, 1.57009245868377517e-16, -1.07328516224002451e-32, 5e-14},
    {"1+2^-54", 1.0, 5.55111512312578270e-17, 1.05367121277235087e-08, -7.69011564893361374e-25, 5e-14},
    {"1.0001", 1.00009999999999999e+00, 1.10134124042815525e-17, 1.41420177752523239e-02,  3.19545991309551628e-19, 1e-13},
    {"1.125", 1.12500000000000000e+00, 0.0, 4.94932923094526911e-01, -5.02253120893777902e-18, 5e-14},
    {"2", 2.0, 0.0, 1.31695789692481680e+00, -8.68225084485202174e-17, 5e-14},
    {"10", 10.0, 0.0, 2.99322284612638079e+00,  1.10644350351014916e-16, 5e-14},
    {"1e40", 1.00000000000000003e+40, -3.03786028427003667e+23, 9.27965509003217761e+01, -3.44178620965229906e-15, 5e-14}
};

static const unary_case atanh_reference_cases[] = {
    {"0", 0.0, 0.0, 0.0, 0.0, 5e-14},
    {"1e-40", 9.99999999999999929e-41, 7.07071206001198558e-57, 9.99999999999999929e-41, 7.07071206001198558e-57, 5e-14},
    {"1e-12", 9.99999999999999980e-13, 2.01133523707443852e-29, 9.99999999999999980e-13, 2.01133527040777193e-29, 5e-14},
    {"1e-8", 1.00000000000000002e-08, -2.09225608301284711e-25, 1.00000000000000002e-08,  1.24107725032048620e-25, 5e-14},
    {"1e-4", 1.00000000000000005e-04, -4.79217360238592994e-21, 1.00000000333333341e-04, -5.97774837353982333e-21, 5e-14},
    {"0.125", 1.25000000000000000e-01, 0.0, 1.25657214140453027e-01,  1.14469746886047311e-17, 5e-14},
    {"0.5", 5.00000000000000000e-01, 0.0, 5.49306144334054891e-01, -4.53564861750076498e-17, 5e-14},
    {"-0.5", -5.00000000000000000e-01, 0.0, -5.49306144334054891e-01,  4.53564861750076498e-17, 5e-14},
    {"1-2^-54", 1.0, -5.55111512312578270e-17, 1.90615474653984975e+01, -1.54107481201913735e-15, 5e-14},
    {"-1+2^-54", -1.0, 5.55111512312578270e-17, -1.90615474653984975e+01,  1.54107481201913735e-15, 5e-14}
};

void test_sinh(CheckIt& test)
{
    for (const auto& sample : sinh_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, sinh(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("sinh ") + sample.label);
        assert_dd_close(test, x.sinh(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member sinh ") + sample.label);
    }
    assert_isnan(test, sinh(DoubleDouble(NAN)));
    assert_true(test, std::isinf(sinh(dd_inf).upper) && sinh(dd_inf).upper > 0, "sinh(+inf) == +inf");
    assert_true(test, std::isinf(sinh(-dd_inf).upper) && sinh(-dd_inf).upper < 0, "sinh(-inf) == -inf");
}

void test_cosh(CheckIt& test)
{
    for (const auto& sample : cosh_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, cosh(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("cosh ") + sample.label);
        assert_dd_close(test, x.cosh(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member cosh ") + sample.label);
    }
    assert_isnan(test, cosh(DoubleDouble(NAN)));
    assert_true(test, std::isinf(cosh(dd_inf).upper) && cosh(dd_inf).upper > 0, "cosh(+inf) == +inf");
    assert_true(test, std::isinf(cosh(-dd_inf).upper) && cosh(-dd_inf).upper > 0, "cosh(-inf) == +inf");
}

void test_sincosh_hy(CheckIt& test)
{
    for (const auto& sample : sinh_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        DoubleDouble s, c;
        sincosh(x, s, c);
        assert_dd_close(test, s, sample.yhi, sample.ylo, sample.reltol,
                        std::string("sincosh sinh ") + sample.label);
    }
    for (const auto& sample : cosh_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        DoubleDouble s, c;
        sincosh(x, s, c);
        assert_dd_close(test, c, sample.yhi, sample.ylo, sample.reltol,
                        std::string("sincosh cosh ") + sample.label);
    }
    DoubleDouble s, c;
    sincosh(DoubleDouble(NAN), s, c);
    assert_isnan(test, s);
    assert_isnan(test, c);
}

void test_tanh(CheckIt &test)
{
  DoubleDouble x, y;
  x = 439.;
  y = tanh(x);
  assert_equal_fp(test, y.upper, 1.0, "tanh(438) (upper)");
  assert_equal_fp(test, y.lower, 0.0, "tanh(438) (lower)");

  for (const auto& sample : tanh_reference_cases) {
      DoubleDouble tx{sample.xhi, sample.xlo};
      assert_dd_close(test, tanh(tx), sample.yhi, sample.ylo, sample.reltol,
                      std::string("tanh ") + sample.label);
      assert_dd_close(test, tx.tanh(), sample.yhi, sample.ylo, sample.reltol,
                      std::string("member tanh ") + sample.label);
  }
  assert_isnan(test, tanh(DoubleDouble(NAN)));
  assert_equal_fp(test, tanh(dd_inf).upper, 1.0, "tanh(+inf) == 1");
  assert_equal_fp(test, tanh(-dd_inf).upper, -1.0, "tanh(-inf) == -1");
}

void test_asinh(CheckIt& test)
{
    for (const auto& sample : asinh_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, asinh(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("asinh ") + sample.label);
        assert_dd_close(test, x.asinh(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member asinh ") + sample.label);
    }
    assert_isnan(test, asinh(DoubleDouble(NAN)));
    assert_true(test, std::isinf(asinh(dd_inf).upper) && asinh(dd_inf).upper > 0, "asinh(+inf) == +inf");
    assert_true(test, std::isinf(asinh(-dd_inf).upper) && asinh(-dd_inf).upper < 0, "asinh(-inf) == -inf");
}

void test_acosh(CheckIt& test)
{
    for (const auto& sample : acosh_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, acosh(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("acosh ") + sample.label);
        assert_dd_close(test, x.acosh(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member acosh ") + sample.label);
    }
    assert_isnan(test, acosh(DoubleDouble(NAN)));
    assert_true(test, std::isinf(acosh(dd_inf).upper) && acosh(dd_inf).upper > 0, "acosh(+inf) == +inf");
    assert_isnan(test, acosh(DoubleDouble(0.5)));
    assert_isnan(test, acosh(DoubleDouble(-dd_inf)));
}

void test_atanh(CheckIt& test)
{
    for (const auto& sample : atanh_reference_cases) {
        DoubleDouble x{sample.xhi, sample.xlo};
        assert_dd_close(test, atanh(x), sample.yhi, sample.ylo, sample.reltol,
                        std::string("atanh ") + sample.label);
        assert_dd_close(test, x.atanh(), sample.yhi, sample.ylo, sample.reltol,
                        std::string("member atanh ") + sample.label);
    }
    assert_isnan(test, atanh(DoubleDouble(NAN)));
    assert_true(test, std::isinf(atanh(DoubleDouble(1.0)).upper) && atanh(DoubleDouble(1.0)).upper > 0, "atanh(1) == +inf");
    assert_true(test, std::isinf(atanh(DoubleDouble(-1.0)).upper) && atanh(DoubleDouble(-1.0)).upper < 0, "atanh(-1) == -inf");
    assert_isnan(test, atanh(DoubleDouble(1.5)));
    assert_isnan(test, atanh(DoubleDouble(-1.5)));
}

void test_hyperbolic_identities(CheckIt& test)
{
    double test_values[] = {-5.0, -1.0, -0.5, 0.0, 0.5, 1.0, 5.0};
    for (double val : test_values) {
        DoubleDouble x{val};
        DoubleDouble s = sinh(x);
        DoubleDouble c = cosh(x);
        DoubleDouble t = tanh(x);
        
        DoubleDouble one_check = c*c - s*s;
        assert_close_fp(test, one_check.upper, 1.0, 1e-15, std::string("cosh^2 - sinh^2 near 1 for x=") + std::to_string(val));
        assert_true(test, std::fabs(one_check.lower) < 1e-25, std::string("cosh^2 - sinh^2 lower near 0 for x=") + std::to_string(val));

        DoubleDouble t_check = s/c;
        assert_close_fp(test, t.upper, t_check.upper, 1e-15, std::string("tanh vs sinh/cosh upper for x=") + std::to_string(val));
        if (std::fabs(t_check.lower) < 1e-30) {
            assert_true(test, std::fabs(t.lower) < 1e-30, std::string("tanh vs sinh/cosh lower near 0 for x=") + std::to_string(val));
        } else {
            assert_close_fp(test, t.lower, t_check.lower, 1e-14, std::string("tanh vs sinh/cosh lower for x=") + std::to_string(val));
        }
    }

    double inv_test_values[] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    for (double val : inv_test_values) {
        DoubleDouble x{val};
        assert_close_fp(test, asinh(sinh(x)).upper, x.upper, 1e-15, std::string("asinh(sinh(x)) for x=") + std::to_string(val));
        assert_close_fp(test, atanh(tanh(x)).upper, x.upper, 1e-15, std::string("atanh(tanh(x)) for x=") + std::to_string(val));
    }

    double acosh_test_values[] = {1.0, 1.5, 2.0, 5.0};
    for (double val : acosh_test_values) {
        DoubleDouble x{val};
        assert_close_fp(test, acosh(cosh(x)).upper, x.upper, 1e-15, std::string("acosh(cosh(x)) for x=") + std::to_string(val));
    }
}

void test_dsum(CheckIt& test)
{
    double data1[]{1.0, 3.0, 99.0};
    double s1 = dsum(3, data1);
    assert_equal_fp(test, s1, 103.0, "Test dsum([1.0, 2.0, 99.0])");

    double data2[]{1.0, 2.0, 2e-17, -2.0, 10.0, -1.0, -10.0};
    double s2 = dsum(7, data2);
    assert_equal_fp(test, s2, 2e-17, "Test of dsum()");

    std::vector<double> data3{1.0, 2.0, 2e-17, -2.0, 10.0, -1.0, -10.0};
    double s3 = dsum(data3);
    assert_equal_fp(test, s3, 2e-17, "Test of dsum(vector)");

    std::array<double, 7> data4{1.0, 2.0, 2e-17, -2.0, 10.0, -1.0, -10.0};
    double s4 = dsum(data4);
    assert_equal_fp(test, s4, 2e-17, "Test of dsum(vector)");
}

void test_pow(CheckIt& test)
{
    /* The contents of test_pow was generated using Google's LLM: Gemini 3 Pro
       Preview. The prompt contained the, then current, implementation in doubledouble.h,
       cppreference's description of std::pow, and the current test suite.

       The generated test suite uncovered bugs in the implmentation, after
       fixing those manually (fixing typos, re-reading the cppreference spec),
       all generated tests passed without any modifications.
    */
    DoubleDouble y;
    DoubleDouble nan = DoubleDouble(NAN);
    DoubleDouble inf = dd_inf;
    DoubleDouble ninf = -dd_inf;
    DoubleDouble zero = dd_zero;
    DoubleDouble nzero = -dd_zero; // -0.0
    DoubleDouble one = dd_one;
    DoubleDouble mone = -dd_one;

    // ---------------------------------------------------------
    // 1. Exact Special Cases (Standard Compliance)
    // ---------------------------------------------------------

    // pow(+1, exp) returns 1 for any exp, even when exp is NaN
    assert_equal_fp(test, doubledouble::pow(one, DoubleDouble(5.0)).upper, 1.0, "pow(1, 5)");
    assert_equal_fp(test, doubledouble::pow(one, nan).upper, 1.0, "pow(1, NaN)");
    assert_equal_fp(test, doubledouble::pow(one, inf).upper, 1.0, "pow(1, inf)");

    // pow(base, ±0) returns 1 for any base, even when base is NaN
    assert_equal_fp(test, doubledouble::pow(DoubleDouble(5.0), zero).upper, 1.0, "pow(5, 0)");
    assert_equal_fp(test, doubledouble::pow(nan, zero).upper, 1.0, "pow(NaN, 0)");
    assert_equal_fp(test, doubledouble::pow(inf, zero).upper, 1.0, "pow(inf, 0)");
    assert_equal_fp(test, doubledouble::pow(DoubleDouble(5.0), nzero).upper, 1.0, "pow(5, -0)");

    // pow(-1, ±∞) returns 1
    assert_equal_fp(test, doubledouble::pow(mone, inf).upper, 1.0, "pow(-1, +inf)");
    assert_equal_fp(test, doubledouble::pow(mone, ninf).upper, 1.0, "pow(-1, -inf)");

    // ---------------------------------------------------------
    // 2. Base Zero Cases (±0)
    // ---------------------------------------------------------

    // pow(±0, -inf) returns +∞
    y = doubledouble::pow(zero, ninf);
    assert_true(test, std::isinf(y.upper) && y.upper > 0, "pow(+0, -inf) -> +inf");
    y = doubledouble::pow(nzero, ninf);
    assert_true(test, std::isinf(y.upper) && y.upper > 0, "pow(-0, -inf) -> +inf");

    // pow(+0, neg odd int) -> +∞
    y = doubledouble::pow(zero, DoubleDouble(-3.0));
    assert_true(test, std::isinf(y.upper) && y.upper > 0, "pow(+0, -3) -> +inf");

    // pow(-0, neg odd int) -> -∞
    y = doubledouble::pow(nzero, DoubleDouble(-3.0));
    assert_true(test, std::isinf(y.upper) && y.upper < 0, "pow(-0, -3) -> -inf");

    // pow(±0, neg even int/non-int) -> +∞
    y = doubledouble::pow(zero, DoubleDouble(-2.0));
    assert_true(test, std::isinf(y.upper) && y.upper > 0, "pow(+0, -2) -> +inf");
    y = doubledouble::pow(nzero, DoubleDouble(-2.0));
    assert_true(test, std::isinf(y.upper) && y.upper > 0, "pow(-0, -2) -> +inf");
    y = doubledouble::pow(nzero, DoubleDouble(-2.5));
    assert_true(test, std::isinf(y.upper) && y.upper > 0, "pow(-0, -2.5) -> +inf");

    // pow(+0, pos odd int) -> +0
    y = doubledouble::pow(zero, DoubleDouble(3.0));
    assert_equal_fp(test, y.upper, 0.0, "pow(+0, 3) -> +0");
    assert_true(test, !std::signbit(y.upper), "pow(+0, 3) sign is positive");

    // pow(-0, pos odd int) -> -0
    y = doubledouble::pow(nzero, DoubleDouble(3.0));
    assert_equal_fp(test, y.upper, 0.0, "pow(-0, 3) -> -0");
    // Note: checking sign bit of 0.0 relies on compiler/platform, but DD preserves it usually
    assert_true(test, std::signbit(y.upper), "pow(-0, 3) sign is negative");

    // pow(±0, pos even/non-int) -> +0
    y = doubledouble::pow(nzero, DoubleDouble(2.0));
    assert_equal_fp(test, y.upper, 0.0, "pow(-0, 2) -> +0");
    assert_true(test, !std::signbit(y.upper), "pow(-0, 2) sign is positive");
    y = doubledouble::pow(nzero, DoubleDouble(2.5));
    assert_equal_fp(test, y.upper, 0.0, "pow(-0, 2.5) -> +0");

    // ---------------------------------------------------------
    // 3. Exponent Infinite Cases (±∞)
    // ---------------------------------------------------------

    DoubleDouble half(0.5);
    DoubleDouble two(2.0);
    DoubleDouble mhalf(-0.5);
    DoubleDouble mtwo(-2.0);

    // pow(base, -∞) returns +∞ for |base| < 1
    assert_true(test, std::isinf(doubledouble::pow(half, ninf).upper), "pow(0.5, -inf) -> +inf");
    assert_true(test, std::isinf(doubledouble::pow(mhalf, ninf).upper), "pow(-0.5, -inf) -> +inf");

    // pow(base, -∞) returns +0 for |base| > 1
    assert_equal_fp(test, doubledouble::pow(two, ninf).upper, 0.0, "pow(2.0, -inf) -> 0");
    assert_equal_fp(test, doubledouble::pow(mtwo, ninf).upper, 0.0, "pow(-2.0, -inf) -> 0");

    // pow(base, +∞) returns +0 for |base| < 1
    assert_equal_fp(test, doubledouble::pow(half, inf).upper, 0.0, "pow(0.5, inf) -> 0");
    assert_equal_fp(test, doubledouble::pow(mhalf, inf).upper, 0.0, "pow(-0.5, inf) -> 0");

    // pow(base, +∞) returns +∞ for |base| > 1
    assert_true(test, std::isinf(doubledouble::pow(two, inf).upper), "pow(2.0, inf) -> inf");
    assert_true(test, std::isinf(doubledouble::pow(mtwo, inf).upper), "pow(-2.0, inf) -> inf");

    // ---------------------------------------------------------
    // 4. Base Infinite Cases (±∞)
    // ---------------------------------------------------------

    // pow(-∞, neg odd int) -> -0
    y = doubledouble::pow(ninf, DoubleDouble(-3.0));
    assert_equal_fp(test, y.upper, 0.0, "pow(-inf, -3) -> -0");
    assert_true(test, std::signbit(y.upper), "pow(-inf, -3) sign negative");

    // pow(-∞, neg even int) -> +0
    y = doubledouble::pow(ninf, DoubleDouble(-2.0));
    assert_equal_fp(test, y.upper, 0.0, "pow(-inf, -2) -> +0");
    assert_true(test, !std::signbit(y.upper), "pow(-inf, -2) sign positive");

    // pow(-∞, pos odd int) -> -∞
    y = doubledouble::pow(ninf, DoubleDouble(3.0));
    assert_true(test, std::isinf(y.upper) && y.upper < 0, "pow(-inf, 3) -> -inf");

    // pow(-∞, pos even int) -> +∞
    y = doubledouble::pow(ninf, DoubleDouble(2.0));
    assert_true(test, std::isinf(y.upper) && y.upper > 0, "pow(-inf, 2) -> +inf");

    // pow(+∞, neg) -> +0
    assert_equal_fp(test, doubledouble::pow(inf, DoubleDouble(-2.0)).upper, 0.0, "pow(inf, -2)");

    // pow(+∞, pos) -> +∞
    assert_true(test, std::isinf(doubledouble::pow(inf, DoubleDouble(2.0)).upper), "pow(inf, 2)");

    // ---------------------------------------------------------
    // 5. Negative Base Finite Cases
    // ---------------------------------------------------------

    // pow(neg, finite non-int) -> NaN (Domain Error)
    y = doubledouble::pow(mtwo, DoubleDouble(2.5));
    assert_isnan(test, y);

    // pow(neg, integer) -> Correct value
    y = doubledouble::pow(mtwo, DoubleDouble(3.0));
    assert_equal_fp(test, y.upper, -8.0, "pow(-2, 3)");

    y = doubledouble::pow(mtwo, DoubleDouble(2.0));
    assert_equal_fp(test, y.upper, 4.0, "pow(-2, 2)");

    // ---------------------------------------------------------
    // 6. General Arithmetic and Accuracy
    // ---------------------------------------------------------

    // 3^4 = 81
    y = doubledouble::pow(DoubleDouble(3.0), DoubleDouble(4.0));
    assert_equal_fp(test, y.upper, 81.0, "pow(3, 4) approx 81");
    assert_close_fp(test, y.lower, 0.0, 1e-14, "pow(3, 4) precision check");

    // sqrt via pow: 4^0.5 = 2
    y = doubledouble::pow(DoubleDouble(4.0), DoubleDouble(0.5));
    assert_equal_fp(test, y.upper, 2.0, "pow(4, 0.5)");

    // ---------------------------------------------------------
    // 7. Potential Bug Triggers / Edge Cases
    // ---------------------------------------------------------

    // Case: Exponent is a large even integer (Double).
    // The implementation casts to `long long int`.
    // 1e20 > 2^63. Casting 1e20 to long long is Undefined Behavior.
    // However, math says (-2)^(even) = positive.
    DoubleDouble large_even(1e20);
    // We expect +Infinity because 2^1e20 is huge, but we specifically check the SIGN.
    // If the UB results in a negative interpretation of parity, this might fail or crash.
    // Commented out to prevent crash if running in strict env, but valid for bug hunting:
    /*
      y = doubledouble::pow(mtwo, large_even);
      assert_true(test, std::isinf(y.upper), "pow(-2, 1e20) overflow to inf");
      assert_true(test, y.upper > 0, "pow(-2, 1e20) should be positive inf (even power)");
    */

    // Case: Base is 1.0 + epsilon
    // Implementation uses strict (base == dd_one).
    // pow(1+eps, inf) should be inf, NOT 1.
    DoubleDouble one_plus_eps = one + DoubleDouble(1e-18); // still > 1
    y = doubledouble::pow(one_plus_eps, inf);
    assert_true(test, std::isinf(y.upper), "pow(1+eps, inf) -> inf");
    assert_true(test, y.upper > 0, "pow(1+eps, inf) is pos");

    // Case: Base is 1.0 - epsilon
    // pow(1-eps, inf) should be 0.
    DoubleDouble one_minus_eps = one - DoubleDouble(1e-18);
    y = doubledouble::pow(one_minus_eps, inf);
    assert_equal_fp(test, y.upper, 0.0, "pow(1-eps, inf) -> 0");

}

int main(int argc, char *argv[])
{
    auto test = CheckIt(std::cerr);

    test_constructor(test);
    test_add(test);
    test_inplace_add(test);
    test_subtract(test);
    test_inplace_subtract(test);
    test_multiply(test);
    test_inplace_multiply(test);
    test_divide(test);
    test_inplace_divide(test);
    test_expressions(test);
    test_comparisons(test);
    test_abs(test);
    test_powi(test);
    test_sqrt(test);
    test_log(test);
    test_log1p(test);
    test_exp(test);
    test_expm1(test);
    test_sin(test);
    test_cos(test);
    test_sincos(test);
    test_cosm1(test);
    test_tan(test);
    test_atan(test);
    test_atan2(test);
    test_asin(test);
    test_acos(test);
    test_hypot(test);
    test_sinh(test);
    test_cosh(test);
    test_sincosh_hy(test);
    test_tanh(test);
    test_asinh(test);
    test_acosh(test);
    test_atanh(test);
    test_hyperbolic_identities(test);
    test_dsum(test);
    test_pow(test);
    return test.print_summary("Summary: ");
}
