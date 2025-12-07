/**
Abstract over C math functions to get `pure` working.
The issue is that core.stdc.math is missing annotations, so 
this is a workaround, also provide float/double overloads.

Copyright (c) 2014-2025 Godot Engine contributors.
Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur. 
Copyright (c) 2025 Guillaume Piolat
*/
module godotmath.funcs;


import numem;
import numem.core.traits;
static import libc = core.stdc.math;

// See: https://docs.godotengine.org/en/stable/classes/class_%40globalscope.html
// Those are math functions defined at global scope in Godot Math, so here they
// are prefixed with a gm_*** prefix to avoid polluting global namespace.

pure nothrow @nogc @safe:

// Godot math constants
// FUTURE: use the better constants from std.math, like PI with hex float precision**
// Phobos constants have more precision than Godot's ones.
enum double GM_SQRT2  = 1.4142135623730950488016887242;
enum double GM_SQRT3  = 1.7320508075688772935274463415059;
enum double GM_SQRT12 = 0.7071067811865475244008443621048490;
enum double GM_SQRT13 = 0.57735026918962576450914878050196;
enum double GM_LN2    = 0.6931471805599453094172321215;
enum double GM_TAU    = 6.2831853071795864769252867666;
enum double GM_PI     = 0x1.921fb54442d18469898cc51701b84p+1;
enum double GM_E      = 2.7182818284590452353602874714;
enum double GM_INF    = double.infinity;
enum double GM_NaN    = double.nan;

// PRECISE_MATH_CHECKS is considered defined in this translation
enum double GM_CMP_EPSILON = 0.00001;
enum double GM_CMP_EPSILON2 = (GM_CMP_EPSILON * GM_CMP_EPSILON);

// PRECISE_MATH_CHECKS
enum GM_UNIT_EPSILON = 0.00001;

float gm_abs(float x) => libc.fabsf(x); ///
double gm_abs(double x) => libc.fabs(x); ///
int gm_abs(int x) => (x >= 0 ? x : -x); ///
float  gm_acos(float x)  @trusted => (x < -1) ? GM_PI : (x > 1 ? 0 : assumePureNothrowNogc(&libc.acosf, x)); ///
double gm_acos(double x)  => (x < -1) ? GM_PI : (x > 1 ? 0 : assumePureNothrowNogc(&libc.acos, x)); ///
float  gm_acosh(float x)  => assumePureNothrowNogc(&libc.acoshf, x); ///
double gm_acosh(double x) => assumePureNothrowNogc(&libc.acosh, x); ///
float  gm_asin(float x)   => x < -1 ? (-GM_PI / 2) : (x > 1 ? (GM_PI / 2) : assumePureNothrowNogc(&libc.asinf, x)); ///
double gm_asin(double x)  => x < -1 ? (-GM_PI / 2) : (x > 1 ? (GM_PI / 2) : assumePureNothrowNogc(&libc.asin, x)); ///
float  gm_asinh(float x) => assumePureNothrowNogc(&libc.asinhf, x); ///
double gm_asinh(double x) => assumePureNothrowNogc(&libc.asinh, x); ///
float  gm_atan2(float y, float x) => assumePureNothrowNogc(&libc.atan2f, y, x); ///
double gm_atan2(double y, double x) => assumePureNothrowNogc(&libc.atan2, y, x); ///

float gm_bezier_derivative(float p_start, float p_control_1, 
                           float p_control_2, float p_end, float p_t) ///
{
    float omt = (1.0f - p_t);
    float omt2 = omt * omt;
    float t2 = p_t * p_t;
    float d = (p_control_1 - p_start) * 3.0f * omt2 + (p_control_2 - p_control_1) * 6.0f * omt * p_t + (p_end - p_control_2) * 3.0f * t2;
    return d;
}

double gm_bezier_derivative(double p_start, double p_control_1, 
                            double p_control_2, double p_end, double p_t) /// 
{
    double omt = (1.0 - p_t);
    double omt2 = omt * omt;
    double t2 = p_t * p_t;

    double d = (p_control_1 - p_start) * 3.0 * omt2 + (p_control_2 - p_control_1) * 6.0 * omt * p_t + (p_end - p_control_2) * 3.0 * t2;
    return d;
}

float gm_bezier_interpolate(float p_start, float p_control_1, float p_control_2, float p_end, float p_t) ///
{
    float omt = (1.0f - p_t);
    float omt2 = omt * omt;
    float omt3 = omt2 * omt;
    float t2 = p_t * p_t;
    float t3 = t2 * p_t;
    return p_start * omt3 + p_control_1 * omt2 * p_t * 3.0f + p_control_2 * omt * t2 * 3.0f + p_end * t3;
}

double gm_bezier_interpolate(double p_start, double p_control_1, double p_control_2, double p_end, double p_t) ///
{
    double omt = (1.0 - p_t);
    double omt2 = omt * omt;
    double omt3 = omt2 * omt;
    double t2 = p_t * p_t;
    double t3 = t2 * p_t;
    return p_start * omt3 + p_control_1 * omt2 * p_t * 3.0 + p_control_2 * omt * t2 * 3.0 + p_end * t3;
}

float  gm_ceil(float x)  => libc.ceilf(x); ///
double gm_ceil(double x) => libc.ceil(x); ///

double gm_clamp(double value, double min, double max) ///
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

float gm_clamp(float value, float min, float max) ///
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

int gm_clamp(int value, int min, int max) ///
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

float  gm_cos(float x)  => libc.cosf(x); ///
double gm_cos(double x) => libc.cos(x); ///

double gm_cubic_interpolate(double p_from, double p_to, double p_pre, double p_post, double p_weight) ///
{
    return 0.5 *
        ((p_from * 2.0) +
         (-p_pre + p_to) * p_weight +
         (2.0 * p_pre - 5.0 * p_from + 4.0 * p_to - p_post) * (p_weight * p_weight) +
         (-p_pre + 3.0 * p_from - 3.0 * p_to + p_post) * (p_weight * p_weight * p_weight));
}

float gm_cubic_interpolate(float p_from, float p_to, float p_pre, float p_post, float p_weight) ///
{
    return 0.5f *
        ((p_from * 2.0f) +
         (-p_pre + p_to) * p_weight +
         (2.0f * p_pre - 5.0f * p_from + 4.0f * p_to - p_post) * (p_weight * p_weight) +
         (-p_pre + 3.0f * p_from - 3.0f * p_to + p_post) * (p_weight * p_weight * p_weight));
}

double gm_cubic_interpolate_in_time(double p_from, double p_to, double p_pre, double p_post, double p_weight,
                                    double p_to_t, double p_pre_t, double p_post_t) ///
{
    /* Barry-Goldman method */
    double t = gm_lerp(0.0, p_to_t, p_weight);
    double a1 = gm_lerp(p_pre, p_from, p_pre_t == 0 ? 0.0 : (t - p_pre_t) / -p_pre_t);
    double a2 = gm_lerp(p_from, p_to, p_to_t == 0 ? 0.5 : t / p_to_t);
    double a3 = gm_lerp(p_to, p_post, p_post_t - p_to_t == 0 ? 1.0 : (t - p_to_t) / (p_post_t - p_to_t));
    double b1 = gm_lerp(a1, a2, p_to_t - p_pre_t == 0 ? 0.0 : (t - p_pre_t) / (p_to_t - p_pre_t));
    double b2 = gm_lerp(a2, a3, p_post_t == 0 ? 1.0 : t / p_post_t);
    return gm_lerp(b1, b2, p_to_t == 0 ? 0.5 : t / p_to_t);
}
float gm_cubic_interpolate_in_time(float p_from, float p_to, float p_pre, float p_post, float p_weight,
                                   float p_to_t, float p_pre_t, float p_post_t) ///
{
    /* Barry-Goldman method */
    float t = gm_lerp(0.0f, p_to_t, p_weight);
    float a1 = gm_lerp(p_pre, p_from, p_pre_t == 0 ? 0.0f : (t - p_pre_t) / -p_pre_t);
    float a2 = gm_lerp(p_from, p_to, p_to_t == 0 ? 0.5f : t / p_to_t);
    float a3 = gm_lerp(p_to, p_post, p_post_t - p_to_t == 0 ? 1.0f : (t - p_to_t) / (p_post_t - p_to_t));
    float b1 = gm_lerp(a1, a2, p_to_t - p_pre_t == 0 ? 0.0f : (t - p_pre_t) / (p_to_t - p_pre_t));
    float b2 = gm_lerp(a2, a3, p_post_t == 0 ? 1.0f : t / p_post_t);
    return gm_lerp(b1, b2, p_to_t == 0 ? 0.5f : t / p_to_t);
}

float  gm_floor(float x)  => libc.floorf(x); ///
double gm_floor(double x) => libc.floor(x); ///
float gm_fmod(float x, float y) => assumePureNothrowNogc(&libc.fmodf, x, y); ///
double gm_fmod(double x, double y) => assumePureNothrowNogc(&libc.fmod, x, y); ///

double gm_fposmod(double x, double y) ///
{
    double value = gm_fmod(x, y);
    if (((value < 0) && (y > 0)) || ((value > 0) && (y < 0))) 
    {
        value += y;
    }
    value += 0.0;
    return value;
}

float gm_fposmod(float x, float y) ///
{
    float value = gm_fmod(x, y);
    if (((value < 0) && (y > 0)) || ((value > 0) && (y < 0))) 
    {
        value += y;
    }
    value += 0.0f;
    return value;
}

bool gm_is_equal_approx(float p_left, float p_right) ///
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) {
        return true;
    }
    // Then check for approximate equality.
    float tolerance = cast(float)GM_CMP_EPSILON * gm_abs(p_left);
    if (tolerance < cast(float)GM_CMP_EPSILON) 
    {
        tolerance = cast(float)GM_CMP_EPSILON;
    }
    return gm_abs(p_left - p_right) < tolerance;
}

bool gm_is_equal_approx(double p_left, double p_right) ///
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) 
    {
        return true;
    }
    // Then check for approximate equality.
    double tolerance = GM_CMP_EPSILON * gm_abs(p_left);
    if (tolerance < GM_CMP_EPSILON) 
    {
        tolerance = GM_CMP_EPSILON;
    }
    return gm_abs(p_left - p_right) < tolerance;
}

bool gm_is_equal_approx(float p_left, float p_right, float p_tolerance) /// 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right)
        return true;

    // Then check for approximate equality.
    return gm_abs(p_left - p_right) < p_tolerance;
}

bool gm_is_equal_approx(double p_left, double p_right, double p_tolerance) ///
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) 
    {
        return true;
    }
    // Then check for approximate equality.
    return gm_abs(p_left - p_right) < p_tolerance;
}

bool gm_is_finite(float x)  => libc.isfinite(x) != 0; ///
bool gm_is_finite(double x) => libc.isfinite(x) != 0; ///

bool gm_is_zero_approx(float  value) => gm_abs(value) < cast(float)GM_CMP_EPSILON;
bool gm_is_zero_approx(double value) => gm_abs(value) < GM_CMP_EPSILON;

float gm_lerp(float from, float to, float weight) 
{
    return from + (to - from) * weight;
}

double gm_lerp(double from, double to, double weight) ///
{
    return from + (to - from) * weight;
}

float  gm_round(float x)  => libc.roundf(x); ///
double gm_round(double x) => libc.round(x); ///

int    gm_sign(int v)    => (v > 0) ? 1 : (v < 0 ? -1 : 0); ///
float  gm_sign(float v)  => (v > 0) ? 1.0f : (v < 0 ? -1.0f : 0.0f); ///
double gm_sign(double v) => (v > 0) ? 1.0 : (v < 0 ? -1.0 : 0.0); ///
float  gm_sin(float x)   => libc.sinf(x); ///
double gm_sin(double x)  => libc.sin(x); ///

double gm_snapped(double value, double step) ///
{
    // BUG: floor(x + 0.5) not exactly the same as round
    // will hardly be met in practice
	if (step != 0) 
    {
		value = gm_floor(value / step + 0.5) * step;
	}
	return value;
}

int gm_snapped(int value, int step) ///
{
    // Strangely enough, Godot doesn't do anything special for integers,
    // it doesn't use a correct integer rounding.
    return cast(int) gm_snapped(value, step);
}

float  gm_sqrt(float x)  => assumePureNothrowNogc(&libc.sqrtf, x); ///
double gm_sqrt(double x) => assumePureNothrowNogc(&libc.sqrt, x); ///

void gm_swap(T)(ref T a, ref T b)
{
    T tmp = a;
    a = b;
    b = tmp;
}

private:

private auto assumePureNothrowNogc(T, Args...)(T expr, auto ref Args args) pure nothrow @nogc @trusted if (isSomeFunction!T) {
    static if (is(T Fptr : Fptr*) && is(Fptr == function))
        alias ft = pure nothrow @nogc ReturnType!T function(Parameters!T);
    else static if (is(T Fdlg == delegate))
        alias ft = pure nothrow @nogc ReturnType!T delegate(Parameters!T);
    else
        static assert(0);

    return (cast(ft)expr)(args);
}
