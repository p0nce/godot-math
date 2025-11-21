// Godot top-level functions, were separated in a module
// so as to only import libc here.
module godotmath.globals;

import std.math: PI;
static import core.stdc.math;

nothrow @nogc @safe:

// See: https://docs.godotengine.org/en/stable/classes/class_%40globalscope.html
// Those are math functions defined at global scope in Godot Math.
// INTERNAL USE ONLY, unlike Godot where this is part of API.

// godot math constants
enum double GM_SQRT2  = 1.4142135623730950488016887242;
enum double GM_SQRT3  = 1.7320508075688772935274463415059;
enum double GM_SQRT12 = 0.7071067811865475244008443621048490;
enum double GM_SQRT13 = 0.57735026918962576450914878050196;
enum double GM_LN2    = 0.6931471805599453094172321215;
enum double GM_TAU    = 6.2831853071795864769252867666;
enum double GM_PI     = 3.1415926535897932384626433833;
enum double GM_E      = 2.7182818284590452353602874714;
enum double GM_INF    = double.infinity;
enum double GM_NaN    = double.nan;

enum double GM_CMP_EPSILON = 0.00001;
enum double GM_CMP_EPSILON2 = (GM_CMP_EPSILON * GM_CMP_EPSILON);

// PRECISE_MATH_CHECKS
enum GM_UNIT_EPSILON = 0.00001;


float abs(float x) => core.stdc.math.fabsf(x);
double abs(double x) => core.stdc.math.fabs(x);
int absi(int x) => (x >= 0 ? x : -x);

float  acos(float x)   => (x < -1) ? PI : (x > 1 ? 0 : core.stdc.math.acosf(x));
double acos(double x)  => (x < -1) ? PI : (x > 1 ? 0 : core.stdc.math.acos(x));
float  acosh(float x)  => core.stdc.math.acoshf(x);
double acosh(double x) => core.stdc.math.acosh(x);

// TODO: angle_difference



float  asin(float x)   => x < -1 ? (-PI / 2) : (x > 1 ? (PI / 2) : core.stdc.math.asinf(x));
double asin(double x)  => x < -1 ? (-PI / 2) : (x > 1 ? (PI / 2) : core.stdc.math.asin(x));
float  asinh(float x) => core.stdc.math.asinhf(x);
double asinh(double x) => core.stdc.math.acosh(x);
float  atan2(float y, float x) => core.stdc.math.atan2(y, x);
double atan2(double y, double x) => core.stdc.math.atan2(y, x);

float bezier_derivative(float p_start, float p_control_1, 
                        float p_control_2, float p_end, float p_t) 
{
    float omt = (1.0f - p_t);
    float omt2 = omt * omt;
    float t2 = p_t * p_t;
    float d = (p_control_1 - p_start) * 3.0f * omt2 + (p_control_2 - p_control_1) * 6.0f * omt * p_t + (p_end - p_control_2) * 3.0f * t2;
    return d;
}

double bezier_derivative(double p_start, double p_control_1, 
                         double p_control_2, double p_end, double p_t) 
{
    double omt = (1.0 - p_t);
    double omt2 = omt * omt;
    double t2 = p_t * p_t;

    double d = (p_control_1 - p_start) * 3.0 * omt2 + (p_control_2 - p_control_1) * 6.0 * omt * p_t + (p_end - p_control_2) * 3.0 * t2;
    return d;
}

float bezier_interpolate(float p_start, float p_control_1, float p_control_2, float p_end, float p_t) 
{
    float omt = (1.0f - p_t);
    float omt2 = omt * omt;
    float omt3 = omt2 * omt;
    float t2 = p_t * p_t;
    float t3 = t2 * p_t;
    return p_start * omt3 + p_control_1 * omt2 * p_t * 3.0f + p_control_2 * omt * t2 * 3.0f + p_end * t3;
}

double bezier_interpolate(double p_start, double p_control_1, double p_control_2, double p_end, double p_t) 
{
    double omt = (1.0 - p_t);
    double omt2 = omt * omt;
    double omt3 = omt2 * omt;
    double t2 = p_t * p_t;
    double t3 = t2 * p_t;
    return p_start * omt3 + p_control_1 * omt2 * p_t * 3.0 + p_control_2 * omt * t2 * 3.0 + p_end * t3;
}

float  ceil(float x)  => core.stdc.math.ceilf(x);
double ceil(double x) => core.stdc.math.ceil(x);

float clampf(float value, float min, float max)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

int clampi(int value, int min, int max)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

float  cos(float x)  => core.stdc.math.cosf(x);
double cos(double x) => core.stdc.math.cos(x);

double cubic_interpolate(double p_from, double p_to, double p_pre, double p_post, double p_weight) 
{
    return 0.5 *
        ((p_from * 2.0) +
         (-p_pre + p_to) * p_weight +
         (2.0 * p_pre - 5.0 * p_from + 4.0 * p_to - p_post) * (p_weight * p_weight) +
         (-p_pre + 3.0 * p_from - 3.0 * p_to + p_post) * (p_weight * p_weight * p_weight));
}

float cubic_interpolate(float p_from, float p_to, float p_pre, float p_post, float p_weight) 
{
    return 0.5f *
        ((p_from * 2.0f) +
         (-p_pre + p_to) * p_weight +
         (2.0f * p_pre - 5.0f * p_from + 4.0f * p_to - p_post) * (p_weight * p_weight) +
         (-p_pre + 3.0f * p_from - 3.0f * p_to + p_post) * (p_weight * p_weight * p_weight));
}

double cubic_interpolate_in_time(double p_from, double p_to, double p_pre, double p_post, double p_weight,
                                 double p_to_t, double p_pre_t, double p_post_t) 
{
    /* Barry-Goldman method */
    double t = lerp(0.0, p_to_t, p_weight);
    double a1 = lerp(p_pre, p_from, p_pre_t == 0 ? 0.0 : (t - p_pre_t) / -p_pre_t);
    double a2 = lerp(p_from, p_to, p_to_t == 0 ? 0.5 : t / p_to_t);
    double a3 = lerp(p_to, p_post, p_post_t - p_to_t == 0 ? 1.0 : (t - p_to_t) / (p_post_t - p_to_t));
    double b1 = lerp(a1, a2, p_to_t - p_pre_t == 0 ? 0.0 : (t - p_pre_t) / (p_to_t - p_pre_t));
    double b2 = lerp(a2, a3, p_post_t == 0 ? 1.0 : t / p_post_t);
    return lerp(b1, b2, p_to_t == 0 ? 0.5 : t / p_to_t);
}
float cubic_interpolate_in_time(float p_from, float p_to, float p_pre, float p_post, float p_weight,
                                float p_to_t, float p_pre_t, float p_post_t) 
{
    /* Barry-Goldman method */
    float t = lerp(0.0f, p_to_t, p_weight);
    float a1 = lerp(p_pre, p_from, p_pre_t == 0 ? 0.0f : (t - p_pre_t) / -p_pre_t);
    float a2 = lerp(p_from, p_to, p_to_t == 0 ? 0.5f : t / p_to_t);
    float a3 = lerp(p_to, p_post, p_post_t - p_to_t == 0 ? 1.0f : (t - p_to_t) / (p_post_t - p_to_t));
    float b1 = lerp(a1, a2, p_to_t - p_pre_t == 0 ? 0.0f : (t - p_pre_t) / (p_to_t - p_pre_t));
    float b2 = lerp(a2, a3, p_post_t == 0 ? 1.0f : t / p_post_t);
    return lerp(b1, b2, p_to_t == 0 ? 0.5f : t / p_to_t);
}

float  floor(float x)  => core.stdc.math.floorf(x);
double floor(double x) => core.stdc.math.floor(x);

double fposmod(double x, double y) 
{
    double value = core.stdc.math.fmod(x, y);
    if (((value < 0) && (y > 0)) || ((value > 0) && (y < 0))) 
    {
        value += y;
    }
    value += 0.0;
    return value;
}

float fposmod(float x, float y) 
{
    float value = core.stdc.math.fmodf(x, y);
    if (((value < 0) && (y > 0)) || ((value > 0) && (y < 0))) 
    {
        value += y;
    }
    value += 0.0f;
    return value;
}

bool is_equal_approx(float p_left, float p_right) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) {
        return true;
    }
    // Then check for approximate equality.
    float tolerance = cast(float)GM_CMP_EPSILON * abs(p_left);
    if (tolerance < cast(float)GM_CMP_EPSILON) 
    {
        tolerance = cast(float)GM_CMP_EPSILON;
    }
    return abs(p_left - p_right) < tolerance;
}

bool is_equal_approx(double p_left, double p_right) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) 
    {
        return true;
    }
    // Then check for approximate equality.
    double tolerance = GM_CMP_EPSILON * abs(p_left);
    if (tolerance < GM_CMP_EPSILON) 
    {
        tolerance = GM_CMP_EPSILON;
    }
    return abs(p_left - p_right) < tolerance;
}

bool is_equal_approx(float p_left, float p_right, float p_tolerance) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right)
        return true;

    // Then check for approximate equality.
    return abs(p_left - p_right) < p_tolerance;
}

bool is_equal_approx(double p_left, double p_right, double p_tolerance) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) 
    {
        return true;
    }
    // Then check for approximate equality.
    return abs(p_left - p_right) < p_tolerance;
}

float  is_finite(float x)  => core.stdc.math.isfinite(x);
double is_finite(double x) => core.stdc.math.isfinite(x);
bool is_zero_approx(float  value) => abs(value) < cast(float)GM_CMP_EPSILON;
bool is_zero_approx(double value) => abs(value) < GM_CMP_EPSILON;


float lerp(float from, float to, float weight) 
{
    return from + (to - from) * weight;
}

double lerp(double from, double to, double weight) 
{
    return from + (to - from) * weight;
}

float  round(float x)  => core.stdc.math.roundf(x);
double round(double x) => core.stdc.math.round(x);

float  sign(float v)  => (v > 0) ? 1.0f : (v < 0 ? -1.0f : 0.0f);
double sign(double v) =>(v > 0) ? 1.0f : (v < 0 ? -1.0f : 0.0f);

float  sin(float x)  => core.stdc.math.sinf(x);
double sin(double x) => core.stdc.math.sin(x);

float  sqrt(float x)  => core.stdc.math.sqrtf(x);
double sqrt(double x) => core.stdc.math.sqrt(x);