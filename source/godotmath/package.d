/*
    Copyright (c) 2014-2025 Godot Engine contributors.
    Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur. 
    Copyright (c) 2025 Guillaume Piolat
*/
module godotmath;

import numem;
import numem.core.traits;
static import libc = core.stdc.math;

pure nothrow @nogc @safe:

// Changes vs Godot
// - Both float and double versions exist at once, with a ***d suffix.
// - in global scope, symbols get a gm_ prefix.
// - .clampf/.clampi replaced by overloaded .clamp
//   Likewise .snappedf /.snappedi replaced by overloaded .snapped. 
//   Same for minf/mini/maxf/maxi => replaced by min/max
// - in Godot names "f" would mean "float", and "float" is `double` in Godot.
//   Here when something is "float" it is actually `float` in the interface too,
//   and same for `double`.

// TODO all operators for Vector3
// TODO all operators for Vector4
// TODO all operators for Transform2D
// TODO all operators for Transform3D
// TODO all operators for Basis

// Provide both float and double versions, should the need arise.
alias Vector2  = Vector2Impl!float;  ///
alias Vector2i = Vector2Impl!int;    ///
alias Vector2d = Vector2Impl!double; ///

alias Size2    = Vector2; ///
alias Point2   = Vector2; ///

alias Vector3  = Vector3Impl!float;  ///
alias Vector3i = Vector3Impl!int;    ///
alias Vector3d = Vector3Impl!double; ///

alias Vector4  = Vector4Impl!float;  ///
alias Vector4i = Vector4Impl!int;    ///
alias Vector4d = Vector4Impl!double; ///

alias Quaternion  = QuaternionImpl!float;  ///
alias Quaterniond = QuaternionImpl!double; ///

// 2x3 matrix
alias Transform2D  = Transform2DImpl!float;  ///
alias Transform2Dd = Transform2DImpl!double; ///

// 3x3 matrix
alias Basis    = BasisImpl!float;  ///
alias Basisd   = BasisImpl!double; ///

// 3x4 matrix
alias Transform3D  = Transform3DImpl!float;  ///
alias Transform3Dd = Transform3DImpl!double; ///

// 4x4 matrix
alias Projection  = ProjectionImpl!float;  ///
alias Projectiond = ProjectionImpl!double; ///



// EulerOrder
alias EulerOrder = int; ///
enum : EulerOrder
{
    GM_EULER_ORDER_XYZ = 0, ///
    GM_EULER_ORDER_XZY = 1, ///
    GM_EULER_ORDER_YXZ = 2, ///
    GM_EULER_ORDER_YZX = 3, ///
    GM_EULER_ORDER_ZXY = 4, ///
    GM_EULER_ORDER_ZYX = 5, ///
}

alias Planes = int;
enum : Planes
{
    GM_PLANE_NEAR   = 0,
    GM_PLANE_FAR    = 1,
    GM_PLANE_LEFT   = 2,
    GM_PLANE_TOP    = 3,
    GM_PLANE_RIGHT  = 4,
    GM_PLANE_BOTTOM = 5
}









/* 
    ███╗   ███╗ █████╗ ████████╗██╗  ██╗███████╗██╗   ██╗███╗   ██╗ ██████╗███████╗
    ████╗ ████║██╔══██╗╚══██╔══╝██║  ██║██╔════╝██║   ██║████╗  ██║██╔════╝██╔════╝
    ██╔████╔██║███████║   ██║   ███████║█████╗  ██║   ██║██╔██╗ ██║██║     ███████╗
    ██║╚██╔╝██║██╔══██║   ██║   ██╔══██║██╔══╝  ██║   ██║██║╚██╗██║██║     ╚════██║
    ██║ ╚═╝ ██║██║  ██║   ██║   ██║  ██║██║     ╚██████╔╝██║ ╚████║╚██████╗███████║
*/
/*
Abstract over C math functions to get `pure` working.
The issue is that core.stdc.math is missing annotations, so 
this is a workaround, also provide float/double overloads.
*/

// See: https://docs.godotengine.org/en/stable/classes/class_%40globalscope.html
// Those are math functions defined at global scope in Godot Math, so here they
// are prefixed with a gm_*** prefix to avoid polluting global namespace.



// Godot math constants
// TODO: use the better constants from std.math, like PI with hex float precision**
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

// PRECISE_MATH_CHECKS is considered defined in this D translation
enum double GM_CMP_EPSILON = 0.00001;
enum double GM_CMP_EPSILON2 = (GM_CMP_EPSILON * GM_CMP_EPSILON);

// PRECISE_MATH_CHECKS
enum GM_UNIT_EPSILON = 0.00001;

float gm_abs(float x)   => libc.fabsf(x); ///
double gm_abs(double x) => libc.fabs(x); ///
int gm_abs(int x) => (x >= 0 ? x : -x); ///
float  gm_acos(float x)  @trusted => (x < -1) ? GM_PI : (x > 1 ? 0 : assumePureNothrowNogc(&libc.acosf, x)); ///
double gm_acos(double x)  => (x < -1) ? GM_PI : (x > 1 ? 0 : assumePureNothrowNogc(&libc.acos, x)); ///
float  gm_acosh(float x)  => assumePureNothrowNogc(&libc.acoshf, x); ///
double gm_acosh(double x) => assumePureNothrowNogc(&libc.acosh, x); ///

float gm_angle_difference(float from, float to) 
{
    float difference = gm_fmod(to - from, cast(float)GM_TAU);
    return gm_fmod(2.0f * difference, cast(float)GM_TAU) - difference;
}

double gm_angle_difference(double from, double to) 
{
    double difference = gm_fmod(to - from, GM_TAU);
    return gm_fmod(2.0 * difference, GM_TAU) - difference;
}

float  gm_asin(float x)   => x < -1 ? (-GM_PI / 2) : (x > 1 ? (GM_PI / 2) : assumePureNothrowNogc(&libc.asinf, x)); ///
double gm_asin(double x)  => x < -1 ? (-GM_PI / 2) : (x > 1 ? (GM_PI / 2) : assumePureNothrowNogc(&libc.asin, x)); ///
float  gm_asinh(float x)  => assumePureNothrowNogc(&libc.asinhf, x); ///
double gm_asinh(double x) => assumePureNothrowNogc(&libc.asinh, x); ///
float  gm_atan(float x)   => libc.atanf(x); //assumePureNothrowNogc(&libc.atanf, x); ///
double  gm_atan(double x)   => libc.atan(x); //assumePureNothrowNogc(&libc.atan, x); ///
float  gm_atan2(float y, float x)   => assumePureNothrowNogc(&libc.atan2f, y, x); ///
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

float gm_deg_to_rad(float y) => y * (cast(float)GM_PI / 180.0f); ///
double gm_deg_to_rad(double y) => y * (GM_PI / 180.0); ///

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

bool gm_is_zero_approx(float  value) => gm_abs(value) < cast(float)GM_CMP_EPSILON;  ///
bool gm_is_zero_approx(double value) => gm_abs(value) < GM_CMP_EPSILON;  ///

float gm_lerp(float from, float to, float weight) ///
{
    return from + (to - from) * weight;
}

double gm_lerp(double from, double to, double weight) ///
{
    return from + (to - from) * weight;
}

float gm_lerp_angle(float from, float to, float weight) => from + gm_angle_difference(from, to) * weight; ///
double gm_lerp_angle(double from, double to, double weight) => from + gm_angle_difference(from, to) * weight; ///

float gm_rad_to_deg(float y) => y * (180.0f / cast(float)GM_PI); ///
double gm_rad_to_deg(double y) => y * (180.0 / GM_PI); ///

float  gm_round(float x)  => libc.roundf(x); ///
double gm_round(double x) => libc.round(x); ///

int    gm_sign(int v)    => (v > 0) ? 1 : (v < 0 ? -1 : 0); ///
float  gm_sign(float v)  => (v > 0) ? 1.0f : (v < 0 ? -1.0f : 0.0f); ///
double gm_sign(double v) => (v > 0) ? 1.0 : (v < 0 ? -1.0 : 0.0); ///

bool gm_signbit(float num) => num < 0;  ///
bool gm_signbit(double num) => num < 0;  ///

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
    return cast(int) gm_snapped(cast(double)value, cast(double)step);
}

float  gm_sqrt(float x)  => assumePureNothrowNogc(&libc.sqrtf, x); ///
double gm_sqrt(double x) => assumePureNothrowNogc(&libc.sqrt, x); ///

void gm_swap(T)(ref T a, ref T b)
{
    T tmp = a;
    a = b;
    b = tmp;
}

float gm_tan(float x)   => libc.tanf(x);
double gm_tan(double x) => libc.tan(x);









/* 
    ██╗   ██╗███████╗ ██████╗████████╗ ██████╗ ██████╗ ██████╗ 
    ██║   ██║██╔════╝██╔════╝╚══██╔══╝██╔═══██╗██╔══██╗╚════██╗
    ██║   ██║█████╗  ██║        ██║   ██║   ██║██████╔╝ █████╔╝
    ╚██╗ ██╔╝██╔══╝  ██║        ██║   ██║   ██║██╔══██╗██╔═══╝ 
     ╚████╔╝ ███████╗╚██████╗   ██║   ╚██████╔╝██║  ██║███████╗
      ╚═══╝  ╚══════╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚══════╝
*/
/// See_also: https://docs.godotengine.org/en/stable/classes/class_vector2.html
struct Vector2Impl(T) 
    if (is(T == int) || is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private
    {
        alias V = Vector2Impl!T;
        
        enum bool isInt = is(T == int);
        enum bool isFloat = is(T == float) || is(T == double);
        static if (isFloat)
            alias F = T;
        else 
            alias F = float;
        alias Elem = T;
    }

    union
    {
        struct
        {
            T x = 0;
            T y = 0;
            alias width = x;
            alias height = y;
        }
        T[2] array;
    }

    enum : int
    {
        AXIS_X,
        AXIS_Y,
    }

    enum V ZERO  = V( 0,  0);
    enum V ONE   = V( 1,  1);
    static if (isFloat)
        enum V INF   = V( GM_INF,  GM_INF);
    enum V LEFT  = V(-1,  0);
    enum V RIGHT = V( 1,  0);
    enum V UP    = V( 0, -1);
    enum V DOWN  = V( 0,  1);

    // All functions as defined at: 
    this(T x, T y) { this.x = x; this.y = y; }
    this(T[2] v) { this.x = v[0]; this.y = v[1]; }
    V abs() => V(gm_abs(x), gm_abs(y));

    static if (isFloat)
    {
        T angle() const => gm_atan2(y, x);
        T angle_to(const V to) const => gm_atan2(cross(to), dot(to));
        T angle_to_point(const V to) const => (to - this).angle();
    }
    
    F aspect() const => width / cast(F)height;

    static if (isFloat)
    {
        V bezier_derivative(const V c1, const V c2, const V end, T t) const 
           => V( gm_bezier_derivative(x, c1.x, c2.x, end.x, t), 
                 gm_bezier_derivative(y, c1.y, c2.y, end.y, t) );
        V bezier_interpolate(const V c1, const V c2, const V end, T t) const 
            => V( gm_bezier_interpolate(x, c1.x, c2.x, end.x, t), 
                  gm_bezier_interpolate(y, c1.y, c2.y, end.y, t) );
        V bounce(const V normal) const => -reflect(normal);
        V ceil() => V(gm_ceil(x), gm_ceil(y));
    }

    V clamp(const V min, const V max) const => V(gm_clamp(x, min.x, max.x), gm_clamp(y, min.y, max.y));
    V clamp(T min, T max) const => V(gm_clamp(x, min, max), gm_clamp(y, min, max));
    
    static if (isFloat)
    {
        T cross(const V other) const => x * other.y - y * other.x;
        V cubic_interpolate(const V b, const V pre_a, const V post_b, T weight) const
            => V( gm_cubic_interpolate(x, b.x, pre_a.x, post_b.x, weight),
                  gm_cubic_interpolate(y, b.y, pre_a.y, post_b.y, weight) );
        V cubic_interpolate_in_time(const V b, const V pre_a, const V post_b, T weight, T b_t, T pre_a_t, T post_b_t) const
            => V( gm_cubic_interpolate_in_time(x, b.x, pre_a.x, post_b.x, weight, b_t, pre_a_t, post_b_t),
                  gm_cubic_interpolate_in_time(y, b.y, pre_a.y, post_b.y, weight, b_t, pre_a_t, post_b_t) );
        V direction_to(const V to) const => (to - this).normalized();
    }

    T distance_squared_to(const V v) const => (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y);
    F distance_to(const V v) const => gm_sqrt(cast(F) distance_squared_to(v));

    static if (isFloat)
    {
        T dot(const V other) const => x * other.x + y * other.y;
        V floor() => V(gm_floor(x), gm_floor(y));
        static V from_angle(T angle) => V( gm_cos(angle), gm_sin(angle) );
        bool is_equal_approx(const V other) const => gm_is_equal_approx(x, other.x) && gm_is_equal_approx(y, other.y);
        bool is_finite() const => gm_is_finite(x) && gm_is_finite(y);
        bool is_normalized() const => gm_is_equal_approx(length_squared(), cast(T)1, cast(T)GM_UNIT_EPSILON);
        bool is_zero_approx() const => gm_is_zero_approx(x) && gm_is_zero_approx(y);
    }
    F length() const => gm_sqrt(cast(F) length_squared());
    T length_squared() const => x * x + y * y;
    
    static if (isFloat)
    {
        V lerp(const V to, T weight) const => V( gm_lerp(x, to.x, weight), gm_lerp(y, to.y, weight) );
        V limit_length(T len) const 
        {
            T l = length();
            V v = this;
            if (l > 0 && len < l) 
            {
                v /= l;
                v *= len;
            }
            return v;
        }
    }

    V max(const V other) const => V( x > other.x ? x : other.x, y > other.y ? y : other.y );
    int max_axis_index() const => x < y ? AXIS_Y : AXIS_X;
    V max(T v) const => V( x > v ? x : v, y > v ? y : v );    
    V min(const V other) const => V( x < other.x ? x : other.x, y < other.y ? y : other.y );
    int min_axis_index() const => x < y ? AXIS_X : AXIS_Y;
    V min(T v) const => V( x < v ? x : v, y < v ? y : v );

    static if (isFloat)
    {
        V move_toward(const V to, T delta) const 
        {
            V v = this;
            V vd = to - v;
            T len = vd.length();
            return len <= delta || len < cast(T)GM_CMP_EPSILON ? to : v + vd / len * delta;
        }
        void normalize() // #BONUS
        {
            T l = x * x + y * y;
            if (l != 0) 
            {
                l = gm_sqrt(l);
                x /= l;
                y /= l;
            }
        }
        V normalized() const 
        {
            V v = this;
            v.normalize();
            return v;
        }
        V orthogonal() const => V(y, -x);
        V posmod(T mod) const => V(gm_fposmod(x, mod), gm_fposmod(y, mod));
        V posmodv(const V modv) const => V(gm_fposmod(x, modv.x), gm_fposmod(y, modv.y));
        V project(const V to) const => to * (dot(to) / to.length_squared());
        V reflect(const V normal) const
        {
            assert(normal.is_normalized());
            return normal * 2 * dot(normal) - this;
        }
        V rotated(T by_radians) const 
        {
            T sine = gm_sin(by_radians);
            T cosi = gm_cos(by_radians);
            return V(x * cosi - y * sine,
                     x * sine + y * cosi);
        }
        V round() const => V(gm_round(x), gm_round(y));
    }

    V sign() const  => V(gm_sign(x), gm_sign(y));

    static if (isFloat)
    {
        V slerp(const V to, T weight) const 
        {
            T start_length_sq = length_squared();
            T end_length_sq = to.length_squared();
            if (start_length_sq == 0 || end_length_sq == 0) 
            {
                // Zero length vectors have no angle, so the best we can do is either lerp or throw an error.
                return lerp(to, weight);
            }
            T start_length = gm_sqrt(start_length_sq);
            T result_length = gm_lerp(start_length, gm_sqrt(end_length_sq), weight);
            T angle = angle_to(to);
            return rotated(angle * weight) * (result_length / start_length);
        }
    }

    static if (isFloat)
        V slide(const V normal) const => this - normal * dot(normal);
    V snapped(const V step) const => V(gm_snapped(x, step.x), gm_snapped(y, step.y));
    V snapped(T step) const => V(gm_snapped(x, step), gm_snapped(y, step));

    V withX(T newX) const => V(newX, y);
    V withY(T newY) const => V(x, newY);

    // operators
    ref inout(T) opIndex(size_t n) inout return => array[n];
    size_t opDollar() => 2;
    inout(T)[] opSlice() inout return => array[];
    inout(T)[] opSlice(size_t a, size_t b) inout return => array[a..b];

    inout(T)* ptr() inout return => array.ptr;

    bool opEquals(V v) const => (x == v.x) && (y == v.y);

    int opCmp(V other) const 
    {
        if (x == other.x)
        {
            if (y == other.y) 
                return 0;
            else 
                return y < other.y ? -1 : 1;
        }
        else
        {
            return x < other.x ? -1 : 1;
        }
    }

    U opCast(U)() const if (isVector2Impl!U)
    {    
        static if (is(U.Elem == float))
            return U(cast(float)x, cast(float)y);
        else static if (is(U.Elem == double))
            return U(cast(double)x, cast(double)y);
        else static if (is(U.Elem == int))
            return U(cast(int)x, cast(int)y);
        else
            static assert(false);
    }

    V opBinary(string op)(const V v) const if (op == "*") => V(x * v.x  , y * v.y  );
    V opBinary(string op)(T scale)   const if (op == "*") => V(x * scale, y * scale);
    V opBinary(string op)(const V v) const if (op == "+") => V(x + v.x  , y + v.y  );
    V opBinary(string op)(T add)     const if (op == "+") => V(x + add  , y + add  );
    V opBinary(string op)(const V v) const if (op == "-") => V(x - v.x  , y - v.y  );
    V opBinary(string op)(T sub)     const if (op == "-") => V(x - sub  , y - sub  );
    V opBinary(string op)(const V v) const if (op == "/") => V(x / v.x  , y / v.y  );
    V opBinary(string op)(T scale)   const if (op == "/") => V(x / scale, y / scale);
    V opBinary(string op)(const V v) const if (op == "%") => V(x % v.x  , y % v.y  );
    V opBinary(string op)(T mod)     const if (op == "%") => V(x % mod  , y % mod  );

    V opBinaryRight(string op)(T scale) const if (op == "*") => V(scale * x, scale * y);
    V opBinaryRight(string op)(T add)   const if (op == "+") => V(add + x  , add + y  );
    V opBinaryRight(string op)(T sub)   const if (op == "-") => V(sub - x  , sub - y  );
    V opBinaryRight(string op)(T scale) const if (op == "/") => V(scale / x, scale / y);
    V opBinaryRight(string op)(T mod)   const if (op == "%") => V(mod % x  , mod % y  );

    static if (isFloat)
        V opBinary(string op)(const Transform2DImpl!T transform) const if (op == "*") => transform.inverse() * this;

    V opOpAssign(string op)(const V v) if (op == "*") { x *= v.x;   y *= v.y;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "*") { x *= scale; y *= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "+") { x += v.x;   y += v.y;   return this; }
    V opOpAssign(string op)(T add)     if (op == "+") { x += add;   y += add;   return this; }
    V opOpAssign(string op)(const V v) if (op == "-") { x -= v.x;   y -= v.y;   return this; }
    V opOpAssign(string op)(T sub)     if (op == "-") { x -= sub;   y -= sub;   return this; }
    V opOpAssign(string op)(const V v) if (op == "/") { x /= v.x;   y /= v.y;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "/") { x /= scale; y /= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "%") { x %= v.x;   y %= v.y;   return this; }
    V opOpAssign(string op)(T mod)     if (op == "%") { x %= mod;   y %= mod;   return this; }
    
    V opUnary(string op)() const if (op == "+") => this;    
    V opUnary(string op)() const if (op == "-") => V(-x, -y);
}










/*
    ██╗   ██╗███████╗ ██████╗████████╗ ██████╗ ██████╗ ██████╗ 
    ██║   ██║██╔════╝██╔════╝╚══██╔══╝██╔═══██╗██╔══██╗╚════██╗
    ██║   ██║█████╗  ██║        ██║   ██║   ██║██████╔╝ █████╔╝
    ╚██╗ ██╔╝██╔══╝  ██║        ██║   ██║   ██║██╔══██╗ ╚═══██╗
     ╚████╔╝ ███████╗╚██████╗   ██║   ╚██████╔╝██║  ██║██████╔╝
      ╚═══╝  ╚══════╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═════╝ 
*/
/// See_also: https://docs.godotengine.org/en/stable/classes/class_vector3.html
struct Vector3Impl(T) 
    if (is(T == int) || is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private 
    {
        alias V = Vector3Impl!T,
              V2 = Vector2Impl!T;

        enum bool isInt = is(T == int);
        enum bool isFloat = is(T == float) || is(T == double);
        static if (isFloat)        
            alias F = T;
        else 
            alias F = float;
        alias Elem = T;
    }

    union
    {
        struct
        {
            T x = 0;
            T y = 0;
            T z = 0;
            alias width  = x;
            alias height = y;
            alias depth  = z;
        }
        T[3] array;
    }

    enum : int
    {
        AXIS_X,
        AXIS_Y,
        AXIS_Z,
    }

    // Common 3D direction vectors (assuming Y-up and Z-back/forward convention)
    enum V ZERO    = V( 0,  0,  0);
    enum V ONE     = V( 1,  1,  1);
    static if (isFloat)
        enum V INF = V( GM_INF, GM_INF, GM_INF);
    enum V LEFT    = V(-1,  0,  0);
    enum V RIGHT   = V( 1,  0,  0);
    enum V UP      = V( 0,  1,  0);
    enum V DOWN    = V( 0, -1,  0);
    enum V FORWARD = V( 0,  0, -1);
    enum V BACK    = V( 0,  0,  1);

    this(T x, T y, T z) { this.x = x; this.y = y; this.z = z; }
    this(T[3] v) { this.x = v[0]; this.y = v[1]; this.z = v[2];}
    V abs() const => V(gm_abs(x), gm_abs(y), gm_abs(z));

    static if (isFloat)
    {
        T angle_to(const V to) const => gm_atan2(cross(to).length(), dot(to));
        V bezier_derivative(const V c1, const V c2, const V end, T t) const
            => V( gm_bezier_derivative(x, c1.x, c2.x, end.x, t),
                  gm_bezier_derivative(y, c1.y, c2.y, end.y, t),
                  gm_bezier_derivative(z, c1.z, c2.z, end.z, t) );
        V bezier_interpolate(const V c1, const V c2, const V end, T t) const
            => V( gm_bezier_interpolate(x, c1.x, c2.x, end.x, t),
                  gm_bezier_interpolate(y, c1.y, c2.y, end.y, t),
                  gm_bezier_interpolate(z, c1.z, c2.z, end.z, t) );
        V bounce(const V normal) const => -reflect(normal);
        V ceil() const => V(gm_ceil(x), gm_ceil(y), gm_ceil(z));
    }
    V clamp(const V min, const V max) const 
        => V(gm_clamp(x, min.x, max.x), gm_clamp(y, min.y, max.y), gm_clamp(z, min.z, max.z));
    V clamp(T min, T max) const 
        => V(gm_clamp(x, min, max), gm_clamp(y, min, max), gm_clamp(z, min, max));
    static if (isFloat)
    {
        V cross(const V other) const
            => V( y * other.z - z * other.y,
                  z * other.x - x * other.z,
                  x * other.y - y * other.x );
        V cubic_interpolate(const V b, const V pre_a, const V post_b, T weight) const
            => V( gm_cubic_interpolate(x, b.x, pre_a.x, post_b.x, weight),
                  gm_cubic_interpolate(y, b.y, pre_a.y, post_b.y, weight),
                  gm_cubic_interpolate(z, b.z, pre_a.z, post_b.z, weight) );
        V cubic_interpolate_in_time(const V b, const V pre_a, const V post_b, T weight, T b_t, T pre_a_t, T post_b_t) const
            => V( gm_cubic_interpolate_in_time(x, b.x, pre_a.x, post_b.x, weight, b_t, pre_a_t, post_b_t),
                  gm_cubic_interpolate_in_time(y, b.y, pre_a.y, post_b.y, weight, b_t, pre_a_t, post_b_t),
                  gm_cubic_interpolate_in_time(z, b.z, pre_a.z, post_b.z, weight, b_t, pre_a_t, post_b_t) );
        V direction_to(const V to) const => (to - this).normalized();
    }
    T distance_squared_to(const V v) const 
        => (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z);
    F distance_to(const V v) const => gm_sqrt(cast(F) distance_squared_to(v));

    static if (isFloat)
    {
        T dot(const V other) const => x * other.x + y * other.y + z * other.z;
        V floor() const => V(gm_floor(x), gm_floor(y), gm_floor(z));
        V get_any_perpendicular()
        {
            assert(! is_zero_approx() );
            return cross((gm_abs(x) <= gm_abs(y) && gm_abs(x) <= gm_abs(z)) ? V.RIGHT : V.UP).normalized();
        }
        V inverse() const => V(1 / x, 1 / y, 1 / z);
        bool is_equal_approx(const V other) const
           => gm_is_equal_approx(x, other.x) && gm_is_equal_approx(y, other.y) && gm_is_equal_approx(z, other.z);
        bool is_finite() const => gm_is_finite(x) && gm_is_finite(y) && gm_is_finite(z);
        bool is_normalized() const 
            => gm_is_equal_approx(length_squared(), cast(T)1, cast(T)GM_UNIT_EPSILON);
        bool is_zero_approx() const 
            => gm_is_zero_approx(x) && gm_is_zero_approx(y) && gm_is_zero_approx(z);
    }

    F length() const => gm_sqrt(cast(F) length_squared());
    T length_squared() const => x * x + y * y + z * z;

    static if (isFloat)
    {
        V lerp(const V to, T weight) const 
            => V( gm_lerp(x, to.x, weight), gm_lerp(y, to.y, weight), gm_lerp(z, to.z, weight) );
        V limit_length(T len) const
        {
            T l = length();
            V v = this;
            if (l > 0 && len < l)
            {
                v /= l;
                v *= len;
            }
            return v;
        }
    }

    V max(const V other) const 
        => V( x > other.x ? x : other.x, y > other.y ? y : other.y, z > other.z ? z : other.z );
    int max_axis_index() const
    {
        if (x > y && x > z) return AXIS_X;
        if (y > z) return AXIS_Y;
        return AXIS_Z;
    }
    V max(T v) const 
        => V( x > v ? x : v, y > v ? y : v, z > v ? z : v );
    V min(const V other) const 
        => V( x < other.x ? x : other.x, y < other.y ? y : other.y, z < other.z ? z : other.z );
    int min_axis_index() const
    {
        if (x < y && x < z) return AXIS_X;
        if (y < z) return AXIS_Y;
        return AXIS_Z;
    }
    V min(T v) const 
        => V( x < v ? x : v, y < v ? y : v, z < v ? z : v );

    static if (isFloat)
    {
        V move_toward(const V to, T delta) const
        {
            V v = this;
            V vd = to - v;
            T len = vd.length();
            return len <= delta || len < cast(T)GM_CMP_EPSILON ? to : v + vd / len * delta;
        }
        void normalize() // #BONUS
        {
            T l = x * x + y * y + z * z;
            if (l != 0)
            {
                l = gm_sqrt(l);
                x /= l;
                y /= l;
                z /= l;
            }
        }
        V normalized() const
        {
            V v = this;
            v.normalize();
            return v;
        }

        static V octahedron_decode(const V2 oct) 
        {
            V2 f = V2(oct.x * 2 - 1, oct.y * 2 - 1);
            V n = V(f.x, f.y, 1 - gm_abs(f.x) - gm_abs(f.y));
            T t = gm_clamp(-n.z, cast(T)0, cast(T)1);
            n.x += n.x >= 0 ? -t : t;
            n.y += n.y >= 0 ? -t : t;
            return n.normalized();
        }

        V2 octahedron_encode() const
        {
            V n = this;
            n /= gm_abs(n.x) + gm_abs(n.y) + gm_abs(n.z);
            V2 o;
            if (n.z >= 0) 
            {
                o.x = n.x;
                o.y = n.y;
            } 
            else 
            {
                o.x = (1 - gm_abs(n.y)) * (n.x >= 0 ? 1 : -1);
                o.y = (1 - gm_abs(n.x)) * (n.y >= 0 ? 1 : -1);
            }
            o.x = o.x * 0.5f + 0.5f;
            o.y = o.y * 0.5f + 0.5f;
            return o;
        }
    
        BasisImpl!T outer(const V with_) const 
        {
            BasisImpl!T basis;
            basis.rows[0] = V(x * with_.x, x * with_.y, x * with_.z);
            basis.rows[1] = V(y * with_.x, y * with_.y, y * with_.z);
            basis.rows[2] = V(z * with_.x, z * with_.y, z * with_.z);
            return basis;
        }

        V posmod(T mod) const => V(gm_fposmod(x, mod), gm_fposmod(y, mod), gm_fposmod(z, mod));
        V posmodv(const V modv) const => V(gm_fposmod(x, modv.x), gm_fposmod(y, modv.y), gm_fposmod(z, modv.z));
        V project(const V to) const => to * (dot(to) / to.length_squared());

        V reflect(const V normal) const
        {
            assert(normal.is_normalized());
            return normal * 2 * dot(normal) - this;
        }

        void rotate(const V axis, T angle)
        {
            T axis_length_sq = axis.length_squared();
            if (gm_is_zero_approx(axis_length_sq))
                return;

            V norm_axis = axis / gm_sqrt(axis_length_sq);

            T sin_angle = gm_sin(angle);
            T cos_angle = gm_cos(angle);

            V term1 = this * cos_angle;
            V term2 = norm_axis.cross(this) * sin_angle;
            V term3 = norm_axis * (norm_axis.dot(this) * (1.0f - cos_angle));

            this = term1 + term2 + term3;
        }

        V rotated(const V axis, T angle) const
        {
            V v = this;
            v.rotate(axis, angle);
            return v;
        }

        V round() const => V(gm_round(x), gm_round(y), gm_round(z));
    }

    V sign() const => V(gm_sign(x), gm_sign(y), gm_sign(z));

    static if (isFloat)
    {
        T signed_angle_to(const V to, const V axis) const
        {
            V cross_to = cross(to);
            T unsigned_angle = gm_atan2(cross_to.length(), dot(to));
            T sign = cross_to.dot(axis);
            return (sign < 0) ? -unsigned_angle : unsigned_angle;
        }

        V slerp(const V to, T weight) const 
        {
            T start_length_sq = length_squared();
            T end_length_sq = to.length_squared();
            if (start_length_sq == 0 || end_length_sq == 0) 
                return lerp(to, weight);

            V axis = cross(to);
            T axis_length_sq = axis.length_squared();
            if (axis_length_sq == 0) 
                return lerp(to, weight);

            axis /= gm_sqrt(axis_length_sq);
            T start_length = gm_sqrt(start_length_sq);
            T result_length = gm_lerp(start_length, gm_sqrt(end_length_sq), weight);
            T angle = angle_to(to);
            return rotated(axis, angle * weight) * (result_length / start_length);
        }

        V slide(const V normal) const => this - normal * dot(normal);
    }

    V snapped(const(V) step) const => V(gm_snapped(x, step.x), gm_snapped(y, step.y),  gm_snapped(z, step.z));
    V snapped(T step) const => V(gm_snapped(x, step), gm_snapped(y, step), gm_snapped(z, step));

    V withX(T newX) const => V(newX, y, z);
    V withY(T newY) const => V(x, newY, z);
    V withZ(T newZ) const => V(x, y, newZ);

    // operators
    ref inout(T) opIndex(size_t n) inout return => array[n];
    size_t opDollar() => 3;
    inout(T)[] opSlice() inout return => array[];
    inout(T)[] opSlice(size_t a, size_t b) inout return => array[a..b];

    inout(T)* ptr() inout return => array.ptr;

    bool opEquals(V v) const => (x == v.x) && (y == v.y) && (z == v.z);

    int opCmp(V other) const 
    {
        if (x == other.x)
        {
            if (y == other.y)
            {
                if (z == other.z) 
                    return 0;
                else
                    return z < other.z ? -1 : 1;
            }            
            else 
                return y < other.y ? -1 : 1;
        }
        else
            return x < other.x ? -1 : 1;
    }

    U opCast(U)() const if (isVector3Impl!U)
    {
        static if (is(U.Elem == float))
            return U(cast(float)x, cast(float)y, cast(float)z);
        else static if (is(U.Elem == double))
            return U(cast(double)x, cast(double)y, cast(double)z);
        else static if (is(U.Elem == int))
            return U(cast(int)x, cast(int)y, cast(int)z );
        else
            static assert(false);
    }

    V opBinary(string op)(const V v) const if (op == "*") => V(x * v.x  , y * v.y  , z * v.z  );
    V opBinary(string op)(T scale)   const if (op == "*") => V(x * scale, y * scale, z * scale);
    V opBinary(string op)(const V v) const if (op == "+") => V(x + v.x  , y + v.y  , z + v.z  );
    V opBinary(string op)(T add)     const if (op == "+") => V(x + add  , y + add  , z + add  );
    V opBinary(string op)(const V v) const if (op == "-") => V(x - v.x  , y - v.y  , z - v.z  );
    V opBinary(string op)(T sub)     const if (op == "-") => V(x - sub  , y - sub  , z - sub  );
    V opBinary(string op)(const V v) const if (op == "/") => V(x / v.x  , y / v.y  , z / v.z  );
    V opBinary(string op)(T scale)   const if (op == "/") => V(x / scale, y / scale, z / scale);
    V opBinary(string op)(const V v) const if (op == "%") => V(x % v.x  , y % v.y  , z % v.z  );
    V opBinary(string op)(T mod)     const if (op == "%") => V(x % mod  , y % mod  , z % mod  );

    V opBinaryRight(string op)(T scale) const if (op == "*") => V(scale * x, scale * y, scale * z);
    V opBinaryRight(string op)(T add)   const if (op == "+") => V(add + x  , add + y  , add + z  );
    V opBinaryRight(string op)(T sub)   const if (op == "-") => V(sub - x  , sub - y  , sub - z  );
    V opBinaryRight(string op)(T scale) const if (op == "/") => V(scale / x, scale / y, scale / z);
    V opBinaryRight(string op)(T mod)   const if (op == "%") => V(mod % x  , mod % y  , mod % z  );

    static if (isFloat)
    {
        V opBinary(string op)(const BasisImpl!T basis) const if (op == "*") => basis.transposed() * this;
        V opBinary(string op)(const QuaternionImpl!T quaternion) const if (op == "*") => quaternion.inverse() * this;
    }

    V opOpAssign(string op)(const V v) if (op == "*") { x *= v.x;   y *= v.y;   z *= v.z;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "*") { x *= scale; y *= scale; z *= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "+") { x += v.x;   y += v.y;   z += v.z;   return this; }
    V opOpAssign(string op)(T add)     if (op == "+") { x += add;   y += add;   z += add;   return this; }
    V opOpAssign(string op)(const V v) if (op == "-") { x -= v.x;   y -= v.y;   z -= v.z;   return this; }
    V opOpAssign(string op)(T sub)     if (op == "-") { x -= sub;   y -= sub;   z -= sub;   return this; }
    V opOpAssign(string op)(const V v) if (op == "/") { x /= v.x;   y /= v.y;   z /= v.z;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "/") { x /= scale; y /= scale; z /= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "%") { x %= v.x;   y %= v.y;   z %= v.z;   return this; }
    V opOpAssign(string op)(T mod)     if (op == "%") { x %= mod;   y %= mod;   z %= mod;   return this; }
    
    V opUnary(string op)() const if (op == "+") => this;
    V opUnary(string op)() const if (op == "-") => V(-x, -y, -z);
}










/*
██╗   ██╗███████╗ ██████╗████████╗ ██████╗ ██████╗ ██╗  ██╗
██║   ██║██╔════╝██╔════╝╚══██╔══╝██╔═══██╗██╔══██╗██║  ██║
██║   ██║█████╗  ██║        ██║   ██║   ██║██████╔╝███████║
╚██╗ ██╔╝██╔══╝  ██║        ██║   ██║   ██║██╔══██╗╚════██║
 ╚████╔╝ ███████╗╚██████╗   ██║   ╚██████╔╝██║  ██║     ██║
  ╚═══╝  ╚══════╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝     ╚═╝
*/
/// See also: 
struct Vector4Impl(T) 
    if (is(T == int) || is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private 
    {
        alias V = Vector4Impl!T,
              V3 = Vector3Impl!T;

        enum bool isFloat = (is(T == float) || is(T == double));
        static if (isFloat)
            alias F = T;
        else 
            alias F = float;
        alias Elem = T;
    }

    union
    {
        struct
        {
            T x = 0;
            T y = 0;
            T z = 0;
            T w = 0;
            alias width  = x;
            alias height = y;
            alias depth  = z;
        }
        T[4] array;
    }

    enum : int
    {
        AXIS_X,
        AXIS_Y,
        AXIS_Z,
        AXIS_W,
    }

    enum V ZERO    = V( 0,  0,  0,  0);
    enum V ONE     = V( 1,  1,  1,  1);
    static if (isFloat)
        enum V INF = V( GM_INF, GM_INF, GM_INF, GM_INF);

    this(T x, T y, T z, T w) { this.x = x; this.y = y; this.z = z; this.w = w; }
    this(T[4] v) { this.x = v[0]; this.y = v[1]; this.z = v[2]; this.w = v[3]; }
    V abs() const => V(gm_abs(x), gm_abs(y), gm_abs(z), gm_abs(w));

    static if (isFloat)
        V ceil() const => V(gm_ceil(x), gm_ceil(y), gm_ceil(z), gm_ceil(w));
    
    V clamp(const V min, const V max) const 
        => V(gm_clamp(x, min.x, max.x), gm_clamp(y, min.y, max.y), gm_clamp(z, min.z, max.z), gm_clamp(w, min.w, max.w));
    V clamp(T min, T max) const 
        => V(gm_clamp(x, min, max), gm_clamp(y, min, max), gm_clamp(z, min, max), gm_clamp(w, min, max));
    static if (isFloat)
    {
        V cubic_interpolate(const V b, const V pre_a, const V post_b, T weight) const
            => V( gm_cubic_interpolate(x, b.x, pre_a.x, post_b.x, weight),
                  gm_cubic_interpolate(y, b.y, pre_a.y, post_b.y, weight),
                  gm_cubic_interpolate(z, b.z, pre_a.z, post_b.z, weight),
                  gm_cubic_interpolate(w, b.w, pre_a.w, post_b.w, weight) );
        V cubic_interpolate_in_time(const V b, const V pre_a, const V post_b, T weight, T b_t, T pre_a_t, T post_b_t) const
            => V( gm_cubic_interpolate_in_time(x, b.x, pre_a.x, post_b.x, weight, b_t, pre_a_t, post_b_t),
                  gm_cubic_interpolate_in_time(y, b.y, pre_a.y, post_b.y, weight, b_t, pre_a_t, post_b_t),
                  gm_cubic_interpolate_in_time(z, b.z, pre_a.z, post_b.z, weight, b_t, pre_a_t, post_b_t),
                  gm_cubic_interpolate_in_time(w, b.w, pre_a.w, post_b.w, weight, b_t, pre_a_t, post_b_t) );

        V direction_to (const V v) const => (v - this).normalized();
    }
    T distance_squared_to(const V v) const 
        => (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z) + (w - v.w) * (w - v.w);
    F distance_to(const V v) const => gm_sqrt(cast(F) distance_squared_to(v));

    static if (isFloat)
    {
        T dot(const V other) const => x * other.x + y * other.y + z * other.z + w * other.w;
        V floor() const => V(gm_floor(x), gm_floor(y), gm_floor(z), gm_floor(w));
        V inverse() const => V(1 / x, 1 / y, 1 / z, 1 / w);
        bool is_equal_approx(const V other) const
           => gm_is_equal_approx(x, other.x) && gm_is_equal_approx(y, other.y) && gm_is_equal_approx(z, other.z) && gm_is_equal_approx(w, other.w);
        bool is_finite() const => gm_is_finite(x) && gm_is_finite(y) && gm_is_finite(z) && gm_is_finite(w);
        bool is_normalized() const 
            => gm_is_equal_approx(length_squared(), cast(T)1, cast(T)GM_UNIT_EPSILON);
        bool is_zero_approx() const 
            => gm_is_zero_approx(x) && gm_is_zero_approx(y) && gm_is_zero_approx(z) && gm_is_zero_approx(w);
    }

    F length() const => gm_sqrt(cast(F) length_squared());
    T length_squared() const => x * x + y * y + z * z + w * w;

    static if (isFloat)
    {
        V lerp(const V to, T weight) const 
            => V( gm_lerp(x, to.x, weight), gm_lerp(y, to.y, weight), gm_lerp(z, to.z, weight), gm_lerp(w, to.w, weight) );
        V limit_length(T len) const
        {
            T l = length();
            V v = this;
            if (l > 0 && len < l)
            {
                v /= l;
                v *= len;
            }
            return v;
        }
    }

    V max(const V other) const 
        => V( x > other.x ? x : other.x, y > other.y ? y : other.y, z > other.z ? z : other.z, w > other.w ? w : other.w );
    int max_axis_index() const
    {
        if (x > y && x > z && x > w) return AXIS_X;
        if (y > z && y > w) return AXIS_Y;
        if (z > w) return AXIS_Z;
        return AXIS_W;
    }
    V max(T v) const 
        => V( x > v ? x : v, y > v ? y : v, z > v ? z : v, w > v ? w : v );
    V min(const V other) const 
        => V( x < other.x ? x : other.x, y < other.y ? y : other.y, z < other.z ? z : other.z, w < other.w ? w : other.w );
    int min_axis_index() const
    {
        if (x < y && x < z && x < w) return AXIS_X;
        if (y < z && y < w) return AXIS_Y;
        if (z < w) return AXIS_Z;
        return AXIS_W;
    }
    V min(T v) const 
        => V( x < v ? x : v, y < v ? y : v, z < v ? z : v, w < v ? w : v );

    static if (isFloat)
    {
        void normalize() // #BONUS
        {
            T l = x * x + y * y + z * z + w * w;
            if (l != 0)
            {
                l = gm_sqrt(l);
                x /= l;
                y /= l;
                z /= l;
                w /= l;
            }
        }
        V normalized() const
        {
            V v = this;
            v.normalize();
            return v;
        }

        V posmod(T mod) const => V(gm_fposmod(x, mod), gm_fposmod(y, mod), gm_fposmod(z, mod), gm_fposmod(w, mod));
        V posmodv(const V modv) const => V(gm_fposmod(x, modv.x), gm_fposmod(y, modv.y), gm_fposmod(z, modv.z), gm_fposmod(w, modv.w));
        V round() const => V(gm_round(x), gm_round(y), gm_round(z), gm_round(w));
        V slide(const V normal) const => this - normal * dot(normal);
    }

    V sign() const => V(gm_sign(x), gm_sign(y), gm_sign(z), gm_sign(w));

    V snapped(const(V) step) const => V(gm_snapped(x, step.x), gm_snapped(y, step.y), gm_snapped(z, step.z), gm_snapped(w, step.w));
    V snapped(T step) const => V(gm_snapped(x, step), gm_snapped(y, step), gm_snapped(z, step), gm_snapped(w, step));

    V withX(T newX) const => V(newX, y, z, w);
    V withY(T newY) const => V(x, newY, z, w);
    V withZ(T newZ) const => V(x, y, newZ, w);
    V withW(T newW) const => V(x, y, z, newW);
 
    // operators
    ref inout(T) opIndex(size_t n) inout return => array[n];
    size_t opDollar() => 4;
    inout(T)[] opSlice() inout return => array[];
    inout(T)[] opSlice(size_t a, size_t b) inout return => array[a..b];

    inout(T)* ptr() inout return => array.ptr;

    bool opEquals(V v) const => (x == v.x) && (y == v.y) && (z == v.z) && (w == v.w);

    int opCmp(V other) const 
    {
        if (x == other.x)
        {
            if (y == other.y)
            {
                if (z == other.z)
                {
                    if (w == other.w)
                        return 0;
                    else
                        return w < other.w ? -1 : 1;
                }
                else
                    return z < other.z ? -1 : 1;
            }            
            else 
                return y < other.y ? -1 : 1;
        }
        else
            return x < other.x ? -1 : 1;
    }

    U opCast(U)() const if (isVector4Impl!U)
    {    
        static if (is(U.Elem == float))
            return U(cast(float)x, cast(float)y, cast(float)z, cast(float)w);
        else static if (is(U.Elem == double))
            return U(cast(double)x, cast(double)y, cast(double)z, cast(double)w);
        else static if (is(U.Elem == int))
            return U(cast(int)x, cast(int)y, cast(int)z, cast(int)w);
        else
            static assert(false);
    }

    V opBinary(string op)(const V v) const if (op == "*") => V(x * v.x  , y * v.y  , z * v.z  , w * v.w  );
    V opBinary(string op)(T scale)   const if (op == "*") => V(x * scale, y * scale, z * scale, w * scale);
    V opBinary(string op)(const V v) const if (op == "+") => V(x + v.x  , y + v.y  , z + v.z  , w + v.w  );
    V opBinary(string op)(T add)     const if (op == "+") => V(x + add  , y + add  , z + add  , w + add  );
    V opBinary(string op)(const V v) const if (op == "-") => V(x - v.x  , y - v.y  , z - v.z  , w - v.w  );
    V opBinary(string op)(T sub)     const if (op == "-") => V(x - sub  , y - sub  , z - sub  , w - sub  );
    V opBinary(string op)(const V v) const if (op == "/") => V(x / v.x  , y / v.y  , z / v.z  , w / v.w  );
    V opBinary(string op)(T scale)   const if (op == "/") => V(x / scale, y / scale, z / scale, w / scale);
    V opBinary(string op)(const V v) const if (op == "%") => V(x % v.x  , y % v.y  , z % v.z  , w % v.w  );
    V opBinary(string op)(T mod)     const if (op == "%") => V(x % mod  , y % mod  , z % mod  , w % mod  );

    V opBinaryRight(string op)(T scale) const if (op == "*") => V(scale * x, scale * y, scale * z, scale * w);
    V opBinaryRight(string op)(T add)   const if (op == "+") => V(add + x  , add + y  , add + z  , add + w  );
    V opBinaryRight(string op)(T sub)   const if (op == "-") => V(sub - x  , sub - y  , sub - z  , sub - w  );
    V opBinaryRight(string op)(T scale) const if (op == "/") => V(scale / x, scale / y, scale / z, scale / w);
    V opBinaryRight(string op)(T mod)   const if (op == "%") => V(mod % x  , mod % y  , mod % z  , mod % w  );

    V opOpAssign(string op)(const V v) if (op == "*") { x *= v.x;   y *= v.y;   z *= v.z;   w *= v.w;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "*") { x *= scale; y *= scale; z *= scale; w *= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "+") { x += v.x;   y += v.y;   z += v.z;   w += v.w;   return this; }
    V opOpAssign(string op)(T add)     if (op == "+") { x += add;   y += add;   z += add;   w += add;   return this; }
    V opOpAssign(string op)(const V v) if (op == "-") { x -= v.x;   y -= v.y;   z -= v.z;   w -= v.w;   return this; }
    V opOpAssign(string op)(T sub)     if (op == "-") { x -= sub;   y -= sub;   z -= sub;   w -= sub;   return this; }
    V opOpAssign(string op)(const V v) if (op == "/") { x /= v.x;   y /= v.y;   z /= v.z;   w /= v.w;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "/") { x /= scale; y /= scale; z /= scale; w /= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "%") { x %= v.x;   y %= v.y;   z %= v.z;   w %= v.w;   return this; }
    V opOpAssign(string op)(T mod)     if (op == "%") { x %= mod;   y %= mod;   z %= mod;   w %= mod;   return this; }

    V opUnary(string op)() const if (op == "+") => this;
    V opUnary(string op)() const if (op == "-") => V(-x, -y, -z, -w);
}









/*
      ██████╗ ██╗   ██╗ █████╗ ████████╗███████╗██████╗ ███╗   ██╗██╗ ██████╗ ███╗   ██╗
     ██╔═══██╗██║   ██║██╔══██╗╚══██╔══╝██╔════╝██╔══██╗████╗  ██║██║██╔═══██╗████╗  ██║
     ██║   ██║██║   ██║███████║   ██║   █████╗  ██████╔╝██╔██╗ ██║██║██║   ██║██╔██╗ ██║
     ██║▄▄ ██║██║   ██║██╔══██║   ██║   ██╔══╝  ██╔══██╗██║╚██╗██║██║██║   ██║██║╚██╗██║
     ╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████╗██║  ██║██║ ╚████║██║╚██████╔╝██║ ╚████║
      ╚══▀▀═╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝
*/
struct QuaternionImpl(T)
    if (is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private
    {
        alias Q = QuaternionImpl!T;
        alias V3 = Vector3Impl!T;
        alias Elem = T;
    }

    enum Q IDENTITY = Q.init;

    union 
    {
        struct 
        {
            T x = 0;
            T y = 0;
            T z = 0;
            T w = 1;
        }
        T[4] array;
    }

    this(V3 arc_from, V3 arc_to) // shortest arc
    {
        alias v0 = arc_from;
        alias v1 = arc_to;
        assert(!v0.is_zero_approx() && !v1.is_zero_approx());

        static if (is(T == double))
            enum T ALMOST_ONE = 0.999999999999999;
        else
            enum T ALMOST_ONE = 0.99999975f;

        V3 n0 = v0.normalized();
        V3 n1 = v1.normalized();
        T d = n0.dot(n1);
        if (gm_abs(d) > ALMOST_ONE) 
        {
            if (d >= 0)
                return; // Vectors are same.

            V3 axis = n0.get_any_perpendicular();
            x = axis.x;
            y = axis.y;
            z = axis.z;
            w = 0;
        } 
        else 
        {
            V3 c = n0.cross(n1);
            T s = gm_sqrt((1.0f + d) * 2.0f);
            T rs = 1.0f / s;
            x = c.x * rs;
            y = c.y * rs;
            z = c.z * rs;
            w = s * 0.5f;
        }
        normalize();
    }

    this(const V3 axis, T angle)
    {
        assert(axis.is_normalized, "Quaternion axis must be normalized when constructed");

        T d = axis.length();
        if (d == 0) 
        {
            x = 0;
            y = 0;
            z = 0;
            w = 0;
        } 
        else 
        {
            T sin_angle = gm_sin(angle * 0.5f);
            T cos_angle = gm_cos(angle * 0.5f);
            T s = sin_angle / d;
            x = axis.x * s;
            y = axis.y * s;
            z = axis.z * s;
            w = cos_angle;
        }
    }

    this(BasisImpl!T from)
    {
        // Note: original source has implicit conversion from Basis to Quaternion
        this = from.get_quaternion();
    }

    this(T x, T y, T z, T w)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    T angle_to(const Q to) const 
    {
        T d = dot(to);
        // acos does clamping.
        return gm_acos(d * d * 2 - 1);
    }
    T dot(const Q q) const => x * q.x + y * q.y + z * q.z + w * q.w;
    Q exp() const 
    {
        Q src = this;
        V3 src_v = V3(src.x, src.y, src.z);
        T theta = src_v.length();
        src_v = src_v.normalized();
        if (theta < GM_CMP_EPSILON || !src_v.is_normalized()) {
            return Q(0, 0, 0, 1);
        }
        return Q(src_v, theta);
    }

    // This implementation uses YXZ convention (Z is the first rotation).
    static Q from_euler(const V3 euler) 
    {
        T half_a1 = euler.y * 0.5f;
        T half_a2 = euler.x * 0.5f;
        T half_a3 = euler.z * 0.5f;

        T cos_a1 = gm_cos(half_a1);
        T sin_a1 = gm_sin(half_a1);
        T cos_a2 = gm_cos(half_a2);
        T sin_a2 = gm_sin(half_a2);
        T cos_a3 = gm_cos(half_a3);
        T sin_a3 = gm_sin(half_a3);

        return Q(sin_a1 * cos_a2 * sin_a3 + cos_a1 * sin_a2 * cos_a3,
                 sin_a1 * cos_a2 * cos_a3 - cos_a1 * sin_a2 * sin_a3,
                -sin_a1 * sin_a2 * cos_a3 + cos_a1 * cos_a2 * sin_a3,
                 sin_a1 * sin_a2 * sin_a3 + cos_a1 * cos_a2 * cos_a3);
    }

    T get_angle() const => 2 * gm_acos(w);
    V3 get_axis() const 
    {
        if (gm_abs(w) > 1 - GM_CMP_EPSILON) 
        {
            return V3(x, y, z);
        }
        T r = (cast(T)1) / gm_sqrt(1 - w * w);
        return V3(x * r, y * r, z * r);
    }

    V3 get_euler(EulerOrder order = GM_EULER_ORDER_YXZ) const 
    {
        assert(is_normalized());
        return BasisImpl!T(this).get_euler(order);
    }

    Q inverse() const => Q(-x, -y, -z, w);

    bool is_equal_approx(const Q to) const => gm_is_equal_approx(x, to.x) && gm_is_equal_approx(y, to.y) && gm_is_equal_approx(z, to.z) && gm_is_equal_approx(w, to.w);
    bool is_finite() const => gm_is_finite(x) && gm_is_finite(y) && gm_is_finite(z) && gm_is_finite(w);
    bool is_normalized() const => gm_is_equal_approx(length_squared(), cast(T)1, cast(T)GM_UNIT_EPSILON);

    T length() const => gm_sqrt(length_squared());
    T length_squared() const => dot(this);
    Q log() const
    {
        V3 src_v = get_axis() * get_angle();
        return Q(src_v.x, src_v.y, src_v.z, 0);
    }

    void normalize() // #BONUS
    {
        this /= length();
    }

    Q normalized() const
    {
        Q r = this;
        r.normalize();
        return r;
    }

    inout(T*) ptr() inout return => array.ptr;

    Q slerp(const Q to, T weight) const 
    {
        assert(is_normalized);
        assert(to.is_normalized);
        Q to1;

        // calc cosine
        T cosom = dot(to);

        // adjust signs (if necessary)
        if (cosom < 0) 
        {
            cosom = -cosom;
            to1 = -to;
        } 
        else 
        {
            to1 = to;
        }

        // calculate coefficients

        T scale0, scale1;
        if ((1 - cosom) > cast(T)GM_CMP_EPSILON) 
        {
            // standard case (slerp)
            T omega = gm_acos(cosom);
            T sinom = gm_sin(omega);
            scale0 = gm_sin((1 - weight) * omega) / sinom;
            scale1 = gm_sin(weight * omega) / sinom;
        } 
        else 
        {
            // "from" and "to" quaternions are very close
            //  ... so we can do a linear interpolation
            scale0 = 1 - weight;
            scale1 = weight;
        }
        // calculate final values
        return Q(
            scale0 * x + scale1 * to1.x,
            scale0 * y + scale1 * to1.y,
            scale0 * z + scale1 * to1.z,
            scale0 * w + scale1 * to1.w);
    }

    Q slerpni(const Q to, T weight) const 
    {
        assert(is_normalized() && to.is_normalized());
        const Q from = this;

        T dot = from.dot(to);

        if (gm_abs(dot) > 0.9999f) 
            return from;
        
        T theta = gm_acos(dot),
           sinT = 1 / gm_sin(theta),
           newFactor = gm_sin(weight * theta) * sinT,
           invFactor = gm_sin((1 - weight) * theta) * sinT;

        return Q(invFactor * from.x + newFactor * to.x,
                 invFactor * from.y + newFactor * to.y,
                 invFactor * from.z + newFactor * to.z,
                 invFactor * from.w + newFactor * to.w);
    }

    Q spherical_cubic_interpolate(const Q b, const Q pre_a, const Q post_b, T weight) const
    {
        assert(is_normalized());
        assert(b.is_normalized());

        Q from_q = this;
        Q pre_q = pre_a;
        Q to_q = b;
        Q post_q = post_b;

        // Align flip phases.
        from_q = BasisImpl!T(from_q).get_rotation_quaternion();
        pre_q = BasisImpl!T(pre_q).get_rotation_quaternion();
        to_q = BasisImpl!T(to_q).get_rotation_quaternion();
        post_q = BasisImpl!T(post_q).get_rotation_quaternion();

        // Flip quaternions to shortest path if necessary.
        bool flip1 = gm_signbit(from_q.dot(pre_q));
        pre_q = flip1 ? -pre_q : pre_q;
        bool flip2 = gm_signbit(from_q.dot(to_q));
        to_q = flip2 ? -to_q : to_q;
        bool flip3 = flip2 ? to_q.dot(post_q) <= 0 : gm_signbit(to_q.dot(post_q));
        post_q = flip3 ? -post_q : post_q;

        // Calc by Expmap in from_q space.
        Q ln_from = Q(0, 0, 0, 0);
        Q ln_to = (from_q.inverse() * to_q).log();
        Q ln_pre = (from_q.inverse() * pre_q).log();
        Q ln_post = (from_q.inverse() * post_q).log();
        Q ln = Q(0, 0, 0, 0);
        ln.x = gm_cubic_interpolate(ln_from.x, ln_to.x, ln_pre.x, ln_post.x, weight);
        ln.y = gm_cubic_interpolate(ln_from.y, ln_to.y, ln_pre.y, ln_post.y, weight);
        ln.z = gm_cubic_interpolate(ln_from.z, ln_to.z, ln_pre.z, ln_post.z, weight);
        Q q1 = from_q * ln.exp();

        // Calc by Expmap in to_q space.
        ln_from = (to_q.inverse() * from_q).log();
        ln_to = Q(0, 0, 0, 0);
        ln_pre = (to_q.inverse() * pre_q).log();
        ln_post = (to_q.inverse() * post_q).log();
        ln = Q(0, 0, 0, 0);
        ln.x = gm_cubic_interpolate(ln_from.x, ln_to.x, ln_pre.x, ln_post.x, weight);
        ln.y = gm_cubic_interpolate(ln_from.y, ln_to.y, ln_pre.y, ln_post.y, weight);
        ln.z = gm_cubic_interpolate(ln_from.z, ln_to.z, ln_pre.z, ln_post.z, weight);
        Q q2 = to_q * ln.exp();

        // To cancel error made by Expmap ambiguity, do blending.
        return q1.slerp(q2, weight);
    }

    Q spherical_cubic_interpolate_in_time(const Q b, const Q pre_a, const Q post_b, T weight, T b_t, T pre_a_t, T post_b_t) const 
    {
        assert(is_normalized() && b.is_normalized());

        Q from_q = this;
        Q pre_q = pre_a;
        Q to_q = b;
        Q post_q = post_b;

        // Align flip phases.
        from_q = BasisImpl!T(from_q).get_rotation_quaternion();
        pre_q = BasisImpl!T(pre_q).get_rotation_quaternion();
        to_q = BasisImpl!T(to_q).get_rotation_quaternion();
        post_q = BasisImpl!T(post_q).get_rotation_quaternion();

        // Flip quaternions to shortest path if necessary.
        bool flip1 = gm_signbit(from_q.dot(pre_q));
        pre_q = flip1 ? -pre_q : pre_q;
        bool flip2 = gm_signbit(from_q.dot(to_q));
        to_q = flip2 ? -to_q : to_q;
        bool flip3 = flip2 ? to_q.dot(post_q) <= 0 : gm_signbit(to_q.dot(post_q));
        post_q = flip3 ? -post_q : post_q;

        // Calc by Expmap in from_q space.
        Q ln_from = Q(0, 0, 0, 0);
        Q ln_to = (from_q.inverse() * to_q).log();
        Q ln_pre = (from_q.inverse() * pre_q).log();
        Q ln_post = (from_q.inverse() * post_q).log();
        Q ln = Q(0, 0, 0, 0);
        ln.x = gm_cubic_interpolate_in_time(ln_from.x, ln_to.x, ln_pre.x, ln_post.x, weight, b_t, pre_a_t, post_b_t);
        ln.y = gm_cubic_interpolate_in_time(ln_from.y, ln_to.y, ln_pre.y, ln_post.y, weight, b_t, pre_a_t, post_b_t);
        ln.z = gm_cubic_interpolate_in_time(ln_from.z, ln_to.z, ln_pre.z, ln_post.z, weight, b_t, pre_a_t, post_b_t);
        Q q1 = from_q * ln.exp();

        // Calc by Expmap in to_q space.
        ln_from = (to_q.inverse() * from_q).log();
        ln_to = Q(0, 0, 0, 0);
        ln_pre = (to_q.inverse() * pre_q).log();
        ln_post = (to_q.inverse() * post_q).log();
        ln = Q(0, 0, 0, 0);
        ln.x = gm_cubic_interpolate_in_time(ln_from.x, ln_to.x, ln_pre.x, ln_post.x, weight, b_t, pre_a_t, post_b_t);
        ln.y = gm_cubic_interpolate_in_time(ln_from.y, ln_to.y, ln_pre.y, ln_post.y, weight, b_t, pre_a_t, post_b_t);
        ln.z = gm_cubic_interpolate_in_time(ln_from.z, ln_to.z, ln_pre.z, ln_post.z, weight, b_t, pre_a_t, post_b_t);
        Q q2 = to_q * ln.exp();

        // To cancel error made by Expmap ambiguity, do blending.
        return q1.slerp(q2, weight);
    }

    private V3 xform(const V3 v) const 
    {
        assert(is_normalized());
        V3 u = V3(x, y, z);
        V3 uv = u.cross(v);
        return v + ((uv * w) + u.cross(uv)) * 2;
    }

    private V3 xform_inv(const V3 v) const => inverse().xform(v);


    // operators

    Q opBinary(string op)(const Q v) const if (op == "*")
    {
        Q r = this;
        r *= v;
        return r;
    }
    V3 opBinary(string op)(const V3 v) const if (op == "*") => xform(v);
    Q opBinary(string op)(T s) if (op == "*")             => Q(x * s, y * s, z * s, w * s);
    Q opBinary(string op)(const Q v) const if (op == "+") => Q(x + v.x  , y + v.y  , z + v.z  , w + v.w  );
    Q opBinary(string op)(const Q v) const if (op == "-") => Q(x - v.x  , y - v.y  , z - v.z  , w - v.w  );
    Q opBinary(string op)(T s) const if (op == "/")       => Q(x / s, y / s, z / s, w / s);

    Q opOpAssign(string op)(const Q v) if (op == "+") { x += v.x;   y += v.y;   z += v.z;   w += v.w;   return this; }
    Q opOpAssign(string op)(const Q v) if (op == "-") { x -= v.x;   y -= v.y;   z -= v.z;   w -= v.w;   return this; }
    Q opOpAssign(string op)(T s) if (op == "*")       { x *= s;     y *= s;     z *= s;     w *= s;     return this; }
    Q opOpAssign(string op)(T s) if (op == "/")           => this *= (1 / s);
    Q opOpAssign(string op)(const Q v) if (op == "*") 
    {
        T xx = w * v.x + x * v.w + y * v.z - z * v.y;
        T yy = w * v.y + y * v.w + z * v.x - x * v.z;
        T zz = w * v.z + z * v.w + x * v.y - y * v.x;
        w = w * v.w - x * v.x - y * v.y - z * v.z;
        x = xx;
        y = yy;
        z = zz;
        return this;
    }

    ref inout(T) opIndex(size_t n) inout => array[n];
    Q opUnary(string op)() const if (op == "+") => this;
    Q opUnary(string op)() const if (op == "-") => Q(-x, -y, -z, -w);
}










/*
    ████████╗██████╗  █████╗ ███╗   ██╗███████╗███████╗ ██████╗ ██████╗ ███╗   ███╗██████╗ ██████╗ 
    ╚══██╔══╝██╔══██╗██╔══██╗████╗  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗████╗ ████║╚════██╗██╔══██╗
       ██║   ██████╔╝███████║██╔██╗ ██║███████╗█████╗  ██║   ██║██████╔╝██╔████╔██║ █████╔╝██║  ██║
       ██║   ██╔══██╗██╔══██║██║╚██╗██║╚════██║██╔══╝  ██║   ██║██╔══██╗██║╚██╔╝██║██╔═══╝ ██║  ██║
       ██║   ██║  ██║██║  ██║██║ ╚████║███████║██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║███████╗██████╔╝
       ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚══════╝╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝╚═════╝ 
*/
struct Transform2DImpl(T)
    if (is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private
    {
        alias V2 = Vector2Impl!T;
        alias Elem = T;
        alias T2D = Transform2DImpl!T;
    }

    union
    {
        V2[3] columns = 
        [
            [ 1, 0 ],
            [ 0, 1 ],
            [ 0, 0 ],
        ];

        struct
        {
            V2 x;
            V2 y;
            V2 origin;
        }
    }

    enum T2D IDENTITY = T2D(V2( 1, 0), V2(0,  1), V2(0, 0));
    enum T2D FLIP_X   = T2D(V2(-1, 0), V2(0,  1), V2(0, 0));
    enum T2D FLIP_Y   = T2D(V2( 1, 0), V2(0, -1), V2(0, 0));

    this(T rotation, V2 position)
    {
        T cr = gm_cos(rotation);
        T sr = gm_sin(rotation);
        columns[0][0] = cr;
        columns[0][1] = sr;
        columns[1][0] = -sr;
        columns[1][1] = cr;
        columns[2] = position;
    }

    this(T rotation, const V2 scale, T skew, const V2 position)
    {
        columns[0][0] = gm_cos(rotation) * scale.x;
        columns[1][1] = gm_cos(rotation + skew) * scale.y;
        columns[1][0] = -gm_sin(rotation + skew) * scale.y;
        columns[0][1] = gm_sin(rotation) * scale.x;
        columns[2] = position;
    }

    this(V2 xAxis, V2 yAxis, V2 originPos)
    {
        x = xAxis;
        y = yAxis;
        origin = originPos;
    }

    private void affine_invert()
    {
        T det = determinant();
        assert(det != 0);
        T idet = 1.0f / det;

        gm_swap(columns[0][0], columns[1][1]);
        columns[0] *= V2(idet, -idet);
        columns[1] *= V2(-idet, idet);
        columns[2] = basis_xform(-columns[2]);
    }

    T2D affine_inverse() const
    {
        T2D r = this;
        r.affine_invert();
        return r;
    }

    V2 basis_xform(V2 v) const => V2(tdotx(v), tdoty(v));
    V2 basis_xform_inv(V2 v) const => V2(x.dot(v), y.dot(v));

    T determinant() const => x.x * y.y - x.y * y.x;
    private T tdotx(const V2 v) const => x.x * v.x + y.x * v.y;
    private T tdoty(const V2 v) const => x.y * v.x + y.y * v.y;
    V2 get_origin() const => origin;
    T get_rotation() const => gm_atan2(x.y, x.x); 
    V2 get_scale() const
    {
        T det_sign = gm_sign(determinant());
        return V2(x.length(), det_sign * y.length());
    }
    T get_skew() const => gm_acos(x.normalized().dot(gm_sign(determinant()) * y.normalized())) - cast(T)GM_PI * 0.5f;

    T2D interpolate_with(T2D xform, T weight) const
    {
        return T2D(gm_lerp_angle(get_rotation(), xform.get_rotation(), weight),
                   get_scale().lerp(xform.get_scale(), weight),
                   gm_lerp_angle(get_skew(), xform.get_skew(), weight),
                   get_origin().lerp(xform.get_origin(), weight));
    }

    private void invert()
    {
        // FIXME: this function assumes the basis is a rotation matrix, with no scaling.
        // Transform2D::affine_inverse can handle matrices with scaling, so GDScript should eventually use that.
        gm_swap(x.y, y.x);
        origin = basis_xform(-origin);
    }

    T2D inverse() const
    {
        T2D r = this;
        r.invert();
        return r;
    }

    bool is_conformal() const
    {
        if (gm_is_equal_approx(x.x, y.y) && gm_is_equal_approx(x.y, -y.x))
            return true;

        if (gm_is_equal_approx(x.x, -y.y) && gm_is_equal_approx(x.y, y.x))
            return true;

        return false;
    }

    bool is_equal_approx(const T2D transform) const 
        => x.is_equal_approx(transform.x) && y.is_equal_approx(transform.y) && origin.is_equal_approx(transform.origin);

    bool is_finite() const
        => x.is_finite() && y.is_finite() && origin.is_finite();

    T2D looking_at(V2 target) const
    {
        T2D r = T2D(get_rotation(), get_origin());
        V2 tpos = affine_inverse().xform(target);
        r.set_rotation(r.get_rotation() + (tpos * get_scale()).angle());
        return r;
    }

    private void orthonormalize()
    {
        // Gram-Schmidt Process
        x.normalize();
        y = y - x * x.dot(y);
        y.normalize();
    }

    T2D orthonormalized() const
    {
        T2D r = this;
        r.orthonormalize();
        return r;
    }

    T2D rotated(float angle) const => T2D(angle, V2.ZERO) * this; /// Equivalent to left multiplication
    T2D rotated_local(float angle) const => this * T2D(angle, V2.ZERO); /// Equivalent to right multiplication

    private void set_rotation(T rotation) 
    {
        V2 scale = get_scale();
        T cr = gm_cos(rotation);
        T sr = gm_sin(rotation);
        columns[0][0] = cr;
        columns[0][1] = sr;
        columns[1][0] = -sr;
        columns[1][1] = cr;
        set_scale(scale);
    }

    private void set_scale(const V2 scale) 
    {
        columns[0].normalize();
        columns[1].normalize();
        columns[0] *= scale.x;
        columns[1] *= scale.y;
    }
    
    private void scale(const V2 v)
    {
        scale_basis(v);
        origin *= v;
    }

    private void scale_basis(const V2 scale) 
    {
        x.x *= scale.x;
        x.y *= scale.y;
        y.x *= scale.x;
        y.y *= scale.y;
    }

    T2D scaled(const V2 scale) const
    {
        T2D copy = this;
        copy.scale(scale);
        return copy;
    }

    T2D scaled_local(const V2 scale) const => T2D(x * scale.x, y * scale.y, origin);

    T2D translated(const V2 offset) const => T2D(x, y, origin + offset);
    T2D translated_local(const V2 offset) const => T2D(x, y, origin + basis_xform(offset));

    V2 xform(const V2 vec) const => V2(tdotx(vec), tdoty(vec)) + origin;
    V2 xform_inv(const V2 vec) const
    {
        V2 v = vec - origin;
        return V2(x.dot(v), y.dot(v));
    }

    // operators
    V2 opBinary(string op)(const V2 v) const if (op == "*") => xform(v);

    T2D opBinary(string op)(const T2D transform) const if (op == "*")
    {
        T2D t = this;
        t *= transform;
        return t;
    }

    T2D opOpAssign(string op)(const T2D transform) if (op == "*") 
    {
        origin = xform(transform.origin);
        T x0 = tdotx(transform.x);
        T x1 = tdoty(transform.x);
        T y0 = tdotx(transform.y);
        T y1 = tdoty(transform.y);
        x.x = x0;
        x.y = x1;
        y.x = y0;
        y.y = y1;
        return this;
    }
}










/*
    ██████╗  █████╗ ███████╗██╗███████╗
    ██╔══██╗██╔══██╗██╔════╝██║██╔════╝
    ██████╔╝███████║███████╗██║███████╗
    ██╔══██╗██╔══██║╚════██║██║╚════██║
    ██████╔╝██║  ██║███████║██║███████║
    ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝╚══════╝
*/
struct BasisImpl(T)
    if (is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private
    {
        alias V3 = Vector3Impl!T;
        alias B = BasisImpl!T;
        alias Q = QuaternionImpl!T;
        alias Elem = T;
    }

    union
    {
        V3[3] rows =
        [
            V3(1, 0, 0),
            V3(0, 1, 0),
            V3(0, 0, 1)
        ];
    }

    V3 x() const => get_column(0);
    V3 y() const=> get_column(1);
    V3 z() const=> get_column(2);

    enum B IDENTITY = B(V3( 1, 0, 0), V3(0,  1, 0), V3(0, 0,  1));
    enum B FLIP_X   = B(V3(-1, 0, 0), V3(0,  1, 0), V3(0, 0,  1));
    enum B FLIP_Y   = B(V3( 1, 0, 0), V3(0, -1, 0), V3(0, 0,  1));
    enum B FLIP_Z   = B(V3( 1, 0, 0), V3(0,  1, 0), V3(0, 0, -1));

    this(V3 axis, T angle)
    {
        set_axis_angle(axis, angle);
    }

    this(V3 x, V3 y, V3 z)
    {
        rows[0] = V3(x.x, y.x, z.x);
        rows[1] = V3(x.y, y.y, z.y);
        rows[2] = V3(x.z, y.z, z.z);
    }

    
    private this(T xx, T xy, T xz,
                 T yx, T yy, T yz,
                 T zx, T zy, T zz)
    {
        rows[0] = V3(xx, xy, xz);
        rows[1] = V3(yx, yy, yz);
        rows[2] = V3(zx, zy, zz);
    }

    this(const Q quaternion) { set_quaternion(quaternion); }

    T determinant() const
    {
        return rows[0][0] * (rows[1][1] * rows[2][2] - rows[2][1] * rows[1][2]) -
               rows[1][0] * (rows[0][1] * rows[2][2] - rows[2][1] * rows[0][2]) +
               rows[2][0] * (rows[0][1] * rows[1][2] - rows[1][1] * rows[0][2]);
    }

    static B from_euler(V3 euler, EulerOrder order = GM_EULER_ORDER_YXZ)
    {
        B b;
        b.set_euler(euler, order);
        return b;
    }

    static B from_scale(const V3 scale)
         => B(scale.x, 0, 0, 0, scale.y, 0, 0, 0, scale.z);


    V3 get_euler(EulerOrder order = GM_EULER_ORDER_YXZ) const 
    {
        // This epsilon value results in angles within a +/- 0.04 degree range being simplified/truncated.
        // Based on testing, this is the largest the epsilon can be without the angle truncation becoming
        // visually noticeable.
        const T epsilon = 0.00000025;

        switch (order) 
        {
            case GM_EULER_ORDER_XYZ: 
            {
                // Euler angles in XYZ convention.
                // See https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
                //
                // rot =  cy*cz          -cy*sz           sy
                //        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
                //       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy

                V3 euler;
                T sy = rows[0][2];
                if (sy < (1.0f - epsilon)) 
                {
                    if (sy > -(1.0f - epsilon)) 
                    {
                        // is this a pure Y rotation?
                        if (rows[1][0] == 0 && rows[0][1] == 0 && rows[1][2] == 0 && rows[2][1] == 0 && rows[1][1] == 1) 
                        {
                            // return the simplest form (human friendlier in editor and scripts)
                            euler.x = 0;
                            euler.y = gm_atan2(rows[0][2], rows[0][0]);
                            euler.z = 0;
                        } 
                        else 
                        {
                            euler.x = gm_atan2(-rows[1][2], rows[2][2]);
                            euler.y = gm_asin(sy);
                            euler.z = gm_atan2(-rows[0][1], rows[0][0]);
                        }
                    } 
                    else 
                    {
                        euler.x = gm_atan2(rows[2][1], rows[1][1]);
                        euler.y = -GM_PI / 2.0f;
                        euler.z = 0.0f;
                    }
                } 
                else 
                {
                    euler.x = gm_atan2(rows[2][1], rows[1][1]);
                    euler.y = GM_PI / 2.0f;
                    euler.z = 0.0f;
                }
                return euler;
            }

            case GM_EULER_ORDER_XZY: 
            {
                // Euler angles in XZY convention.
                // See https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
                //
                // rot =  cz*cy             -sz             cz*sy
                //        sx*sy+cx*cy*sz    cx*cz           cx*sz*sy-cy*sx
                //        cy*sx*sz          cz*sx           cx*cy+sx*sz*sy

                V3 euler;
                T sz = rows[0][1];
                if (sz < (1.0f - epsilon)) 
                {
                    if (sz > -(1.0f - epsilon)) 
                    {
                        euler.x = gm_atan2(rows[2][1], rows[1][1]);
                        euler.y = gm_atan2(rows[0][2], rows[0][0]);
                        euler.z = gm_asin(-sz);
                    } 
                    else 
                    {
                        // It's -1
                        euler.x = -gm_atan2(rows[1][2], rows[2][2]);
                        euler.y = 0.0f;
                        euler.z = GM_PI / 2.0f;
                    }
                } 
                else 
                {
                    // It's 1
                    euler.x = -gm_atan2(rows[1][2], rows[2][2]);
                    euler.y = 0.0f;
                    euler.z = -GM_PI / 2.0f;
                }
                return euler;
            }

            case GM_EULER_ORDER_YXZ: 
            {
                // Euler angles in YXZ convention.
                // See https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
                //
                // rot =  cy*cz+sy*sx*sz    cz*sy*sx-cy*sz        cx*sy
                //        cx*sz             cx*cz                 -sx
                //        cy*sx*sz-cz*sy    cy*cz*sx+sy*sz        cy*cx

                V3 euler;

                T m12 = rows[1][2];

                if (m12 < (1 - epsilon)) 
                {
                    if (m12 > -(1 - epsilon)) 
                    {
                        // is this a pure X rotation?
                        if (rows[1][0] == 0 && rows[0][1] == 0 && rows[0][2] == 0 && rows[2][0] == 0 && rows[0][0] == 1) {
                            // return the simplest form (human friendlier in editor and scripts)
                            euler.x = gm_atan2(-m12, rows[1][1]);
                            euler.y = 0;
                            euler.z = 0;
                        } 
                        else 
                        {
                            euler.x = gm_asin(-m12);
                            euler.y = gm_atan2(rows[0][2], rows[2][2]);
                            euler.z = gm_atan2(rows[1][0], rows[1][1]);
                        }
                    } 
                    else 
                    { // m12 == -1
                        euler.x = GM_PI * 0.5f;
                        euler.y = gm_atan2(rows[0][1], rows[0][0]);
                        euler.z = 0;
                    }
                } 
                else 
                { // m12 == 1
                    euler.x = -GM_PI * 0.5f;
                    euler.y = -gm_atan2(rows[0][1], rows[0][0]);
                    euler.z = 0;
                }

                return euler;
            }

            case GM_EULER_ORDER_YZX: 
            {
                // Euler angles in YZX convention.
                // See https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
                //
                // rot =  cy*cz             sy*sx-cy*cx*sz     cx*sy+cy*sz*sx
                //        sz                cz*cx              -cz*sx
                //        -cz*sy            cy*sx+cx*sy*sz     cy*cx-sy*sz*sx

                V3 euler;
                T sz = rows[1][0];
                if (sz < (1.0f - epsilon)) 
                {
                    if (sz > -(1.0f - epsilon)) 
                    {
                        euler.x = gm_atan2(-rows[1][2], rows[1][1]);
                        euler.y = gm_atan2(-rows[2][0], rows[0][0]);
                        euler.z = gm_asin(sz);
                    } 
                    else 
                    {
                        // It's -1
                        euler.x = gm_atan2(rows[2][1], rows[2][2]);
                        euler.y = 0.0f;
                        euler.z = -GM_PI / 2.0f;
                    }
                } 
                else 
                {
                    // It's 1
                    euler.x = gm_atan2(rows[2][1], rows[2][2]);
                    euler.y = 0.0f;
                    euler.z = GM_PI / 2.0f;
                }
                return euler;
            } 

            case GM_EULER_ORDER_ZXY: 
            {
                // Euler angles in ZXY convention.
                // See https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
                //
                // rot =  cz*cy-sz*sx*sy    -cx*sz                cz*sy+cy*sz*sx
                //        cy*sz+cz*sx*sy    cz*cx                 sz*sy-cz*cy*sx
                //        -cx*sy            sx                    cx*cy
                V3 euler;
                T sx = rows[2][1];
                if (sx < (1.0f - epsilon)) 
                {
                    if (sx > -(1.0f - epsilon)) 
                    {
                        euler.x = gm_asin(sx);
                        euler.y = gm_atan2(-rows[2][0], rows[2][2]);
                        euler.z = gm_atan2(-rows[0][1], rows[1][1]);
                    } else 
                    {
                        // It's -1
                        euler.x = -GM_PI / 2.0f;
                        euler.y = gm_atan2(rows[0][2], rows[0][0]);
                        euler.z = 0;
                    }
                } 
                else 
                {
                    // It's 1
                    euler.x = GM_PI / 2.0f;
                    euler.y = gm_atan2(rows[0][2], rows[0][0]);
                    euler.z = 0;
                }
                return euler;
            }

            case GM_EULER_ORDER_ZYX: 
            {
                // Euler angles in ZYX convention.
                // See https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
                //
                // rot =  cz*cy             cz*sy*sx-cx*sz        sz*sx+cz*cx*cy
                //        cy*sz             cz*cx+sz*sy*sx        cx*sz*sy-cz*sx
                //        -sy               cy*sx                 cy*cx
                V3 euler;
                T sy = rows[2][0];
                if (sy < (1.0f - epsilon)) 
                {
                    if (sy > -(1.0f - epsilon)) 
                    {
                        euler.x = gm_atan2(rows[2][1], rows[2][2]);
                        euler.y = gm_asin(-sy);
                        euler.z = gm_atan2(rows[1][0], rows[0][0]);
                    } 
                    else 
                    {
                        // It's -1
                        euler.x = 0;
                        euler.y = GM_PI / 2.0f;
                        euler.z = -gm_atan2(rows[0][1], rows[1][1]);
                    }
                } 
                else 
                {
                    // It's 1
                    euler.x = 0;
                    euler.y = -GM_PI / 2.0f;
                    euler.z = -gm_atan2(rows[0][1], rows[1][1]);
                }
                return euler;
        }

        default:
            assert(0); // bad Euler order
        }
    }

    Q get_quaternion() const 
    {
        assert(is_rotation());

        /* Allow getting a quaternion from an unnormalized transform */
        B m = this;
        T trace = m.rows[0][0] + m.rows[1][1] + m.rows[2][2];
        T[4] temp;

        if (trace > 0)
        {
            T s = gm_sqrt(trace + 1);
            temp[3] = (s * 0.5f);
            s = 0.5f / s;
            temp[0] = ((m.rows[2][1] - m.rows[1][2]) * s);
            temp[1] = ((m.rows[0][2] - m.rows[2][0]) * s);
            temp[2] = ((m.rows[1][0] - m.rows[0][1]) * s);
        } 
        else
        {
            int i = m.rows[0][0] < m.rows[1][1]
                 ? (m.rows[1][1] < m.rows[2][2] ? 2 : 1)
                 : (m.rows[0][0] < m.rows[2][2] ? 2 : 0);
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;

            T s = gm_sqrt(m.rows[i][i] - m.rows[j][j] - m.rows[k][k] + 1);
            temp[i] = s * 0.5f;
            s = 0.5f / s;

            temp[3] = (m.rows[k][j] - m.rows[j][k]) * s;
            temp[j] = (m.rows[j][i] + m.rows[i][j]) * s;
            temp[k] = (m.rows[k][i] + m.rows[i][k]) * s;
        }
        return Q(temp[0], temp[1], temp[2], temp[3]);
    }

    Q get_rotation_quaternion() const
    {
        // Assumes that the matrix can be decomposed into a proper rotation and scaling matrix as M = R.S,
        // and returns the Euler angles corresponding to the rotation part, complementing get_scale().
        // See the comment in get_scale() for further information.
        B m = orthonormalized();
        T det = m.determinant();
        if (det < 0) 
        {
            // Ensure that the determinant is 1, such that result is a proper rotation matrix which can be represented by Euler angles.
            m.scale(V3(-1, -1, -1));
        }
        return m.get_quaternion();
    }

    private V3 get_column(int index) const => V3(rows[0][index], rows[1][index], rows[2][index]);

    // get_scale works with get_rotation, use get_scale_abs if you need to enforce positive signature.
    V3 get_scale() const 
    {
        // FIXME: We are assuming M = R.S (R is rotation and S is scaling), and use polar decomposition to extract R and S.
        // A polar decomposition is M = O.P, where O is an orthogonal matrix (meaning rotation and reflection) and
        // P is a positive semi-definite matrix (meaning it contains absolute values of scaling along its diagonal).
        //
        // Despite being different from what we want to achieve, we can nevertheless make use of polar decomposition
        // here as follows. We can split O into a rotation and a reflection as O = R.Q, and obtain M = R.S where
        // we defined S = Q.P. Now, R is a proper rotation matrix and S is a (signed) scaling matrix,
        // which can involve negative scalings. However, there is a catch: unlike the polar decomposition of M = O.P,
        // the decomposition of O into a rotation and reflection matrix as O = R.Q is not unique.
        // Therefore, we are going to do this decomposition by sticking to a particular convention.
        // This may lead to confusion for some users though.
        //
        // The convention we use here is to absorb the sign flip into the scaling matrix.
        // The same convention is also used in other similar functions such as get_rotation_axis_angle, get_rotation, ...
        //
        // A proper way to get rid of this issue would be to store the scaling values (or at least their signs)
        // as a part of Basis. However, if we go that path, we need to disable direct (write) access to the
        // matrix elements.
        //
        // The rotation part of this decomposition is returned by get_rotation* functions.
        T det_sign = gm_sign(determinant());
        return get_scale_abs() * det_sign;
    }

    private V3 get_scale_abs() const
    {
        return V3( V3(rows[0][0], rows[1][0], rows[2][0]).length(),
                   V3(rows[0][1], rows[1][1], rows[2][1]).length(),
                   V3(rows[0][2], rows[1][2], rows[2][2]).length() );
    }

    B inverse() const
    {
        B inv;
        inv = this;
        inv.invert();
        return inv;
    }

    private void invert()
    {
        T cofac(int row1, int col1, int row2, int col2)
            => (rows[row1][col1] * rows[row2][col2] - rows[row1][col2] * rows[row2][col1]);

        T[3] co = 
        [
            cofac(1, 1, 2, 2), cofac(1, 2, 2, 0), cofac(1, 0, 2, 1)
        ];

        T det = rows[0][0] * co[0] +
                rows[0][1] * co[1] +
                rows[0][2] * co[2];

        assert(det != 0);

        T s = 1.0 / det;

        set(co[0] * s, cofac(0, 2, 2, 1) * s, cofac(0, 1, 1, 2) * s,
            co[1] * s, cofac(0, 0, 2, 2) * s, cofac(0, 2, 1, 0) * s,
            co[2] * s, cofac(0, 1, 2, 0) * s, cofac(0, 0, 1, 1) * s);
    }

    bool is_conformal() const
    {
        const V3 x = get_column(0);
        const V3 y = get_column(1);
        const V3 z = get_column(2);
        const T x_len_sq = x.length_squared();
        return gm_is_equal_approx(x_len_sq, y.length_squared())
            && gm_is_equal_approx(x_len_sq, z.length_squared())
            && gm_is_zero_approx(x.dot(y))
            && gm_is_zero_approx(x.dot(z))
            && gm_is_zero_approx(y.dot(z));
    }

    bool is_rotation() const => is_conformal() && gm_is_equal_approx(determinant(), cast(T)1, cast(T)GM_UNIT_EPSILON);
    bool is_equal_approx(const B b) const => rows[0].is_equal_approx(b.rows[0]) && rows[1].is_equal_approx(b.rows[1]) && rows[2].is_equal_approx(b.rows[2]);
    bool is_finite() const => rows[0].is_finite() && rows[1].is_finite() && rows[2].is_finite();

    static B looking_at(V3 target, V3 up = V3(0, 1, 0), bool use_model_front = false)
    {
        assert(!target.is_zero_approx());
        assert(!up.is_zero_approx());
        V3 z = target.normalized();
        if (!use_model_front)
        {
            z = -z;
        }
        V3 x = up.cross(z);
        if (x.is_zero_approx())
        {
            //WARN_PRINT("Target and up vectors are colinear. This is not advised as it may cause unwanted rotation around local Z axis.");
            x = up.get_any_perpendicular(); // Vectors are almost parallel.
        }
        x.normalize();
        V3 y = z.cross(x);

        B basis;
        basis.set_columns(x, y, z);
        return basis;
    }

    B orthonormalized() const
    {
        B m = this;
        m.orthonormalize();
        return m;
    }

    private void orthonormalize()
    {
        // Gram-Schmidt Process
        V3 x = get_column(0);
        V3 y = get_column(1);
        V3 z = get_column(2);
        x.normalize();
        y = (y - x * (x.dot(y)));
        y.normalize();
        z = (z - x * (x.dot(z)) - y * (y.dot(z)));
        z.normalize();
        set_column(0, x);
        set_column(1, y);
        set_column(2, z);
    }

    private void rotate(const V3 axis, T angle) { this = rotated(axis, angle); }
    B rotated(const V3 axis, T angle) const => B(axis, angle) * this;

    private void scale(const V3 scale)
    {
        rows[0] *= scale.x;
        rows[1] *= scale.y;
        rows[2] *= scale.z;
    }

    B scaled(const V3 scale) const
    {
        B m = this;
        m.scale(scale);
        return m;
    }

    private void scale_local(const V3 scale)
    {
        rows[0] *= scale;
        rows[1] *= scale;
        rows[2] *= scale;
    }

    B scaled_local(const V3 scale) const
    {
        B m = this;
        m.scale_local(scale);
        return m;
    }

    private void _set_diagonal(const V3 diag) 
    {
        rows[0][0] = diag.x;
        rows[0][1] = 0;
        rows[0][2] = 0;

        rows[1][0] = 0;
        rows[1][1] = diag.y;
        rows[1][2] = 0;

        rows[2][0] = 0;
        rows[2][1] = 0;
        rows[2][2] = diag.z;
    }

    private void set_quaternion_scale(const Q quaternion, const V3 scale) 
    {
        _set_diagonal(scale);
        rotate(quaternion);
    }

    private void rotate(const Q quaternion)
    {
        this = B(quaternion) * this;
    }

    B slerp(B to, T weight) const
    {
        Q from = Q(this);
        Q qto = Q(to);
        B b = B(from.slerp(qto, weight));
        b.rows[0] *= gm_lerp(rows[0].length(), to.rows[0].length(), weight);
        b.rows[1] *= gm_lerp(rows[1].length(), to.rows[1].length(), weight);
        b.rows[2] *= gm_lerp(rows[2].length(), to.rows[2].length(), weight);
        return b;
    }

    private void set(T xx, T xy, T xz,
                     T yx, T yy, T yz,
                     T zx, T zy, T zz)
    {
        rows[0] = V3(xx, xy, xz);
        rows[1] = V3(yx, yy, yz);
        rows[2] = V3(zx, zy, zz);
    }

    private void set_axis_angle(const V3 axis, T angle)
    {
        // Rotation matrix from axis and angle, see https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_angle
        assert(axis.is_normalized());
        V3 axis_sq = V3(axis.x * axis.x, axis.y * axis.y, axis.z * axis.z);
        T cosine = gm_cos(angle);
        rows[0][0] = axis_sq.x + cosine * (1 - axis_sq.x);
        rows[1][1] = axis_sq.y + cosine * (1 - axis_sq.y);
        rows[2][2] = axis_sq.z + cosine * (1 - axis_sq.z);

        T sine = gm_sin(angle);
        T t = 1 - cosine;

        T xyzt = axis.x * axis.y * t;
        T zyxs = axis.z * sine;
        rows[0][1] = xyzt - zyxs;
        rows[1][0] = xyzt + zyxs;

        xyzt = axis.x * axis.z * t;
        zyxs = axis.y * sine;
        rows[0][2] = xyzt + zyxs;
        rows[2][0] = xyzt - zyxs;

        xyzt = axis.y * axis.z * t;
        zyxs = axis.x * sine;
        rows[1][2] = xyzt - zyxs;
        rows[2][1] = xyzt + zyxs;
    }

    private void set_column(int index, const V3 c)
    {
        rows[0][index] = c.x;
        rows[1][index] = c.y;
        rows[2][index] = c.z;
    }

    private void set_columns(const V3 x, const V3 y, const V3 z)
    {
        set_column(0, x);
        set_column(1, y);
        set_column(2, z);
    }

    private void set_euler(const V3 euler, EulerOrder order)
    {
        T c, s;

        c = gm_cos(euler.x);
        s = gm_sin(euler.x);
        B xmat = B(1, 0, 0, 0, c, -s, 0, s, c);

        c = gm_cos(euler.y);
        s = gm_sin(euler.y);
        B ymat = B(c, 0, s, 0, 1, 0, -s, 0, c);

        c = gm_cos(euler.z);
        s = gm_sin(euler.z);
        B zmat = B(c, -s, 0, s, c, 0, 0, 0, 1);

        switch (order)
        {
            case GM_EULER_ORDER_XYZ: this = xmat * ymat * zmat; break;
            case GM_EULER_ORDER_XZY: this = xmat * zmat * ymat; break;
            case GM_EULER_ORDER_YXZ: this = ymat * xmat * zmat; break;
            case GM_EULER_ORDER_YZX: this = ymat * zmat * xmat; break;
            case GM_EULER_ORDER_ZXY: this = zmat * xmat * ymat; break;
            case GM_EULER_ORDER_ZYX: this = zmat * ymat * xmat; break;
            default:
                assert(false);
        }
    }

    private void set_quaternion(const Q q)
    {
        T d = q.length_squared();
        T s = 2.0f / d;
        T xs = q.x * s, ys = q.y * s, zs = q.z * s;
        T wx = q.w * xs, wy = q.w * ys, wz = q.w * zs;
        T xx = q.x * xs, xy = q.x * ys, xz = q.x * zs;
        T yy = q.y * ys, yz = q.y * zs, zz = q.z * zs;
        set(1.0f - (yy + zz), xy - wz, xz + wy,
            xy + wz, 1.0f - (xx + zz), yz - wx,
            xz - wy, yz + wx, 1.0f - (xx + yy));
    }

    // transposed dot products
    T tdotx(const V3 v) const =>
        rows[0][0] * v.x + rows[1][0] * v.y + rows[2][0] * v.z;

    T tdoty(const V3 v) const =>
        rows[0][1] * v.x + rows[1][1] * v.y + rows[2][1] * v.z;

    T tdotz(const V3 v) const =>
        rows[0][2] * v.x + rows[1][2] * v.y + rows[2][2] * v.z;

    private void transpose() 
    {
        gm_swap!T(rows[0][1], rows[1][0]);
        gm_swap!T(rows[0][2], rows[2][0]);
        gm_swap!T(rows[1][2], rows[2][1]);
    }

    B transposed() const 
    {
        B tr = this;
        tr.transpose();
        return tr;
    }

    V3 xform(const V3 v) const => V3(rows[0].dot(v), rows[1].dot(v), rows[2].dot(v));
    
    V3 xform_inv(const V3 v) const
    {
        assert(is_conformal());
        return V3((rows[0][0] * v.x) + (rows[1][0] * v.y) + (rows[2][0] * v.z),
                  (rows[0][1] * v.x) + (rows[1][1] * v.y) + (rows[2][1] * v.z),
                  (rows[0][2] * v.x) + (rows[1][2] * v.y) + (rows[2][2] * v.z));
    }

    U opCast(U)() const if (isBasisImpl!U)
    {
        static if (is(U.Elem == float))
        {
            alias VD = Vector3Impl!float;
            return U(cast(VD)rows[0], cast(VD)rows[1], cast(VD)rows[2]);
        }
        else static if (is(U.Elem == double))
        {
            alias VD = Vector3Impl!double;
            return U(cast(VD)rows[0], cast(VD)rows[1], cast(VD)rows[2]);
        }
        else
            static assert(false);
    }

    // operators
    V3 opIndex(size_t n) const => rows[n];

    B opBinary(string op)(const B m) const if (op == "*")
        => B(m.tdotx(rows[0]), m.tdoty(rows[0]), m.tdotz(rows[0]),
             m.tdotx(rows[1]), m.tdoty(rows[1]), m.tdotz(rows[1]),
             m.tdotx(rows[2]), m.tdoty(rows[2]), m.tdotz(rows[2]));

    V3 opBinary(string op)(const V3 v) const if (op == "*") => xform(v);

    B opOpAssign(string op)(const B m) if (op == "*")
    {
        this = this * m;
        return this;
    }
}









/*
    ████████╗██████╗  █████╗ ███╗   ██╗███████╗███████╗ ██████╗ ██████╗ ███╗   ███╗██████╗ ██████╗ 
    ╚══██╔══╝██╔══██╗██╔══██╗████╗  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗████╗ ████║╚════██╗██╔══██╗
       ██║   ██████╔╝███████║██╔██╗ ██║███████╗█████╗  ██║   ██║██████╔╝██╔████╔██║ █████╔╝██║  ██║
       ██║   ██╔══██╗██╔══██║██║╚██╗██║╚════██║██╔══╝  ██║   ██║██╔══██╗██║╚██╔╝██║ ╚═══██╗██║  ██║
       ██║   ██║  ██║██║  ██║██║ ╚████║███████║██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║██████╔╝██████╔╝
       ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚══════╝╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═════╝ 
*/
struct Transform3DImpl(T)
    if (is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private
    {
        alias V3   = Vector3Impl!T;
        alias B    = BasisImpl!T;
        alias Q    = QuaternionImpl!T;
        alias T3D  = Transform3DImpl!T;
        alias Elem = T;
    }

    // Identity by default
    B basis;
    V3 origin;

    enum T3D IDENTITY = T3D(V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1), V3(0, 0, 0));
    enum T3D FLIP_X = T3D(V3(-1, 0, 0), V3(0, 1, 0), V3(0, 0, 1), V3(0, 0, 0));
    enum T3D FLIP_Y = T3D(V3(1, 0, 0), V3(0, -1, 0), V3(0, 0, 1), V3(0, 0, 0));
    enum T3D FLIP_Z = T3D(V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, -1), V3(0, 0, 0));

    this(B basis, V3 origin)
    {
        this.basis = basis;
        this.origin = origin;
    }

    // TODO this(Projection projection)

    this(V3 x_axis, V3 y_axis, V3 z_axis, V3 origin)
    {
        basis = B(x_axis, y_axis, z_axis);
        this.origin = origin;
    }

    private void affine_invert() 
    {
        basis.invert();        
        origin = basis.xform(-origin);
    }

    T3D affine_inverse() const
    {
        T3D r = this;
        r.affine_invert();
        return r;
    }

    T3D interpolate_with(const T3D xform, T weight) const 
    {
        T3D interp;
        V3 src_scale = basis.get_scale();
        Q src_rot = basis.get_rotation_quaternion();
        V3 src_loc = origin;
        V3 dst_scale = xform.basis.get_scale();
        Q dst_rot = xform.basis.get_rotation_quaternion();
        V3 dst_loc = xform.origin;
        interp.basis.set_quaternion_scale(src_rot.slerp(dst_rot, weight).normalized(), src_scale.lerp(dst_scale, weight));
        interp.origin = src_loc.lerp(dst_loc, weight);
        return interp;
    }

    private void invert()
    {
        // If the assertion fail:
        // Transform3D.inverse() should be used on a rotation. Maybe use affine_inverse()?
        assert(basis.is_rotation());

        basis.transpose();
        origin = basis.xform(-origin);
    }

    T3D inverse() const
    {
        T3D r = this;
        r.invert();
        return r;
    }
    bool is_equal_approx(T3D xform) const => basis.is_equal_approx(xform.basis) && origin.is_equal_approx(xform.origin);
    bool is_finite() const => basis.is_finite() && origin.is_finite();
    
    T3D looking_at(V3 target, V3 up = V3(0, 1, 0), bool use_model_front) const
    {
        assert(!origin.is_equal_approx(target));
        T3D t = this;
        t.basis = B.looking_at(target - origin, up, use_model_front);
        return t;
    }

    private void orthonormalize()
    {
        basis.orthonormalize();
    }
    
    T3D orthonormalized() const
    {
        T3D t = this;
        t.orthonormalize();
        return t;
    }

    private void rotate(const V3 axis, T angle) { this = rotated(axis, angle); }
    T3D rotated(const V3 axis, T angle) const 
    {
        // Equivalent to left multiplication
        B nbasis = B(axis, angle);
        return T3D(nbasis * basis, nbasis.xform(origin));
    }

    T3D rotated_local(const V3 axis, T angle) const 
    {
        // Equivalent to right multiplication
        B nbasis = B(axis, angle);
        return T3D(basis * nbasis, origin);
    }

    private void scale(const V3 scale) 
    {
        basis.scale(scale);
        origin *= scale;
    }

    T3D scaled(const V3 scale) const 
    {
        // Equivalent to left multiplication
        return T3D(basis.scaled(scale), origin * scale);
    }

    T3D scaled_local(const V3 scale) const 
    {
        // Equivalent to right multiplication
        return T3D(basis.scaled_local(scale), origin);
    }

    T3D translated(V3 offset) const => T3D(basis, origin + offset);
    T3D translated_local(V3 offset) const => T3D(basis, origin + basis.xform(offset));

    V3 xform(const V3 v) const => V3(basis[0].dot(v) + origin.x,
                                     basis[1].dot(v) + origin.y,
                                     basis[2].dot(v) + origin.z);

    // operators
    T3D opBinary(string op)(const T3D transform) const if (op == "*")
    {
        T3D r = this;
        r *= transform;
        return r;
    }
    T3D opOpAssign(string op)(const T3D transform) if (op == "*") 
    {
        origin = xform(transform.origin);
        basis *= transform.basis;
        return this;
    }
    V3 opBinary(string op)(const V3 v) const if (op == "*") => xform(v);
    T3D opBinary(string op)(const T fact) const if (op == "*")
    {
        T3D r = this;
        r.origin *= fact;
        r.basis *= fact;
        return r;
    }

    
}










/*
    ██████╗ ██████╗  ██████╗      ██╗███████╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗
    ██╔══██╗██╔══██╗██╔═══██╗     ██║██╔════╝██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║
    ██████╔╝██████╔╝██║   ██║     ██║█████╗  ██║        ██║   ██║██║   ██║██╔██╗ ██║
    ██╔═══╝ ██╔══██╗██║   ██║██   ██║██╔══╝  ██║        ██║   ██║██║   ██║██║╚██╗██║
    ██║     ██║  ██║╚██████╔╝╚█████╔╝███████╗╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║
*/
struct ProjectionImpl(T)
    if (is(T == float) || is(T == double))
{
pure nothrow @nogc @safe:

    private
    {
        alias V2   = Vector2Impl!T;
        alias V4   = Vector4Impl!T;
        alias P    = ProjectionImpl!T;
        alias T3D  = Transform3DImpl!T;
        alias Elem = T;
    }

    // Identity by default
    union
    {
        V4[4] columns = [V4(1, 0, 0, 0),
                         V4(0, 1, 0, 0),
                         V4(0, 0, 1, 0),
                         V4(0, 0, 0, 1)];
        struct
        {
            V4 x, y, z, w;
        }

        T[16] m; // added to avoid @trusted
    }

    enum P IDENTITY = P.init;
    enum P ZERO = P(V4.ZERO, V4.ZERO, V4.ZERO, V4.ZERO);

    this(V4 x, V4 y, V4 z, V4 w)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    this(Transform3D tr)
    {
        m[0] = tr.basis.rows[0][0];
        m[1] = tr.basis.rows[1][0];
        m[2] = tr.basis.rows[2][0];
        m[3] = 0.0;
        m[4] = tr.basis.rows[0][1];
        m[5] = tr.basis.rows[1][1];
        m[6] = tr.basis.rows[2][1];
        m[7] = 0.0;
        m[8] = tr.basis.rows[0][2];
        m[9] = tr.basis.rows[1][2];
        m[10] = tr.basis.rows[2][2];
        m[11] = 0.0;
        m[12] = tr.origin.x;
        m[13] = tr.origin.y;
        m[14] = tr.origin.z;
        m[15] = 1.0;
    }

    static P create_depth_correction(bool flip_y)
    {
        P p;
        p.set_depth_correction(flip_y);
        return p;
    }

    //TODO Projection create_fit_aabb(aabb: AABB) static

    static P create_for_hmd(int eye, T aspect, T intraocular_dist, T display_width, T display_to_lens, T oversample, T z_near, T z_far)
    {
        P proj;
        proj.set_for_hmd(eye, aspect, intraocular_dist, display_width, display_to_lens, oversample, z_near, z_far);
        return proj;
    }

    static P create_frustum(T left, T right, T bottom, T top, T z_near, T z_far)
    {
        P proj;
        proj.set_frustum(left, right, bottom, top, z_near, z_far);
        return proj;
    }

    static P create_frustum_aspect(T size, T aspect, V2 offset, T z_near, T z_far, bool flip_fov = false)
    {
        P proj;
        proj.set_frustum(size, aspect, offset, z_near, z_far, flip_fov);
        return proj;
    }

    //TODO static P create_light_atlas_rect(Rect2 rect)

    static P create_orthogonal(T left, T right, T bottom, T top, T z_near, T z_far)
    {
        P proj;
        proj.set_orthogonal(left, right, bottom, top, z_near, z_far);
        return proj;
    }

    static P create_orthogonal_aspect(T size, T aspect, T z_near, T z_far, bool flip_fov) 
    {
        P proj;
        proj.set_orthogonal(size, aspect, z_near, z_far, flip_fov);
        return proj;
    }

    static P create_perspective(T fovy, T aspect, T z_near, T z_far, bool flip_fov = false)
    {
        P proj;
        proj.set_perspective(fovy, aspect, z_near, z_far, flip_fov);
        return proj;
    }

    static T get_fovy(T fovx, T aspect) 
    {
        return gm_rad_to_deg(gm_atan(aspect * gm_tan(gm_deg_to_rad(fovx) * 0.5f)) * 2.0f);
    }

    private void set_depth_correction(bool flip_y = true, bool reverse_z  =true, bool remap_z = true)
    {
        m[0] = 1;
        m[1] = 0.0;
        m[2] = 0.0;
        m[3] = 0.0;
        m[4] = 0.0;
        m[5] = flip_y ? -1 : 1;
        m[6] = 0.0;
        m[7] = 0.0;
        m[8] = 0.0;
        m[9] = 0.0;
        m[10] = remap_z ? (reverse_z ? -0.5 : 0.5) : (reverse_z ? -1.0 : 1.0);
        m[11] = 0.0;
        m[12] = 0.0;
        m[13] = 0.0;
        m[14] = remap_z ? 0.5 : 0.0;
        m[15] = 1.0;
    }

    private void set_for_hmd(int eye, T aspect, T intraocular_dist, T display_width, T display_to_lens, T oversample, T z_near, T z_far) 
    {
        // we first calculate our base frustum on our values without taking our lens magnification into account.
        T f1 = (intraocular_dist * 0.5) / display_to_lens;
        T f2 = ((display_width - intraocular_dist) * 0.5) / display_to_lens;
        T f3 = (display_width / 4.0) / display_to_lens;

        // now we apply our oversample factor to increase our FOV. how much we oversample is always a balance we strike between performance and how much
        // we're willing to sacrifice in FOV.
        T add = ((f1 + f2) * (oversample - 1.0)) / 2.0;
        f1 += add;
        f2 += add;
        f3 *= oversample;

        // always apply KEEP_WIDTH aspect ratio
        f3 /= aspect;

        switch (eye) 
        {
        case 1:  // left eye
            set_frustum(-f2 * z_near, f1 * z_near, -f3 * z_near, f3 * z_near, z_near, z_far);
            break;
        case 2:  // right eye
            set_frustum(-f1 * z_near, f2 * z_near, -f3 * z_near, f3 * z_near, z_near, z_far);
            break;
        default:  // mono, does not apply here!
            break;
        }
    }


    private void set_frustum(T left, T right, T bottom, T top, T near, T far)
    {
        assert(right > left);
        assert(top > bottom);
        assert(far > near);
        T x = 2 * near / (right - left);
        T y = 2 * near / (top - bottom);
        T a = (right + left) / (right - left);
        T b = (top + bottom) / (top - bottom);
        T c = -(far + near) / (far - near);
        T d = -2 * far * near / (far - near);
        m[0] = x;
        m[1] = 0;
        m[2] = 0;
        m[3] = 0;
        m[4] = 0;
        m[5] = y;
        m[6] = 0;
        m[7] = 0;
        m[8] = a;
        m[9] = b;
        m[10] = c;
        m[11] = -1;
        m[12] = 0;
        m[13] = 0;
        m[14] = d;
        m[15] = 0;
    }

    private void set_frustum(T size, T aspect, V2 offset, T near, T far, bool flip_fov)
    {
        if (!flip_fov)
            size *= aspect;
        set_frustum(-size / 2 + offset.x, +size / 2 + offset.x, -size / aspect / 2 + offset.y, size / aspect / 2 + offset.y, near, far);
    }

    private void set_identity() 
    {
        this = P.init;
    }

    private void set_orthogonal(T left, T right, T bottom, T top, T znear, T zfar) 
    {
        set_identity();
        columns[0][0] = 2.0 / (right - left);
        columns[3][0] = -((right + left) / (right - left));
        columns[1][1] = 2.0 / (top - bottom);
        columns[3][1] = -((top + bottom) / (top - bottom));
        columns[2][2] = -2.0 / (zfar - znear);
        columns[3][2] = -((zfar + znear) / (zfar - znear));
        columns[3][3] = 1.0;
    }

    private void set_orthogonal(T size, T aspect, T znear, T zfar, bool flip_fov) 
    {
        if (!flip_fov)
            size *= aspect;

        set_orthogonal(-size / 2, +size / 2, -size / aspect / 2, +size / aspect / 2, znear, zfar);
    }

    private void set_perspective(T fovy_degrees, T aspect, T z_near, T z_far, bool flip_fov) 
    {
        if (flip_fov) 
        {
            fovy_degrees = get_fovy(fovy_degrees, 1.0 / aspect);
        }

        T sine, cotangent, deltaZ;
        T radians = gm_deg_to_rad(fovy_degrees / 2.0);

        deltaZ = z_far - z_near;
        sine = gm_sin(radians);

        if ((deltaZ == 0) || (sine == 0) || (aspect == 0)) 
        {
            return;
        }
        cotangent = gm_cos(radians) / sine;

        set_identity();

        columns[0][0] = cotangent / aspect;
        columns[1][1] = cotangent;
        columns[2][2] = -(z_far + z_near) / deltaZ;
        columns[2][3] = -1;
        columns[3][2] = -2 * z_near * z_far / deltaZ;
        columns[3][3] = 0;
    }

    private void set_perspective(T p_fovy_degrees, T p_aspect, T p_z_near, T p_z_far, bool p_flip_fov, int p_eye, T p_intraocular_dist, T p_convergence_dist) 
    {
        if (p_flip_fov) 
        {
            p_fovy_degrees = get_fovy(p_fovy_degrees, 1.0 / p_aspect);
        }

        T left, right, modeltranslation, ymax, xmax, frustumshift;

        ymax = p_z_near * gm_tan(gm_deg_to_rad(p_fovy_degrees / 2.0));
        xmax = ymax * p_aspect;
        frustumshift = (p_intraocular_dist / 2.0) * p_z_near / p_convergence_dist;

        switch (p_eye) 
        {
            case 1: { // left eye
                left = -xmax + frustumshift;
                right = xmax + frustumshift;
                modeltranslation = p_intraocular_dist / 2.0;
            } break;
            case 2: 
            { // right eye
                left = -xmax - frustumshift;
                right = xmax - frustumshift;
                modeltranslation = -p_intraocular_dist / 2.0;
            } break;
            default: 
            { // mono, should give the same result as set_perspective(p_fovy_degrees,p_aspect,p_z_near,p_z_far,p_flip_fov)
                left = -xmax;
                right = xmax;
                modeltranslation = 0.0;
            } break;
        }

        set_frustum(left, right, -ymax, ymax, p_z_near, p_z_far);

        // translate matrix by (modeltranslation, 0.0, 0.0)
        P cm;
        cm.set_identity();
        cm.columns[3][0] = modeltranslation;
        this = this * cm;
    }

    // operators

    P opBinary(string op)(const P matrix) const if (op == "*")
    {
        P r;
        for (int j = 0; j < 4; j++) 
        {
            for (int i = 0; i < 4; i++) 
            {
                T ab = 0;
                for (int k = 0; k < 4; k++) 
                    ab += columns[k][i] * matrix.columns[j][k];
                r.columns[j][i] = ab;
            }
        }
        return r;
    }
}


// internal

private:

enum bool isVector2Impl(T)     = is(T : Vector2Impl!U, U...);
enum bool isVector3Impl(T)     = is(T : Vector3Impl!U, U...);
enum bool isVector4Impl(T)     = is(T : Vector4Impl!U, U...);
enum bool isQuaternionImpl(T)  = is(T : QuaternionImpl!U, U...);
enum bool isTransform2DImpl(T) = is(T : Transform2DImpl!U, U...);
enum bool isTransform3DImpl(T) = is(T : Transform3DImpl!U, U...);
enum bool isBasisImpl(T)       = is(T : BasisImpl!U, U...);
enum bool isProjectionImpl(T)  = is(T : ProjectionImpl!U, U...);

auto assumePureNothrowNogc(T, Args...)(T expr, auto ref Args args) pure nothrow @nogc @trusted if (isSomeFunction!T) {
    static if (is(T Fptr : Fptr*) && is(Fptr == function))
        alias ft = pure nothrow @nogc ReturnType!T function(Parameters!T);
    else static if (is(T Fdlg == delegate))
        alias ft = pure nothrow @nogc ReturnType!T delegate(Parameters!T);
    else
        static assert(0);

    return (cast(ft)expr)(args);
}
