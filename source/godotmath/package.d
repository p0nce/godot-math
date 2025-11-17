module godotmath;

import core.stdc.math: sinf, cosf, sqrtf, atan2f, ceilf, floorf, fabsf, fmodf,
                       sin, cos, sqrt, atan2,  ceil,  floor, fabs, fmod,
                       isfinite, isinf;

nothrow @nogc @safe:


// TODO: lack of pure because eg. atan2f isn't pure.

/// See_also: https://docs.godotengine.org/en/stable/classes/class_vector2.html
struct Vector2
{
nothrow @nogc @safe:

    union { float x = 0; float width;  }
    union { float y = 0; float height; }

    enum : int
    {
        AXIS_X,
        AXIS_Y,
    };

    enum Vector2 LEFT  = Vector2(-1.0f,  0.0f);
    enum Vector2 RIGHT = Vector2( 1.0f,  0.0f);
    enum Vector2 UP    = Vector2( 0.0f, -1.0f);
    enum Vector2 DOWN  = Vector2( 0.0f,  1.0f);

    // All functions as defined at: 
    this(float x, float y) { this.x = x; this.y = y; }
    Vector2 abs() => Vector2(x < 0 ? -x : x, y < 0 ? -y : y); 
    float angle() const => atan2f(y, x);
    float angle_to(const Vector2 to) const => atan2f(cross(to), dot(to));
    float angle_to_point(const Vector2 to) const => (to - this).angle();
    float aspect() const => width / height;
    Vector2 bezier_derivative(const Vector2 c1, const Vector2 c2, const Vector2 end, float p_t) const 
       => Vector2( .bezier_derivative(x, c1.x, c2.x, end.x, p_t), 
                   .bezier_derivative(y, c1.y, c2.y, end.y, p_t) );
    Vector2 bezier_interpolate(const Vector2 c1, const Vector2 c2, const Vector2 end, float p_t) const 
        => Vector2( .bezier_interpolate(x, c1.x, c2.x, end.x, p_t), 
                    .bezier_interpolate(y, c1.y, c2.y, end.y, p_t) );
    Vector2 bounce(const Vector2 normal) const => -reflect(normal);
    Vector2 ceil() => Vector2(.ceilf(x), .ceilf(y));
    Vector2 clamp(const Vector2 min, const Vector2 max) const => Vector2(.clampf(x, min.x, max.x), .clampf(y, min.y, max.y));
    Vector2 clampf(float min, float max) const => Vector2(.clampf(x, min, max), .clampf(y, min, max));
    float cross(const Vector2 p_other) const => x * p_other.y - y * p_other.x;
    Vector2 cubic_interpolate(const Vector2 b, const Vector2 pre_a, const Vector2 post_b, float weight) const
        => Vector2( .cubic_interpolate(x, b.x, pre_a.x, post_b.x, weight),
                    .cubic_interpolate(y, b.y, pre_a.y, post_b.y, weight) );
    Vector2 cubic_interpolate_in_time(const Vector2 b, const Vector2 pre_a, const Vector2 post_b, float weight, float b_t, float pre_a_t, float post_b_t) const
        => Vector2( .cubic_interpolate_in_time(x, b.x, pre_a.x, post_b.x, weight, b_t, pre_a_t, post_b_t),
                    .cubic_interpolate_in_time(y, b.y, pre_a.y, post_b.y, weight, b_t, pre_a_t, post_b_t) );
    Vector2 direction_to(const Vector2 to) const => (to - this).normalized();
    float distance_squared_to(const Vector2 v) const => (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y);
    float distance_to(const Vector2 v) const => sqrtf(distance_squared_to(v));
    float dot(const Vector2 other) const => x * other.x + y * other.y;
    Vector2 floor() => Vector2(.floorf(x), .floorf(y));
    static Vector2 from_angle(float angle) => Vector2( .cos(angle), .sin(angle) );
    bool is_equal_approx(const Vector2 other) => .is_equal_approx(x, other.x) && .is_equal_approx(y, other.y);

    bool is_finite() const => .isfinite(x) && .isfinite(y);
    bool is_normalized() const => .is_equal_approx(length_squared(), 1.0f, cast(float)GM_UNIT_EPSILON);
    bool is_zero_approx() const => .is_zero_approx(x) && .is_zero_approx(y);
    float length() const => sqrtf(x * x + y * y);
    float length_squared() const => x * x + y * y;
    Vector2 lerp(const Vector2 to, float weight) const => Vector2( .lerp(x, to.x, weight), .lerp(y, to.y, weight) );
    Vector2 limit_length(float len) const 
    {
        float l = length();
        Vector2 v = this;
        if (l > 0 && len < l) 
        {
            v /= l;
            v *= len;
        }
        return v;
    }
    Vector2 max(const Vector2 other) const => Vector2( x > other.x ? x : other.x, y > other.y ? y : other.y );
    int max_axis_index() const => x < y ? AXIS_Y : AXIS_X;
    Vector2 maxf(float v) const => Vector2( x > v ? x : v, y > v ? y : v );    
    Vector2 min(const Vector2 other) const => Vector2( x < other.x ? x : other.x, y < other.y ? y : other.y );
    int min_axis_index() const => x < y ? AXIS_X : AXIS_Y;
    Vector2 minf(float v) const => Vector2( x < v ? x : v, y < v ? y : v );
    Vector2 move_toward(const Vector2 to, float delta) const 
    {
        Vector2 v = this;
        Vector2 vd = to - v;
        float len = vd.length();
        return len <= delta || len < cast(float)GM_CMP_EPSILON ? to : v + vd / len * delta;
    }
    void normalize()
    {
        float l = x * x + y * y;
        if (l != 0) 
        {
            l = (l);
            x /= l;
            y /= l;
        }
    }
    Vector2 normalized() const 
    {
        Vector2 v = this;
        v.normalize();
        return v;
    }
    Vector2 orthogonal() const => Vector2(y, -x);
    

    Vector2 reflect(const Vector2 normal) const
    {
        assert(normal.is_normalized());
        return normal * 2 * dot(normal) - this;
    }

    // operators
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "*") => Vector2(x * v.x, y * v.y);
    Vector2 opBinary(string op)(float scale) const if (op == "*") => Vector2(x * scale, y * scale);
    Vector2 opBinary(string op)(int scale) const if (op == "*") => Vector2(x * scale, y * scale);
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "+") => Vector2(x + v.x, y + v.y);
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "-") => Vector2(x - v.x, y - v.y);
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "/") => Vector2(x / v.x, y / v.y);
    Vector2 opBinary(string op)(float scale) const if (op == "/") => Vector2(x / scale, y / scale);
    Vector2 opBinary(string op)(int scale) const if (op == "/") => Vector2(x / scale, y / scale);
    Vector2 opOpAssign(string op)(const Vector2 v) if (op == "*") { x *= v.x; y *= v.y; return this; }
    Vector2 opOpAssign(string op)(float scale) if (op == "*") { x *= scale; y *= scale; return this; }
    Vector2 opOpAssign(string op)(int scale) if (op == "*") { x *= scale; y *= scale; return this; }
    Vector2 opOpAssign(string op)(const Vector2 v) if (op == "+") { x += v.x; y += v.y; return this; }
    Vector2 opOpAssign(string op)(const Vector2 v) if (op == "-") { x -= v.x; y -= v.y; return this; }
    Vector2 opOpAssign(string op)(const Vector2 v) if (op == "/") { x /= v.x; y /= v.y; return this; }
    Vector2 opOpAssign(string op)(float scale) if (op == "/") { x /= scale; y /= scale; return this; }
    Vector2 opOpAssign(string op)(int scale) if (op == "/") { x /= scale; y /= scale; return this; }
    float opIndex(size_t n) const { assert(n < 2); return n ? y : x; }
    Vector2 opUnary(string op)() const if (op == "+") => this;    
    Vector2 opUnary(string op)() const if (op == "-") => Vector2(-x, -y);

}


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

// math funcs

double bezier_interpolate(double p_start, double p_control_1, double p_control_2, double p_end, double p_t) 
{
    double omt = (1.0 - p_t);
    double omt2 = omt * omt;
    double omt3 = omt2 * omt;
    double t2 = p_t * p_t;
    double t3 = t2 * p_t;
    return p_start * omt3 + p_control_1 * omt2 * p_t * 3.0 + p_control_2 * omt * t2 * 3.0 + p_end * t3;
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

double bezier_derivative(double p_start, double p_control_1, 
    double p_control_2, double p_end, double p_t) 
{
    double omt = (1.0 - p_t);
    double omt2 = omt * omt;
    double t2 = p_t * p_t;

    double d = (p_control_1 - p_start) * 3.0 * omt2 + (p_control_2 - p_control_1) * 6.0 * omt * p_t + (p_end - p_control_2) * 3.0 * t2;
    return d;
}

float bezier_derivative(float p_start, float p_control_1, 
    float p_control_2, float p_end, float p_t) 
{
    float omt = (1.0f - p_t);
    float omt2 = omt * omt;
    float t2 = p_t * p_t;
    float d = (p_control_1 - p_start) * 3.0f * omt2 + (p_control_2 - p_control_1) * 6.0f * omt * p_t + (p_end - p_control_2) * 3.0f * t2;
    return d;
}

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

double fposmod(double x, double y) 
{
    double value = fmod(x, y);
    if (((value < 0) && (y > 0)) || ((value > 0) && (y < 0))) 
    {
        value += y;
    }
    value += 0.0;
    return value;
}

float fposmod(float x, float y) 
{
    float value = fmodf(x, y);
    if (((value < 0) && (y > 0)) || ((value > 0) && (y < 0))) 
    {
        value += y;
    }
    value += 0.0f;
    return value;
}

double lerp(double from, double to, double weight) 
{
    return from + (to - from) * weight;
}

float lerp(float from, float to, float weight) 
{
    return from + (to - from) * weight;
}

bool is_equal_approx(double p_left, double p_right, double p_tolerance) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) 
    {
        return true;
    }
    // Then check for approximate equality.
    return fabs(p_left - p_right) < p_tolerance;
}

bool is_equal_approx(float p_left, float p_right, float p_tolerance) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right)
        return true;

    // Then check for approximate equality.
    return fabsf(p_left - p_right) < p_tolerance;
}

bool is_equal_approx(double p_left, double p_right) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) 
    {
        return true;
    }
    // Then check for approximate equality.
    double tolerance = GM_CMP_EPSILON * fabs(p_left);
    if (tolerance < GM_CMP_EPSILON) 
    {
        tolerance = GM_CMP_EPSILON;
    }
    return fabs(p_left - p_right) < tolerance;
}

bool is_equal_approx(float p_left, float p_right) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right) {
        return true;
    }
    // Then check for approximate equality.
    float tolerance = cast(float)GM_CMP_EPSILON * fabsf(p_left);
    if (tolerance < cast(float)GM_CMP_EPSILON) 
    {
        tolerance = cast(float)GM_CMP_EPSILON;
    }
    return fabsf(p_left - p_right) < tolerance;
}

bool is_zero_approx(double p_value) 
{
    return fabs(p_value) < GM_CMP_EPSILON;
}
bool is_zero_approx(float p_value) 
{
    return fabsf(p_value) < cast(float)GM_CMP_EPSILON;
}


