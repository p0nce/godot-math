module godotmath;

import core.stdc.math: atan2f;

nothrow @nogc @safe:


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

struct Vector2
{
nothrow @nogc @safe:

    union { float x = 0; float width;  }
    union { float y = 0; float height; }
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
    

    float dot(const(Vector2) other) const => x * other.x + y * other.y;
    float cross(const(Vector2) p_other) const => x * p_other.y - y * p_other.x;

    // operators
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "*") => Vector2(x * v.x, y * v.y);
    Vector2 opBinary(string op)(float scale) const if (op == "*") => Vector2(x * scale, y * scale);
    Vector2 opBinary(string op)(int scale) const if (op == "*") => Vector2(x * scale, y * scale);
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "+") => Vector2(x + v.x, y + v.y);
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "-") => Vector2(x - v.x, y - v.y);
    Vector2 opBinary(string op)(const Vector2 v) const if (op == "/") => Vector2(x / v.x, y / v.y);
    Vector2 opBinary(string op)(float scale) const if (op == "/") => Vector2(x / scale, y / scale);
    Vector2 opBinary(string op)(int scale) const if (op == "/") => Vector2(x / scale, y / scale);
    float opIndex(size_t n) const { assert(n < 2); return n ? y : x; }
    Vector2 opUnary(string op)() const if (op == "+") => this;    
    Vector2 opUnary(string op)() const if (op == "-") => Vector2(-x, -y);

}