/*
    Copyright (c) 2014-2025 Godot Engine contributors.
    Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur. 
    Copyright (c) 2025 Guillaume Piolat
*/
module godotmath;

import std.math;

nothrow @nogc @safe:

// Implementation done for: 
// Vector2 => https://docs.godotengine.org/en/stable/classes/class_vector2.html
// ...TBD
//
// Note: In case of semantic conflict, the way Godot does it is favoured.

// TODO: lack of pure because eg. atan2f isn't pure.

/* 
    ██╗   ██╗███████╗ ██████╗████████╗ ██████╗ ██████╗ ██████╗ 
    ██║   ██║██╔════╝██╔════╝╚══██╔══╝██╔═══██╗██╔══██╗╚════██╗
    ██║   ██║█████╗  ██║        ██║   ██║   ██║██████╔╝ █████╔╝
    ╚██╗ ██╔╝██╔══╝  ██║        ██║   ██║   ██║██╔══██╗██╔═══╝ 
     ╚████╔╝ ███████╗╚██████╗   ██║   ╚██████╔╝██║  ██║███████╗
      ╚═══╝  ╚══════╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚══════╝
*/
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

    enum Vector2 ZERO  = Vector2( 0.0f,  0.0f);
    enum Vector2 ONE   = Vector2( 1.0f,  1.0f);
    enum Vector2 INF   = Vector2( GM_INF,  GM_INF);
    enum Vector2 LEFT  = Vector2(-1.0f,  0.0f);
    enum Vector2 RIGHT = Vector2( 1.0f,  0.0f);
    enum Vector2 UP    = Vector2( 0.0f, -1.0f);
    enum Vector2 DOWN  = Vector2( 0.0f,  1.0f);

    // All functions as defined at: 
    this(float x, float y) { this.x = x; this.y = y; }
    Vector2 abs() => Vector2(x < 0 ? -x : x, y < 0 ? -y : y); 
    float angle() const => atan2(y, x);
    float angle_to(in Vector2 to) const => atan2(cross(to), dot(to));
    float angle_to_point(const Vector2 to) const => (to - this).angle();
    float aspect() const => width / height;
    Vector2 bezier_derivative(const Vector2 c1, const Vector2 c2, const Vector2 end, float p_t) const 
       => Vector2( .bezier_derivative(x, c1.x, c2.x, end.x, p_t), 
                   .bezier_derivative(y, c1.y, c2.y, end.y, p_t) );
    Vector2 bezier_interpolate(const Vector2 c1, const Vector2 c2, const Vector2 end, float p_t) const 
        => Vector2( .bezier_interpolate(x, c1.x, c2.x, end.x, p_t), 
                    .bezier_interpolate(y, c1.y, c2.y, end.y, p_t) );
    Vector2 bounce(const Vector2 normal) const => -reflect(normal);
    Vector2 ceil() => Vector2(.ceil(x), .ceil(y));
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
    float distance_to(const Vector2 v) const => sqrt(distance_squared_to(v));
    float dot(const Vector2 other) const => x * other.x + y * other.y;
    Vector2 floor() => Vector2(.floor(x), .floor(y));
    static Vector2 from_angle(float angle) => Vector2( .cos(angle), .sin(angle) );
    bool is_equal_approx(const Vector2 other) => .is_equal_approx(x, other.x) && .is_equal_approx(y, other.y);

    bool is_finite() const => .isFinite(x) && .isFinite(y);
    bool is_normalized() const => .is_equal_approx(length_squared(), 1.0f, cast(float)GM_UNIT_EPSILON);
    bool is_zero_approx() const => .is_zero_approx(x) && .is_zero_approx(y);
    float length() const => sqrt(x * x + y * y);
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
            l = sqrt(l);
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
    Vector2 posmod(float mod) const => Vector2(.fposmod(x, mod), .fposmod(y, mod));
    Vector2 posmodv(const Vector2 modv) const => Vector2(.fposmod(x, modv.x), .fposmod(y, modv.y));
    Vector2 project(const Vector2 to) const => to * (dot(to) / to.length_squared());
    Vector2 reflect(const Vector2 normal) const
    {
        assert(normal.is_normalized());
        return normal * 2 * dot(normal) - this;
    }
    Vector2 rotated(float by_radians) const 
    {
        float sine = .sin(by_radians);
        float cosi = .cos(by_radians);
        return Vector2(x * cosi - y * sine,
                       x * sine + y * cosi);
    }
    Vector2 round() const => Vector2(.round(x), .round(y));
    Vector2 sign() const => Vector2(.signf(x), .signf(y));

    Vector2 slerp(const Vector2 to, float weight) const 
    {
        float start_length_sq = length_squared();
	    float end_length_sq = to.length_squared();
	    if (start_length_sq == 0.0f || end_length_sq == 0.0f) 
        {
		    // Zero length vectors have no angle, so the best we can do is either lerp or throw an error.
		    return lerp(to, weight);
	    }
	    float start_length = .sqrt(start_length_sq);
	    float result_length = .lerp(start_length, .sqrt(end_length_sq), weight);
	    float angle = angle_to(to);
	    return rotated(angle * weight) * (result_length / start_length);
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
    ref inout(float) opIndex(size_t n) inout return { assert(n < 2); return n ? y : x; }
    Vector2 opUnary(string op)() const if (op == "+") => this;    
    Vector2 opUnary(string op)() const if (op == "-") => Vector2(-x, -y);

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
struct Vector3
{
nothrow @nogc @safe:

    union { float x = 0; float width; }
    union { float y = 0; float height; }
    union { float z = 0; float depth; }

    enum : int
    {
        AXIS_X,
        AXIS_Y,
        AXIS_Z,
    }

    // Common 3D direction vectors (assuming Y-up and Z-back/forward convention)
    enum Vector3 ZERO    = Vector3( 0.0f,  0.0f,  0.0f);
    enum Vector3 ONE     = Vector3( 1.0f,  1.0f,  1.0f);
    enum Vector3 INF     = Vector3( GM_INF, GM_INF, GM_INF);
    enum Vector3 LEFT    = Vector3(-1.0f,  0.0f,  0.0f);
    enum Vector3 RIGHT   = Vector3( 1.0f,  0.0f,  0.0f);
    enum Vector3 UP      = Vector3( 0.0f,  1.0f,  0.0f);
    enum Vector3 DOWN    = Vector3( 0.0f, -1.0f,  0.0f);
    enum Vector3 FORWARD = Vector3( 0.0f,  0.0f, -1.0f);
    enum Vector3 BACK    = Vector3( 0.0f,  0.0f,  1.0f);

    this(float x, float y, float z) { this.x = x; this.y = y; this.z = z; }
    Vector3 abs() const => Vector3(x < 0 ? -x : x, y < 0 ? -y : y, z < 0 ? -z : z);
    float angle_to(in Vector3 to) => atan2(cross(to).length(), dot(to));
    Vector3 bezier_derivative(const Vector3 c1, const Vector3 c2, const Vector3 end, float p_t) const
        => Vector3( .bezier_derivative(x, c1.x, c2.x, end.x, p_t),
                    .bezier_derivative(y, c1.y, c2.y, end.y, p_t),
                   .bezier_derivative(z, c1.z, c2.z, end.z, p_t) );
    Vector3 bezier_interpolate(const Vector3 c1, const Vector3 c2, const Vector3 end, float p_t) const
        => Vector3( .bezier_interpolate(x, c1.x, c2.x, end.x, p_t),
                    .bezier_interpolate(y, c1.y, c2.y, end.y, p_t),
                   .bezier_interpolate(z, c1.z, c2.z, end.z, p_t) );
    Vector3 bounce(const Vector3 normal) const => -reflect(normal);
    Vector3 ceil() const => Vector3(.ceil(x), .ceil(y), .ceil(z));
    Vector3 clamp(const Vector3 min, const Vector3 max) const 
        => Vector3(.clampf(x, min.x, max.x), .clampf(y, min.y, max.y), .clampf(z, min.z, max.z));
    Vector3 clampf(float min, float max) const 
        => Vector3(.clampf(x, min, max), .clampf(y, min, max), .clampf(z, min, max));
    Vector3 cross(const Vector3 other) const
        => Vector3( y * other.z - z * other.y,
                    z * other.x - x * other.z,
                    x * other.y - y * other.x );
    Vector3 cubic_interpolate(const Vector3 b, const Vector3 pre_a, const Vector3 post_b, float weight) const
        => Vector3( .cubic_interpolate(x, b.x, pre_a.x, post_b.x, weight),
                    .cubic_interpolate(y, b.y, pre_a.y, post_b.y, weight),
                    .cubic_interpolate(z, b.z, pre_a.z, post_b.z, weight) );
    Vector3 cubic_interpolate_in_time(const Vector3 b, const Vector3 pre_a, const Vector3 post_b, float weight, float b_t, float pre_a_t, float post_b_t) const
        => Vector3( .cubic_interpolate_in_time(x, b.x, pre_a.x, post_b.x, weight, b_t, pre_a_t, post_b_t),
                    .cubic_interpolate_in_time(y, b.y, pre_a.y, post_b.y, weight, b_t, pre_a_t, post_b_t),
                    .cubic_interpolate_in_time(z, b.z, pre_a.z, post_b.z, weight, b_t, pre_a_t, post_b_t) );
    Vector3 direction_to(const Vector3 to) const => (to - this).normalized();
    float distance_squared_to(const Vector3 v) const 
        => (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z);
    float distance_to(const Vector3 v) const => sqrt(distance_squared_to(v));
    float dot(const Vector3 other) const => x * other.x + y * other.y + z * other.z;
    Vector3 floor() const => Vector3(.floor(x), .floor(y), .floor(z));
    Vector3 inverse() const => Vector3(1.0f / x, 1.0f / y, 1.0f / z);
    bool is_equal_approx(const Vector3 other) const
        => .is_equal_approx(x, other.x) && .is_equal_approx(y, other.y) && .is_equal_approx(z, other.z);
    bool is_finite() const => .isFinite(x) && .isFinite(y) && .isFinite(z);
    bool is_normalized() const 
        => .is_equal_approx(length_squared(), 1.0f, cast(float)GM_UNIT_EPSILON);
    bool is_zero_approx() const 
        => .is_zero_approx(x) && .is_zero_approx(y) && .is_zero_approx(z);
    float length() const => sqrt(x * x + y * y + z * z);
    float length_squared() const => x * x + y * y + z * z;
    Vector3 lerp(const Vector3 to, float weight) const 
        => Vector3( .lerp(x, to.x, weight), .lerp(y, to.y, weight), .lerp(z, to.z, weight) );
    Vector3 limit_length(float len) const
    {
        float l = length();
        Vector3 v = this;
        if (l > 0 && len < l)
        {
            v /= l;
            v *= len;
        }
        return v;
    }

    Vector3 max(const Vector3 other) const 
        => Vector3( x > other.x ? x : other.x, y > other.y ? y : other.y, z > other.z ? z : other.z );
    int max_axis_index() const
    {
        if (x > y && x > z) return AXIS_X;
        if (y > z) return AXIS_Y;
        return AXIS_Z;
    }
    Vector3 maxf(float v) const 
        => Vector3( x > v ? x : v, y > v ? y : v, z > v ? z : v );
    Vector3 min(const Vector3 other) const 
        => Vector3( x < other.x ? x : other.x, y < other.y ? y : other.y, z < other.z ? z : other.z );
    int min_axis_index() const
    {
        if (x < y && x < z) return AXIS_X;
        if (y < z) return AXIS_Y;
        return AXIS_Z;
    }
    Vector3 minf(float v) const 
        => Vector3( x < v ? x : v, y < v ? y : v, z < v ? z : v );
    Vector3 move_toward(const Vector3 to, float delta) const
    {
        Vector3 v = this;
        Vector3 vd = to - v;
        float len = vd.length();
        return len <= delta || len < cast(float)GM_CMP_EPSILON ? to : v + vd / len * delta;
    }
    void normalize()
    {
        float l = x * x + y * y + z * z;
        if (l != 0)
        {
            l = sqrt(l);
            x /= l;
            y /= l;
            z /= l;
        }
    }
    Vector3 normalized() const
    {
        Vector3 v = this;
        v.normalize();
        return v;
    }
    //TODO octahedron_decode
    //TODO octahedron_encode
    //TODO outer
    Vector3 posmod(float mod) const => Vector3(.fposmod(x, mod), .fposmod(y, mod), .fposmod(z, mod));
    Vector3 posmodv(const Vector3 modv) const => Vector3(.fposmod(x, modv.x), .fposmod(y, modv.y), .fposmod(z, modv.z));
    Vector3 project(const Vector3 to) const => to * (dot(to) / to.length_squared());

    Vector3 reflect(const Vector3 normal) const
    {
        assert(normal.is_normalized());
        return normal * 2 * dot(normal) - this;
    }

    void rotate(const Vector3 axis, float angle)
    {
        float axis_length_sq = axis.length_squared();
        if (.is_zero_approx(axis_length_sq))
            return;

        Vector3 norm_axis = axis / sqrt(axis_length_sq);

        float sin_angle = sin(angle);
        float cos_angle = cos(angle);

        Vector3 term1 = this * cos_angle;
        Vector3 term2 = norm_axis.cross(this) * sin_angle;
        Vector3 term3 = norm_axis * (norm_axis.dot(this) * (1.0f - cos_angle));

        this = term1 + term2 + term3;
    }

    Vector3 rotated(const Vector3 axis, float angle) const
    {
        Vector3 v = this;
        v.rotate(axis, angle);
        return v;
    }

    Vector3 round() const => Vector3(.round(x), .round(y), .round(z));
    Vector3 sign() const => Vector3(.signf(x), .signf(y), .signf(z));

    // float signed_angle_to(to: Vector3, axis: Vector3
    // Vector3 slerp(to: Vector3, weight: float) const 

// Vector3 slide(n: Vector3) const
// snapped
// Vector3 snappedf(step: float) const 

    // operators
    Vector3 opBinary(string op)(const Vector3 v) const if (op == "*") => Vector3(x * v.x, y * v.y, z * v.z);
    Vector3 opBinary(string op)(float scale) const if (op == "*") => Vector3(x * scale, y * scale, z * scale);
    Vector3 opBinary(string op)(int scale) const if (op == "*") => Vector3(x * scale, y * scale, z * scale);
    Vector3 opBinary(string op)(const Vector3 v) const if (op == "+") => Vector3(x + v.x, y + v.y, z + v.z);
    Vector3 opBinary(string op)(const Vector3 v) const if (op == "-") => Vector3(x - v.x, y - v.y, z - v.z);
    Vector3 opBinary(string op)(const Vector3 v) const if (op == "/") => Vector3(x / v.x, y / v.y, z / v.z);
    Vector3 opBinary(string op)(float scale) const if (op == "/") => Vector3(x / scale, y / scale, z / scale);
    Vector3 opBinary(string op)(int scale) const if (op == "/") => Vector3(x / scale, y / scale, z / scale);

    Vector3 opOpAssign(string op)(const Vector3 v) if (op == "*") { x *= v.x; y *= v.y; z *= v.z; return this; }
    Vector3 opOpAssign(string op)(float scale) if (op == "*") { x *= scale; y *= scale; z *= scale; return this; }
    Vector3 opOpAssign(string op)(int scale) if (op == "*") { x *= scale; y *= scale; z *= scale; return this; }
    Vector3 opOpAssign(string op)(const Vector3 v) if (op == "+") { x += v.x; y += v.y; z += v.z; return this; }
    Vector3 opOpAssign(string op)(const Vector3 v) if (op == "-") { x -= v.x; y -= v.y; z -= v.z; return this; }
    Vector3 opOpAssign(string op)(const Vector3 v) if (op == "/") { x /= v.x; y /= v.y; z /= v.z; return this; }
    Vector3 opOpAssign(string op)(float scale) if (op == "/") { x /= scale; y /= scale; z /= scale; return this; }
    Vector3 opOpAssign(string op)(int scale) if (op == "/") { x /= scale; y /= scale; z /= scale; return this; }

    ref inout(float) opIndex(size_t n) inout return { assert(n < 3); switch (n) { case 0: return x; case 1: return y; default: return z; } }

    Vector3 opUnary(string op)() const if (op == "+") => this;
    Vector3 opUnary(string op)() const if (op == "-") => Vector3(-x, -y, -z);
}








/*
    ██████╗  █████╗ ███████╗██╗███████╗
    ██╔══██╗██╔══██╗██╔════╝██║██╔════╝
    ██████╔╝███████║███████╗██║███████╗
    ██╔══██╗██╔══██║╚════██║██║╚════██║
    ██████╔╝██║  ██║███████║██║███████║
    ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝╚══════╝
*/
struct Basis
{
    alias real_t = float;

nothrow @nogc @safe:

    Vector3[3] rows = 
    [
        Vector3(1, 0, 0),
        Vector3(0, 1, 0),
        Vector3(0, 0, 1)
    ];

    this(Vector3 axis, float angle)
    {
        set_axis_angle(axis, angle);
    }

    this(Vector3 x, Vector3 y, Vector3 z)
    {
        rows[0] = x;
        rows[1] = y;
        rows[2] = z;
    }

    this(float xx, float xy, float xz, 
         float yx, float yy, float yz, 
         float zx, float zy, float zz)
    {
        rows[0] = Vector3(xx, xy, xz);
        rows[1] = Vector3(yx, yy, yz);
        rows[2] = Vector3(zx, zy, zz);
    }

    float determinant() const
    {
        return rows[0][0] * (rows[1][1] * rows[2][2] - rows[2][1] * rows[1][2]) -
               rows[1][0] * (rows[0][1] * rows[2][2] - rows[2][1] * rows[0][2]) +
               rows[2][0] * (rows[0][1] * rows[1][2] - rows[1][1] * rows[0][2]);
    }

    static Basis from_euler(Vector3 euler, int order = EULER_ORDER_YXZ)
    {
        Basis b;
        b.set_euler(euler, order);
        return b;
    }

    // operators
    Vector3 opIndex(size_t n) const => rows[n];
    Basis opBinary(string op)(const Basis m) const if (op == "*")
        => Basis(m.tdotx(rows[0]), m.tdoty(rows[0]), m.tdotz(rows[0]),
                 m.tdotx(rows[1]), m.tdoty(rows[1]), m.tdotz(rows[1]),
                 m.tdotx(rows[2]), m.tdoty(rows[2]), m.tdotz(rows[2]));

    void set_axis_angle(const Vector3 axis, float angle)
    {
       // Rotation matrix from axis and angle, see https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_angle
       assert(axis.is_normalized);
       Vector3 axis_sq = Vector3(axis.x * axis.x, axis.y * axis.y, axis.z * axis.z);
        float cosine = cos(angle);
        rows[0][0] = axis_sq.x + cosine * (1.0f - axis_sq.x);
        rows[1][1] = axis_sq.y + cosine * (1.0f - axis_sq.y);
        rows[2][2] = axis_sq.z + cosine * (1.0f - axis_sq.z);

        float sine = sin(angle);
        float t = 1 - cosine;

        float xyzt = axis.x * axis.y * t;
        float zyxs = axis.z * sine;
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

    void set_euler(const Vector3 euler, EulerOrder order) 
    {
        float c, s;

        c = cos(euler.x);
        s = sin(euler.x);
        Basis xmat = Basis(1, 0, 0, 0, c, -s, 0, s, c);

        c = cos(euler.y);
        s = sin(euler.y);
        Basis ymat = Basis(c, 0, s, 0, 1, 0, -s, 0, c);

        c = cos(euler.z);
        s = sin(euler.z);
        Basis zmat = Basis(c, -s, 0, s, c, 0, 0, 0, 1);

        switch (order) 
        {
            case EULER_ORDER_XYZ: this = xmat * ymat * zmat; break;
            case EULER_ORDER_XZY: this = xmat * zmat * ymat; break;
            case EULER_ORDER_YXZ: this = ymat * xmat * zmat; break;
            case EULER_ORDER_YZX: this = ymat * zmat * xmat; break;
            case EULER_ORDER_ZXY: this = zmat * xmat * ymat; break;
            case EULER_ORDER_ZYX: this = zmat * ymat * xmat; break;
            default: 
                assert(false);
        }
    }

    // transposed dot products
	real_t tdotx(const Vector3 v) const =>
		rows[0][0] * v[0] + rows[1][0] * v[1] + rows[2][0] * v[2];

	real_t tdoty(const Vector3 v) const =>
		rows[0][1] * v[0] + rows[1][1] * v[1] + rows[2][1] * v[2];

	real_t tdotz(const Vector3 v) const =>
		rows[0][2] * v[0] + rows[1][2] * v[1] + rows[2][2] * v[2];
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


// EulerOrder
alias EulerOrder = int;
enum : EulerOrder
{
    EULER_ORDER_XYZ = 0,
    EULER_ORDER_XZY = 1,
    EULER_ORDER_YXZ = 2,
    EULER_ORDER_YZX = 3,
    EULER_ORDER_ZXY = 4,
    EULER_ORDER_ZYX = 5,
}

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
    float value = fmod(x, y);
    if (((value < 0) && (y > 0)) || ((value > 0) && (y < 0))) 
    {
        value += y;
    }
    value += 0.0f;
    return value;
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

bool is_equal_approx(float p_left, float p_right, float p_tolerance) 
{
    // Check for exact equality first, required to handle "infinity" values.
    if (p_left == p_right)
        return true;

    // Then check for approximate equality.
    return abs(p_left - p_right) < p_tolerance;
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
    float tolerance = cast(float)GM_CMP_EPSILON * abs(p_left);
    if (tolerance < cast(float)GM_CMP_EPSILON) 
    {
        tolerance = cast(float)GM_CMP_EPSILON;
    }
    return abs(p_left - p_right) < tolerance;
}

bool is_zero_approx(double p_value) 
{
    return abs(p_value) < GM_CMP_EPSILON;
}
bool is_zero_approx(float p_value) 
{
    return abs(p_value) < cast(float)GM_CMP_EPSILON;
}

double lerp(double from, double to, double weight) 
{
    return from + (to - from) * weight;
}

float lerp(float from, float to, float weight) 
{
    return from + (to - from) * weight;
}

float signf(float v) 
{
    return v > 0 ? +1.0f : (v < 0 ? -1.0f : 0.0f);
}

double sign(double v) 
{
    return v > 0 ? +1.0f : (v < 0 ? -1.0f : 0.0f);
}

