/*
    Copyright (c) 2014-2025 Godot Engine contributors.
    Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur. 
    Copyright (c) 2025 Guillaume Piolat
*/
module godotmath;

import godotmath.globals;

nothrow @nogc @safe:


// Implementation done for: 
// Vector2 => https://docs.godotengine.org/en/stable/classes/class_vector2.html
// ...TBD
//
// Note: In case of semantic conflict, the way Godot does it is favoured.

// TODO: lack of pure because eg. atan2f isn't pure.

// TODO: != for vectors

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
{
nothrow @nogc @safe:

    alias V = Vector2Impl!T;

    union { T x = 0; T width;  }
    union { T y = 0; T height; }

    enum : int
    {
        AXIS_X,
        AXIS_Y,
    };

    enum V ZERO  = V( 0.0f,  0.0f);
    enum V ONE   = V( 1.0f,  1.0f);
    enum V INF   = V( GM_INF,  GM_INF);
    enum V LEFT  = V(-1.0f,  0.0f);
    enum V RIGHT = V( 1.0f,  0.0f);
    enum V UP    = V( 0.0f, -1.0f);
    enum V DOWN  = V( 0.0f,  1.0f);

    // All functions as defined at: 
    this(T x, T y) { this.x = x; this.y = y; }
    V abs() => V(x < 0 ? -x : x, y < 0 ? -y : y); 
    T angle() const => .atan2(y, x);
    T angle_to(const V to) const => .atan2(cross(to), dot(to));
    T angle_to_point(const V to) const => (to - this).angle();
    T aspect() const => width / height;
    V bezier_derivative(const V c1, const V c2, const V end, T p_t) const 
       => V( .bezier_derivative(x, c1.x, c2.x, end.x, p_t), 
             .bezier_derivative(y, c1.y, c2.y, end.y, p_t) );
    V bezier_interpolate(const V c1, const V c2, const V end, T p_t) const 
        => V( .bezier_interpolate(x, c1.x, c2.x, end.x, p_t), 
              .bezier_interpolate(y, c1.y, c2.y, end.y, p_t) );
    V bounce(const V normal) const => -reflect(normal);
    V ceil() => V(.ceil(x), .ceil(y));
    V clamp(const V min, const V max) const => V(.clampf(x, min.x, max.x), .clampf(y, min.y, max.y));
    V clampf(T min, T max) const => V(.clampf(x, min, max), .clampf(y, min, max));
    T cross(const V p_other) const => x * p_other.y - y * p_other.x;
    V cubic_interpolate(const V b, const V pre_a, const V post_b, T weight) const
        => V( .cubic_interpolate(x, b.x, pre_a.x, post_b.x, weight),
              .cubic_interpolate(y, b.y, pre_a.y, post_b.y, weight) );
    V cubic_interpolate_in_time(const V b, const V pre_a, const V post_b, T weight, T b_t, T pre_a_t, T post_b_t) const
        => V( .cubic_interpolate_in_time(x, b.x, pre_a.x, post_b.x, weight, b_t, pre_a_t, post_b_t),
              .cubic_interpolate_in_time(y, b.y, pre_a.y, post_b.y, weight, b_t, pre_a_t, post_b_t) );
    V direction_to(const V to) const => (to - this).normalized();
    T distance_squared_to(const V v) const => (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y);
    T distance_to(const V v) const => sqrt(distance_squared_to(v));
    T dot(const V other) const => x * other.x + y * other.y;
    V floor() => V(.floor(x), .floor(y));
    static V from_angle(T angle) => V( .cos(angle), .sin(angle) );
    bool is_equal_approx(const V other) => .is_equal_approx(x, other.x) && .is_equal_approx(y, other.y);
    bool is_finite() const => .is_finite(x) && .is_finite(y);
    bool is_normalized() const => .is_equal_approx(length_squared(), cast(T)1, cast(T)GM_UNIT_EPSILON);
    bool is_zero_approx() const => .is_zero_approx(x) && .is_zero_approx(y);
    T length() const => sqrt(x * x + y * y);
    T length_squared() const => x * x + y * y;
    V lerp(const V to, T weight) const => V( .lerp(x, to.x, weight), .lerp(y, to.y, weight) );
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
    V max(const V other) const => V( x > other.x ? x : other.x, y > other.y ? y : other.y );
    int max_axis_index() const => x < y ? AXIS_Y : AXIS_X;
    V maxf(T v) const => V( x > v ? x : v, y > v ? y : v );    
    V min(const V other) const => V( x < other.x ? x : other.x, y < other.y ? y : other.y );
    int min_axis_index() const => x < y ? AXIS_X : AXIS_Y;
    V minf(T v) const => V( x < v ? x : v, y < v ? y : v );
    V move_toward(const V to, T delta) const 
    {
        V v = this;
        V vd = to - v;
        T len = vd.length();
        return len <= delta || len < cast(T)GM_CMP_EPSILON ? to : v + vd / len * delta;
    }
    void normalize()
    {
        T l = x * x + y * y;
        if (l != 0) 
        {
            l = sqrt(l);
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
    V posmod(T mod) const => V(.fposmod(x, mod), .fposmod(y, mod));
    V posmodv(const V modv) const => V(.fposmod(x, modv.x), .fposmod(y, modv.y));
    V project(const V to) const => to * (dot(to) / to.length_squared());
    V reflect(const V normal) const
    {
        assert(normal.is_normalized());
        return normal * 2 * dot(normal) - this;
    }
    V rotated(T by_radians) const 
    {
        T sine = .sin(by_radians);
        T cosi = .cos(by_radians);
        return V(x * cosi - y * sine,
                 x * sine + y * cosi);
    }
    V round() const => V(.round(x), .round(y));
    V sign() const  => V(.sign(x), .sign(y));

    V slerp(const V to, T weight) const 
    {
        T start_length_sq = length_squared();
	    T end_length_sq = to.length_squared();
	    if (start_length_sq == 0 || end_length_sq == 0) 
        {
		    // Zero length vectors have no angle, so the best we can do is either lerp or throw an error.
		    return lerp(to, weight);
	    }
	    T start_length = .sqrt(start_length_sq);
	    T result_length = .lerp(start_length, .sqrt(end_length_sq), weight);
	    T angle = angle_to(to);
	    return rotated(angle * weight) * (result_length / start_length);
    }

    // operators
    ref inout(T) opIndex(size_t n) inout return { assert(n < 2); return n ? y : x; }

    V opBinary(string op)(const V v)   const if (op == "*") => V(x * v.x  , y * v.y  );
    V opBinary(string op)(float scale) const if (op == "*") => V(x * scale, y * scale);
    V opBinary(string op)(int scale)   const if (op == "*") => V(x * scale, y * scale);
    V opBinary(string op)(const V v)   const if (op == "+") => V(x + v.x  , y + v.y  );
    V opBinary(string op)(float add)   const if (op == "+") => V(x + add  , y + add  );
    V opBinary(string op)(const V v)   const if (op == "-") => V(x - v.x  , y - v.y  );
    V opBinary(string op)(const V v)   const if (op == "/") => V(x / v.x  , y / v.y  );
    V opBinary(string op)(float scale) const if (op == "/") => V(x / scale, y / scale);
    V opBinary(string op)(int scale)   const if (op == "/") => V(x / scale, y / scale);

    V opOpAssign(string op)(const V v)   if (op == "*") { x *= v.x;   y *= v.y;   return this; }
    V opOpAssign(string op)(float scale) if (op == "*") { x *= scale; y *= scale; return this; }
    V opOpAssign(string op)(int scale)   if (op == "*") { x *= scale; y *= scale; return this; }
    V opOpAssign(string op)(const V v)   if (op == "+") { x += v.x;   y += v.y;   return this; }
    V opOpAssign(string op)(float add)   if (op == "+") { x += add;   y += add;   return this; }
    V opOpAssign(string op)(const V v)   if (op == "-") { x -= v.x;   y -= v.y;   return this; }
    V opOpAssign(string op)(const V v)   if (op == "/") { x /= v.x;   y /= v.y;   return this; }
    V opOpAssign(string op)(float scale) if (op == "/") { x /= scale; y /= scale; return this; }
    V opOpAssign(string op)(int scale)   if (op == "/") { x /= scale; y /= scale; return this; }
    
    V opUnary(string op)() const if (op == "+") => this;    
    V opUnary(string op)() const if (op == "-") => V(-x, -y);
}

alias Vector2 = Vector2Impl!float;
alias Vector2d = Vector2Impl!double;

unittest
{
    Vector2 a = Vector2(1.2f, 4);
    a += Vector2.ONE;
    assert(a.is_equal_approx( Vector2(2.2f, 5)));
    a = a * Vector2.ZERO;
    assert(a.is_zero_approx);

    Vector2d b = Vector2d(2.5, 4.0) + 4.0;
    b /= 2.0;
    assert(b.is_equal_approx( Vector2d(3.25, 4.0)));
    b = b * Vector2d.ZERO;
    assert(b.is_zero_approx);
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
    bool is_finite() const => .is_finite(x) && .is_finite(y) && .is_finite(z);
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
    Vector3 sign() const => Vector3(.sign(x), .sign(y), .sign(z));

    float signed_angle_to(const Vector3 to, const Vector3 axis)
    {
        Vector3 cross_to = cross(to);
        float unsigned_angle = atan2(cross_to.length(), dot(to));
        float sign = cross_to.dot(axis);
        return (sign < 0) ? -unsigned_angle : unsigned_angle;
    }

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



