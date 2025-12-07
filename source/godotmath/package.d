/*
    Copyright (c) 2014-2025 Godot Engine contributors.
    Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur. 
    Copyright (c) 2025 Guillaume Piolat
*/
module godotmath;

public import godotmath.funcs;

//import numem.core.traits;

pure nothrow @nogc @safe:

// Changes vs Godot
// - in global scope, symbols get a gm_ prefix.
// - .clampf/.clampi replaced by overloaded .clamp
// - .snappedf /.snappedi replaced by overloaded .snapped
// - same for minf/mini/maxf/maxi => replaced by min/max

// Provide both float and double versions should the need arise.
alias Vector2  = Vector2Impl!float;
alias Vector2d = Vector2Impl!double;
alias Vector2i = Vector2Impl!int;

alias Vector3  = Vector3Impl!float;
alias Vector3d = Vector3Impl!double;
alias Vector3i = Vector3Impl!int;

alias Vector4  = Vector4Impl!float;
alias Vector4d = Vector4Impl!double;
alias Vector4i = Vector4Impl!int;

alias Basis    = BasisImpl!float;
alias Basisd   = BasisImpl!double;

alias Quaternion  = QuaternionImpl!float;
alias Quaterniond = QuaternionImpl!double;


// EulerOrder
alias EulerOrder = int;
enum : EulerOrder
{
    GM_EULER_ORDER_XYZ = 0,
    GM_EULER_ORDER_XZY = 1,
    GM_EULER_ORDER_YXZ = 2,
    GM_EULER_ORDER_YZX = 3,
    GM_EULER_ORDER_ZXY = 4,
    GM_EULER_ORDER_ZYX = 5,
}


// Implementation done for: 
// - Vector2 => https://docs.godotengine.org/en/stable/classes/class_vector2.html
// - Vector3 => https://docs.godotengine.org/en/stable/classes/class_vector3.html
// - Vector4 => https://docs.godotengine.org/en/stable/classes/class_vector3.html
// - Vector2i => https://docs.godotengine.org/en/stable/classes/class_vector2i.html
// - Vector3i => https://docs.godotengine.org/en/stable/classes/class_vector3i.html
// - Vector4i => https://docs.godotengine.org/en/stable/classes/class_vector4i.html
//
// Note: In case of semantic conflict, the way Godot does it is favoured.

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
        bool is_equal_approx(const V other) => gm_is_equal_approx(x, other.x) && gm_is_equal_approx(y, other.y);
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
        void normalize()
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
        void normalize()
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
    if (is(T : float) || is(T : double))
{
pure nothrow @nogc @safe:

    private 
    {
        alias V = Vector4Impl!T,
              V3 = Vector3Impl!T;

        enum bool isFloat = (is(T : float) || is(T : double));
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
        void normalize()
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
    if (is(T : float) || is(T : double))
{
pure nothrow @nogc @safe:

    alias Q = QuaternionImpl!T;
    alias V3 = Vector3Impl!T;
    alias Elem = T;

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

    this(const Vector3 axis, T angle)
    {
        assert(axis.is_normalized);

        T d = p_axis.length();
        if (d == 0) 
        {
            x = 0;
            y = 0;
            z = 0;
            w = 0;
        } 
        else 
        {
            T sin_angle = gm_sin(p_angle * 0.5f);
            T cos_angle = gm_sin(p_angle * 0.5f);
            T s = sin_angle / d;
            x = p_axis.x * s;
            y = p_axis.y * s;
            z = p_axis.z * s;
            w = cos_angle;
        }
    }

    this(BasisImpl!T from);

    this(T x, T y, T z, T w)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    T dot(const Q q) const => x * q.x + y * q.y + z * q.z + w * q.w;

    bool is_normalized() const => gm_is_equal_approx(length_squared(), cast(T)1, cast(T)GM_UNIT_EPSILON);

    T length() const => gm_sqrt(length_squared());
    T length_squared() const => dot(this);

    void normalize()
    {
        this /= length();
    }

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


    // operators
    Q opBinary(string op)(const Q v) const if (op == "+") => V(x + v.x  , y + v.y  , z + v.z  , w + v.w  );
    Q opBinary(string op)(const Q v) const if (op == "-") => V(x - v.x  , y - v.y  , z - v.z  , w - v.w  );

    Q opOpAssign(string op)(const Q v) if (op == "+") { x += v.x;   y += v.y;   z += v.z;   w += v.w;   return this; }
    Q opOpAssign(string op)(const Q v) if (op == "-") { x -= v.x;   y -= v.y;   z -= v.z;   w -= v.w;   return this; }
    Q opOpAssign(string op)(T s) if (op == "*")       { x *= s;     y *= s;     z *= s;     w *= s;     return this; }
    Q opOpAssign(string op)(T s) if (op == "/") => this *= (1 / s);

    Q opUnary(string op)() const if (op == "+") => this;    
    Q opUnary(string op)() const if (op == "-") => Q(-x, -y, -z, -w);
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
{
pure nothrow @nogc @safe:

    alias V3 = Vector3Impl!T;
    alias B = BasisImpl!T;
    alias Q = QuaternionImpl!T;
    alias Elem = T;

    V3[3] rows = 
    [
        V3(1, 0, 0),
        V3(0, 1, 0),
        V3(0, 0, 1)
    ];

    this(V3 axis, T angle)
    {
        set_axis_angle(axis, angle);
    }

    this(V3 x, V3 y, V3 z)
    {
        rows[0] = x;
        rows[1] = y;
        rows[2] = z;
    }

    this(T xx, T xy, T xz, 
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

    static B from_scale(const Vector3 scale)
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


// TODO Quaternion get_rotation_quaternion() const

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

    void invert() 
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

    bool is_equal_approx(const B b) const => rows[0].is_equal_approx(b.rows[0]) && rows[1].is_equal_approx(b.rows[1]) && rows[2].is_equal_approx(b.rows[2]);
    bool is_finite() const => rows[0].is_finite() && rows[1].is_finite() && rows[2].is_finite();

    static B looking_at(V3 target, V3 up = V3(0, 1, 0), bool use_model_front)
    {
        assert(!target.is_zero_approx());
        assert(!up.is_zero_approx());        
        V3 z = target.normalized();
        if (!use_model_front) {
            z = -z;
        }
        V3 x = up.cross(z);
        if (x.is_zero_approx()) {
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

    void orthonormalize() 
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

    void rotate(const V3 axis, T angle) { this = rotated(axis, angle); }
    B rotated(const V3 axis, T angle) const => B(axis, angle) * this;

    void scale(const V3 scale) 
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

    void scale_local(const V3 scale) 
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

    void set_axis_angle(const V3 axis, T angle)
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

    void set_column(int index, const V3 c)
    {
        rows[0][index] = c.x;
        rows[1][index] = c.x;
        rows[2][index] = c.x;
    }

    void set_columns(const V3 x, const V3 y, const V3 z)
    {
        set_column(0, x);
        set_column(1, y);
        set_column(2, z);
    }

    void set_euler(const V3 euler, EulerOrder order) 
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
        rows[0][0] * v[0] + rows[1][0] * v[1] + rows[2][0] * v[2];

    T tdoty(const V3 v) const =>
        rows[0][1] * v[0] + rows[1][1] * v[1] + rows[2][1] * v[2];

    T tdotz(const V3 v) const =>
        rows[0][2] * v[0] + rows[1][2] * v[1] + rows[2][2] * v[2];

    void transpose() 
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

    U opCast(U)() const if (isBasisImpl!U)
    {    
        static if (is(U.Elem == float))
        {
            alias VD = Vector3Impl!float;
            return U(cast(VD)rows[0], cast(VD)rows[1], cast(VD)rows[1]);
        }
        else static if (is(U.Elem == double))
        {
            alias VD = Vector3Impl!double;
            return U(cast(VD)rows[0], cast(VD)rows[1], cast(VD)rows[1]);
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
}






// internal

private:

enum bool isVector2Impl(T) = is(T : Vector2Impl!U, U...);
enum bool isVector3Impl(T) = is(T : Vector3Impl!U, U...);
enum bool isVector4Impl(T) = is(T : Vector4Impl!U, U...);
enum bool isQuaternionImpl(T) = is(T : QuaternionImpl!U, U...);
enum bool isBasisImpl(T)   = is(T : BasisImpl!U, U...);
