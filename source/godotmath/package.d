/*
    Copyright (c) 2014-2025 Godot Engine contributors.
    Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur. 
    Copyright (c) 2025 Guillaume Piolat
*/
module godotmath;

public import godotmath.funcs;

import numem.core.traits;

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

// TODO direction_to

// Implementation done for: 
// - Vector2 => https://docs.godotengine.org/en/stable/classes/class_vector2.html
// - Vector3 => https://docs.godotengine.org/en/stable/classes/class_vector3.html
// - Vector2i => https://docs.godotengine.org/en/stable/classes/class_vector2i.html
// - Vector3i => https://docs.godotengine.org/en/stable/classes/class_vector3i.html
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
        enum bool isFloat = isFloatingPoint!T;
        static if (isFloat)        
            alias F = T;
        else 
            alias F = float;
        alias Elem = T;
    }
    
    union { T x = 0; T width;  }
    union { T y = 0; T height; }

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
        V bezier_derivative(const V c1, const V c2, const V end, T p_t) const 
           => V( gm_bezier_derivative(x, c1.x, c2.x, end.x, p_t), 
                 gm_bezier_derivative(y, c1.y, c2.y, end.y, p_t) );
        V bezier_interpolate(const V c1, const V c2, const V end, T p_t) const 
            => V( gm_bezier_interpolate(x, c1.x, c2.x, end.x, p_t), 
                  gm_bezier_interpolate(y, c1.y, c2.y, end.y, p_t) );
        V bounce(const V normal) const => -reflect(normal);
        V ceil() => V(gm_ceil(x), gm_ceil(y));
    }

    V clamp(const V min, const V max) const => V(gm_clamp(x, min.x, max.x), gm_clamp(y, min.y, max.y));
    V clamp(T min, T max) const => V(gm_clamp(x, min, max), gm_clamp(y, min, max));
    
    static if (isFloat)
    {
        T cross(const V p_other) const => x * p_other.y - y * p_other.x;
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

    // operators
    ref inout(T) opIndex(size_t n) inout return { assert(n < 2); return n ? y : x; }

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
    V opBinary(string op)(T scale)   const if (op == "%") => V(x % scale, y % scale);

    V opOpAssign(string op)(const V v) if (op == "*") { x *= v.x;   y *= v.y;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "*") { x *= scale; y *= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "+") { x += v.x;   y += v.y;   return this; }
    V opOpAssign(string op)(T add)     if (op == "+") { x += add;   y += add;   return this; }
    V opOpAssign(string op)(const V v) if (op == "-") { x -= v.x;   y -= v.y;   return this; }
    V opOpAssign(string op)(T sub)     if (op == "-") { x -= sub;   y -= sub;   return this; }
    V opOpAssign(string op)(const V v) if (op == "/") { x /= v.x;   y /= v.y;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "/") { x /= scale; y /= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "%") { x %= v.x;   y %= v.y;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "%") { x %= scale; y %= scale; return this; }
    
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
        enum bool isFloat = isFloatingPoint!T;
        static if (isFloat)        
            alias F = T;
        else 
            alias F = float;
        alias Elem = T;
    }

    union { T x = 0; T width;  }
    union { T y = 0; T height; }
    union { T z = 0; T depth;  }

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
    
        BasisImpl!T outer(const V p_with) const 
        {
            BasisImpl!T basis;
            basis.rows[0] = V(x * p_with.x, x * p_with.y, x * p_with.z);
            basis.rows[1] = V(y * p_with.x, y * p_with.y, y * p_with.z);
            basis.rows[2] = V(z * p_with.x, z * p_with.y, z * p_with.z);
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
 
    // operators
    ref inout(T) opIndex(size_t n) inout return { assert(n < 3); switch (n) { case 0: return x; case 1: return y; default: return z; } }

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
    V opBinary(string op)(T scale)   const if (op == "%") => V(x % scale, y % scale, z % scale);

    V opOpAssign(string op)(const V v) if (op == "*") { x *= v.x;   y *= v.y;   z *= v.z;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "*") { x *= scale; y *= scale; z *= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "+") { x += v.x;   y += v.y;   z += v.z;   return this; }
    V opOpAssign(string op)(T add)     if (op == "+") { x += add;   y += add;   z += add;   return this; }
    V opOpAssign(string op)(const V v) if (op == "-") { x -= v.x;   y -= v.y;   z -= v.z;   return this; }
    V opOpAssign(string op)(T sub)     if (op == "-") { x -= sub;   y -= sub;   z -= sub;   return this; }
    V opOpAssign(string op)(const V v) if (op == "/") { x /= v.x;   y /= v.y;   z /= v.z;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "/") { x /= scale; y /= scale; z /= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "%") { x %= v.x;   y %= v.y;   z %= v.z;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "%") { x %= scale; y %= scale; z %= scale; return this; }
    
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

        enum bool isFloat = isFloatingPoint!T;
        static if (isFloat)
            alias F = T;
        else 
            alias F = float;
        alias Elem = T;
    }

    union { T x = 0; T width;  }
    union { T y = 0; T height; }
    union { T z = 0; T depth;  }
    union { T w = 0; }

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
 
    // operators
    ref inout(T) opIndex(size_t n) inout return 
    { 
        assert(n < 4); 
        switch (n) 
        { 
            case 0: return x; 
            case 1: return y; 
            case 2: return z; 
            default: return w; 
        } 
    }

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
    V opBinary(string op)(T scale)   const if (op == "%") => V(x % scale, y % scale, z % scale, w % scale);

    V opOpAssign(string op)(const V v) if (op == "*") { x *= v.x;   y *= v.y;   z *= v.z;   w *= v.w;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "*") { x *= scale; y *= scale; z *= scale; w *= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "+") { x += v.x;   y += v.y;   z += v.z;   w += v.w;   return this; }
    V opOpAssign(string op)(T add)     if (op == "+") { x += add;   y += add;   z += add;   w += add;   return this; }
    V opOpAssign(string op)(const V v) if (op == "-") { x -= v.x;   y -= v.y;   z -= v.z;   w -= v.w;   return this; }
    V opOpAssign(string op)(T sub)     if (op == "-") { x -= sub;   y -= sub;   z -= sub;   w -= sub;   return this; }
    V opOpAssign(string op)(const V v) if (op == "/") { x /= v.x;   y /= v.y;   z /= v.z;   w /= v.w;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "/") { x /= scale; y /= scale; z /= scale; w /= scale; return this; }
    V opOpAssign(string op)(const V v) if (op == "%") { x %= v.x;   y %= v.y;   z %= v.z;   w %= v.w;   return this; }
    V opOpAssign(string op)(T scale)   if (op == "%") { x %= scale; y %= scale; z %= scale; w %= scale; return this; }

    V opUnary(string op)() const if (op == "+") => this;
    V opUnary(string op)() const if (op == "-") => V(-x, -y, -z, -w);
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

    T determinant() const
    {
        return rows[0][0] * (rows[1][1] * rows[2][2] - rows[2][1] * rows[1][2]) -
               rows[1][0] * (rows[0][1] * rows[2][2] - rows[2][1] * rows[0][2]) +
               rows[2][0] * (rows[0][1] * rows[1][2] - rows[1][1] * rows[0][2]);
    }

    static B from_euler(V3 euler, int order = EULER_ORDER_YXZ)
    {
        B b;
        b.set_euler(euler, order);
        return b;
    }

    // operators
    V3 opIndex(size_t n) const => rows[n];
    B opBinary(string op)(const B m) const if (op == "*")
        => B(m.tdotx(rows[0]), m.tdoty(rows[0]), m.tdotz(rows[0]),
             m.tdotx(rows[1]), m.tdoty(rows[1]), m.tdotz(rows[1]),
             m.tdotx(rows[2]), m.tdoty(rows[2]), m.tdotz(rows[2]));

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
	T tdotx(const V3 v) const =>
		rows[0][0] * v[0] + rows[1][0] * v[1] + rows[2][0] * v[2];

	T tdoty(const V3 v) const =>
		rows[0][1] * v[0] + rows[1][1] * v[1] + rows[2][1] * v[2];

	T tdotz(const V3 v) const =>
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



// internal

private:

enum bool isVector2Impl(T) = is(T : Vector2Impl!U, U...);
enum bool isVector3Impl(T) = is(T : Vector3Impl!U, U...);
enum bool isVector4Impl(T) = is(T : Vector4Impl!U, U...);
