# godot-math

`godot-math` is a D port of Godot's linear algebra with unchanged semantics.


```json
"dependencies":
{
    "godot-math": "~>1.0"
}
```

```
dependency "godot-math" version="~>1.0"
```

## Status

All types are fully `pure nothrow @nogc @safe`.

| Type  | Status | Description | Documentation                  |
|-------|--------|-------------|--------------------------------|
| Vector2   | ✅ | 2D `float` vector, or size, or point | [Vector2](https://docs.godotengine.org/en/stable/classes/class_vector2.html) |
| Vector3   | ✅ | 3D `float` vector, or size, or point | [Vector3](https://docs.godotengine.org/en/stable/classes/class_vector3.html) |
| Vector4   | ✅ | 4D `float` vector or point | [Vector4](https://docs.godotengine.org/en/stable/classes/class_vector4.html) |
| Vector2i  | ✅ | 2D `int` vector, or size, or point | [Vector2i](https://docs.godotengine.org/en/stable/classes/class_vector2i.html) |
| Vector3i  | ✅ | 3D `int` vector, or size, or point | [Vector3i](https://docs.godotengine.org/en/stable/classes/class_vector3i.html) |
| Vector4i  | ✅ | 4D `int` vector or point | [Vector4i](https://docs.godotengine.org/en/stable/classes/class_vector4i.html) |



### Global function implemented

Both `float` and `double` overload exist, with identical names.
```d
    float gm_abs(float x);
    int gm_abs(int x);
    float  gm_acos(float x);
    float  gm_acosh(float x);
    float gm_angle_difference(float from, float to);
    float  gm_asin(float x);
    float  gm_asinh(float x);
    float  gm_atan(float x);
    float  gm_atan2(float y, float x);
    float gm_bezier_derivative(float start, float control_1, float control_2, float end, float t);
    float gm_bezier_interpolate(float start, float control_1, float control_2, float end, float t);
    float  gm_ceil(float x);
    float gm_clamp(float value, float min, float max);
    int gm_clamp(int value, int min, int max);
    float  gm_cos(float x);
    float gm_cubic_interpolate(float from, float to, float pre, float post, float weight);
    float gm_cubic_interpolate_in_time(float from, float to, float pre, float post, float weight, float to_t, float pre_t, float post_t);
    float gm_deg_to_rad(float y);
    float  gm_floor(float x);
    float gm_fmod(float x, float y);
    float gm_fposmod(float x, float y);
    bool gm_is_equal_approx(float left, float right);
    bool gm_is_equal_approx(float p_left, float p_right, float p_tolerance);
    bool gm_is_finite(float x);
    bool gm_is_zero_approx(float  value);
    float gm_lerp(float from, float to, float weight);
    float gm_lerp_angle(float from, float to, float weight) => from + gm_angle_difference(from, to) * weight; ///
    bool gm_likely(bool b);
    bool gm_unlikely(bool b);
    float gm_rad_to_deg(float y);
    float  gm_round(float x);
    int    gm_sign(int v);
    float  gm_sign(float v);
    bool gm_signbit(float num);
    bool gm_signbit(double num);
    float  gm_sin(float x);
    float gm_snapped(float value, float step);
    int gm_snapped(int value, int step);
    float  gm_sqrt(float x);
    void gm_swap(T)(ref T a, ref T b);
    float gm_tan(float x);
```


// Implementation of full public API done for: 
// - Vector2 => https://docs.godotengine.org/en/stable/classes/class_vector2.html
// - Vector3 => https://docs.godotengine.org/en/stable/classes/class_vector3.html
// - Vector4 => https://docs.godotengine.org/en/stable/classes/class_vector3.html
// - Vector2i => https://docs.godotengine.org/en/stable/classes/class_vector2i.html
// - Vector3i => https://docs.godotengine.org/en/stable/classes/class_vector3i.html
// - Vector4i => https://docs.godotengine.org/en/stable/classes/class_vector4i.html
// - Quaternion => https://docs.godotengine.org/en/stable/classes/class_quaternion.html
// - Basis => https://docs.godotengine.org/en/stable/classes/class_basis.html
//
// Note: In case of semantic conflict, the way Godot does it is favoured.



## Differing semantics

- `Rect2`/`Rect2i`/`Rect2d` have a `.merge_non_empty` method that, in case of a union with an empty rectangle, return the other rectangle.
- Bonus methods for rectangles, such has `.left`, `.top`, `.right`, `.bottom`.
- `Projection` inverse is less precise than in original Godot. You can always go double to get more precise `.inverse()`