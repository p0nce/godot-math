# godot-math

[`godot-math`]() is a D port of Godot's linear algebra with unchanged semantics. 

 ➡️ [See Godot Math documentation](https://docs.godotengine.org/en/stable/tutorials/math/index.html)


```json
"dependencies":
{
    "godot-math": "~>1.0"
}
```

```perl
dependency "godot-math" version="~>1.0"
```

## Features

- one .d file
- Everything is fully `pure nothrow @nogc @safe`.
- LDC, GDC, DMD

| Type  | Status | Description | Documentation                  |
|-------|--------|-------------|--------------------------------|
| Vector2   | ✅ | 2D `float` vector, or size, or point | [Vector2](https://docs.godotengine.org/en/stable/classes/class_vector2.html) |
| Vector3   | ✅ | 3D `float` vector, or size, or point | [Vector3](https://docs.godotengine.org/en/stable/classes/class_vector3.html) |
| Vector4   | ✅ | 4D `float` vector or point | [Vector4](https://docs.godotengine.org/en/stable/classes/class_vector4.html) |
| Vector2i  | ✅ | 2D `int` vector, or size, or point | [Vector2i](https://docs.godotengine.org/en/stable/classes/class_vector2i.html) |
| Vector3i  | ✅ | 3D `int` vector, or size, or point | [Vector3i](https://docs.godotengine.org/en/stable/classes/class_vector3i.html) |
| Vector4i  | ✅ | 4D `int` vector or point | [Vector4i](https://docs.godotengine.org/en/stable/classes/class_vector4i.html) |
| Rect2  | ✅ | 2D `float` rectangle | [Rect2](https://docs.godotengine.org/en/stable/classes/class_rect2.html) |
| Rect2i  | ✅ | 2D `int` rectangle | [Rect2i](https://docs.godotengine.org/en/stable/classes/class_rect2i.html) |
| Transform2D | ✅ | 2x3 matrix | [Transform2D](https://docs.godotengine.org/en/stable/classes/class_transform2d.html) |
| Basis | ✅ | 3x3 matrix | [Basis](https://docs.godotengine.org/en/stable/classes/class_basis.html)
| Transform3D | ✅ | 3x4 matrix | [Transform3D](https://docs.godotengine.org/en/stable/classes/class_transform3d.html) |
| Projection | ✅ | 4x4 matrix | [Projection](https://docs.godotengine.org/en/stable/classes/class_projection.html) |
| Plane | ❌ | 3D Plane matrix | [Plane](https://docs.godotengine.org/en/stable/classes/class_plane.html) |
| AABB | ❌ | 3D bounding box | [AABB](https://docs.godotengine.org/en/stable/classes/class_aabb.html) |



### Math functions

Note: When a `float` function exist, a `double` overload also exist with identical name.

```d
float gm_abs(float x);
int   gm_abs(int x);
float gm_acos(float x);
float gm_acosh(float x);
float gm_angle_difference(float from, float to);
float gm_asin(float x);
float gm_asinh(float x);
float gm_atan(float x);
float gm_atan2(float y, float x);
float gm_bezier_derivative(float start, float control_1, float control_2, float end, float t);
float gm_bezier_interpolate(float start, float control_1, float control_2, float end, float t);
float gm_ceil(float x);
float gm_clamp(float value, float min, float max);
int   gm_clamp(int value, int min, int max);
float gm_cos(float x);
float gm_cubic_interpolate(float from, float to, float pre, float post, float weight);
float gm_cubic_interpolate_in_time(float from, float to, float pre, float post, float weight, float to_t, float pre_t, float post_t);
float gm_deg_to_rad(float y);
float gm_floor(float x);
float gm_fmod(float x, float y);
float gm_fposmod(float x, float y);
bool  gm_is_equal_approx(float left, float right);
bool  gm_is_equal_approx(float left, float right, float tolerance);
bool  gm_is_finite(float x);
bool  gm_is_zero_approx(float value);
float gm_lerp(float from, float to, float weight);
float gm_lerp_angle(float from, float to, float weight);
bool  gm_likely(bool b);
bool  gm_unlikely(bool b);
float gm_rad_to_deg(float y);
float gm_round(float x);
int   gm_sign(int v);
float gm_sign(float v);
bool  gm_signbit(float num);
bool  gm_signbit(double num);
float gm_sin(float x);
float gm_snapped(float value, float step);
int   gm_snapped(int value, int step);
float gm_sqrt(float x);
void  gm_swap(T)(ref T a, ref T b);
float gm_tan(float x);
```



## Differing names

- Both `float` and `double` versions exist at once, with a `***d` suffix (eg: `Projectiond`, `Vector2d`, `Rect2d`).
- In global scope, function symbols get a `gm_` prefix.
- `.clampf`/`.clampi` replaced by overloaded `.clamp`
  Likewise `.snappedf` / `.snappedi` replaced by overloaded `.snapped`.
  Same for `minf`/`mini`, `maxf`/`maxi` => replaced by `min`/`max`
// - in Godot names "f" would mean "float", and "float" is `double` in Godot.
//   Here when something is "float" it is actually `float` in the interface too,
//   and same for `double`.


## Small semantic differences

- `Rect2`/`Rect2i`/`Rect2d` have a `.merge_non_empty` method that, in case of a union with an empty rectangle, return the other rectangle.
- Bonus methods for rectangles, such has `.left`, `.top`, `.right`, `.bottom`.
- `Projection` inverse is less precise than in original Godot. You can always go double to get more precise `.inverse()`