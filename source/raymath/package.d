module raymath;

import godotmath;


pure nothrow @nogc @safe:

/// Clamp float value
float Clamp(float value, float min, float max)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

/// Calculate linear interpolation between two floats
float Lerp(float start, float end, float amount)
{
   return start + (end - start) * amount;
}

/// Normalize input value within input range
float Normalize(float value, float start, float end)
{
   return (value - start) / (end - start);
}

/// Remap input value within input range to output range
float Remap(float value, float inputStart, float inputEnd, float outputStart, float outputEnd)
{
    float t = Normalize(value, inputStart, inputEnd);
    return Lerp(outputStart, outputEnd, t);
}

/// Wrap input value from min to max
float Wrap(float value, float min, float max)
{
    return value - (max - min)*gm_floor((value - min)/(max - min));
}

/// Check whether two given floats are almost equal
int FloatEquals(float x, float y)
{
    return gm_is_equal_approx(x, y);
}

/* 
██╗   ██╗███████╗ ██████╗████████╗ ██████╗ ██████╗ ██████╗ 
██║   ██║██╔════╝██╔════╝╚══██╔══╝██╔═══██╗██╔══██╗╚════██╗
██║   ██║█████╗  ██║        ██║   ██║   ██║██████╔╝ █████╔╝
╚██╗ ██╔╝██╔══╝  ██║        ██║   ██║   ██║██╔══██╗██╔═══╝ 
╚████╔╝ ███████╗╚██████╗   ██║   ╚██████╔╝██║  ██║███████╗
╚═══╝  ╚══════╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚══════╝
*/

/// Vector with components value 0.0f
Vector2 Vector2Zero() => Vector2(0, 0);

/// Vector with components value 1.0f
Vector2 Vector2One() => Vector2(1, 1);                                                   

/// Add two vectors (v1 + v2)
Vector2 Vector2Add(Vector2 v1, Vector2 v2) => v1 + v2;

/// Add vector and float value
Vector2 Vector2AddValue(Vector2 v, float add) => Vector2(v.x + add, v.y + add);

/// Subtract two vectors (v1 - v2)
Vector2 Vector2Subtract(Vector2 v1, Vector2 v2) => Vector2(v1.x - v2.x, v1.y - v2.y);

/// Subtract vector by float value
Vector2 Vector2SubtractValue(Vector2 v, float sub) => Vector2(v.x - sub, v.y - sub);

/// Calculate vector length
float Vector2Length(Vector2 v) => v.length();

/// Calculate vector square length
float Vector2LengthSqr(Vector2 v) => v.length_squared();

/// Calculate two vectors dot product
float Vector2DotProduct(Vector2 v1, Vector2 v2) => v1.dot(v2);

/// Calculate distance between two vectors
float Vector2Distance(Vector2 v1, Vector2 v2) => v1.distance_to(v2);                              

/// Calculate square distance between two vectors
float Vector2DistanceSqr(Vector2 v1, Vector2 v2) => v1.distance_squared_to(v2);

/// Calculate angle from two vectors
/// Fortunately raylib's semantic is identical with Godot here.
float Vector2Angle(Vector2 v1, Vector2 v2) => v1.angle_to(v2);

/// Scale vector (multiply by value)
Vector2 Vector2Scale(Vector2 v, float scale) => Vector2(v.x * scale, v.y * scale);

/// Multiply vector by vector
Vector2 Vector2Multiply(Vector2 v1, Vector2 v2) => v1 * v2;                            

/// Negate vector
Vector2 Vector2Negate(Vector2 v) => -v;                                           

/// Divide vector by vector
Vector2 Vector2Divide(Vector2 v1, Vector2 v2) => v1 / v2;                              

/// Normalize provided vector
Vector2 Vector2Normalize(Vector2 v) => v.normalized();

/// Transforms a Vector2 by a given Matrix
//TODO
//Vector2 Vector2Transform(Vector2 v, Matrix mat);                            

/// Calculate linear interpolation between two vectors
Vector2 Vector2Lerp(Vector2 v1, Vector2 v2, float amount) => v1 + (v2 * v1) * amount;                  

/// Calculate reflected vector to normal
/// raylib semantic: doesn't verify normal is normalized, but Godot does.
/// Warning: raylib's reflect = Godot's bounce = minus Godot's reflect.
Vector2 Vector2Reflect(Vector2 v, Vector2 normal) 
    => v.bounce(normal);

/// Rotate vector by angle
Vector2 Vector2Rotate(Vector2 v, float angle)
{
    float c = gm_cos(angle);
    float s = gm_sin(angle);
    return Vector2(c * v.x - s * v.y, s * v.x + c * v.y);
}

/// Move Vector towards target
Vector2 Vector2MoveTowards(Vector2 v, Vector2 target, float maxDistance)
{
    return v.move_toward(target, maxDistance);
}

/// Invert the given vector
Vector2 Vector2Invert(Vector2 v) => Vector2(1.0f / v.x, 1.0f / v.y);

/// Clamp the components of the vector between min and max values specified by the given vectors
Vector2 Vector2Clamp(Vector2 v, Vector2 min, Vector2 max) => v.clamp(min, max);

/// Clamp the magnitude of the vector between two min and max values
Vector2 Vector2ClampValue(Vector2 v, float min, float max) => v.clampf(min, max);

/// Check whether two given vectors are almost equal
int Vector2Equals(Vector2 p, Vector2 q) => p.is_equal_approx(q);
