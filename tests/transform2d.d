module transform2d;

import godotmath;

pure nothrow @nogc @safe:

@("Transform2D * Vector2")
unittest
{
    Vector2 A = Vector2(4, 5);
    float angle = 0;
    Transform2D T = Transform2D(angle, Vector2(4, 1));
    A = T * A;
    assert(A.is_equal_approx(Vector2(8, 6)));
}

@("Vector2 * Transform2D")
unittest
{
    Vector2 A = Vector2(4, 5);
    Transform2D T = Transform2D(GM_PI, Vector2(4, 1));    
    A = A * T;
    assert(A.is_equal_approx(Vector2(0, -4)));

}

@("Transform2D interpolate_with")
unittest
{
    Transform2D T1 = Transform2D(0, Vector2(0, 0));
    Transform2D T2 = Transform2D(GM_PI / 2.0f, Vector2(10, 10));
    
    Transform2D interpolated = T1.interpolate_with(T2, 0.5f);
    
    // Interpolated rotation should be halfway
    assert(gm_is_equal_approx(interpolated.get_rotation(), cast(float)(GM_PI / 4.0f)));
    // Interpolated origin should be halfway
    assert(interpolated.get_origin().is_equal_approx(Vector2(5, 5)));
}