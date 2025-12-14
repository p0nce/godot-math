module transform2d;

import godotmath;

pure nothrow @nogc @safe:

@("Transform2D * Vector2")
unittest
{
    Vector2 A = Vector2(4, 5);
    Transform2D T = Transform2D(0, Vector2(4, 1)); 
    A = T * A;
    assert(A.is_equal_approx(Vector2(8, 6)));

    T.set_scale(Vector2(-1, 1));
    T.set_rotation(-GM_PI);
    T.origin = Vector2.ZERO;
    A = Vector2(4, 5) * T;
    assert(A.is_equal_approx(Vector2(-4, 5)));
}

@("Vector2 * Transform2D")
unittest
{
    Vector2 A = Vector2(4, 5);
    Transform2D T = Transform2D(GM_PI, Vector2(4, 1));    
    A = A * T;
    assert(A.is_equal_approx(Vector2(0, -4)));

}