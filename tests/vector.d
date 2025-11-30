module vector;

import godotmath;

@("Vector basics")
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

    Vector2i a2 = Vector2i(0, 1);
    immutable Vector2i b2 = Vector2i(0, 0);
    assert(b2[0] == 0 && b2[1] == 0);
    Vector2i c = [0, 1];
    float[2] arr2 = [4.0f, 1];
    Vector2 d = arr2;
    assert(a2 == c);
    assert(cast(Vector2d)a2 != Vector2d.ZERO);
}

@("Octahedron encoding")
unittest
{
    Vector3 up = Vector3.UP;
    Vector3 down = Vector3.DOWN;
    Vector3 A = Vector3.octahedron_decode(up.octahedron_encode);
    Vector3 B = Vector3.octahedron_decode(down.octahedron_encode);
    assert(up.is_equal_approx(A));
    assert(down.is_equal_approx(B));
}

@("Vector slerp")
unittest
{
    Vector3 s = Vector3.UP.slerp(Vector3.LEFT, 0.5);
    Vector3d sd = Vector3d.UP.slerp(Vector3d.LEFT, 0.5);
    assert(s.is_equal_approx(Vector3(-1, 1, 0).normalized));
    assert(sd.is_equal_approx(Vector3d(-1, 1, 0).normalized));
}



    /*
    Vector4i x = [4, 5, 6, 7];
    assert(x == x);
    --x[0];
    assert(x[0] == 3);
    ++x[0];
    assert(x[0] == 4);
    x[1] &= 1;
    x[2] = 77 + x[2];
    x[3] += 3;
    assert(x == [4, 1, 83, 10]);
    */
/*
    assert(Vector2i(4, 5) + 1 == Vector2i(5,6));
    assert(Vector2i(4, 5) - 1 == Vector2i(3,4));
    assert(1 + Vector2i(4, 5) == Vector2i(5,6));
    assert(Vector3(1,1,1) * 0 == 0);
    assert(1.0 * Vector3d(4,5,6) == Vector3(4,5.0f,6.0));

    auto dx = Vector2i(1,2);
    auto dy = Vector2i(4,5);
    auto dp = dx.dot(dy);
    assert(dp == 14 );

    Vector3i h = cast(Vector3i)(Vector3d(0.5, 1.1, -2.2));
    assert(h == [0, 1, -2]);

    //assert(h[] == [0, 1, -2]);
    //assert(h[1..3] == [1, -2]);
    assert(-h[1] == -5);
    assert(++h[0] == 3);
    assert(vec3i(-1, 0, 2).abs == vec3i(1, 0, 2));*/

