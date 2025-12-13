module vector;

import godotmath;

pure nothrow @nogc @safe:

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
    assert(b.is_equal_approx( Vector2d(3.25, 4.00000001)));
    b = b * Vector2d.ZERO;
    assert(b.is_zero_approx);

    Vector2i a2 = Vector2i(0, 1);
    int sqlen = a2.length_squared();
    float len = a2.length();
    immutable Vector2i b2 = Vector2i(0, 0);

    assert(b2[0] == 0 && b2[1] == 0);
    Vector2i c = [0, 1];
    float[2] arr2 = [4.0f, 1];
    Vector2 d = arr2;
    assert(a2 == c);
    assert(cast(Vector2d)a2 != Vector2d.ZERO);

    assert(Vector2i(1, 1) != Vector2i(1, 3));
    assert(Vector2(1, 2) < Vector2(2, 1));
    assert(Vector2(2, 1) > Vector2(2, 0));
    assert(Vector2(2, 1) > Vector2(1, 0));
    assert(Vector2d(1, 2) <= Vector2d(1, 2));

    assert(Vector3i(1, 1, 1) != Vector3i(1, 3, 1));
    assert(Vector3(1, 2, 4) < Vector3(2, 1, 1));
    assert(Vector3(2, 1, 1) > Vector3(2, 1, 0));
    assert(Vector3(2, 1, 0) > Vector3(1, 4, 8));
    assert(Vector3d(1, 2, 4) <= Vector3d(1, 2, 5));

    assert(Vector4i(1, 1, 1, 0) != Vector4i(1, 3, 1, 0));
    assert(Vector4(1, 2, 4, 0) < Vector4(2, 1, 1, 0));
    assert(Vector4(2, 1, 1, 1) > Vector4(2, 1, 1, 0));
    assert(Vector4(2, 1, 0, 0) > Vector4(1, 4, 8, 0));
    assert(Vector4d(1, 2, 4, 0) <= Vector4d(1, 2, 5, 0));
} 

    
// More tests from GFM

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
    assert(Vector3i(-1, 0, 2).abs == Vector3i(1, 0, 2));*/


@("Vector casts")
unittest 
{
    cast(Vector2) cast(Vector2i) cast(Vector2d) Vector2.ONE;
    cast(Vector3) cast(Vector3i) cast(Vector3d) Vector3.ONE;
    cast(Vector4) cast(Vector4i) cast(Vector4d) Vector4.ONE;
}

@("Vector as static arrays")
unittest 
{
    Vector4i A, B;
    assert(A == Vector4i.ZERO);

    A.array[2] = 42;
    *A.ptr = 1 + A[2];
    assert(A == Vector4i(43, 0, 42, 0));

    B.array = [7, 8, 9, 10];
    assert(B.array == [7, 8, 9, 10]);
}

@("Vector as slices")
unittest
{
    Vector3i V;
    int[] sl = V[];
    sl = V[1..$];
    sl[] = 2;
    sl = V[1..2];
    sl[0] = 4;
    assert(V == Vector3i(0, 4, 2));
}

@("Vector arithmetics")
unittest 
{
    Vector2i A = Vector2i(-42, 87);
    A = ((((A + 8) - 8) * 8) / 8) % 100;
    A += 2;
    A -= 2;
    A *= 4;
    A /= 4;
    A %= 100;
    assert(A == Vector2i(-42, 87));

    Vector3i B = Vector3i(1, 2, 3);
    B = ((((B + 1) - 1) * 1) / 1) % 100;;
    B += 1;
    B -= 1;
    B *= 1;
    B /= 1;
    B %= 100;
    assert(B == Vector3i(1, 2, 3));

    Vector4i C = Vector4i(1, 2, 3, 4);
    C = ((((C + 1) - 1) * 1) / 1) % 100;
    C += 1;
    C -= 1;
    C *= 1;
    C /= 1;
    C %= 100;
    assert(C == Vector4i(1, 2, 3, 4));

    Vector4i x = [4, 5, 6, 7];
    assert(x == x);
    --x[0];
    assert(x[0] == 3);
    ++x[0];
    assert(x[0] == 4);
    x[1] &= 1;
    x[2] = 77 + x[2];
    x[3] += 3;
    assert(x == Vector4i(4, 1, 83, 10));
}

@("Vector right operands")
unittest
{
    Vector2i A = Vector2i(1, 2);
    A = -(1 - (1 + A));
    A = 2 * A;
    A = 4 / A;
    assert(A == Vector2i(2, 1));

    Vector3i B = Vector3i(9, 6, 3);
    B = -(1 - (1 + B));
    B = 2 * B;
    B = 36 / B;
    assert(B == Vector3i(2, 3, 6));

    Vector4i C = Vector4i(1, 2, 3, 4);
    C = -(1 - (1 + C));
    C = 2 * C;
    C = 24 / C;
    assert(C == Vector4i(12, 6, 4, 3));
}

@("Vector3 octahedron encoding")
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





