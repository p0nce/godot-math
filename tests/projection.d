module projection;

import godotmath;

pure nothrow @nogc @safe:

@("Projection inverse")
unittest
{
    // Using https://math.stackexchange.com/questions/2574905/invertible-4x4-matrix

    // Test identity multiplication
    Projection A = Projection.IDENTITY;
    Projection B = Projection(Vector4(5, 2, 6, 2),
                              Vector4(6, 2, 6, 3),
                              Vector4(6, 2, 2, 6),
                              Vector4(8, 8, 8, 7));
    
    // test identity multiplication
    assert((B * A) == B);
    assert((A * B) == B);

    assert(B.determinant() == -8);

    Projection C = B.inverse();

    Projection IB = Projection(Vector4(-17, 17, -4, 1),
                               Vector4(-9, 8.75, -2.25, 0.75),
                               Vector4(12, -11.75, 2.75, -0.75),
                               Vector4(16, -16, 4, -1));

    assert(C.is_equal_approx(IB));
    assert( (C*B).is_equal_approx(Projection.IDENTITY));
    assert( (B*C).is_equal_approx(Projection.IDENTITY));
}