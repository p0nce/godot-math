module transform3d;

import godotmath;

pure nothrow @nogc @safe:

@("Transform3D rotated")
unittest
{
    // rotated applies rotation in global space, origin gets rotated
    Transform3D T = Transform3D.IDENTITY;
    T.origin = Vector3(2, 9, 0);
    Transform3D R = T.rotated(Vector3(0, 0, 1), cast(float)(GM_PI / 2));
    assert(R.origin.is_equal_approx(Vector3(-9, 2, 0)));

    // rotated_local applies rotation in local space, origin unchanged
    T = Transform3D.IDENTITY;
    T.origin = Vector3(2, 9, 0);
    R = T.rotated_local(Vector3(0, 0, 1), cast(float)(GM_PI / 2));
    // Origin should be unchanged for local rotation
    assert(R.origin.is_equal_approx(Vector3(2, 9, 0)));
}

@("Transform3D inverse")
unittest
{
    Transform3D T = Transform3D.IDENTITY.rotated(Vector3(0, 1, 0), cast(float)(GM_PI / 5));
    Transform3D invT = T.inverse();
    assert( (T * invT).is_equal_approx(Transform3D.IDENTITY) );
    assert( (invT * T).is_equal_approx(Transform3D.IDENTITY) );
}

@("Transform3D affine_inverse")
unittest
{
    Transform3D T = Transform3D.IDENTITY.rotated(Vector3(0, 1, 2).normalized, 0.7);
    T.origin = Vector3(4, 2, 5);
    Transform3D invT = T.affine_inverse();
    assert( (T * invT).is_equal_approx(Transform3D.IDENTITY) );
    assert( (invT * T).is_equal_approx(Transform3D.IDENTITY) );
}

/+
@("Transform3D interpolate_with")
unittest
{
    import godotmath : GM_PI;
    // Use rotated transforms which have proper rotation bases
    auto T1 = Transform3D.IDENTITY.rotated(Vector3(0, 1, 0), cast(float)(GM_PI / 6));
    T1.origin = Vector3(0, 0, 0);
    auto T2 = Transform3D.IDENTITY.rotated(Vector3(0, 1, 0), cast(float)(GM_PI / 6));
    T2.origin = Vector3(10, 20, 30);
    
    auto mid = T1.interpolate_with(T2, 0.5);
    
    // Origin should be halfway
    assert(mid.origin.is_equal_approx(Vector3(5, 10, 15)));
}

@("Transform3D looking_at")
unittest
{
    auto T = Transform3D.IDENTITY;
    T.origin = Vector3(0, 0, 5);
    
    // Look at origin from (0,0,5)
    auto L = T.looking_at(Vector3(0, 0, 0), Vector3(0, 1, 0), true);
    
    // After looking_at, the Z basis column should point from origin toward target
    // Forward direction should be (0,0,-1) normalized
    auto forward = L.basis.xform(Vector3(0, 0, 1));
    assert(forward.is_equal_approx(Vector3(0, 0, -1)));
}
+/

@("Transform3D orthonormalized")
unittest
{
    Transform3D T = Transform3D.IDENTITY.rotated(Vector3(0, 1, 0), cast(float)(GM_PI / 6));    
    T = T.scaled(Vector3(2, 2, 2));    
    T.origin = Vector3(5, 5, 5);    
    Transform3D O = T.orthonormalized();
    // After orthonormalization, determinant should be 1
    assert(gm_is_equal_approx(O.basis.determinant(), 1.0f));
    // Origin preserved
    assert(O.origin.is_equal_approx(Vector3(5, 5, 5)));
}


@("Transform3D xform idempotence")
unittest
{
    Transform3Dd T = Transform3Dd.IDENTITY.rotated(Vector3d(1, 1, 1).normalized, cast(double)(GM_PI / 3));
    T.origin = Vector3d(7, 8, 9);

    Transform3Dd inv = T.affine_inverse();
    Vector3d v = Vector3d(1, 2, 3);
    
    // Transform then inverse should return original
    Vector3d transformed = T * v;
    Vector3d back = inv * transformed;
    
    assert(back.is_equal_approx(v));
}
