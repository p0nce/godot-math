module basis;

import godotmath;

pure nothrow @nogc @safe:

@("Basis * Vector3")
unittest
{
    Basis B = Basis(Vector3(0, 2, 1),
                    Vector3(2, 0, 1),
                    Vector3(4, 0, -2));
    Vector3 V = Vector3(1, 2, 3);
    V = B * V;    
    assert(V.is_equal_approx(Vector3(16, 2, -3)));
}

@("Vector3 * Basis")
unittest
{
    Basis B = Basis(Vector3(0, 0, 1),
                    Vector3(0, 1, 0),
                    Vector3(-1, 0, 0));
    assert(B.is_rotation());
    assert(B.is_conformal());
    Vector3 V = Vector3(4, 2, -6);
    V = V * B;
    assert(V.is_equal_approx(Vector3(-6, 2, -4)));
    V = B * V;
    assert(V.is_equal_approx(Vector3(4, 2, -6)));
}

@("Basis <=> Euler angles")
unittest
{
    Vector3d angles = [0.1, 0.1, 0.1];
    for (EulerOrder order = GM_EULER_ORDER_XYZ; order <= GM_EULER_ORDER_ZYX; ++order)
    {
        Basisd B = Basisd.from_euler(angles, order);
        Vector3d angles2 = B.get_euler(order);
        assert(angles.is_equal_approx(angles2));
    }

    Vector3 anglesf = [-0.1f, -0.1f, -0.1f];
    for (EulerOrder order = GM_EULER_ORDER_XYZ; order <= GM_EULER_ORDER_ZYX; ++order)
    {
        Basis B = Basis.from_euler(anglesf, order);
        Vector3 angles2 = B.get_euler(order);
        assert(anglesf.is_equal_approx(angles2));
    }
}

@("Basis scale")
unittest
{
    Basis my_basis = Basis(
        Vector3(2, 0, 0),
        Vector3(0, 4, 0),
        Vector3(0, 0, 8)
    );
    // Rotating the Basis in any way preserves its scale.
    my_basis = my_basis.rotated(Vector3.UP, GM_TAU / 2);
    my_basis = my_basis.rotated(Vector3.RIGHT, GM_TAU / 4);
    assert(my_basis.get_scale().is_equal_approx( Vector3(2.0, 4.0, 8.0)));
}

@("Basis scaled")
unittest
{
    Basis ma_base = Basis(
        Vector3(1, 1, 1),
        Vector3(2, 2, 2),
        Vector3(3, 3, 3)
    );
    ma_base = ma_base.scaled(Vector3(0, 2, -2));
    assert(ma_base.x.is_equal_approx(Vector3(0.0, 2.0, -2.0)));
    assert(ma_base.y.is_equal_approx(Vector3(0.0, 4.0, -4.0)));
    assert(ma_base.z.is_equal_approx(Vector3(0.0, 6.0, -6.0)));
}

@("Basis casts")
unittest
{
    Basis B = Basis.from_scale(Vector3(1, 2, 1));
    Basisd Bd = cast(Basisd)B;
    Basis C = cast(Basis)Bd;
}

@("Basis inverse rotation")
unittest
{
    Basis B = Basis.from_euler(Vector3(4, 5, 6));
    assert(B.is_rotation());
    assert(B.is_conformal());
    Basis invB = B.inverse();
    assert(invB.is_rotation());
    assert(invB.is_conformal());
    assert((B * invB).is_equal_approx(Basis.IDENTITY));
    assert((invB * B).is_equal_approx(Basis.IDENTITY));
}

@("Basis inverse conformal")
unittest
{
    Basis B = Basis.from_euler(Vector3(4, 5, 6)).scaled(Vector3(2, 2, 2));
    assert(!B.is_rotation());
    assert(B.is_conformal());
    Basis invB = B.inverse();
    assert(invB.is_conformal());
    assert((B * invB).is_equal_approx(Basis.IDENTITY));
    assert((invB * B).is_equal_approx(Basis.IDENTITY));
}

@("Basis inverse any")
unittest
{
    Basis B = Basis(Vector3(-0.9, -0.25, 0.15),
                    Vector3(0.6, 0.5,-0.1),
                    Vector3(0.1,-0.25,0.15));
    assert(!B.is_rotation());
    assert(!B.is_conformal());
    Basis invB = B.inverse();
    assert(!invB.is_rotation());
    assert(!invB.is_conformal());
    assert((B * invB).is_equal_approx(Basis.IDENTITY));
    assert((invB * B).is_equal_approx(Basis.IDENTITY));
}

@("Basis orthonormalized")
unittest
{
    Basis B = Basis(Vector3(-0.9, -0.25,  0.15),
                    Vector3(-0.6,  -0.5,  0.1),
                    Vector3( 0.1, -0.25,  0.15));
    B = B.orthonormalized();
    assert(B.is_conformal);
    assert(B.is_rotation);
}