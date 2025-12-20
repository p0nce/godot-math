/*
    Copyright (c) 2014-2025 Godot Engine contributors.
    Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur. 
    Copyright (c) 2025 Guillaume Piolat
*/
module quaternion;

import godotmath;

pure nothrow @nogc @safe:

@("Quaternion basics")
unittest
{
    Quaternion A = Quaternion.IDENTITY;
    Quaternion B = Quaternion(Vector3(1, 1, 1).normalized, 2 * GM_PI / 3);
    B = B * B * B;
    
    // B is also the identity, but with anti-rotation 360°
    assert(A.is_equal_approx(-B)); 
    A = (A + B) - (B / 4 * 4);
    assert(A.is_equal_approx(-B));
}

@("Quaternion * Vector3")
unittest
{
    Vector3 axis = [1.0, 1, 1];
    Quaternion Q = Quaternion(axis.normalized, 2 * GM_PI / 3);    
    Vector3 V = Vector3(1, 2, 3);
    V = Q * V;    
    assert(V.is_equal_approx(Vector3(3, 1, 2)));
}

@("Vector3 * Quaternion")
unittest
{
    Vector3 axis = [1.0, 1, 1];
    Quaternion Q = Quaternion(axis.normalized, 2 * GM_PI / 3);
    Vector3 V = Vector3(3, 1, 2);
    V = V * Q;
    assert(V.is_equal_approx(Vector3(1, 2, 3)));
}

@("Quaternion slerp")
unittest
{
    Vector3 axis = Vector3(-2, 5, 8).normalized;
    // Test normal SLERP between two different quaternions
    Quaternion q1 = Quaternion(axis, GM_PI / 4); // 45° around Z
    Quaternion q2 = Quaternion(axis, GM_PI / 2); // 90° around Z

    // At weight 0, should return q1
    Quaternion result0 = q1.slerp(q2, 0.0f);
    assert(result0.is_equal_approx(q1));

    // At weight 1, should return q2
    Quaternion result1 = q1.slerp(q2, 1.0f);
    assert(result1.is_equal_approx(q2));

    // At weight 0.5, should be halfway between
    Quaternion result05 = q1.slerp(q2, 0.5f);
    // The angle should be halfway: (45° + 90°) / 2 = 67.5°
    float expected_angle = (GM_PI / 4 + GM_PI / 2) / 2;
    Quaternion expected = Quaternion(axis, expected_angle);
    assert(result05.is_equal_approx(expected));
}

@("Quaternion log/exp")
unittest
{
    Quaternion Q = Quaternion(Vector3(1.1, 2.4, 3.5).normalized, GM_PI);
    Quaternion LQ = Q.log();
    Quaternion ELG = LQ.exp();
    assert(Q.is_equal_approx(ELG));
    assert(Q.is_equal_approx(Q.exp.log.exp.log));
}

@("Quaternion <=> Euler angles")
unittest
{
    // Test round-trip conversion: euler -> quaternion -> euler
    Vector3 original_euler = Vector3(GM_PI / 6, GM_PI / 4, GM_PI / 3); // 30°, 45°, 60°
    Quaternion q = Quaternion.from_euler(original_euler);
    Vector3 recovered_euler = q.get_euler(GM_EULER_ORDER_YXZ);
    assert(recovered_euler.is_equal_approx(original_euler));

    Vector3 identity_euler = Quaternion.IDENTITY.get_euler();
    assert(identity_euler.is_equal_approx(Vector3.ZERO));
}