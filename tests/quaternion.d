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