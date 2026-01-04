module rect2;

import godotmath;

pure nothrow @nogc @safe:

@("Rect2 merge vs merge_non_empty")
unittest
{
    // union with empty box is NOT identity in Godot
    assert(Rect2i(0, 0, 1, 1).merge(Rect2i(10, 4, 0, 2)) != Rect2i(0, 0, 1, 1));

    // .merge_non_empty perform like dplug:math wrt empty union.
    assert(Rect2i(0, 0, 1, 1).merge_non_empty(Rect2i(10, 4, 0, 2)) == Rect2i(0, 0, 1, 1));
    assert(Rect2i(10, 4, 0, 2).merge_non_empty(Rect2i(0, 0, 1, 1)) == Rect2i(0, 0, 1, 1));
}

@("Rect2 intersection with empty")
unittest
{
    // intersection with empty box is empty
    assert(Rect2i(10, 4, 15, 7).intersection(Rect2i(10, 4, 0, 2)).has_no_area);

    // This behaviour is similar in dplug:math and godot-math
    assert(Rect2i(8, 4, 2, 7).intersection(Rect2i(10, 4, 0, 2)).has_no_area);
}


@("Rect2 basics")
unittest
{
    // construction from position and size
    Rect2i a = Rect2i(1, 2, 2, 2);
    a = Rect2i(Point2i(1, 2), Size2i(2, 2)); // equivalent
    assert(a.width == 2);
    assert(a.height == 2);
    assert(a.has_area);
    assert(a.get_area == 4);

    // rectangle equality (will not cast)
    Rect2i b = Rect2i(Vector2i(1, 2), Vector2i(2, 2));
    assert(a == b);

    // casts between rectangles
    Rect2 bf = cast(Rect2)b;
    assert(bf == Rect2(1, 2, 2, 2));

    // Assigning .left .right .top .bottom
    a = Rect2i.init; // (0, 0, 0, 0)    
    a.right = 100;
    a.left = -10;
    a.bottom = 120;
    a.top = -20;    
    assert(a.position == Vector2i(-10, -20));
    assert(a.size == Vector2i(110, 140));
    
    assert(!Rect2(0, 0, 0, 0).has_area());
    assert(Rect2i(0, 0, 1, 1).has_area());

    Transform2D T;
    Rect2 C = Rect2(0, 0, 4, 5);
    C = T * C;
    C = C * T;
    assert(C.is_equal_approx(Rect2(0, 0, 4, 5)));
}

@("Rect2 abs")
unittest
{
    Rect2d c = Rect2d(25, 25, -100, -50);
    assert(c.abs() == Rect2d(-75, -25, 100, 50));
}

@("Rect2 encloses and has_point")
unittest
{
    Rect2i c = Rect2i(0, 0, 1, 1);
    assert(c.translated(Vector2i(3, 3)) == Rect2i(3, 3, 1, 1));
    assert(c.encloses(Rect2i(0, 0, 1, 1)));
    

    // C is (0, 0)-(1, 1)

    // 0,0 --------- 1,0
    //  |             |
    //  |             |
    //  |             |
    //  |             |
    //  |             |
    // 0,1 --------- 1,1
    // So we consider it contains point 0,0 but not the 
    // others (1,0) (0,1) (1,1)
    assert(c.has_point(Point2i(0, 0)));
    assert(!c.has_point(Point2i(1, 0)));
    assert(!c.has_point(Point2i(0, 1)));
    assert(!c.has_point(Point2i(1, 1)));

    // WARNING:
    // dplug:math and godot-math disagree on whether rectangle 
    // of size 0,0 and with position (1,1) is included in c.
    assert(c.encloses(Rect2i(1, 1, 0, 0)));

    assert(!c.encloses(Rect2i(1, 1, 1, 1))); // has only point (1,1)
    assert(!c.has_point(Vector2i(1, 1)));

    // a rectangle .encloses itself
    Rect2i b = Rect2i(1, 2, 2, 2);
    assert(b.encloses(b));

    //Rect2i d = c.expand(Rect2i(3, 3));
    //assert(d.contains(Rect2i(2, 2)));
}

@("Rect2 expand")
unittest
{
   Rect2 A = Rect2(0, 0, 5, 2);
   A = A.expand(Vector2(10, 0));
   assert(A == Rect2(0, 0, 10, 2));
   A = A.expand(Vector2(-5, 5));
   assert(A == Rect2(-5, 0, 15, 5));
}

@("Rect2 intersection")
unittest
{
    {
        Rect2 rect1 = Rect2(0, 0, 5, 10);
        Rect2 rect2 = Rect2(2, 0, 8, 4);
        Rect2 a = rect1.intersection(rect2);
        assert(a.is_equal_approx(Rect2(2, 0, 3, 4)));
    }
    {
        Rect2i rect1 = Rect2i(0, 0, 5, 10);
        Rect2i rect2 = Rect2i(2, 0, 8, 4);
        Rect2i a = rect1.intersection(rect2);
        assert(a == Rect2i(2, 0, 3, 4));
    }
}

@("Rect2 merge")
unittest
{
    Rect2i A = Rect2i(0, 0, 4, 4);
    Rect2i B = Rect2i(2, 2, 4, 4);
    Rect2i C = Rect2i(0, 0, 2, 2);
    assert(A == A.merge(A));
    assert(!A.encloses(B));
    assert(A.encloses(C));
    assert(Rect2i(260, 100, 100, 100).intersection(Rect2i(100, 100, 100, 100)).has_no_area());
}
    
@("Rect2 scale_by_factor")
unittest
{
    assert(Rect2i(10, 10, 20, 20).scaled(1.5f) == Rect2i(15, 15, 30, 30));
    assert(Rect2i(10, 10, 20, 20).scaled(1.5f, 2.0f) == Rect2i(15, 20, 30, 40));
}