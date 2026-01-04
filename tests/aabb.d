module aabb;

import godotmath;

pure nothrow @nogc @safe:

@("AABB abs")
unittest
{
    AABB box = AABB(Vector3(5, 0, 5), Vector3(-20, -10, -5));
    AABB absolute = box.abs();
    assert(absolute.position == Vector3(-15.0, -10.0, 0.0));
    assert(absolute.size == Vector3(20.0, 10.0, 5.0));
}
  
@("AABB encloses")
unittest
{
    AABB a = AABB(Vector3(0, 0, 0), Vector3(4, 4, 4));
    AABB b = AABB(Vector3(1, 1, 1), Vector3(3, 3, 3));
    AABB c = AABB(Vector3(2, 2, 2), Vector3(8, 8, 8));
    assert(a.encloses(a)); 
    assert(a.encloses(b)); 
    assert(!a.encloses(c));
}

@("AABB expand")
unittest
{
    AABB box = AABB(Vector3(0, 0, 0), Vector3(5, 2, 5));

    box = box.expand(Vector3(10, 0, 0));
    assert(box.position == Vector3(0.0, 0.0, 0.0));
    assert(box.size  == Vector3(10.0, 2.0, 5.0));

    box = box.expand(Vector3(-5, 0, 5));
    assert(box.position == Vector3(-5.0, 0.0, 0.0));
    assert(box.size  == Vector3(15.0, 2.0, 5.0));
}

@("AABB get_longest_axis")
unittest
{
    AABB box = AABB(Vector3(0, 0, 0), Vector3(2, 4, 8));
    assert(box.get_longest_axis() == Vector3(0, 0, 1));
    assert(box.get_longest_axis_index() == 2);
    assert(box.get_longest_axis_size() == 8.0);
}

@("AABB get_shortest_axis")
unittest
{
    AABB box = AABB(Vector3(0, 0, 0), Vector3(2, 4, 8));
    assert(box.get_shortest_axis() == Vector3(1.0, 0.0, 0.0));
    assert(box.get_shortest_axis_index() == 0);
    assert(box.get_shortest_axis_size() == 2.0);
}

@("AABB grow")
unittest
{
    AABB box = AABB(Vector3(4, 4, 4), Vector3(8, 8, 8));
    AABB grown = box.grow(4);
    assert(grown.position == Vector3(0.0, 0.0, 0.0));
    assert(grown.size == Vector3(16.0, 16.0, 16.0));
}

@("AABB intersection")
unittest
{
    AABB box1 = AABB(Vector3(0, 0, 0), Vector3(5, 2, 8));
    AABB box2 = AABB(Vector3(2, 0, 2), Vector3(8, 4, 4));
    AABB intersection = box1.intersection(box2);
    assert(intersection.position == Vector3(2.0, 0.0, 2.0));
    assert(intersection.size == Vector3(3.0, 2.0, 4.0));
}