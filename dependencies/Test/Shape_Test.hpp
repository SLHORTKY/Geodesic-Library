#include "testlib.hpp"
#include "Shapes.hpp"
#include <vector>

using namespace Shapes;

///@brief POLYGON AREA, PERIMETER, FILTERING TESTS

#pragma region POLYGON  // before running planimeter test reconfigure the compute
// reverse parameter to false otherwise the area calculation will result in negative number
static void Planimeter0()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    // Check fix for pole-encircling bug found 2011-03-16
    std::vector<BasePoint *> verticies = {new GeodesicPoint(0, 89), new GeodesicPoint(90, 89), new GeodesicPoint(180, 89), new GeodesicPoint(270, 89)};
    Polygon p1 = Polygon(verticies);
    std::vector<bool> results;

    results.push_back(assertEquals(p1.calculatePerimeter(), 631819.8745, 1e-4));
    results.push_back(assertEquals(p1.calculateArea(), 24952305678.0, 1));
    std::vector<BasePoint *> verticies1 = {new GeodesicPoint(0, -89), new GeodesicPoint(90, -89), new GeodesicPoint(180, -89), new GeodesicPoint(270, -89)};
    Polygon p2 = Polygon(verticies1);

    results.push_back(assertEquals(p2.calculatePerimeter(), 631819.8745, 1e-4));
    results.push_back(assertEquals(p2.calculateArea(), -24952305678.0, 1));

    std::vector<BasePoint *> verticies2 = {new GeodesicPoint(-1, 0), new GeodesicPoint(0, -1), new GeodesicPoint(1, 0), new GeodesicPoint(0, 1)};
    Polygon p3 = Polygon(verticies2);
    results.push_back(assertEquals(p3.calculatePerimeter(), 627598.2731, 1e-4));
    results.push_back(assertEquals(p3.calculateArea(), 24619419146.0, 1));

    std::vector<BasePoint *> verticies3 = {new GeodesicPoint(0, 90), new GeodesicPoint(0, 0), new GeodesicPoint(90, 0)};
    Polygon p4 = Polygon(verticies3);
    results.push_back(assertEquals(p4.calculatePerimeter(), 30022685.0, 1.0));
    results.push_back(assertEquals(p4.calculateArea(), 63758202715511.0, 1));

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon::calculateArea() , Polygon::calculatePerimeter()", success, error_count, failed_cases, "PLAINMETER0 TEST");

    test_counter++;
}

static void Planimeter5()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;

    // Check fix for Planimeter pole crossing bug found 2011-06-24
    std::vector<BasePoint *> verticies = {new GeodesicPoint(0.1, 89), new GeodesicPoint(90.1, 89), new GeodesicPoint(-179.9, 89)};
    Polygon p = Polygon(verticies);
    std::vector<bool> results = {assertEquals(p.calculatePerimeter(), 539297.0, 1),
                                 assertEquals(p.calculateArea(), 12476152838.5, 1)};

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon::calculateArea() , Polygon::calculatePerimeter()", success, error_count, failed_cases, "PLAINMETER5 TEST");

    test_counter++;
}

static void Planimeter6()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;

    std::vector<BasePoint *> verticies1 = {new GeodesicPoint(-0.00000000000001, 9), new GeodesicPoint(180, 9), new GeodesicPoint(0, 9)};
    Polygon p = Polygon(verticies1);
    results.push_back(assertEquals(p.calculatePerimeter(), 36026861.0, 1));
    results.push_back(assertEquals(p.calculateArea(), 0.0, 1));

    std::vector<BasePoint *> verticies2 = {new GeodesicPoint(0.00000000000001, 9), new GeodesicPoint(0, 9), new GeodesicPoint(180, 9)};
    Polygon p1 = Polygon(verticies2);
    results.push_back(assertEquals(p1.calculatePerimeter(), 36026861.0, 1));
    results.push_back(assertEquals(p1.calculateArea(), 0.0, 1));

    std::vector<BasePoint *> verticies3 = {new GeodesicPoint(0.00000000000001, 9), new GeodesicPoint(180, 9), new GeodesicPoint(0, 9)};
    Polygon p2 = Polygon(verticies3);
    results.push_back(assertEquals(p2.calculatePerimeter(), 36026861.0, 1));
    results.push_back(assertEquals(p2.calculateArea(), 0.0, 1));

    std::vector<BasePoint *> verticies4 = {new GeodesicPoint(-0.00000000000001, 9), new GeodesicPoint(0, 9), new GeodesicPoint(180, 9)};
    Polygon p3 = Polygon(verticies4);
    results.push_back(assertEquals(p3.calculatePerimeter(), 36026861.0, 1));
    results.push_back(assertEquals(p3.calculateArea(), 0.0, 1));

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon::calculateArea() , Polygon::calculatePerimeter()", success, error_count, failed_cases, "PLAINMETER6 TEST");

    test_counter++;
}

static void Planimeter12()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    // Area of arctic circle (not really -- adjunct to rhumb-area test)
    std::vector<BasePoint *> vertex = {new GeodesicPoint(0, 66.562222222), new GeodesicPoint(180, 66.562222222), new GeodesicPoint(360, 66.562222222)};
    Polygon p = Polygon(vertex);
    std::vector<bool> results = {assertEquals(p.calculatePerimeter(), 10465729.0, 1),
                                 assertEquals(p.calculateArea(), 0.0, 1)};

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon::calculateArea() , Polygon::calculatePerimeter()", success, error_count, failed_cases, "PLAINMETER 12");

    test_counter++;
}

static void Planimeter12r()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    // Reverse area of arctic circle
    std::vector<BasePoint *> verticies =
        {new GeodesicPoint(-0, 66.562222222), new GeodesicPoint(-180, 66.562222222), new GeodesicPoint(-360, 66.562222222)};
    Polygon p = Polygon(verticies);
    std::vector<bool> results = {assertEquals(p.calculatePerimeter(), 10465729.0, 1),
                                 assertEquals(p.calculateArea(), 0.0, 1)};

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon::calculateArea() , Polygon::calculatePerimeter()", success, error_count, failed_cases, "PLAINMETER12R TEST");

    test_counter++;
}

static void Planimeter13()
{
    // Check encircling pole twice

    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<BasePoint *> verticies =
        {new GeodesicPoint(-360, 89), new GeodesicPoint(-240, 89), new GeodesicPoint(-120, 89), new GeodesicPoint(0, 89), new GeodesicPoint(120, 89), new GeodesicPoint(240, 89)};
    Polygon p = Polygon(verticies);
    std::vector<bool> results = {assertEquals(p.calculatePerimeter(), 1160741.0, 1),
                                 assertEquals(p.calculateArea(), 32415230256.0, 1)};

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon::calculateArea() , Polygon::calculatePerimeter()", success, error_count, failed_cases, "PLAINMETER13 TEST");

    test_counter++;
}

static void POLYGON_TEST1()
{
    std::vector<BasePoint *> verticies = {
        new GeodesicPoint(69.09302, 62.68675),
        new GeodesicPoint(77.96997, 66.76809),
        new GeodesicPoint(90.53833, 65.81403),
        new GeodesicPoint(95.06470, 63.54121),
        new GeodesicPoint(91.41724, 59.41434),
        new GeodesicPoint(81.57349, 57.99937),
    };

    std::vector<BasePoint *> points_inside_Polygon_1 = {
        new GeodesicPoint(74.80591, 63.30528),
        new GeodesicPoint(81.04614, 65.83203),
        new GeodesicPoint(86.71903, 63.06740),
        new GeodesicPoint(81.92505, 59.09985),
        new GeodesicPoint(92.34009, 64.57911),
    };

    std::vector<BasePoint *> points_outside_Polygon_1 = {
        new GeodesicPoint(71.64185, 65.15766),
        new GeodesicPoint(71.24634, 60.16611),
    };

    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;

    Polygon p1 = Polygon(verticies);

    results.push_back(assertEquals(p1.calculateArea(), 821989461394.8, 1e-1));
    results.push_back(assertEquals(p1.calculatePerimeter(), 3478828.485131, 1e-1));

    for (size_t i = 0; i < points_inside_Polygon_1.size(); i++)
    {
        results.push_back(assertTrue(p1.doesContain(points_inside_Polygon_1[i])));
    }
    for (size_t i = 0; i < points_outside_Polygon_1.size(); i++)
    {
        results.push_back(assertFalse(p1.doesContain(points_outside_Polygon_1[i])));
    }

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon->calculateArea(), Polygon->calculatePerimeter(),Polygon->doesContain()", success, error_count, failed_cases, "Polygon 1 TESTS");

    test_counter++;
}

static void POLYGON_TEST2(){
    std::vector<BasePoint *> verticies = {
        new GeodesicPoint(85.29785, 32.76880),
        new GeodesicPoint(93.03223, 41.70573),
        new GeodesicPoint(104.63379, 48.980022),
        new GeodesicPoint(119.92676, 47.33882),
        new GeodesicPoint(112.36816, 45.52174),
        new GeodesicPoint(123.44238, 39.43619),
        new GeodesicPoint(122.03613, 28.99853),
        new GeodesicPoint(105.86426, 11.60919),
        new GeodesicPoint(90.92285, 21.69827),
    };

    std::vector<BasePoint *> points_inside_Polygon_2 = {
        new GeodesicPoint(97.77832, 38.75408),
        new GeodesicPoint(103.75488, 34.52466),
        new GeodesicPoint(115.00488, 37.23033),
        new GeodesicPoint(113.59863, 27.13737),
        new GeodesicPoint(103.93066, 13.49647),
        new GeodesicPoint(93.20801, 28.22697),

    };
    std::vector<BasePoint *> points_outside_Polygon_2 = {
        new GeodesicPoint(90.21973, 13.66734),
        new GeodesicPoint(69.65332, 29.45873),
        new GeodesicPoint(98.65723, 46.37725),
        new GeodesicPoint(127.30957, 42.87596),
        new GeodesicPoint(139.96582, 11.09217),
    };

    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;

    Polygon p2 = Polygon(verticies);

    results.push_back(assertEquals(p2.calculateArea(), 9212630176394.3, 1e-1));
    results.push_back(assertEquals(p2.calculatePerimeter(), 12326512.870055, 1e-1));

    for (size_t i = 0; i < points_inside_Polygon_2.size(); i++)
    {
        results.push_back(assertTrue(p2.doesContain(points_inside_Polygon_2[i])));
    }
    for (size_t i = 0; i < points_outside_Polygon_2.size(); i++)
    {
        results.push_back(assertFalse(p2.doesContain(points_outside_Polygon_2[i])));
    }

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Polygon->calculateArea(), Polygon->calculatePerimeter(),Polygon->doesContain()", success, error_count, failed_cases, "Polygon 2 TESTS");

    test_counter++;
}
#pragma endregion

///@brief CIRCLE AREA, PERIMETER, FILTERING TESTS
#pragma region CIRCLE

// 14 points all around the circle in Geodesic Coordinates
std::vector<BasePoint *> points_inside_circle1 = {
        new GeodesicPoint(31.77383, 51.68959),
        new GeodesicPoint(32.89444, 51.45743),
        new GeodesicPoint(34.38858, 51.80522),
        new GeodesicPoint(31.38931, 50.29987),
        new GeodesicPoint(34.21280, 50.39801),
        new GeodesicPoint(33.99307, 49.64273),
        new GeodesicPoint(32.79556, 50.23667),
        new GeodesicPoint(35.14664, 50.62159)};
std::vector<BasePoint *> points_outside_circle1 = {
    new GeodesicPoint(31.89468, 52.16382),
    new GeodesicPoint(30.67520, 51.20344),
    new GeodesicPoint(35.25650, 49.80609),
    new GeodesicPoint(35.91568, 51.21721)};
/// @brief a test to calculate area perimeter calculations and to check whether a test points lie within the circle
static void CIRCLE1_TEST()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;

    BasePoint *center_point_first_circle = new GeodesicPoint(33.21579, 50.81635);
    Circle c1 = Circle(center_point_first_circle, 149099.6163844825, 0.0);

    results.push_back(assertEquals(c1.calculateArea(), 69839790000.0, 1e-5));
    results.push_back(assertEquals(c1.calculatePerimeter(), 936820.5189730931, 1e-5));

    for (size_t i = 0; i < points_inside_circle1.size(); i++)
    {
        results.push_back(assertTrue(c1.doesContain(points_inside_circle1[i])));
    }
    for (size_t i = 0; i < points_outside_circle1.size(); i++)
    {
        results.push_back(assertFalse(c1.doesContain(points_outside_circle1[i])));
    }

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Circle->calculateArea(), Circle->calculatePerimeter(),Circle->doesContain()", success, error_count, failed_cases, "CIRCLE 1 TESTS");

    test_counter++;
}
// 14 points all around the circle in Geodesic Coordinates
std::vector<BasePoint *> points_inside_circle2 = {
    new GeodesicPoint(121.55411, 73.13610),
    new GeodesicPoint(133.50723, 73.05944),
    new GeodesicPoint(147.65762, 73.28841),
    new GeodesicPoint(122.96036, 69.48645),
    new GeodesicPoint(124.36661, 62.89271),
    new GeodesicPoint(147.48184, 62.93274),
    new GeodesicPoint(122.25723, 56.33786),
    new GeodesicPoint(137.02286, 55.25095),
    new GeodesicPoint(147.65762, 69.30085),
    new GeodesicPoint(111.35880, 55.25095),
    new GeodesicPoint(135.44083, 55.25095),
    new GeodesicPoint(159.69864, 55.25095),
    new GeodesicPoint(148.88809, 55.25095)};
std::vector<BasePoint *> points_outside_circle2 = {
    new GeodesicPoint(99.40567, 74.38886),
    new GeodesicPoint(85.43106, 67.82376),
    new GeodesicPoint(92.19864, 61.83282),
    new GeodesicPoint(104.50333, 56.04443),
    new GeodesicPoint(162.95059, 69.39385),
    new GeodesicPoint(176.13419, 73.83553),
    new GeodesicPoint(160.92911, 61.99833),

};
/// @brief a test to calculate area perimeter calculations and to check whether a test points lie within the circle
static void CIRCLE2_TEST()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;

    BasePoint *center_point_first_circle = new GeodesicPoint(135.44083, 66.10494);
    Circle c1 = Circle(center_point_first_circle, 1248078.601404626, 0.0);

    results.push_back(assertEquals(c1.calculateArea(), 4893659490000.0, 1e-5));
    results.push_back(assertEquals(c1.calculatePerimeter(), 7841909.130550792, 1e-5));

    for (size_t i = 0; i < points_inside_circle2.size(); i++)
    {
        results.push_back(assertTrue(c1.doesContain(points_inside_circle2[i])));
    }
    for (size_t i = 0; i < points_outside_circle2.size(); i++)
    {
        results.push_back(assertFalse(c1.doesContain(points_outside_circle2[i])));
    }

    for (size_t i = 0; i < results.size(); i++)
    {
        if (results[i] == false)
        {
            success = false;
            failed_cases.push_back(i);
            error_count++;
        }
    }

    logResult("Circle->calculateArea(), Circle->calculatePerimeter(),Circle->doesContain()", success, error_count, failed_cases, "CIRCLE 2 TESTS");

    test_counter++;
}
#pragma endregion

///@brief ELLIPSE AREA, PERIMETER, FILTERING TESTS
#pragma region Ellipse

#pragma endregion