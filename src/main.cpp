#include "Point.hpp"
#include "KarneyFF.hpp"
#include "Shapes.hpp"
#include "Test_Main.cpp"
#include <iostream>

int main()
{
    BasePoint* point = new GeodesicPoint(10,0,5000);
    point->repr();

    BasePoint* point1 = point->convert(POINT_TYPE::UTM_59S);
    point1->repr();
}
