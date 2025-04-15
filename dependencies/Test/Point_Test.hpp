#pragma once
#include "testlib.hpp"
#include <vector>
#include <tuple>
#include "Point.hpp"

using namespace Point;

std::vector<std::tuple<BasePoint *, BasePoint *>> Geodesic_Geocentric_Points = {
    {new GeodesicPoint(0, 0, 0), new GeocentricPoint(6378137.0, 0, 0)},
    {new GeodesicPoint(7.523, 11.2423, 5323), new GeocentricPoint(6207868.059, 819816.659, 1236342.758)},
    {new GeodesicPoint(12.123, 15.234, 6534), new GeocentricPoint(6024328.161, 1294033.734, 1666813.083)},
    {new GeodesicPoint(15.456, 22.2341, 7522.432), new GeocentricPoint(5699822.474, 1575987.880, 2401274.180)},
    {new GeodesicPoint(17.435, 26.123, 7742.312), new GeocentricPoint(5473701.661, 1719026.896, 2794712.364)},
    {new GeodesicPoint(21.5323, 29.24234, 8234.1342), new GeocentricPoint(5187740.239, 2046884.307, 3101388.660)},
    {new GeodesicPoint(23.123, 31.2341, 8432.234), new GeocentricPoint(5026688.027, 2146449.502, 3292486.351)},
    {new GeodesicPoint(30.4234, 35.123, 9132.32), new GeocentricPoint(4509914.754, 2648427.744, 3654290.639)},
    {new GeodesicPoint(45.3213, 37.142, 9231.32), new GeocentricPoint(3584454.904, 3624883.443, 3835540.905)},
    {new GeodesicPoint(55.5345, 40.534, 9423.234), new GeocentricPoint(2751195.908, 4008183.126, 4129354.238)},
    {new GeodesicPoint(69.4213, 43.312, 9490.312), new GeocentricPoint(1636255.785, 4358106.590, 4359297.443)},
    {new GeodesicPoint(89.234, 60.2131, 1123.42), new GeocentricPoint(42473.812, 3176795.248, 5513285.034)},
    {new GeodesicPoint(99.234, 80.131, 1232.543), new GeocentricPoint(-176027.088, 1082752.577, 6263280.835)},
    {new GeodesicPoint(123.42, 85.123, 1323.423), new GeocentricPoint(-299715.452, 454197.276, 6334902.093)},
    {new GeodesicPoint(156.234, 87.213, 1523.321), new GeocentricPoint(-284846.501, 125430.416, 6350704.454)},
    {new GeodesicPoint(165.234, 84.213, 1823.321), new GeocentricPoint(-624120.708, 164503.651, 6325953.209)},
    {new GeodesicPoint(175.123, 88.234, 1923.432), new GeocentricPoint(-196564.734, 16772.059, 6355635.193)},
    {new GeodesicPoint(178.2423, 89.234, 2123.123), new GeocentricPoint(-85543.103, 2625.085, 6358303.337)},
    {new GeodesicPoint(-123, -65.42, 6450.234), new GeocentricPoint(-1450443.982, -2233487.871, -5783210.293)},
    {new GeodesicPoint(-113.123, -68.13, 7123.32), new GeocentricPoint(-936761.859, -2193766.784, -5903120.122)},
    {new GeodesicPoint(-45.4234, -32.324, 8234.312), new GeocentricPoint(3791383.377, -3847836.128, -3395249.260)},
};

std::vector<std::tuple<BasePoint *, BasePoint *>> Geocentric_Geodesic_Points = {
    {new GeocentricPoint(50113, 2239432, 2454334), new GeodesicPoint(88.71807482, 47.98109703, -3043562.462)},
    {new GeocentricPoint(10133.123, 323432.876, 446.2344), new GeodesicPoint(88.20551535, 0.09102128, -6054545.074)},
    {new GeocentricPoint(20143.123, 423432.765, 5465.434), new GeodesicPoint(87.27643531, 0.82137939, -5954186.215)},
    {new GeocentricPoint(30153.42, 223432.5435, 6234.234), new GeodesicPoint(82.31406398, 1.95342533, -6152572.663)},
    {new GeocentricPoint(40163.123, 123432.423, 198.7234), new GeodesicPoint(71.97586227, 0.13071609, -6248334.464)},
    {new GeocentricPoint(50173.342, 323432.765, 2876.34), new GeodesicPoint(81.18212888, 0.57903443, -6050821.201)},
    {new GeocentricPoint(60183.534, 223432.86, 5675.634), new GeodesicPoint(74.92468969, 1.72263398, -6146655.254)},
    {new GeocentricPoint(70193.645, 323432.678, 742.34), new GeodesicPoint(77.75515537, 0.14754803, -6047174.045)},
    {new GeocentricPoint(8213.345, 13432.96, 5484), new GeodesicPoint(58.55703938, 4.45598957, -6391404.079)},
    {new GeocentricPoint(93123.756, 523432.976, 422.34), new GeodesicPoint(79.91208023, 0.04948986, -5846484.568)},
    {new GeocentricPoint(14123.756, 22343.9872, 543.344), new GeodesicPoint(57.70279648, 0.44946293, -6404513.844)},
    {new GeocentricPoint(25123.46, 833432.9790, 3451.23634), new GeodesicPoint(88.27336752, 0.24995127, -5544317.911)},
    {new GeocentricPoint(310123.75, 8423432.567, 842.234), new GeodesicPoint(87.89150563, 0.00575410, 2051002.557)},
    {new GeocentricPoint(403123.345, 543432.456, 921453.61234), new GeodesicPoint(53.43162581, 54.72119680, -5220872.634)},
    {new GeocentricPoint(501223.423, 4273432.456, 94612.234), new GeodesicPoint(83.31044149, 1.27228674, -2074360.629)},
    {new GeocentricPoint(601243., 72343.632, 4466.234), new GeodesicPoint(6.8610495, 0.45460750, -5772539.604)},
    {new GeocentricPoint(7012.33, 52743.756432, 14661.34), new GeodesicPoint(82.42688667, 37.33888708, -6319082.360)},
    {new GeocentricPoint(801.213, 628654.456, 2622), new GeodesicPoint(89.92697722, 0.25638118, -5749476.167)},
    {new GeocentricPoint(90122.3, 72643.5742, 3643.34), new GeodesicPoint(38.87074136, 2.85292055, -6262291.684)},
    {new GeocentricPoint(1012.43, 9238.7432, 475.64), new GeodesicPoint(83.74618179, 0.52338901, -6387418.992)},
    {new GeocentricPoint(25012.33, 253632, 546.234), new GeodesicPoint(84.36789475, 0.14751190, -6123273.965)},
};

std::vector<std::tuple<BasePoint *, BasePoint *>> Geodesic_UTM_Points = {
    {new GeodesicPoint(0, 0, 0), new UtmPoint(5223032.968, 19995929.886, 7, HEMISPHERE::NORTH, 0.000)},
    {new GeodesicPoint(7.523, 11.2423, 5323), new UtmPoint(4107917.225, 18543880.841, 7, HEMISPHERE::NORTH, 5323.000)},
    {new GeodesicPoint(12.123, 15.234, 6534), new UtmPoint(12353745.070, 6791533.886, 19, HEMISPHERE::NORTH, 6534.000)},
    {new GeodesicPoint(15.456, 22.2341, 7522.432), new UtmPoint(10648412.118, 8539808.118, 19, HEMISPHERE::NORTH, 7522.432)},
    {new GeodesicPoint(17.435, 26.123, 7742.312), new UtmPoint(5938049.748, 4170181.517, 25, HEMISPHERE::NORTH, 7742.312)},
    {new GeodesicPoint(21.5323, 29.24234, 8234.1342), new UtmPoint(6170377.930, 4880068.000, 25, HEMISPHERE::NORTH, 8234.1342)},
    {new GeodesicPoint(23.123, 31.2341, 8432.234), new UtmPoint(-442771.350, 3497944.465, 36, HEMISPHERE::NORTH, 8432.234)},
    {new GeodesicPoint(30.4234, 35.123, 9132.32), new UtmPoint(265206.394, 3889722.033, 36, HEMISPHERE::NORTH, 9132.32)},
    {new GeodesicPoint(45.3213, 37.142, 9231.32), new UtmPoint(-6369973.822, 9079733.433, 52, HEMISPHERE::NORTH, 9231.32)},
    {new GeodesicPoint(55.5345, 40.534, 9423.234), new UtmPoint(-5406858.554, 7950660.791, 52, HEMISPHERE::NORTH, 9423.234)},
    {new GeodesicPoint(69.4213, 43.312, 9490.312), new UtmPoint(-1968937.751, 24727849.072, 7, HEMISPHERE::SOUTH, 9490.312)},
    {new GeodesicPoint(89.234, 60.2131, 1123.42), new UtmPoint(-2070895.818, 22241098.166, 7, HEMISPHERE::SOUTH, 1123.42)},
    {new GeodesicPoint(99.234, 80.131, 1232.543), new UtmPoint(495522.572, 21099713.404, 17, HEMISPHERE::SOUTH, 1232.543)},
    {new GeodesicPoint(123.42, 85.123, 1323.423), new UtmPoint(275070.480, 204939558.495, 17, HEMISPHERE::SOUTH, 1323.423)},
    {new GeodesicPoint(156.234, 87.213, 1523.321), new UtmPoint(579279.807, 20298868, 29, HEMISPHERE::SOUTH, 1523.321)},
    {new GeodesicPoint(165.234, 84.213, 1823.321), new UtmPoint(564802.165, 20640810.405, 29, HEMISPHERE::SOUTH, 1823.321)},
    {new GeodesicPoint(175.123, 88.234, 1923.432), new UtmPoint(650774.157, 20125052.175, 38, HEMISPHERE::SOUTH, 1923.432)},
    {new GeodesicPoint(178.2423, 89.234, 2123.123), new UtmPoint(562300.687, 20056557.543, 38, HEMISPHERE::SOUTH, 2123.123)},
    {new GeodesicPoint(-123, -65.42, 6450.234), new UtmPoint(3313218.003, -303377.470, 54, HEMISPHERE::SOUTH, 6450.234)},
    {new GeodesicPoint(-113.123, -68.13, 7123.32), new UtmPoint(2897219.244, -697243.517, 54, HEMISPHERE::SOUTH, 7123.32)},
    {new GeodesicPoint(-45.4234, -32.324, 8234.312), new UtmPoint(4633740.851, -5496157.147, 60, HEMISPHERE::SOUTH, 8234.312)},
    {new GeodesicPoint(-15.4234, -72.324, 4165.312), new UtmPoint(918333.460, -1927443.582, 60, HEMISPHERE::SOUTH, 4165.312)}};

std::vector<std::tuple<BasePoint *, BasePoint *>> UTM_Geodesic_Points = {{}};

/// @test
void GeodesictoGeocentric_ConversionTest()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;

    int num_vector = Geodesic_Geocentric_Points.size();
    for (size_t i = 0; i < num_vector; i++)
    {
        BasePoint *point = std::get<0>(Geodesic_Geocentric_Points[i]);

        BasePoint *point1 = std::get<1>(Geodesic_Geocentric_Points[i]);
        BasePoint *point2 = point->convert(POINT_TYPE::GEOCENTRIC);

        results.push_back(point2->compareLocation(point1, 1e-3)); // the cross checked answers are rounded despite this 0.001  error per answer or less
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

    logResult("GeodesicPoint -> convert(POINT_TYPE::GEOCENTRIC)", success, error_count, failed_cases,"GeodesictoGeocentric_ConversionTest()");

    test_counter++;
}

/// @test
void GeocentrictoGeodesic_ConversionTest()
{
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;

    int num_vector = Geocentric_Geodesic_Points.size();
    for (size_t i = 0; i < num_vector; i++)
    {
        BasePoint *point = std::get<0>(Geocentric_Geodesic_Points[i]);

        BasePoint *point1 = std::get<1>(Geocentric_Geodesic_Points[i]);
        BasePoint *point2 = point->convert(POINT_TYPE::GEODESIC);

        results.push_back(point2->compareLocation(point1, 1e-3)); // the cross checked answers are rounded despite this 0.001  error per answer or less
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

    logResult("GeocentricPoint->convert(POINT_TYPE::GEODESIC)", success, error_count, failed_cases,"GeocentrictoGeodesic_ConversionTest()");

    test_counter++;
}


///@test Tests UTM TO GEODESIC and INVERSE CONVERSIONS TEST for both HEMISPHERIC AND ZONE_LETTER types;
void GeodesicToUTMStandard(POINT_TYPE type){
    static int test_counter = 0;
    static int error_count = 0;
    bool success = true;
    std::vector<int> failed_cases;
    std::vector<bool> results;


    int num_vector = Geodesic_Geocentric_Points.size();
    for (size_t i = 0; i < num_vector; i++)
    {
        BasePoint *point = std::get<0>(Geodesic_Geocentric_Points[i]);
        BasePoint* point1 = point->convert(type);
        BasePoint* point2 = point1->convert(POINT_TYPE::GEODESIC);
        results.push_back(point->compareLocation(point2));
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
    logResult("GeodesicPoint->convert(POINT_TYPE::UTM)->convert(POINT_TYPE::GEODESIC)", success, error_count, failed_cases,"GeodesicToUTMStandard()");

    test_counter++;
}



