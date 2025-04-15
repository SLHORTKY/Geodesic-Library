#include "KarneyFF_Test.hpp"
#include "Sign_Test.hpp"
#include "Point_Test.hpp"
#include "Shape_Test.hpp"

void POINT_TESTS();
void KARNEY_FF_TESTS();
void GEOMATH_FUNCTION_TESTS();
void SHAPE_TESTS();

/// @brief MAIN TEST EXECUTION
void TEST_MAIN()
{
   POINT_TESTS();
}


void KARNEY_FF_TESTS()
{
    ///@brief KARNEY_FF FUNCTION TESTS
    DirectCheck();
    InverseCheck();
    ArcDirectCheck();
    GeodSolve0();
    GeodSolve1();
    GeodSolve2();
    GeodSolve4();
    GeodSolve5();
    GeodSolve6();
    GeodSolve9();
    GeodSolve10();
    GeodSolve11();
    GeodSolve12();
    GeodSolve14();
    GeodSolve15();
    GeodSolve17();
    GeodSolve26();
    GeodSolve28();
    GeodSolve29();
    GeodSolve33();
    GeodSolve55();
    GeodSolve59();
    GeodSolve61();
    GeodSolve65();
    GeodSolve69();
    GeodSolve71();
    GeodSolve73();
    GeodSolve76();
    GeodSolve78();
    GeodSolve80();
    GeodSolve84();
    GeodSolve92();
    GeodSolve94();
    GeodSolve96();
}
void GEOMATH_FUNCTION_TESTS()
{
    ///@brief GEOMATH FUNCTION TESTS
    test_sincosd();
    test_atan2d();
    test_sum();
    test_AngNormalize();
    test_AngDiff();
    test_equatorial_coincident();
    test_equatorial_NS();
    test_antipodal();
    test_antipodal_prolate();
    test_azimuth_0_180();
}
void POINT_TESTS(){
    GeodesictoGeocentric_ConversionTest();
    GeocentrictoGeodesic_ConversionTest();
    GeodesicToUTMStandard(POINT_TYPE::UTM_H);
    GeodesicToUTMStandard(POINT_TYPE::UTM_Z);
}
void SHAPE_TESTS(){
    POLYGON_TEST2();
}

