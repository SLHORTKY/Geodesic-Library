#include "KarneyFF_Test.hpp"
#include "testlib.hpp"

/// @test
void test_AngRound()
{
  // Test special cases for AngRound

  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;

  double eps = GeoMath::ulp(1.0);
  std::vector<bool> results = {
      assertTrue(equiv(GeoMath::AngRound(-eps / 32), -eps / 32)),
      assertTrue(equiv(GeoMath::AngRound(-eps / 64), -0.0)),
      assertTrue(equiv(GeoMath::AngRound(-0.0), -0.0)),
      assertTrue(equiv(GeoMath::AngRound(0.0), +0.0)),
      assertTrue(equiv(GeoMath::AngRound(eps / 64), +0.0)),
      assertTrue(equiv(GeoMath::AngRound(eps / 32), +eps / 32)),
      assertTrue(equiv(GeoMath::AngRound((1 - 2 * eps) / 64), (1 - 2 * eps) / 64)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps) / 64), 1.0 / 64)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 2) / 64), 1.0 / 64)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 4) / 64), 1.0 / 64)),
      assertTrue(equiv(GeoMath::AngRound(1.0 / 64), 1.0 / 64)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps / 2) / 64), 1.0 / 64)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps) / 64), 1.0 / 64)),
      assertTrue(equiv(GeoMath::AngRound((1 + 2 * eps) / 64), (1 + 2 * eps) / 64)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps) / 32), (1 - eps) / 32)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 2) / 32), 1.0 / 32)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 4) / 32), 1.0 / 32)),
      assertTrue(equiv(GeoMath::AngRound(1.0 / 32), 1.0 / 32)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps / 2) / 32), 1.0 / 32)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps) / 32), (1 + eps) / 32)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps) / 16), (1 - eps) / 16)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 2) / 16), (1 - eps / 2) / 16)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 4) / 16), 1.0 / 16)),
      assertTrue(equiv(GeoMath::AngRound(1.0 / 16), 1.0 / 16)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps / 4) / 16), 1.0 / 16)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps / 2) / 16), 1.0 / 16)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps) / 16), (1 + eps) / 16)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps) / 8), (1 - eps) / 8)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 2) / 8), (1 - eps / 2) / 8)),
      assertTrue(equiv(GeoMath::AngRound((1 - eps / 4) / 8), 1.0 / 8)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps / 2) / 8), 1.0 / 8)),
      assertTrue(equiv(GeoMath::AngRound((1 + eps) / 8), (1 + eps) / 8)),
      assertTrue(equiv(GeoMath::AngRound(1 - eps), 1 - eps)),
      assertTrue(equiv(GeoMath::AngRound(1 - eps / 2), 1 - eps / 2)),
      assertTrue(equiv(GeoMath::AngRound(1 - eps / 4), 1)),
      assertTrue(equiv(GeoMath::AngRound(1.0), 1)),
      assertTrue(equiv(GeoMath::AngRound(1 + eps / 4), 1)),
      assertTrue(equiv(GeoMath::AngRound(1 + eps / 2), 1)),
      assertTrue(equiv(GeoMath::AngRound(1 + eps), 1 + eps)),
      assertTrue(equiv(GeoMath::AngRound(90.0 - 64 * eps), 90 - 64 * eps)),
      assertTrue(equiv(GeoMath::AngRound(90.0 - 32 * eps), 90)),
      assertTrue(equiv(GeoMath::AngRound(90.0), 90)),
  };

  for (size_t i = 0; i < results.size(); i++)
  {
    if (results[i] == false)
    {
      success = false;
      failed_cases.push_back(i);
      error_count++;
    }
  }

  logResult("Geodesic::GeoMath::AngRound(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
/// @test
void test_sincosd()
{
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;

  std::vector<bool> results;
  // Test special cases for sincosd
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  GeoMath::Pair<double> p = GeoMath::Pair<double>();
  GeoMath::sincosd(p, -inf);
  results.push_back(assertTrue(equiv(p.first, nan) && equiv(p.second, nan)));
  GeoMath::sincosd(p, -810.0);
  results.push_back(assertTrue(equiv(p.first, -1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, -720.0);
  results.push_back(assertTrue(equiv(p.first, -0.0) && equiv(p.second, +1.0)));
  GeoMath::sincosd(p, -630.0);
  results.push_back(assertTrue(equiv(p.first, +1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, -540.0);
  results.push_back(assertTrue(equiv(p.first, -0.0) && equiv(p.second, -1.0)));
  GeoMath::sincosd(p, -450.0);
  results.push_back(assertTrue(equiv(p.first, -1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, -360.0);
  results.push_back(assertTrue(equiv(p.first, -0.0) && equiv(p.second, +1.0)));
  GeoMath::sincosd(p, -270.0);
  results.push_back(assertTrue(equiv(p.first, +1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, -180.0);
  results.push_back(assertTrue(equiv(p.first, -0.0) && equiv(p.second, -1.0)));
  GeoMath::sincosd(p, -90.0);
  results.push_back(assertTrue(equiv(p.first, -1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, -0.0);
  results.push_back(assertTrue(equiv(p.first, -0.0) && equiv(p.second, +1.0)));
  GeoMath::sincosd(p, +0.0);
  results.push_back(assertTrue(equiv(p.first, +0.0) && equiv(p.second, +1.0)));
  GeoMath::sincosd(p, +90.0);
  results.push_back(assertTrue(equiv(p.first, +1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, +180.0);
  results.push_back(assertTrue(equiv(p.first, +0.0) && equiv(p.second, -1.0)));
  GeoMath::sincosd(p, +270.0);
  results.push_back(assertTrue(equiv(p.first, -1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, +360.0);
  results.push_back(assertTrue(equiv(p.first, +0.0) && equiv(p.second, +1.0)));
  GeoMath::sincosd(p, +450.0);
  results.push_back(assertTrue(equiv(p.first, +1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, +540.0);
  results.push_back(assertTrue(equiv(p.first, +0.0) && equiv(p.second, -1.0)));
  GeoMath::sincosd(p, +630.0);
  results.push_back(assertTrue(equiv(p.first, -1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, +720.0);
  results.push_back(assertTrue(equiv(p.first, +0.0) && equiv(p.second, +1.0)));
  GeoMath::sincosd(p, +810.0);
  results.push_back(assertTrue(equiv(p.first, +1.0) && equiv(p.second, +0.0)));
  GeoMath::sincosd(p, +inf);
  results.push_back(assertTrue(equiv(p.first, nan) && equiv(p.second, nan)));
  GeoMath::sincosd(p, nan);
  results.push_back(assertTrue(equiv(p.first, nan) && equiv(p.second, nan)));

  for (size_t i = 0; i < results.size(); i++)
  {
    if (results[i] == false)
    {
      success = false;
      failed_cases.push_back(i);
      error_count++;
    }
  }

  logResult("Geodesic::GeoMath::sincosd(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;

  /// Test accuracy of sincosd
  double s1, c1, s2, c2, s3, c3;
  GeoMath::sincosd(p, 9.0);
  s1 = p.first;
  c1 = p.second;
  GeoMath::sincosd(p, 81.0);
  s2 = p.first;
  c2 = p.second;
  GeoMath::sincosd(p, -123456789.0);
  s3 = p.first;
  c3 = p.second;
  assertTrue(equiv(s1, c2));
  assertTrue(equiv(s1, s3));
  assertTrue(equiv(c1, s2));
  assertTrue(equiv(c1, -c3));
}
/// @test
void test_atan2d()
{
  // Test special cases for atan2d
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;

  double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<bool> results = {
      assertTrue(equiv(GeoMath::atan2d(+0.0, -0.0), +180)),
      assertTrue(equiv(GeoMath::atan2d(-0.0, -0.0), -180)),
      assertTrue(equiv(GeoMath::atan2d(+0.0, +0.0), +0.0)),
      assertTrue(equiv(GeoMath::atan2d(-0.0, +0.0), -0.0)),
      assertTrue(equiv(GeoMath::atan2d(+0.0, -1.0), +180)),
      assertTrue(equiv(GeoMath::atan2d(-0.0, -1.0), -180)),
      assertTrue(equiv(GeoMath::atan2d(+0.0, +1.0), +0.0)),
      assertTrue(equiv(GeoMath::atan2d(-0.0, +1.0), -0.0)),
      assertTrue(equiv(GeoMath::atan2d(-1.0, +0.0), -90)),
      assertTrue(equiv(GeoMath::atan2d(-1.0, -0.0), -90)),
      assertTrue(equiv(GeoMath::atan2d(+1.0, +0.0), +90)),
      assertTrue(equiv(GeoMath::atan2d(+1.0, -0.0), +90)),
      assertTrue(equiv(GeoMath::atan2d(+1.0, -inf), +180)),
      assertTrue(equiv(GeoMath::atan2d(-1.0, -inf), -180)),
      assertTrue(equiv(GeoMath::atan2d(+1.0, +inf), +0.0)),
      assertTrue(equiv(GeoMath::atan2d(-1.0, +inf), -0.0)),
      assertTrue(equiv(GeoMath::atan2d(+inf, +1.0), +90)),
      assertTrue(equiv(GeoMath::atan2d(+inf, -1.0), +90)),
      assertTrue(equiv(GeoMath::atan2d(-inf, +1.0), -90)),
      assertTrue(equiv(GeoMath::atan2d(-inf, -1.0), -90)),
      assertTrue(equiv(GeoMath::atan2d(+inf, -inf), +135)),
      assertTrue(equiv(GeoMath::atan2d(-inf, -inf), -135)),
      assertTrue(equiv(GeoMath::atan2d(+inf, +inf), +45)),
      assertTrue(equiv(GeoMath::atan2d(-inf, +inf), -45)),
      assertTrue(equiv(GeoMath::atan2d(nan, +1.0), nan)),
      assertTrue(equiv(GeoMath::atan2d(+1.0, nan), nan)),
  };
  double s = 7e-16;
  results.push_back(assertEquals(GeoMath::atan2d(s, -1.0), 180 - GeoMath::atan2d(s, 1.0), 0));

  for (size_t i = 0; i < results.size(); i++)
  {
    if (results[i] == false)
    {
      success = false;
      failed_cases.push_back(i);
      error_count++;
    }
  }

  logResult("Geodesic::GeoMath::atan2d(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;

  // Test accuracy of atan2d
}
/// @test
void test_sum()
{
  // Test special cases of sum

  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;

  std::vector<bool> results;

  GeoMath::Pair<double> p = GeoMath::Pair<double>();
  GeoMath::sum(p, +9.0, -9.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::sum(p, -9.0, +9.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::sum(p, -0.0, +0.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::sum(p, +0.0, -0.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::sum(p, -0.0, -0.0);
  results.push_back(assertTrue(equiv(p.first, -0.0)));
  GeoMath::sum(p, +0.0, +0.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));

  for (size_t i = 0; i < results.size(); i++)
  {
    if (results[i] == false)
    {
      success = false;
      failed_cases.push_back(i);
      error_count++;
    }
  }

  logResult("Geodesic::GeoMath::sum(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
/// @test
void test_AngNormalize()
{
  // Test special cases of AngNormalize
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;

  std::vector<bool> results = {
      assertTrue(equiv(GeoMath::AngNormalize(-900.0), -180)),
      assertTrue(equiv(GeoMath::AngNormalize(-720.0), -0.0)),
      assertTrue(equiv(GeoMath::AngNormalize(-540.0), -180)),
      assertTrue(equiv(GeoMath::AngNormalize(-360.0), -0.0)),
      assertTrue(equiv(GeoMath::AngNormalize(-180.0), -180)),
      assertTrue(equiv(GeoMath::AngNormalize(-0.0), -0.0)),
      assertTrue(equiv(GeoMath::AngNormalize(+0.0), +0.0)),
      assertTrue(equiv(GeoMath::AngNormalize(180.0), +180)),
      assertTrue(equiv(GeoMath::AngNormalize(360.0), +0.0)),
      assertTrue(equiv(GeoMath::AngNormalize(540.0), +180)),
      assertTrue(equiv(GeoMath::AngNormalize(720.0), +0.0)),
      assertTrue(equiv(GeoMath::AngNormalize(900.0), +180)),
  };

  for (size_t i = 0; i < results.size(); i++)
  {
    if (results[i] == false)
    {
      success = false;
      failed_cases.push_back(i);
      error_count++;
    }
  }

  logResult("Geodesic::GeoMath::AngNormalize(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
/// @test
void test_AngDiff()
{
  // Test special cases of AngDiff
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;

  double eps = GeoMath::ulp(1.0);
  std::vector<bool> results;
  GeoMath::Pair<double> p = GeoMath::Pair<double>();
  GeoMath::AngDiff(p, +0.0, +0.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::AngDiff(p, +0.0, -0.0);
  results.push_back(assertTrue(equiv(p.first, -0.0)));
  GeoMath::AngDiff(p, -0.0, +0.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::AngDiff(p, -0.0, -0.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::AngDiff(p, +5.0, +365.0);
  results.push_back(assertTrue(equiv(p.first, +0.0)));
  GeoMath::AngDiff(p, +365.0, +5.0);
  results.push_back(assertTrue(equiv(p.first, -0.0)));
  GeoMath::AngDiff(p, +5.0, +185.0);
  results.push_back(assertTrue(equiv(p.first, +180.0)));
  GeoMath::AngDiff(p, +185.0, +5.0);
  results.push_back(assertTrue(equiv(p.first, -180.0)));
  GeoMath::AngDiff(p, +eps, +180.0);
  results.push_back(assertTrue(equiv(p.first, +180.0)));
  GeoMath::AngDiff(p, -eps, +180.0);
  results.push_back(assertTrue(equiv(p.first, -180.0)));
  GeoMath::AngDiff(p, +eps, -180.0);
  results.push_back(assertTrue(equiv(p.first, +180.0)));
  GeoMath::AngDiff(p, -eps, -180.0);
  results.push_back(assertTrue(equiv(p.first, -180.0)));

  // Test accuracy of AngDiff
  double x = 138 + 128 * eps, y = -164;
  GeoMath::AngDiff(p, x, y);
  results.push_back(assertEquals(p.first, 58 - 128 * eps, 0));

  for (size_t i = 0; i < results.size(); i++)
  {
    if (results[i] == false)
    {
      success = false;
      failed_cases.push_back(i);
      error_count++;
    }
  }

  logResult("Geodesic::GeoMath::AngDiff(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
/// @test
void test_equatorial_coincident()
{
  // azimuth with coincident point on equator
  //  lat1 lat2 azi1/2
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;
  std::vector<bool> results;

  double C[2][3] = {
      {+0.0, -0.0, 180},
      {-0.0, +0.0, 0}};
  for (int i = 0; i < 2; ++i)
  {
    Geodesic::GeodesicData inv = kObj.Inverse(C[i][0], 0.0, C[i][1], 0.0);
    results.push_back(assertTrue(equiv(inv.azi1, C[i][2])));
    results.push_back(assertTrue(equiv(inv.azi2, C[i][2])));
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

  logResult("Geodesic::KarneyFF::Inverse(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
/// @test
void test_equatorial_NS()
{
  // Does the nearly antipodal equatorial solution go north or south?
  //  lat1 lat2 azi1 azi2
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;
  std::vector<bool> results;

  double C[2][4] = {
      {+0.0, +0.0, 56, 124},
      {-0.0, -0.0, 124, 56}};
  for (int i = 0; i < 2; ++i)
  {
    Geodesic::GeodesicData inv = kObj.Inverse(C[i][0], 0.0, C[i][1], 179.5);
    results.push_back(assertEquals(inv.azi1, C[i][2], 1));
    results.push_back(assertEquals(inv.azi2, C[i][3], 1));
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

  logResult("Geodesic::KarneyFF::Inverse(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
/// @test
void test_antipodal()
{
  // How does the exact antipodal equatorial path go N/S + E/W"""
  //  lat1 lat2 lon2 azi1 azi2
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;
  std::vector<bool> results;

  double C[4][5] = {
      {+0.0, +0.0, +180, +0.0, +180},
      {-0.0, -0.0, +180, +180, +0.0},
      {+0.0, +0.0, -180, -0.0, -180},
      {-0.0, -0.0, -180, -180, -0.0}};
  for (int i = 0; i < 4; ++i)
  {
    Geodesic::GeodesicData inv = kObj.Inverse(C[i][0], 0.0, C[i][1], C[i][2]);
    results.push_back(assertTrue(equiv(inv.azi1, C[i][3])));
    results.push_back(assertTrue(equiv(inv.azi2, C[i][4])));
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

  logResult("Geodesic::KarneyFF::Inverse(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
/// @test
void test_antipodal_prolate()
{
  // Antipodal points on the equator with prolate ellipsoid
  //  lon2 azi1/2
  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;
  std::vector<bool> results;

  double C[2][2] = {
      {+180, +90},
      {-180, -90}};

  configureGeod(6.4e6, -1 / 300.0);
  for (int i = 0; i < 2; ++i)
  {
    Geodesic::GeodesicData inv = kObj.Inverse(0.0, 0.0, 0.0, C[i][0]);
    results.push_back(assertTrue(equiv(inv.azi1, C[i][1])));
    results.push_back(assertTrue(equiv(inv.azi2, C[i][1])));
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

  logResult("Geodesic::KarneyFF::Inverse(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
  configureGeod(GEOD_TYPE::WGS_84);
}
/// @test
void test_azimuth_0_180()
{
  // azimuths = +/-0 and +/-180 for the direct problem
  //  azi1, lon2, azi2

  static int test_counter = 0;
  static int error_count = 0;
  bool success = true;
  std::vector<int> failed_cases;
  std::vector<bool> results;

  double C[4][3] = {
      {+0.0, +180, +180},
      {-0.0, -180, -180},
      {+180, +180, +0.0},
      {-180, -180, -0.0}};
  for (int i = 0; i < 4; ++i)
  {
    Geodesic::GeodesicData dir =
        kObj.Direct(0.0, 0.0, C[i][0], 15e6,
                    GeodesicMask::STANDARD | GeodesicMask::LONG_UNROLL);
    results.push_back(assertTrue(equiv(dir.lon2, C[i][1])));
    results.push_back(assertTrue(equiv(dir.azi2, C[i][2])));
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

  logResult("Geodesic::KarneyFF::Direct(double lat1, double lon1, double azi1, double s12, int outmask)", success, error_count, failed_cases);

  test_counter++;
}
