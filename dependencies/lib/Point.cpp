#include "Point.hpp"
#include "KarneyFF.hpp"
#include <cmath>

using namespace Point;

#pragma region PointConversion

double cgb[6], cbg[6], utg[6], gtu[6];
double Qn, Zb;

static void update_parameters(double *cgb, double *cbg, double *utg, double *gtu, double &Qn, double &Zb)
{
    double lat0 = 0;
    double k0 = 0.9996;
    double es = pow(Geodesic::geodParameters.eccentricity, 2.0);
    double f = es / (1.0 + sqrt(1 - es));
    double n = f / (2.0 - f);
    double np = n;

    cgb[0] = n * (2.0 + n * (-2.0 / 3.0 + n * (-2.0 + n * (116.0 / 45.0 + n * (26.0 / 45.0 + n * (-2854.0 / 675.0))))));
    cbg[0] = n * (-2.0 + n * (2.0 / 3.0 + n * (4.0 / 3.0 + n * (-82.0 / 45.0 + n * (32.0 / 45.0 + n * (4642.0 / 4725.0))))));

    np = np * n;
    cgb[1] = np * (7.0 / 3.0 + n * (-8.0 / 5.0 + n * (-227.0 / 45.0 + n * (2704.0 / 315.0 + n * (2323.0 / 945.0)))));
    cbg[1] = np * (5.0 / 3.0 + n * (-16.0 / 15.0 + n * (-13.0 / 9.0 + n * (904.0 / 315.0 + n * (-1522.0 / 945.0)))));

    np = np * n;
    cgb[2] = np * (56.0 / 15.0 + n * (-136.0 / 35.0 + n * (-1262.0 / 105.0 + n * (73814.0 / 2835.0))));
    cbg[2] = np * (-26.0 / 15.0 + n * (34.0 / 21.0 + n * (8.0 / 5.0 + n * (-12686.0 / 2835.0))));

    np = np * n;
    cgb[3] = np * (4279.0 / 630.0 + n * (-332.0 / 35.0 + n * (-399572.0 / 14175.0)));
    cbg[3] = np * (1237.0 / 630.0 + n * (-12.0 / 5.0 + n * (-24832.0 / 14175.0)));

    np = np * n;
    cgb[4] = np * (4174.0 / 315.0 + n * (-144838.0 / 6237.0));
    cbg[4] = np * (-734.0 / 315.0 + n * (109598.0 / 31185.0));

    np = np * n;
    cgb[5] = np * (601676.0 / 22275.0);
    cbg[5] = np * (444337.0 / 155925.0);

    np = pow(n, 2.0);
    Qn = k0 / (1 + n) * (1 + np * (1.0 / 4.0 + np * (1.0 / 64.0 + np / 256.0)));

    utg[0] = n * (-0.5 + n * (2.0 / 3.0 + n * (-37.0 / 96.0 + n * (1.0 / 360.0 + n * (81.0 / 512.0 + n * (-96199.0 / 604800.0))))));
    gtu[0] = n * (0.5 + n * (-2.0 / 3.0 + n * (5.0 / 16.0 + n * (41.0 / 180.0 + n * (-127.0 / 288.0 + n * (7891.0 / 37800.0))))));

    utg[1] = np * (-1.0 / 48.0 + n * (-1.0 / 15.0 + n * (437.0 / 1440.0 + n * (-46.0 / 105.0 + n * (1118711.0 / 3870720.0)))));
    gtu[1] = np * (13.0 / 48.0 + n * (-3.0 / 5.0 + n * (557.0 / 1440.0 + n * (281.0 / 630.0 + n * (-1983433.0 / 1935360.0)))));

    np = np * n;
    utg[2] = np * (-17.0 / 480.0 + n * (37.0 / 840.0 + n * (209.0 / 4480.0 + n * (-5569.0 / 90720.0))));
    gtu[2] = np * (61.0 / 240.0 + n * (-103.0 / 140.0 + n * (15061.0 / 26880.0 + n * (167603.0 / 181440.0))));

    np = np * n;
    utg[3] = np * (-4397.0 / 161280.0 + n * (11.0 / 504.0 + n * (830251.0 / 7257600.0)));
    gtu[3] = np * (49561.0 / 161280.0 + n * (-179.0 / 168.0 + n * (6601661.0 / 7257600.0)));

    np = np * n;
    utg[4] = np * (-4583.0 / 161280.0 + n * (108847.0 / 3991680.0));
    gtu[4] = np * (34729.0 / 80640.0 + n * (-3418889.0 / 1995840.0));

    np = np * n;
    utg[5] = np * (-20648693.0 / 638668800.0);
    gtu[5] = np * (212378941.0 / 319334400.0);

    double Z = Geodesic::GeoMath::gatg(cbg, 6, lat0);
    Zb = -Qn * (Z + Geodesic::GeoMath::clens(gtu, 6, 2 * Z));
}

UtmPoint *Point::PointConversion::geocentricToUtm(double X, double Y, double Z, POINT_TYPE type)
{
    BasePoint *temp = new GeocentricPoint(X, Y, Z);
    BasePoint *geodesic_point = temp->convert(POINT_TYPE::GEODESIC);
    UtmPoint *point = dynamic_cast<UtmPoint *>(geodesic_point->convert(type));

    delete temp, geodesic_point;
    return point;
}
UtmPoint *Point::PointConversion::geodesicToUtm(double lat, double lon, double alt, POINT_TYPE type)
{
    const double a = Geodesic::geodParameters.semi_major_axis;
    const double b = Geodesic::geodParameters.semi_minor_axis;
    const double f = Geodesic::geodParameters.flattening;
    const double e = sqrt(1 - ((b * b) / (a * a)));
    const double e2 = e * e;
    const double K0 = 0.9996;

    double latRad = Geodesic::Units::toRad(lat);
    double lonRad = Geodesic::Units::toRad(lon);

    double zone = static_cast<short>((lon + 180) / 6) + 1;

    if (type == POINT_TYPE::UTM_Z)
    {
        if (56.0 <= lat && lat < 64.0 && 3.0 <= lon && lon < 12.0)
            zone = 32;
        if (72.0 <= lat && lat < 84.0)
        {
            if (0.0 <= lon && lon < 9.0)
                zone = 31;
            else if (9.0 <= lon && lon < 21.0)
                zone = 33;
            else if (33.0 <= lon && lon < 42.0)
                zone = 37;
        }
    }

    std::string zone_letters = "CDEFGHJKLMNPQRSTUVWXYZ";
    std::string ZoneLetter = std::string(1, zone_letters[static_cast<int>((lat + 80) / 8)]);

    HEMISPHERE hemisphere = (lat >= 0) ? HEMISPHERE::NORTH : HEMISPHERE::SOUTH;

    double lonOrigin = (zone - 1) * 6 - 180 + 3;
    double lonOriginRad = Geodesic::Units::toRad(lonOrigin);

    double N = a / std::sqrt(1 - e2 * std::sin(latRad) * std::sin(latRad));
    double T = std::tan(latRad) * std::tan(latRad);
    double C = e2 / (1 - e2) * std::cos(latRad) * std::cos(latRad);
    double A = std::cos(latRad) * (lonRad - lonOriginRad);

    double M = a * ((1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256) *
                        latRad -
                    (3 * e2 / 8 + 3 * e2 * e2 / 32 + 45 * e2 * e2 * e2 / 1024) *
                        std::sin(2 * latRad) +
                    (15 * e2 * e2 / 256 + 45 * e2 * e2 * e2 / 1024) *
                        std::sin(4 * latRad) -
                    (35 * e2 * e2 * e2 / 3072) * std::sin(6 * latRad));

    double easting = K0 * N * (A + (1 - T + C) * A * A * A / 6 + (5 - 18 * T + T * T + 72 * C - 58 * e2 / (1 - e2)) * A * A * A * A * A / 120) + 500000.0;

    double northing = K0 * (M + N * std::tan(latRad) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A * A * A * A / 24 + (61 - 58 * T + T * T + 600 * C - 330 * e2 / (1 - e2)) * A * A * A * A * A * A / 720));

    if (hemisphere == HEMISPHERE::SOUTH)
    {
        northing += 10000000.0;
    }
    return new UtmPoint(easting, northing, zone, type, hemisphere, ZoneLetter, alt);
}
UtmPoint *Point::PointConversion::geodesicToUtm(double lat, double lon, double alt, short zone, HEMISPHERE hemisphere)
{
    double D2R = 0.01745329251994329577;
    double lat0 = 0;
    double long0 = 0;
    double k0 = 0.9996;
    double x0 = 500000.0;
    double y0 = 0;

    update_parameters(cgb, cbg, utg, gtu, Qn, Zb);

    lon = Geodesic::Units::toRad(lon);
    lat = Geodesic::Units::toRad(lat);

    y0 = hemisphere == HEMISPHERE::SOUTH ? 10000000.0 : 0;
    long0 = ((6 * std::abs(zone)) - 183) * D2R;
    double Ce = Geodesic::GeoMath::adjust_lon(lon - long0);
    double Cn = lat;

    Cn = Geodesic::GeoMath::gatg(cbg, 6, Cn);
    double sin_Cn = sin(Cn);
    double cos_Cn = cos(Cn);
    double sin_Ce = sin(Ce);
    double cos_Ce = cos(Ce);

    Cn = atan2(sin_Cn, cos_Ce * cos_Cn);
    Ce = atan2(sin_Ce * cos_Cn, Geodesic::GeoMath::Hypot(sin_Cn, cos_Cn * cos_Ce));
    Ce = Geodesic::GeoMath::asinHy(tan(Ce));

    double *tmp = Geodesic::GeoMath::clens_cmplx(gtu, 6, 2 * Cn, 2 * Ce);

    Cn = Cn + tmp[0];
    Ce = Ce + tmp[1];

    double x;
    double y;

    if (std::abs(Ce) <= 2.623395162778)
    {
        x = Geodesic::geodParameters.semi_major_axis * (Qn * Ce) + x0;
        y = Geodesic::geodParameters.semi_major_axis * (Qn * Cn + Zb) + y0;
    }
    else
    {
        x = std::numeric_limits<double>::infinity();
        y = std::numeric_limits<double>::infinity();
    }

    int type = (hemisphere == HEMISPHERE::SOUTH) ? zone + 3 + 60 : zone + 3;

    return new UtmPoint(x, y, zone, (POINT_TYPE)type, hemisphere, "_", alt);
}

GeodesicPoint *Point::PointConversion::geocentricToGeodesic(double X, double Y, double Z)
{
    Geodesic::kObj.EquatorialRadius();
    const double a = Geodesic::geodParameters.semi_major_axis;
    const double b = Geodesic::geodParameters.semi_minor_axis;
    const double f = Geodesic::geodParameters.flattening;
    const double e2 = (a * a - b * b) / (a * a);
    const double ep2 = (a * a - b * b) / (b * b);

    const double AD_C = 1.002600;
    const double COS_67P5 = 0.38268343236508977;

    const double genua = 1.E-12;
    const double genua2 = genua * genua;
    short maxiter = 30;

    double W, W2, T0, T1, S0, S1, Sin_B0, Sin3_B0, Cos_B0, Sin_p1, Cos_p1, Rn, Sum;
    int At_Pole;
    double longitude, latitude, altitude;

    At_Pole = 0;
    if (X != 0.0)
    {
        longitude = atan2(Y, X);
    }
    else
    {
        if (Y > 0)
        {
            longitude = M_PI_2;
        }
        else if (Y < 0)
        {
            longitude = -M_PI_2;
        }
        else
        {
            At_Pole = 1;
            longitude = 0.0;
            if (Z > 0.0)
            {
                latitude = M_PI_2;
            }
            else if (Z < 0.0)
            {
                latitude < -M_PI_2;
            }
            else
            {
                latitude = M_PI_2;
                altitude = b;
                longitude = Geodesic::Units::toDeg(longitude);
                latitude = Geodesic::Units::toDeg(latitude);
                return new GeodesicPoint(longitude, latitude, altitude);
            }
        }
    }
    W2 = X * X + Y * Y;
    W = sqrt(W2);
    T0 = Z * AD_C;
    S0 = sqrt(T0 * T0 + W2);
    Sin_B0 = T0 / S0;
    Cos_B0 = W / S0;
    Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0;
    T1 = Z + b * ep2 * Sin3_B0;
    Sum = W - a * e2 * Cos_B0 * Cos_B0 * Cos_B0;
    S1 = sqrt(T1 * T1 + Sum * Sum);
    Sin_p1 = T1 / S1;
    Cos_p1 = Sum / S1;
    Rn = a / sqrt(1.0 - e2 * Sin_p1 * Sin_p1);
    if (Cos_p1 >= COS_67P5)
    {
        altitude = W / Cos_p1 - Rn;
    }
    else if (Cos_p1 <= -COS_67P5)
    {
        altitude = W / -Cos_p1 - Rn;
    }
    else
    {
        altitude = Z / Sin_p1 + Rn * (e2 - 1.0);
    }
    if (At_Pole == 0)
    {
        latitude = atan(Sin_p1 / Cos_p1);
    }

    double P, RR, CT, ST, RX, RK, RN, CPHI0, SPHI0, CPHI, SPHI, SDPHI;
    int iter = 0;

    At_Pole = 0;
    P = sqrt(X * X + Y * Y);
    RR = sqrt(X * X + Y * Y + Z * Z);
    if (P / a < genua)
    {
        At_Pole = 1;
        longitude = 0.;

        if (RR / a < genua)
        {
            latitude = M_PI_2;
            altitude = -b;
            longitude = Geodesic::Units::toDeg(longitude);
            latitude = Geodesic::Units::toDeg(latitude);
            return new GeodesicPoint(longitude, latitude, altitude);
        }
    }
    else
    {
        longitude = atan2(Y, X);
    }
    CT = Z / RR;
    ST = P / RR;
    RX = 1.0 / sqrt(1.0 - e2 * (2.0 - e2) * ST * ST);
    CPHI0 = ST * (1.0 - e2) * RX;
    SPHI0 = CT * RX;
    iter = 0;

    do
    {
        iter++;
        RN = a / sqrt(1.0 - e2 * SPHI0 * SPHI0);
        altitude = P * CPHI0 + Z * SPHI0 - RN * (1.0 - e2 * SPHI0 * SPHI0);
        RK = e2 * RN / (RN + altitude);
        RX = 1.0 / sqrt(1.0 - RK * (2.0 - RK) * ST * ST);
        CPHI = ST * (1.0 - RK) * RX;
        SPHI = CT * RX;
        SDPHI = SPHI * CPHI0 - CPHI * SPHI0;
        CPHI0 = CPHI;
        SPHI0 = SPHI;

    } while (SDPHI * SDPHI > genua2 && iter < maxiter);

    latitude = atan(SPHI / fabs(CPHI));
    longitude = Geodesic::Units::toDeg(longitude);
    latitude = Geodesic::Units::toDeg(latitude);
    return new GeodesicPoint(longitude, latitude, altitude);
}
GeodesicPoint *Point::PointConversion::utmToGeodesic(double easting, double northing, short zone, HEMISPHERE hemisphere, std::string zoneLetter, double altitude)
{
    const double a = Geodesic::geodParameters.semi_major_axis;
    const double b = Geodesic::geodParameters.semi_minor_axis;
    const double f = Geodesic::geodParameters.flattening;
    const double e = sqrt(1 - ((b * b) / (a * a)));

    const double e2 = e * e;  // Eccentricity squared
    const double K0 = 0.9996; // Scale factor
    double e1 = (1 - std::sqrt(1 - e2)) / (1 + std::sqrt(1 - e2));

    double x = easting - 500000.0;
    double y;

    if (hemisphere != HEMISPHERE::UNKNOWN)
        y = (hemisphere == HEMISPHERE::NORTH) ? northing : northing - 10000000.0;
    else
        y = (zoneLetter[0] >= 'N') ? northing : northing - 10000000.0;

    double m = y / K0;
    double mu = m / (a * (1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256));

    double phi1Rad = mu + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * std::sin(2 * mu) + (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) * std::sin(4 * mu) + (151 * e1 * e1 * e1 / 96) * std::sin(6 * mu) + (1097 * e1 * e1 * e1 * e1 / 512) * std::sin(8 * mu);

    double N1 = a / std::sqrt(1 - e2 * std::sin(phi1Rad) * std::sin(phi1Rad));
    double T1 = std::tan(phi1Rad) * std::tan(phi1Rad);
    double C1 = e2 / (1 - e2) * std::cos(phi1Rad) * std::cos(phi1Rad);
    double R1 = a * (1 - e2) / std::pow(1 - e2 * std::sin(phi1Rad) * std::sin(phi1Rad), 1.5);
    double D = x / (N1 * K0);

    double lat = phi1Rad - (N1 * std::tan(phi1Rad) / R1) * (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * e2 / (1 - e2)) * D * D * D * D / 24 + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * e2 / (1 - e2) - 3 * C1 * C1) * D * D * D * D * D * D / 720);
    lat = Geodesic::Units::toDeg(lat);

    double lonOrigin = (zone - 1) * 6 - 180 + 3; // +3 puts origin in middle of zone
    double lon = (D - (1 + 2 * T1 + C1) * D * D * D / 6 + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * e2 / (1 - e2) + 24 * T1 * T1) * D * D * D * D * D / 120) / std::cos(phi1Rad);
    lon = lonOrigin + Geodesic::Units::toDeg(lon);

    return new GeodesicPoint(lon, lat, altitude);
}
GeodesicPoint *Point::PointConversion::utmToGeodesic(double easting, double northing, double zone, HEMISPHERE hemisphere, double altitude)
{
    double D2R = 0.01745329251994329577;
    double lat0 = 0;
    double long0 = 0;
    double k0 = 0.9996;
    double x0 = 500000.0;
    double y0 = 0;

    update_parameters(cgb, cbg, utg, gtu, Qn, Zb);
    y0 = hemisphere == HEMISPHERE::SOUTH ? 10000000.0 : 0;
    long0 = ((6 * std::abs(zone)) - 183) * D2R;

    double Ce = (easting - x0) * (1 / Geodesic::geodParameters.semi_major_axis);
    double Cn = (northing - y0) * (1 / Geodesic::geodParameters.semi_major_axis);

    Cn = (Cn - Zb) / Qn;
    Ce = Ce / Qn;

    double lon;
    double lat;

    if (std::abs(Ce) <= 2.623395162778)
    {
        double *tmp = Geodesic::GeoMath::clens_cmplx(utg, 6, 2 * Cn, 2 * Ce);

        Cn = Cn + tmp[0];
        Ce = Ce + tmp[1];
        Ce = atan(sinh(Ce));

        double sin_Cn = sin(Cn);
        double cos_Cn = cos(Cn);
        double sin_Ce = sin(Ce);
        double cos_Ce = cos(Ce);

        Cn = atan2(sin_Cn * cos_Ce, Geodesic::GeoMath::Hypot(sin_Ce, cos_Ce * cos_Cn));
        Ce = atan2(sin_Ce, cos_Ce * cos_Cn);

        lon = Geodesic::Units::toDeg(Geodesic::GeoMath::adjust_lon(Ce + long0));
        lat = Geodesic::Units::toDeg(Geodesic::GeoMath::gatg(cgb, 6, Cn));
    }
    else
    {
        lon = std::numeric_limits<double>::infinity();
        lat = std::numeric_limits<double>::infinity();
    }

    return new GeodesicPoint(lon, lat, altitude);
}

GeocentricPoint *Point::PointConversion::utmToGeocentric(double easting, double northing, short zone, HEMISPHERE hemisphere, std::string zoneLetter, double altitude, POINT_TYPE type)
{
    BasePoint *temp = new UtmPoint(easting, northing, zone, type, hemisphere, zoneLetter, altitude);
    BasePoint *geodesic_point = temp->convert(POINT_TYPE::GEODESIC);
    GeocentricPoint *point = dynamic_cast<GeocentricPoint *>(geodesic_point->convert(POINT_TYPE::GEOCENTRIC));

    delete temp, geodesic_point;
    return point;
}
GeocentricPoint *Point::PointConversion::geodesicToGeocentric(double lat, double lon, double alt)
{
    const double a = Geodesic::geodParameters.semi_major_axis;
    const double b = Geodesic::geodParameters.semi_minor_axis;
    const double f = Geodesic::geodParameters.flattening;
    const double e2 = (a * a - b * b) / (a * a);
    const double ep2 = (a * a - b * b) / (b * b);

    lat = Geodesic::Units::toRad(lat);
    lon = Geodesic::Units::toRad(lon);

    double Rn, Sin_Lat, Sin2_Lat, Cos_Lat, X, Y, Z;
    bool error_flag = false;

    if (lat < -M_PI_2 && lat > -1.001 * M_PI_2)
    {
        lat = -M_PI_2;
    }
    else if (lat > M_PI_2 && lat < 1.001 * M_PI_2)
    {
        lat = M_PI_2;
    }
    else if ((lat < -M_PI_2) || (lat > M_PI_2))
    {
        error_flag = true;
    }

    if (!error_flag)
    {
        if (lon > M_PI)
            lon -= (2 * M_PI);
        Sin_Lat = sin(lat);
        Cos_Lat = cos(lat);
        Sin2_Lat = Sin_Lat * Sin_Lat;
        Rn = a / sqrt(1.0e0 - e2 * Sin2_Lat);

        X = (Rn + alt) * Cos_Lat * cos(lon);
        Y = (Rn + alt) * Cos_Lat * sin(lon);
        Z = ((Rn * (1 - e2)) + alt) * Sin_Lat;

        return new GeocentricPoint(X, Y, Z);
    }
    return nullptr;
}
#pragma endregion

#pragma region BasePoint
bool Point::BasePoint::primaryTypeValidation(POINT_TYPE p1_type, POINT_TYPE p2_type)
{
    if (p1_type == p2_type)
        return true;
    return false;
}
POINT_TYPE Point::BasePoint::getPointType()
{
    return this->type;
}
#pragma endregion

#pragma region Geodesic Point
GeodesicPoint::GeodesicPoint() : longitude(0), latitude(0), altitude(0)
{
    this->type = POINT_TYPE::GEODESIC;
}
GeodesicPoint::GeodesicPoint(double longitudeA, double latitudeA, double alt) : longitude(longitudeA), latitude(latitudeA), altitude(alt)
{
    this->type = POINT_TYPE::GEODESIC;
}
GeodesicPoint::~GeodesicPoint() {}

double Point::GeodesicPoint::getLongitude()
{
    return this->longitude;
}
double Point::GeodesicPoint::getLatitude()
{
    return this->latitude;
}
double Point::GeodesicPoint::getAltitude()
{
    return this->altitude;
}

void Point::GeodesicPoint::setLongitude(double value)
{
    this->longitude = value;
}
void Point::GeodesicPoint::setLatitude(double value)
{
    this->latitude = value;
}
void Point::GeodesicPoint::setAltitude(double value)
{
    this->altitude = value;
}

double Point::GeodesicPoint::calculateDistance(BasePoint *point, bool arcmode) // returns the distance in meters
{
    if (!primaryTypeValidation(this->type, point->getPointType()))
        return 0.0;

    double lat1 = this->getLatitude(), lon1 = this->getLongitude(), alt1 = this->getAltitude();
    GeodesicPoint *p2 = dynamic_cast<GeodesicPoint *>(point);
    double lat2 = p2->getLatitude(), lon2 = p2->getLongitude(), alt2 = p2->getAltitude();

    Geodesic::GeodesicData data = Geodesic::kObj.Inverse(lat1, lon1, lat2, lon2, Geodesic::GeodesicMask::DISTANCE);
    if (!arcmode)
        return sqrt(pow(data.s12, 2.0) + pow((alt2 - alt1), 2.0));
    return data.a12;
}
BasePoint *Point::GeodesicPoint::convert(POINT_TYPE type)
{
    if (primaryTypeValidation(this->type, type) && this->isValid())
        return nullptr;

    switch (type)
    {
    case POINT_TYPE::GEOCENTRIC:
        return Point::PointConversion::geodesicToGeocentric(this->latitude, this->longitude, this->altitude);
        break;
    case POINT_TYPE::UTM_H:
        return Point::PointConversion::geodesicToUtm(this->latitude, this->longitude, this->altitude, POINT_TYPE::UTM_H);
        break;
    case POINT_TYPE::UTM_Z:
        return Point::PointConversion::geodesicToUtm(this->latitude, this->longitude, this->altitude, POINT_TYPE::UTM_Z);
        break;

    default:
        HEMISPHERE hem = (int)type >= 64 ? HEMISPHERE::SOUTH : HEMISPHERE::NORTH;
        int index = (int)type - 3;
        int zone = index % 60;
        zone = zone == 0 ? 60 : zone;
        return Point::PointConversion::geodesicToUtm(this->latitude, this->longitude, this->altitude, zone, hem);
        break;
    }
    return nullptr;
}


bool Point::GeodesicPoint::compareLocation(BasePoint *point, double tol)
{
    if (primaryTypeValidation(this->type, type))
    {
        GeodesicPoint *c_p = dynamic_cast<GeodesicPoint *>(point);
        double dist = sqrt((pow(c_p->latitude - this->latitude, 2) + pow(c_p->longitude - this->longitude, 2) + pow(c_p->altitude - this->altitude, 2)));
        if (dist <= tol)
            return true;
    }
    return false;
}
bool Point::GeodesicPoint::isEqual(BasePoint *point)
{
    return this == point;
}
void Point::GeodesicPoint::repr()
{
    std::cout << std::fixed << "longitude:" << this->longitude << std::endl;
    std::cout << std::fixed << "latitude:" << this->latitude << std::endl;
    std::cout << std::fixed << "altitude:" << this->altitude << std::endl;
    std::cout << "Type:" << (int)this->type << std::endl;
}
bool Point::GeodesicPoint::isValid()
{
    bool condition_lat = this->latitude <= 90 && this->latitude >= -90;
    bool condition_lon = this->longitude <= 180 && this->longitude >= -180;
    bool condition_alt = this->altitude <= 100000 && this->altitude >= 0; // these values can be changed;

    return condition_alt && condition_lat && condition_alt;
}
#pragma endregion

#pragma region Cartesian Point
GeocentricPoint::GeocentricPoint() : x(0), y(0), z(0)
{
    this->type = POINT_TYPE::GEOCENTRIC;
}
GeocentricPoint::GeocentricPoint(double X, double Y, double Z)
    : x(X), y(Y), z(Z)
{
    this->type = POINT_TYPE::GEOCENTRIC;
}
GeocentricPoint::~GeocentricPoint() {}

double Point::GeocentricPoint::getX()
{
    return this->x;
}
double Point::GeocentricPoint::getY()
{
    return this->y;
}
double Point::GeocentricPoint::getZ()
{
    return this->z;
}

void Point::GeocentricPoint::setX(double value)
{
    this->x = value;
}
void Point::GeocentricPoint::setY(double value)
{
    this->y = value;
}
void Point::GeocentricPoint::setZ(double value)
{
    this->z = value;
}

double Point::GeocentricPoint::calculateDistance(BasePoint *point, bool arcmode)
{
    BasePoint *pointfrom = this->convert(POINT_TYPE::GEODESIC);
    BasePoint *pointTo = point->convert(POINT_TYPE::GEODESIC);

    double distance = pointfrom->calculateDistance(pointTo, arcmode);

    delete pointfrom, pointTo;
    return distance;
}
BasePoint *Point::GeocentricPoint::convert(POINT_TYPE type)
{
    if (primaryTypeValidation(this->type, type) && this->isValid())
        return nullptr;

    switch (type)
    {
    case POINT_TYPE::GEODESIC:
        return Point::PointConversion::geocentricToGeodesic(this->x, this->y, this->z);
        break;
    case POINT_TYPE::UTM_H:
        return Point::PointConversion::geocentricToUtm(this->x, this->y, this->z, POINT_TYPE::UTM_H);
        break;
    case POINT_TYPE::UTM_Z:
        return Point::PointConversion::geocentricToUtm(this->x, this->y, this->z, POINT_TYPE::UTM_Z);
        break;

    default:
        // determine hemisphere zone from entered type
        BasePoint *temp_point = this->convert(POINT_TYPE::GEODESIC);
        BasePoint *utm_zone_specific_point = temp_point->convert(type);
        delete temp_point;
        return utm_zone_specific_point;
        break;
    };
    return nullptr;
}
bool Point::GeocentricPoint::compareLocation(BasePoint *point, double tol)
{
    if (primaryTypeValidation(this->type, type))
    {
        GeocentricPoint *c_p = dynamic_cast<GeocentricPoint *>(point);
        double dist = sqrt((pow(c_p->x - this->x, 2) + pow(c_p->y - this->y, 2) + pow(c_p->z - this->z, 2)));
        if (dist < tol)
            return true;
    }
    return false;
}
bool Point::GeocentricPoint::isEqual(BasePoint *point)
{
    return this == point;
}
void Point::GeocentricPoint::repr()
{
    std::cout << std::fixed << "X:" << this->x << std::endl;
    std::cout << std::fixed << "Y:" << this->y << std::endl;
    std::cout << std::fixed << "Z:" << this->z << std::endl;
    std::cout << "Type:" << (int)this->type << std::endl;
}
bool Point::GeocentricPoint::isValid()
{
    bool xy_condition = (x <= Geodesic::geodParameters.semi_major_axis && x >= -Geodesic::geodParameters.semi_minor_axis) && (y <= Geodesic::geodParameters.semi_major_axis && y >= -Geodesic::geodParameters.semi_minor_axis);
    bool z_condition = z <= Geodesic::geodParameters.semi_major_axis && z >= -Geodesic::geodParameters.semi_minor_axis;

    return xy_condition && z_condition;
}
#pragma endregion

#pragma region UTM POINT
Point::UtmPoint::UtmPoint() : easting(0), northing(0), zone(0), hemisphere(HEMISPHERE::UNKNOWN), zoneLetter("-") { this->type = POINT_TYPE::UTM_H; }
Point::UtmPoint::UtmPoint(double east, double north, short zone, HEMISPHERE hem, double alt) : easting(east), northing(north), zone(zone), hemisphere(hem), altitude(alt) { this->type = POINT_TYPE::UTM_H; }
Point::UtmPoint::UtmPoint(double east, double north, short zone, std::string zoneL, double alt) : easting(east), northing(north), zone(zone), zoneLetter(zoneL), altitude(alt) { this->type = POINT_TYPE::UTM_Z; }
Point::UtmPoint::UtmPoint(double east, double north, short zone, POINT_TYPE type, HEMISPHERE hem, std::string zoneL, double alt) : easting(east), northing(north), zone(zone), hemisphere(hem), zoneLetter(zoneL), altitude(alt) { this->type = type; }
Point::UtmPoint::~UtmPoint() {}

double Point::UtmPoint::getEasting()
{
    return this->easting;
}
double Point::UtmPoint::getNorthing()
{
    return this->northing;
}
short Point::UtmPoint::getZone()
{
    return this->zone;
}
HEMISPHERE Point::UtmPoint::getHemisphere()
{
    return this->hemisphere;
}
double Point::UtmPoint::getAltitude()
{
    return this->altitude;
}

void Point::UtmPoint::setEasting(double value)
{
    this->easting = value;
}
void Point::UtmPoint::setNorthing(double value)
{
    this->northing = value;
}
void Point::UtmPoint::setZone(short value)
{
    this->zone = value;
}
void Point::UtmPoint::setHemisphere(HEMISPHERE value)
{
    this->hemisphere = value;
}
void Point::UtmPoint::setAltitude(double value)
{
    this->altitude = value;
}

double Point::UtmPoint::calculateDistance(BasePoint *point, bool arcmode)
{
    BasePoint *pointfrom = this->convert(POINT_TYPE::GEODESIC);
    BasePoint *pointTo = point->convert(POINT_TYPE::GEODESIC);

    double distance = pointfrom->calculateDistance(pointTo, arcmode);

    delete pointfrom, pointTo;
    return distance;
}
BasePoint *Point::UtmPoint::convert(POINT_TYPE type)
{
    if (primaryTypeValidation(this->type, type) && this->isValid())
        return nullptr;

    switch (type)
    {
    case POINT_TYPE::GEODESIC:
        if (this->type == POINT_TYPE::UTM_H)
            return Point::PointConversion::utmToGeodesic(this->easting, this->northing, this->zone, this->hemisphere, "_", this->altitude);
        else if (this->type == POINT_TYPE::UTM_Z)
            return Point::PointConversion::utmToGeodesic(this->easting, this->northing, this->zone, HEMISPHERE::UNKNOWN, this->zoneLetter, this->altitude);
        else
            return Point::PointConversion::utmToGeodesic(this->easting, this->northing, this->zone, this->hemisphere, this->altitude);
        break;

    case POINT_TYPE::GEOCENTRIC:
        return Point::PointConversion::utmToGeocentric(this->easting, this->northing, this->zone, this->hemisphere, this->zoneLetter, this->altitude, this->getPointType());
        break;

    case POINT_TYPE::UTM_Z:
        this->type = POINT_TYPE::UTM_Z;
        return this;
        break;

    case POINT_TYPE::UTM_H:
        this->type = POINT_TYPE::UTM_H;
        return this;
        break;

    default:
        BasePoint *temp_point = this->convert(POINT_TYPE::GEODESIC);
        BasePoint *utm_zone_specific_point = temp_point->convert(type);
        delete temp_point;
        return utm_zone_specific_point;
        break;
    }
    return nullptr;
}
bool Point::UtmPoint::compareLocation(BasePoint *point, double tol)
{
    if (primaryTypeValidation(this->type, point->getPointType()))
    {
        UtmPoint *c_p = dynamic_cast<UtmPoint *>(point);
        double dist = sqrt((pow(c_p->northing - this->northing, 2) + pow(c_p->easting - this->easting, 2) + pow(c_p->zone - this->zone, 2)));
        if (dist <= tol)
            return true;
    }
    return false;
}
bool Point::UtmPoint::isEqual(BasePoint *point)
{
    return this == point;
}
void Point::UtmPoint::repr()
{
    std::cout << std::fixed << "easting: " << this->easting << std::endl;
    std::cout << std::fixed << "northing: " << this->northing << std::endl;
    std::cout << std::fixed << "zone: " << this->zone;
    if (this->type == POINT_TYPE::UTM_Z)
        std::cout << this->zoneLetter << std::endl;
    else
        this->hemisphere == HEMISPHERE::NORTH ? std::cout << "N" << std::endl : std::cout << "S" << std::endl;
    std::cout << std::fixed << "altitude: " << this->altitude << std::endl;
    std::cout << std::fixed << "Type : " << (int)this->type << std::endl;
}
bool Point::UtmPoint::isValid()
{
    bool condition_Zone = this->zone <= 60 && this->zone >= 0;
    bool condition_Hemisphere = this->hemisphere != HEMISPHERE::UNKNOWN;
    bool condition_Easting = this->easting <= 900000 && this->easting >= 100000;
    bool condition_Northing = this->northing <= 10000000 && this->northing >= 0;

    return condition_Zone && condition_Hemisphere && condition_Easting && condition_Northing;
}
#pragma endregion
