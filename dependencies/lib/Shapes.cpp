#include "Geodesic.hpp"
#include "Shapes.hpp"
#include "KarneyFF.hpp"
#include <cmath>

using namespace Shapes;
using namespace Geodesic;

#pragma region FILTER

Shapes::Filter::Filter()
{
    this->shapes = {};
    this->type = Shapes::SHAPES::FILTER;
}
Shapes::Filter::Filter(std::vector<Filter *> shapes)
{
    this->shapes = shapes;
    this->type = Shapes::SHAPES::FILTER;
}
Shapes::Filter::~Filter() {}
SHAPES Shapes::Filter::getShapeType()
{
    return this->type;
}
#pragma endregion

#pragma region Circle
Circle::Circle()
{
    this->radius = 0;
    this->type = Shapes::SHAPES::CIRCLE;
    this->centerPoint = nullptr;
}
Circle::Circle(BasePoint *center_point, double radius, double height)
{
    this->centerPoint = center_point;
    this->radius = radius;
    this->height = height;
    this->type = Shapes::SHAPES::CIRCLE;
}
Circle::~Circle()
{
    delete centerPoint;
}
BasePoint *Shapes::Circle::getCenter()
{
    return this->centerPoint;
}
double Shapes::Circle::getRadius()
{
    return this->radius;
}
double Shapes::Circle::getHeight()
{
    return this->height;
}
void Shapes::Circle::setCenter(BasePoint *value)
{
    delete this->centerPoint;
    this->centerPoint = value;
}
void Shapes::Circle::setRadius(double value)
{
    this->radius = value;
}
void Shapes::Circle::setHeight(double value)
{
    this->height = value;
}
double Shapes::Circle::calculatePerimeter()
{
    return 2 * M_PI * radius;
}
double Circle::calculateArea()
{
    return M_PI * radius * radius;
}
void Shapes::Circle::move(double azimuth, double distance, bool isArcLength)
{
    GeodesicPoint *point = dynamic_cast<GeodesicPoint *>(this->centerPoint);
    GeodesicData data = kObj.Direct(point->getLatitude(), point->getLongitude(), azimuth, distance);

    point->setLongitude(data.lon2);
    point->setLatitude(data.lat2);
}
bool Circle::doesContain(BasePoint *point)
{
    bool converted_flag = false;
    GeodesicPoint *pointfrom, *pointTo;
    if (!(this->centerPoint->getPointType() == POINT_TYPE::GEODESIC))
    {
        pointfrom = (GeodesicPoint *)this->centerPoint->convert(POINT_TYPE::GEODESIC);
        pointTo = (GeodesicPoint *)point->convert(POINT_TYPE::GEODESIC);
        converted_flag = true;
    }
    else
    {
        pointfrom = (GeodesicPoint *)this->centerPoint;
        pointTo = (GeodesicPoint *)point;
    }

    if (pointTo->getAltitude() < pointfrom->getAltitude() - this->height || pointTo->getAltitude() > pointfrom->getAltitude() + this->height)
        return false;

    double distance = pointfrom->calculateDistance(pointTo);
    double inside = distance <= radius;
    return inside;
}
#pragma endregion

#pragma region Polygon
int Shapes::Polygon::transit(double lon1, double lon2)
{
    GeoMath::Pair<double> p = GeoMath::Pair<double>();
    GeoMath::AngDiff(p, lon1, lon2);
    lon1 = GeoMath::AngNormalize(lon1);
    lon2 = GeoMath::AngNormalize(lon2);
    double lon12 = p.first;
    return lon12 > 0 && ((lon1 < 0 && lon2 >= 0) ||
                         (lon1 > 0 && lon2 == 0))
               ? 1
               : (lon12 < 0 && lon1 >= 0 && lon2 < 0 ? -1 : 0);
}
double Shapes::Polygon::AreaReduceA(double area, double area0, int crossings, bool reverse, bool sign)
{
    area = GeoMath::IeeeRemainder(area, area0);
    if ((crossings & 1) != 0)
        area += ((area < 0 ? 1 : -1) * area0 / 2);
    if (!reverse)
        area *= -1;
    // If sign put area in (-area0/2, area0/2], else put area in [0, area0)
    if (sign)
    {
        if (area > area0 / 2)
            area += (-area0);
        else if (area <= -area0 / 2)
            area += (+area0);
    }
    else
    {
        if (area >= area0)
            area += (-area0);
        else if (area < 0)
            area += (+area0);
    }
    return 0 + area;
}
void Shapes::Polygon::clear()
{
    _num = 0;
    _crossings = 0;
    _perimetersum = 0.0;
    _areasum = 0.0;
    _areasum0 = 0.0;
    _perimetersum0 = 0.0;
    _lat0 = _lon0 = _lat1 = _lon1 = std::numeric_limits<double>::quiet_NaN();
    _mask = GeodesicMask::LATITUDE | GeodesicMask::LONGITUDE |
            GeodesicMask::DISTANCE | GeodesicMask::AREA | GeodesicMask::LONG_UNROLL;
}
void Shapes::Polygon::compute(bool reverse, bool sign)
{
    GeodesicData g = kObj.Inverse(_lat1, _lon1, _lat0, _lon0, _mask);
    double x = _areasum0 + g.S12;
    _areasum = AreaReduceA(x, _area0, _crossings + transit(_lon1, _lon0), reverse, sign); // this is to prevent cumulative summation
    _perimetersum = _perimetersum0 + g.s12;
}
Polygon::Polygon()
{
    this->vertices = std::vector<BasePoint *>();
    this->type = Shapes::SHAPES::POLYGON;
    this->height = 0;
    clear();
}
Polygon::Polygon(std::vector<BasePoint *> vertices, double height)
{
    this->type = Shapes::SHAPES::POLYGON;
    this->height = height;
    clear();

    for (int i = 0; i < vertices.size(); i++)
    {
        addPoint(vertices[i]);
    }
}
Polygon::~Polygon()
{
    for (size_t i = 0; i < this->vertices.size(); i++)
    {
        delete vertices[i];
    }
}
std::vector<BasePoint *> Shapes::Polygon::getVertices()
{
    return this->vertices;
}
double Shapes::Polygon::getHeight()
{
    return this->height;
}
void Shapes::Polygon::setHeight(double value)
{
    this->height = height;
}
void Shapes::Polygon::addPoint(BasePoint *point)
{
    this->vertices.push_back(point);
    // type check safety is needed before casting
    double lat = dynamic_cast<GeodesicPoint *>(point)->getLatitude();
    double lon = dynamic_cast<GeodesicPoint *>(point)->getLongitude();
    if (_num == 0)
    {
        _lat0 = _lat1 = lat;
        _lon0 = _lon1 = lon;
    }
    else
    {
        GeodesicData g = kObj.Inverse(_lat1, _lon1, lat, lon, _mask);
        _perimetersum0 += g.s12;
        _areasum0 += g.S12;
        _crossings += transit(_lon1, lon);
        _lat1 = lat;
        _lon1 = lon;
    }
    ++_num;
}
double Shapes::Polygon::calculateArea()
{
    compute(true, true);
    return this->_areasum;
}
void Shapes::Polygon::move(double azimuth, double distance, bool isArcLength)
{
    int numVertices = this->vertices.size();
    for (size_t i = 0; i < numVertices; i++)
    {
        GeodesicPoint *c_point = dynamic_cast<GeodesicPoint *>(this->vertices[i]);
        GeodesicData data = kObj.Direct(c_point->getLatitude(), c_point->getLongitude(), azimuth, distance);

        c_point->setLongitude(data.lon2);
        c_point->setLatitude(data.lat2);
    }
}
double Shapes::Polygon::calculatePerimeter()
{
    compute(true, true);
    return this->_perimetersum;
}
bool Polygon::doesContain(BasePoint *point)
{
    int n = this->vertices.size();
    GeodesicPoint *c_point = dynamic_cast<GeodesicPoint *>(point);

    if (n < 3)
        return false; // A polygon must have at least 3 vertices
    // altitude comparison is made relative to the altitude of the point at 0th index
    double refAlt = dynamic_cast<GeodesicPoint *>(this->vertices[0])->getAltitude();

    if (c_point->getAltitude() < (refAlt - this->height / 2) || c_point->getAltitude() > (refAlt + this->height / 2))
        return false;

    bool inside = false;
    double x = c_point->getLongitude(), y = c_point->getLatitude();
    double xi, yi, xj, yj;
    GeodesicPoint *x_i, *x_j, *y_i, *y_j;
    bool converted_flag = false;

    for (int i = 0, j = n - 1; i < n; j = i++)
    {
        if (point->getPointType() == POINT_TYPE::GEODESIC)
        {
            x_i = dynamic_cast<GeodesicPoint *>(vertices[i]);
            xi = x_i->getLongitude();
            y_i = dynamic_cast<GeodesicPoint *>(vertices[i]);
            yi = y_i->getLatitude();
            x_j = dynamic_cast<GeodesicPoint *>(vertices[j]);
            xj = x_j->getLongitude();
            y_j = dynamic_cast<GeodesicPoint *>(vertices[j]);
            yj = y_j->getLatitude();
        }
        else
        {
            x_i = dynamic_cast<GeodesicPoint *>(vertices[i]->convert(POINT_TYPE::GEODESIC));
            xi = x_i->getLongitude();
            y_i = dynamic_cast<GeodesicPoint *>(vertices[i]->convert(POINT_TYPE::GEODESIC));
            yi = y_i->getLatitude();
            x_j = dynamic_cast<GeodesicPoint *>(vertices[j]->convert(POINT_TYPE::GEODESIC));
            xj = x_j->getLongitude();
            y_j = dynamic_cast<GeodesicPoint *>(vertices[j]->convert(POINT_TYPE::GEODESIC));
            yj = y_j->getLatitude();
            converted_flag = true;
        }

        bool intersect = ((yi > y) != (yj > y)) &&
                         (x < (xj - xi) * (y - yi) / (yj - yi) + xi);

        if (converted_flag)
            delete x_i,x_j,y_i,y_j;

        if (intersect)
            inside = !inside;
    }

    return inside;
}
#pragma endregion

#pragma region Ellipsoid
Shapes::Ellipse::Ellipse() : centerPoint(nullptr), semiMajorAxis(0), semiMinorAxis(0), height(0), orientation(0) {}
Shapes::Ellipse::Ellipse(BasePoint *centerPoint, double semiMajorAxis, double semiMinorAxis, double height, double orientation)   
{
    this->centerPoint = centerPoint;
    this->semiMajorAxis = semiMajorAxis;
    this->semiMinorAxis = semiMinorAxis;
    this->orientation = orientation;
    this->height = height;
}
Shapes::Ellipse::~Ellipse()
{
    delete centerPoint;
}
BasePoint *Shapes::Ellipse::getCenter()
{
    return this->centerPoint;
}
double Shapes::Ellipse::getMajorAxis()
{
    return this->semiMajorAxis;
}
double Shapes::Ellipse::getMinorAxis()
{
    return this->semiMinorAxis;
}
double Shapes::Ellipse::getHeight()
{
    return this->height;
}
void Shapes::Ellipse::setCenter(BasePoint *value)
{
    delete this->centerPoint;
    this->centerPoint = value;
}
void Shapes::Ellipse::setMajorAxis(double value)
{
    this->semiMajorAxis = value;
}
void Shapes::Ellipse::setMinorAxis(double value)
{
    this->semiMinorAxis = value;
}
void Shapes::Ellipse::setHeight(double value)
{
    this->height = value;
}
double Shapes::Ellipse::calculatePerimeter()
{
    // approximation
    return M_PI * (3 * (semiMajorAxis + semiMinorAxis) - sqrt((3 * semiMajorAxis + semiMinorAxis) * (semiMajorAxis + 3 * semiMinorAxis)));
}
double Shapes::Ellipse::calculateArea()
{
    return M_PI * semiMajorAxis * semiMinorAxis;
}
void Shapes::Ellipse::move(double azimuth, double distance, bool isArcLength)
{
    GeodesicPoint *point = dynamic_cast<GeodesicPoint *>(this->centerPoint);
    GeodesicData data = kObj.Direct(point->getLatitude(), point->getLongitude(), azimuth, distance);

    point->setLongitude(data.lon2);
    point->setLatitude(data.lat2);
}
bool Shapes::Ellipse::doesContain(BasePoint *point)
{
    bool converted_flag = false;
    GeodesicPoint *pointfrom, *pointTo;
    if (!(this->centerPoint->getPointType() == POINT_TYPE::GEODESIC))
    {
        pointfrom = (GeodesicPoint *)this->centerPoint->convert(POINT_TYPE::GEODESIC);
        pointTo = (GeodesicPoint *)point->convert(POINT_TYPE::GEODESIC);
        converted_flag = true;
    }

    if (pointTo->getAltitude() < pointfrom->getAltitude() - this->height || pointTo->getAltitude() > pointfrom->getAltitude() + this->height)
        return false;

    double r_a = this->semiMajorAxis;
    double r_b = this->semiMinorAxis;
    double thetaoffset = this->orientation;

    pointfrom = dynamic_cast<GeodesicPoint *>(this->centerPoint);
    pointTo = dynamic_cast<GeodesicPoint *>(point);

    GeodesicData data = kObj.Inverse(pointfrom->getLatitude(), pointfrom->getLongitude(), pointTo->getLatitude(), pointTo->getLongitude());
    double theta_base = data.azi1;
    double distance = data.s12;

    double dif_theta = theta_base - thetaoffset;
    double d_major = distance * cos(Geodesic::Units::toRad(dif_theta));
    double d_minor = distance * sin(Geodesic::Units::toRad(dif_theta));

    bool inside = (pow((d_major / r_a), 2.0) + pow((d_minor / r_b), 2.0)) <= 1;

    if (converted_flag)
    {
        delete pointfrom, pointTo;
    }

    return inside;
}
#pragma endregion
