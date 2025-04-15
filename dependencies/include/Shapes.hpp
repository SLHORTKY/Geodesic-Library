#pragma once
#include <vector>
#include "Point.hpp"
#include "KarneyFF.hpp"

using namespace Point;

namespace Shapes
{
    enum class SHAPES
    {
        CIRCLE,
        POLYGON,
        FILTER
    };

    /// @brief non abstract parent class for all of the subsequent shapes this class can be used to represent instances of its child classes
    /// @brief also shapes vector can be used to create a larger air filter by combining individual shapes.
    class Filter
    {
    public:
        std::vector<Filter *> shapes;

    protected:
        SHAPES type; // to differentiate between shapes for type casting

    public:
        Filter();
        Filter(std::vector<Filter *> shapes);
        ~Filter();


        virtual void move(double azimuth, double distance, bool isArcLength) = 0;
        SHAPES getShapeType();
        
        virtual bool doesContain(BasePoint *point) = 0;

        virtual double calculatePerimeter() = 0;
      
        virtual double calculateArea() = 0;
    };
    /// @brief  simple representation of a Circle even though it can be initialized with multiple types of points
    /// be carful while doing so mostly should be used with GeodesicPoints
    class Circle : public Filter
    {
    private:
        BasePoint *centerPoint;
        double radius;
        double height; // this is a referential quantity that
                       // determines the half height of the shape in 3D relative to the altitude of a point defining the shape

    public:
        Circle();
        Circle(BasePoint *centerPoint, double radius, double height = 0);
        ~Circle();

        BasePoint *getCenter();
        double getRadius();
        double getHeight();
        
        /// @attention this function clears the pointer assigned to it before reassigning
        /// implicit objects that are fed through the constructor may be deleted;
        void setCenter(BasePoint *value);
        void setRadius(double value);
        void setHeight(double value);

        double calculatePerimeter() override;
        double calculateArea() override;
        void move(double azimuth, double distance, bool isArcLength = 0) override;
        bool doesContain(BasePoint *point) override;
    };

    /// @brief simple representation of a polygon defined by n BasePoints*
    class Polygon : public Filter
    {
    private:
        std::vector<BasePoint *> vertices;
        double height;     // this is a referential quantity that
                           // determines the height of the shape in 3D space relative to the altitude of a point defining the shape
        double _area0 = Geodesic::kObj.EllipsoidArea();
        double _lat0, _lon0, _lat1, _lon1;
        int _num, _crossings, _mask;
        double _perimetersum, _areasum, _areasum0,_perimetersum0;
        //_areasum0 and _perimetersum0 are the calculations just befor the calculation is done for the last edge;
        //they are added to prevent cumulative increase of variables _areasum and _perimetersum with consecutive calls of the 
        //calculate Area and calculate perimeter  

    private: 
        int transit(double lon1, double lon2);
        double AreaReduceA(double area, double area0, int crossings, bool reverse, bool sign);
        void compute(bool reverse, bool sign);
        void clear();

    public:
        Polygon();
        Polygon(std::vector<BasePoint *> vertices, double height = 0);
        ~Polygon();

        std::vector<BasePoint *> getVertices();

        double getHeight();
        void setHeight(double value);

        void addPoint(BasePoint *point);
        double calculatePerimeter() override;
        double calculateArea() override;
        void move(double azimuth, double distance, bool isArcLength = 0) override;
        bool doesContain(BasePoint *point) override;
    };

    
    /// @brief simple representation of an ellipse
    class Ellipse : public Filter
    {
    private:
        BasePoint *centerPoint;
        double semiMajorAxis;
        double semiMinorAxis;
        double height;      // this is a referential quantity that
                            // determines the positive / negative half height of the shape in 3D relative to the altitude of a point defining the shape
        double orientation; // theta relative to semi major axis of the Earth(Geoid)
    public:
        Ellipse();
        Ellipse(BasePoint *centerPoint, double semiMajorAxis, double semiMinorAxis, double height, double orientation = 0.0);
        ~Ellipse();

        BasePoint *getCenter();
        double getMajorAxis();
        double getMinorAxis();
        double getHeight();

        /// @attention this function clears the pointer assigned to it before reassigning
        /// implicit objects that are fed through the constructor may be deleted;
        void setCenter(BasePoint *value);
        void setMajorAxis(double value);
        void setMinorAxis(double value);
        void setHeight(double value);

        double calculatePerimeter() override;
        double calculateArea();
        void move(double azimuth, double distance, bool isArcLength = 0) override;
        bool doesContain(BasePoint *point) override;
    };

} // namespace GeodesicShapes
