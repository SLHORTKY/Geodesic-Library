#pragma once
#include <string>

namespace Point
{
    // point type used in conversion between points by using convert(POINT_TYPE) and may provide useful in dynamic_cast
    enum class POINT_TYPE // for utm zone specific enums one hemisphere bit and the rest is to represent 60 different types
    {
        GEODESIC,
        GEOCENTRIC,
        UTM_H,  // Hemispheric definition default
        UTM_Z,  // ZoneLetter definition
        UTM_1N, // do not change the locations of utm zone specific conversion a mathematical formula determines the type
        UTM_2N, // if you want to add other types add to the end of the list
        UTM_3N,
        UTM_4N,
        UTM_5N,
        UTM_6N,
        UTM_7N,
        UTM_8N,
        UTM_9N,
        UTM_10N,
        UTM_11N,
        UTM_12N,
        UTM_13N,
        UTM_14N,
        UTM_15N,
        UTM_16N,
        UTM_17N,
        UTM_18N,
        UTM_19N,
        UTM_20N,
        UTM_21N,
        UTM_22N,
        UTM_23N,
        UTM_24N,
        UTM_25N,
        UTM_26N,
        UTM_27N,
        UTM_28N,
        UTM_29N,
        UTM_30N,
        UTM_31N,
        UTM_32N,
        UTM_33N,
        UTM_34N,
        UTM_35N,
        UTM_36N,
        UTM_37N,
        UTM_38N,
        UTM_39N,
        UTM_40N,
        UTM_41N,
        UTM_42N,
        UTM_43N,
        UTM_44N,
        UTM_45N,
        UTM_46N,
        UTM_47N,
        UTM_48N,
        UTM_49N,
        UTM_50N,
        UTM_51N,
        UTM_52N,
        UTM_53N,
        UTM_54N,
        UTM_55N,
        UTM_56N,
        UTM_57N,
        UTM_58N,
        UTM_59N,
        UTM_60N,
        UTM_1S,
        UTM_2S,
        UTM_3S,
        UTM_4S,
        UTM_5S,
        UTM_6S,
        UTM_7S,
        UTM_8S,
        UTM_9S,
        UTM_10S,
        UTM_11S,
        UTM_12S,
        UTM_13S,
        UTM_14S,
        UTM_15S,
        UTM_16S,
        UTM_17S,
        UTM_18S,
        UTM_19S,
        UTM_20S,
        UTM_21S,
        UTM_22S,
        UTM_23S,
        UTM_24S,
        UTM_25S,
        UTM_26S,
        UTM_27S,
        UTM_28S,
        UTM_29S,
        UTM_30S,
        UTM_31S,
        UTM_32S,
        UTM_33S,
        UTM_34S,
        UTM_35S,
        UTM_36S,
        UTM_37S,
        UTM_38S,
        UTM_39S,
        UTM_40S,
        UTM_41S,
        UTM_42S,
        UTM_43S,
        UTM_44S,
        UTM_45S,
        UTM_46S,
        UTM_47S,
        UTM_48S,
        UTM_49S,
        UTM_50S,
        UTM_51S,
        UTM_52S,
        UTM_53S,
        UTM_54S,
        UTM_55S,
        UTM_56S,
        UTM_57S,
        UTM_58S,
        UTM_59S,
        UTM_60S,
    };
    enum class HEMISPHERE
    {
        SOUTH,
        NORTH,
        UNKNOWN,
    };

    /// @brief abstract class to polymorph all the points
    /// UTM Geodesic and Geocentric ... points it can also be used to convert between these types
    class BasePoint
    {
    protected:
        POINT_TYPE type;
    public:
        
        static bool primaryTypeValidation(POINT_TYPE p1_type, POINT_TYPE p2_type);
        POINT_TYPE getPointType();
        
        virtual double calculateDistance(BasePoint *point, bool arcmode = 0) = 0;
       
        virtual BasePoint *convert(POINT_TYPE type) = 0;
        virtual bool compareLocation(BasePoint *point, double tol = 1e-5) = 0;
       
        virtual bool isEqual(BasePoint *point) = 0;
        virtual ~BasePoint() = default;
        virtual void repr() = 0; // this is a simple representation function for testing and debugging can be deleted
        virtual bool isValid() = 0;
    };

    /// @brief a representation for a position in lon,lat,alt has the ability to calculate distances and type conversion
    class GeodesicPoint : public BasePoint
    {
    private:
        double longitude;
        double latitude;
        double altitude;

    public:
        GeodesicPoint();
        GeodesicPoint(double longitude, double latitude, double altitude = 0);
        ~GeodesicPoint();

        double getLongitude();
        double getLatitude(); // getters for Geodesic Point
        double getAltitude();

        void setLongitude(double value);
        void setLatitude(double value); // setters for Geodesic Point
        void setAltitude(double value);

        double calculateDistance(BasePoint *point, bool arcmode = 0) override;
        BasePoint *convert(POINT_TYPE type) override;

    
        bool compareLocation(BasePoint *point, double tol = 1e-5) override;

        virtual bool isEqual(BasePoint *point) override;
        void repr(); // this is a simple representation function for testing and debugging can be deleted
        bool isValid() override;
    };

    /// @brief a representation for a cartesian point x,y,z geocentric because 0,0,0 is the center of the globe
    class GeocentricPoint : public BasePoint
    {
    private:
        // Rectangular
        double x;
        double y;
        double z;

    public:
        GeocentricPoint(); // empty constructor defines initial values most likely sets them to 0 or to a defined default
        GeocentricPoint(double X, double Y, double Z);
        ~GeocentricPoint();

        double getX();
        double getY(); // getter
        double getZ();

        void setX(double value);
        void setY(double value); // setter can be altered to restrict the user no predefined restrictions
        void setZ(double value);

        double calculateDistance(BasePoint *point, bool arcmode = 0) override;
        BasePoint *convert(POINT_TYPE type) override;


        bool compareLocation(BasePoint *point, double tol = 1e-5) override;
        virtual bool isEqual(BasePoint *point) override;
        void repr() override; // this is a simple representation function for testing and debugging can be deleted
        bool isValid() override;
    };

    /// @brief an alterantive representation for Universal Transverse Mercator
    class UtmPoint : public BasePoint
    {
    private:
        double easting;
        double northing;
        short zone;
        HEMISPHERE hemisphere;
        double altitude;
        std::string zoneLetter;

    public:
        UtmPoint(); // empty constructor defines initial values most likely sets them to 0 or to a defined default
        UtmPoint(double east, double north, short zone, HEMISPHERE hem, double alt = 0.0);
        UtmPoint(double east, double north, short zone, std::string zoneL, double alt = 0.0);
        UtmPoint(double east, double north, short zone, POINT_TYPE type, HEMISPHERE hem = HEMISPHERE::UNKNOWN, std::string zoneL = "_", double alt = 0.0);

        ~UtmPoint();

        double getEasting();
        double getNorthing(); // getters
        short getZone();
        HEMISPHERE getHemisphere();
        double getAltitude();

        void setEasting(double value);
        void setNorthing(double value); // setter can be altered to restrict the user no predefined restrictions
        void setZone(short value);
        void setHemisphere(HEMISPHERE value);
        void setAltitude(double value);

        double calculateDistance(BasePoint *point, bool arcmode = 0) override;
        BasePoint *convert(POINT_TYPE type) override;

        bool compareLocation(BasePoint *point, double tol = 1e-12) override;

        virtual bool isEqual(BasePoint *point) override;
        void repr(); // this is a simple representation function for testing and debugging can be deleted
        bool isValid() override;
    };

    /// @brief a class to manage conversions between 3 or more types of points  Geodesic Geocentric and UTM
    /// keep in mind some functions may require the conversion to implicit GeodesicPoints
    class PointConversion
    {
    public:
        static GeodesicPoint *geocentricToGeodesic(double X, double Y, double Z);
        static UtmPoint *geocentricToUtm(double X, double Y, double Z, POINT_TYPE type);

        static GeocentricPoint *geodesicToGeocentric(double lat, double lon, double alt);
        ///@brief conversion to UTM directly we have no say in zone and hemisphere or zoneLetter it configures it
        /// depending on lat lon coordinates
        static UtmPoint *geodesicToUtm(double lat, double lon, double alt, POINT_TYPE type);
        ///@brief conversion to a specific zone and hemisphere of UTM
        static UtmPoint *geodesicToUtm(double lat, double lon, double alt, short zone, HEMISPHERE hemisphere);

        static GeodesicPoint *utmToGeodesic(double easting, double northing, short zone, HEMISPHERE hemisphere = HEMISPHERE::UNKNOWN, std::string zoneLetter = "-", double altitude = 0.0);
        static GeocentricPoint *utmToGeocentric(double easting, double northing, short zone, HEMISPHERE hemisphere, std::string zoneLetter, double altitude, POINT_TYPE type);
        static GeodesicPoint *utmToGeodesic(double easting, double northing, double zone, HEMISPHERE hemisphere, double altitude);
    };
}
