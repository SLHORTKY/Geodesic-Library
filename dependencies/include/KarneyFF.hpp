#pragma once
#include "GeoMath.hpp"
#include <limits>

/// @brief  this file includes the implementation of the GeographicLib Geodesic library Geodesic class in the original implementation has been modified to be Karney FF
/// the required libraries for this implementation to work are GeoMath Units GeodesicLine if you wish to separate this file from this project
/// all you need to do is transfer GeodesicParams and GeodesicData structs to this file finally the complete documentation can be found at
/// https://geographiclib.sourceforge.io/doc/library.html

namespace Geodesic
{
    // Kobj is KarneyFF composition object in Geodesic namespace executes -> karneyFF() constructor and is an object of that class;
    // it is our main gateway to using KarneyFF Geodesic Functions with predefined geodesic reference frames.
    class KarneyFF;
    extern KarneyFF kObj;

    typedef struct GeodesicData
    {
        double lat1;
        double lon1;
        double azi1;
        double lat2;
        double lon2;
        double azi2;
        double s12;
        double a12;
        double m12;
        double M12;
        double M21;
        double S12;
        GeodesicData()
        {
            lat1 = lon1 = azi1 = lat2 = lon2 = azi2 =
                s12 = a12 = m12 = M12 = M21 = S12 = std::numeric_limits<double>::quiet_NaN();
        }
    } GeodesicData;

    enum class GEOD_TYPE
    {
        WGS_84,
        GRS_80,
        AIRY_1830,
        INTL_1924,
        CLARKE_1880,
        GRS_67,
        CUSTOM,
    };

    typedef struct GeodParams
    {
        double semi_major_axis;
        double flattening;
        double semi_minor_axis;
        double eccentricity;
        GEOD_TYPE geod_type;

    } GeodParams;

    extern GeodParams geodParameters;
    // These are functions to instantiate newly created GeodParams structs with values
    // this can be done by either providing GEOD::TYPE, or  of semi_major_axis and flattening parameters
    // if geodParams parameter is left empty it modifies and uses the extern global referenced geodParameters struct
    // defined within Geodesic namespace this function automatically updates kObj object's parameters by updateGeodesicProperties()
    extern void configureGeod(GEOD_TYPE geod_type, GeodParams &geodParams = geodParameters);
    extern void configureGeod(double semi_major_axis, double flattening, GeodParams &geodParams = geodParameters);

    class GeodesicMask
    {
    public:
        static constexpr int CAP_NONE = 0;
        static constexpr int CAP_C1 = 1 << 0;
        static constexpr int CAP_C1p = 1 << 1;
        static constexpr int CAP_C2 = 1 << 2;
        static constexpr int CAP_C3 = 1 << 3;
        static constexpr int CAP_C4 = 1 << 4;
        static constexpr int CAP_ALL = 0x1F;
        static constexpr int CAP_MASK = CAP_ALL;
        static constexpr int OUT_ALL = 0x7F80;
        static constexpr int OUT_MASK = 0xFF80; // Include LONG_UNROLL
        /**
         * No capabilities, no output.
         **********************************************************************/
        static constexpr int NONE = 0;
        /**
         * Calculate latitude <i>lat2</i>.  (It's not necessary to include this as a
         * capability to {@link GeodesicLine} because this is included by default.)
         **********************************************************************/
        static constexpr int LATITUDE = 1 << 7 | CAP_NONE;
        /**
         * Calculate longitude <i>lon2</i>.
         **********************************************************************/
        static constexpr int LONGITUDE = 1 << 8 | CAP_C3;
        /**
         * Calculate azimuths <i>azi1</i> and <i>azi2</i>.  (It's not necessary to
         * include this as a capability to {@link GeodesicLine} because this is
         * included by default.)
         **********************************************************************/
        static constexpr int AZIMUTH = 1 << 9 | CAP_NONE;
        /**
         * Calculate distance <i>s12</i>.
         **********************************************************************/
        static constexpr int DISTANCE = 1 << 10 | CAP_C1;
        /**
         * All of the above, the "standard" output and capabilities.
         **********************************************************************/
        static constexpr int STANDARD = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE;
        /**
         * Allow distance <i>s12</i> to be used as <i>input</i> in the direct
         * geodesic problem.
         **********************************************************************/
        static constexpr int DISTANCE_IN = 1 << 11 | CAP_C1 | CAP_C1p;
        /**
         * Calculate reduced length <i>m12</i>.
         **********************************************************************/
        static constexpr int REDUCEDLENGTH = 1 << 12 | CAP_C1 | CAP_C2;
        /**
         * Calculate geodesic scales <i>M12</i> and <i>M21</i>.
         **********************************************************************/
        static constexpr int GEODESICSCALE = 1 << 13 | CAP_C1 | CAP_C2;
        /**
         * Calculate area <i>S12</i>.
         **********************************************************************/
        static constexpr int AREA = 1 << 14 | CAP_C4;
        /**
         * All capabilities, calculate everything.  (LONG_UNROLL is not included in
         * this mask.)
         **********************************************************************/
        static constexpr int ALL = OUT_ALL | CAP_ALL;
        /**
         * Unroll <i>lon2</i>.
         **********************************************************************/
        static constexpr int LONG_UNROLL = 1 << 15;
    };

    class GeodesicLine;
    class KarneyFF
    {
        friend GeodesicLine;

    protected:
        static double tiny_;
        static double tol0_;
        static double tol1_;
        static double tol2_;
        static double tolb_;
        static double xthresh_;

        static constexpr int GEOGRAPHICLIB_GEODESIC_ORDER = 6;
        static constexpr int nA1_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nC1_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nC1p_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nA2_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nC2_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nA3_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nA3x_ = nA3_;
        static constexpr int nC3_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nC3x_ = (nC3_ * (nC3_ - 1)) / 2;
        static constexpr int nC4_ = GEOGRAPHICLIB_GEODESIC_ORDER;
        static constexpr int nC4x_ = (nC4_ * (nC4_ + 1)) / 2;
        static constexpr int maxit1_ = 20;
        static constexpr int maxit2_ = maxit1_ + GeoMath::digits + 10;

        double _a, _f, _f1, _e2, _ep2, _b, _c2;
        double _n, _etol2;
        std::vector<double> _A3x, _C3x, _C4x;

    private:
        /// @brief Dont worry about it.
        typedef struct Lambda12V
        {
            double lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12;
            Lambda12V()
            {
                lam12 = salp2 = calp2 = sig12 = ssig1 = csig1 = ssig2 = csig2 = eps = domg12 = dlam12 = std::numeric_limits<double>::quiet_NaN();
            }
        } Lambda12V;
        typedef struct LengthsV
        {
            double s12b, m12b, m0, M12, M21;
            LengthsV()
            {
                s12b = m12b = m0 = M12 = M21 = std::numeric_limits<double>::quiet_NaN();
            }
        } LengthsV;
        typedef struct InverseStartV
        {
            double sig12, salp1, calp1, salp2, calp2, dnm;
            InverseStartV()
            {
                sig12 = salp1 = calp1 = salp2 = calp2 = dnm = std::numeric_limits<double>::quiet_NaN();
            }
        } InverseStartV;
        typedef struct InverseData
        {
            Geodesic::GeodesicData g;
            double salp1, calp1, salp2, calp2;
            InverseData()
            {
                g = Geodesic::GeodesicData();
                salp1 = calp1 = salp2 = calp2 = std::numeric_limits<double>::quiet_NaN();
            }
        } InverseData;

        void Lambda12(Lambda12V &w, double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double salp1, double calp1, double slam120, double clam120, bool diffp, std::vector<double> &C1a, std::vector<double> &C2a, std::vector<double> &C3a, GeoMath::Pair<double> &p, LengthsV &v);
        void Lengths(LengthsV &v, double eps, double sig12, double ssig1, double csig1, double dn1, double ssig2, double csig2, double dn2, double cbet1, double cbet2, int outmask, std::vector<double> &C1a, std::vector<double> &C2a);
        InverseData InverseInt(double lat1, double lon1, double lat2, double lon2, int outmask);
        InverseStartV InverseStart(double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double lam12, double slam12, double clam12, std::vector<double> &C1a, std::vector<double> &C2a, GeoMath::Pair<double> &p, LengthsV &v);
        static double Astroid(double x, double y);

    public:
        KarneyFF();
        /// @brief  takes no parameter updates the Geod parameters(_a,_b,_f,_e,.....) according to the current state of GeodParams struct.
        /// automatically called within configure geod to update Kobj
        void updateKarneyFFProperties(KarneyFF &obj = kObj, GeodParams &params = Geodesic::geodParameters);
        Geodesic::GeodesicData Inverse(double lat1, double lon1, double lat2, double lon2, int outmask);
        Geodesic::GeodesicData Inverse(double lat1, double lon1, double lat2, double lon2);
        Geodesic::GeodesicData Direct(double lat1, double lon1, double azi1, double s12);
        Geodesic::GeodesicData Direct(double lat1, double lon1, double azi1, double s12, int outmask);
        Geodesic::GeodesicData ArcDirect(double lat1, double lon1, double azi1, double a12);
        Geodesic::GeodesicData ArcDirect(double lat1, double lon1, double azi1, double a12, int outmask);
        Geodesic::GeodesicData Direct(double lat1, double lon1, double azi1, bool arcmode, double s12_a12, int outmask);
        GeodesicLine InverseLine(double lat1, double lon1, double lat2, double lon2);
        GeodesicLine InverseLine(double lat1, double lon1, double lat2, double lon2, int caps);
        GeodesicLine Line(double lat1, double lon1, double azi1);
        GeodesicLine Line(double lat1, double lon1, double azi1, int caps);
        GeodesicLine DirectLine(double lat1, double lon1, double azi1, double s12);
        GeodesicLine DirectLine(double lat1, double lon1, double azi1, double s12, int caps);
        GeodesicLine ArcDirectLine(double lat1, double lon1, double azi1, double a12);
        GeodesicLine ArcDirectLine(double lat1, double lon1, double azi1, double a12, int caps);
        GeodesicLine GenDirectLine(double lat1, double lon1, double azi1, bool arcmode, double s12_a12, int caps);
        double EllipsoidArea();
        double EquatorialRadius();
        double Flattening();

    protected:
        double A3f(double eps)
        {
            return GeoMath::polyval(nA3_ - 1, _A3x, 0, eps);
        }
        void C3f(double eps, std::vector<double> &c)
        {
            double mult = 1;
            int o = 0;
            for (int l = 1; l < nC3_; ++l)
            {
                int m = nC3_ - l - 1;
                mult *= eps;
                c[l] = mult * GeoMath::polyval(m, _C3x, o, eps);
                o += m + 1;
            }
        }
        void C4f(double eps, std::vector<double> &c)
        {
            double mult = 1;
            int o = 0;
            for (int l = 0; l < nC4_; ++l)
            {                         // l is index of C4[l]
                int m = nC4_ - l - 1; // order of polynomial in eps
                c[l] = mult * GeoMath::polyval(m, _C4x, o, eps);
                o += m + 1;
                mult *= eps;
            }
        }
        static double A1m1f(double eps)
        {
            const std::vector<double> coeff = {
                // (1-eps)*A1-1, polynomial in eps2 of order 3
                1,
                4,
                64,
                0,
                256,
            };
            int m = nA1_ / 2;
            double t = GeoMath::polyval(m, coeff, 0, GeoMath::sq(eps)) / coeff[m + 1];
            return (t + eps) / (1 - eps);
        }
        static void C1f(double eps, std::vector<double> &c)
        {
            const std::vector<double> coeff = {
                // C1[1]/eps^1, polynomial in eps2 of order 2
                -1,
                6,
                -16,
                32,
                // C1[2]/eps^2, polynomial in eps2 of order 2
                -9,
                64,
                -128,
                2048,
                // C1[3]/eps^3, polynomial in eps2 of order 1
                9,
                -16,
                768,
                // C1[4]/eps^4, polynomial in eps2 of order 1
                3,
                -5,
                512,
                // C1[5]/eps^5, polynomial in eps2 of order 0
                -7,
                1280,
                // C1[6]/eps^6, polynomial in eps2 of order 0
                -7,
                2048,
            };
            double
                eps2 = GeoMath::sq(eps),
                d = eps;
            int o = 0;
            for (int l = 1; l <= nC1_; ++l)
            {                           // l is index of C1p[l]
                int m = (nC1_ - l) / 2; // order of polynomial in eps^2
                c[l] = d * GeoMath::polyval(m, coeff, o, eps2) / coeff[o + m + 1];
                o += m + 2;
                d *= eps;
            }
        }
        static void C1pf(double eps, std::vector<double> &c)
        {
            const std::vector<double> coeff = {
                // C1p[1]/eps^1, polynomial in eps2 of order 2
                205,
                -432,
                768,
                1536,
                // C1p[2]/eps^2, polynomial in eps2 of order 2
                4005,
                -4736,
                3840,
                12288,
                // C1p[3]/eps^3, polynomial in eps2 of order 1
                -225,
                116,
                384,
                // C1p[4]/eps^4, polynomial in eps2 of order 1
                -7173,
                2695,
                7680,
                // C1p[5]/eps^5, polynomial in eps2 of order 0
                3467,
                7680,
                // C1p[6]/eps^6, polynomial in eps2 of order 0
                38081,
                61440,
            };
            double
                eps2 = GeoMath::sq(eps),
                d = eps;
            int o = 0;
            for (int l = 1; l <= nC1p_; ++l)
            {                            // l is index of C1p[l]
                int m = (nC1p_ - l) / 2; // order of polynomial in eps^2
                c[l] = d * GeoMath::polyval(m, coeff, o, eps2) / coeff[o + m + 1];
                o += m + 2;
                d *= eps;
            }
        }
        static double A2m1f(double eps)
        {
            const std::vector<double> coeff = {
                // (eps+1)*A2-1, polynomial in eps2 of order 3
                -11,
                -28,
                -192,
                0,
                256,
            };
            int m = nA2_ / 2;
            double t = GeoMath::polyval(m, coeff, 0, GeoMath::sq(eps)) / coeff[m + 1];
            return (t - eps) / (1 + eps);
        }
        static void C2f(double eps, std::vector<double> &c)
        {
            const std::vector<double> coeff = {
                // C2[1]/eps^1, polynomial in eps2 of order 2
                1,
                2,
                16,
                32,
                // C2[2]/eps^2, polynomial in eps2 of order 2
                35,
                64,
                384,
                2048,
                // C2[3]/eps^3, polynomial in eps2 of order 1
                15,
                80,
                768,
                // C2[4]/eps^4, polynomial in eps2 of order 1
                7,
                35,
                512,
                // C2[5]/eps^5, polynomial in eps2 of order 0
                63,
                1280,
                // C2[6]/eps^6, polynomial in eps2 of order 0
                77,
                2048,
            };
            double
                eps2 = GeoMath::sq(eps),
                d = eps;
            int o = 0;
            for (int l = 1; l <= nC2_; ++l)
            {                           // l is index of C2[l]
                int m = (nC2_ - l) / 2; // order of polynomial in eps^2
                c[l] = d * GeoMath::polyval(m, coeff, o, eps2) / coeff[o + m + 1];
                o += m + 2;
                d *= eps;
            }
        }

        void A3coeff()
        {
            const std::vector<double> coeff = {
                // A3, coeff of eps^5, polynomial in n of order 0
                -3,
                128,
                // A3, coeff of eps^4, polynomial in n of order 1
                -2,
                -3,
                64,
                // A3, coeff of eps^3, polynomial in n of order 2
                -1,
                -3,
                -1,
                16,
                // A3, coeff of eps^2, polynomial in n of order 2
                3,
                -1,
                -2,
                8,
                // A3, coeff of eps^1, polynomial in n of order 1
                1,
                -1,
                2,
                // A3, coeff of eps^0, polynomial in n of order 0
                1,
                1,
            };
            int o = 0, k = 0;
            for (int j = nA3_ - 1; j >= 0; --j)
            {                                              // coeff of eps^j
                int m = std::min<double>(nA3_ - j - 1, j); // order of polynomial in n
                _A3x[k++] = GeoMath::polyval(m, coeff, o, _n) / coeff[o + m + 1];
                o += m + 2;
            }
        }
        void C3coeff()
        {
            const std::vector<double> coeff = {
                // C3[1], coeff of eps^5, polynomial in n of order 0
                3,
                128,
                // C3[1], coeff of eps^4, polynomial in n of order 1
                2,
                5,
                128,
                // C3[1], coeff of eps^3, polynomial in n of order 2
                -1,
                3,
                3,
                64,
                // C3[1], coeff of eps^2, polynomial in n of order 2
                -1,
                0,
                1,
                8,
                // C3[1], coeff of eps^1, polynomial in n of order 1
                -1,
                1,
                4,
                // C3[2], coeff of eps^5, polynomial in n of order 0
                5,
                256,
                // C3[2], coeff of eps^4, polynomial in n of order 1
                1,
                3,
                128,
                // C3[2], coeff of eps^3, polynomial in n of order 2
                -3,
                -2,
                3,
                64,
                // C3[2], coeff of eps^2, polynomial in n of order 2
                1,
                -3,
                2,
                32,
                // C3[3], coeff of eps^5, polynomial in n of order 0
                7,
                512,
                // C3[3], coeff of eps^4, polynomial in n of order 1
                -10,
                9,
                384,
                // C3[3], coeff of eps^3, polynomial in n of order 2
                5,
                -9,
                5,
                192,
                // C3[4], coeff of eps^5, polynomial in n of order 0
                7,
                512,
                // C3[4], coeff of eps^4, polynomial in n of order 1
                -14,
                7,
                512,
                // C3[5], coeff of eps^5, polynomial in n of order 0
                21,
                2560,
            };
            int o = 0, k = 0;
            for (int l = 1; l < nC3_; ++l)
            { // l is index of C3[l]
                for (int j = nC3_ - 1; j >= l; --j)
                {                                              // coeff of eps^j
                    int m = std::min<double>(nC3_ - j - 1, j); // order of polynomial in n
                    _C3x[k++] = GeoMath::polyval(m, coeff, o, _n) / coeff[o + m + 1];
                    o += m + 2;
                }
            }
        }
        void C4coeff()
        {
            const std::vector<double> coeff = {
                // C4[0], coeff of eps^5, polynomial in n of order 0
                97,
                15015,
                // C4[0], coeff of eps^4, polynomial in n of order 1
                1088,
                156,
                45045,
                // C4[0], coeff of eps^3, polynomial in n of order 2
                -224,
                -4784,
                1573,
                45045,
                // C4[0], coeff of eps^2, polynomial in n of order 3
                -10656,
                14144,
                -4576,
                -858,
                45045,
                // C4[0], coeff of eps^1, polynomial in n of order 4
                64,
                624,
                -4576,
                6864,
                -3003,
                15015,
                // C4[0], coeff of eps^0, polynomial in n of order 5
                100,
                208,
                572,
                3432,
                -12012,
                30030,
                45045,
                // C4[1], coeff of eps^5, polynomial in n of order 0
                1,
                9009,
                // C4[1], coeff of eps^4, polynomial in n of order 1
                -2944,
                468,
                135135,
                // C4[1], coeff of eps^3, polynomial in n of order 2
                5792,
                1040,
                -1287,
                135135,
                // C4[1], coeff of eps^2, polynomial in n of order 3
                5952,
                -11648,
                9152,
                -2574,
                135135,
                // C4[1], coeff of eps^1, polynomial in n of order 4
                -64,
                -624,
                4576,
                -6864,
                3003,
                135135,
                // C4[2], coeff of eps^5, polynomial in n of order 0
                8,
                10725,
                // C4[2], coeff of eps^4, polynomial in n of order 1
                1856,
                -936,
                225225,
                // C4[2], coeff of eps^3, polynomial in n of order 2
                -8448,
                4992,
                -1144,
                225225,
                // C4[2], coeff of eps^2, polynomial in n of order 3
                -1440,
                4160,
                -4576,
                1716,
                225225,
                // C4[3], coeff of eps^5, polynomial in n of order 0
                -136,
                63063,
                // C4[3], coeff of eps^4, polynomial in n of order 1
                1024,
                -208,
                105105,
                // C4[3], coeff of eps^3, polynomial in n of order 2
                3584,
                -3328,
                1144,
                315315,
                // C4[4], coeff of eps^5, polynomial in n of order 0
                -128,
                135135,
                // C4[4], coeff of eps^4, polynomial in n of order 1
                -2560,
                832,
                405405,
                // C4[5], coeff of eps^5, polynomial in n of order 0
                128,
                99099,
            };
            int o = 0, k = 0;
            for (int l = 0; l < nC4_; ++l)
            { // l is index of C4[l]
                for (int j = nC4_ - 1; j >= l; --j)
                {                         // coeff of eps^j
                    int m = nC4_ - j - 1; // order of polynomial in n
                    _C4x[k++] = GeoMath::polyval(m, coeff, o, _n) / coeff[o + m + 1];
                    o += m + 2;
                }
            }
        };
    };

    /// @brief a class to represent GeodesicLines and to solve inverse Direct ArcDirect Geodesic and Inverse Problems.
    class GeodesicLine
    {
        friend KarneyFF;
        
    private:
        static const int nC1_;
        static const int nC1p_;
        static const int nC2_;
        static const int nC3_;
        static const int nC4_;

        double _lat1, _lon1, _azi1;
        double _a, _f, _b, _c2, _f1, _salp0, _calp0, _k2,
            _salp1, _calp1, _ssig1, _csig1, _dn1, _stau1, _ctau1, _somg1, _comg1,
            _A1m1, _A2m1, _A3c, _B11, _B21, _B31, _A4, _B41;

        double _a13, _s13;
        std::vector<double> _C1a, _C1pa, _C2a, _C3a, _C4a;
        int _caps;

    public:
        GeodesicLine();
        GeodesicLine(KarneyFF &g, double lat1, double lon1, double azi1);
        GeodesicLine(KarneyFF &g, double lat1, double lon1, double azi1, int caps);

    private:
        void LineInit(KarneyFF &g, double lat1, double lon1, double azi1, double salp1, double calp1, int caps, GeoMath::Pair<double> &p);
    protected:
        GeodesicLine(KarneyFF &g, double lat1, double lon1, double azi1, double salp1, double calp1, int caps, bool arcmode, double s13_a13);

    public:
        Geodesic::GeodesicData Position(double s12);
        Geodesic::GeodesicData Position(double s12, int outmask);
        Geodesic::GeodesicData ArcPosition(double a12);
        Geodesic::GeodesicData ArcPosition(double a12, int outmask);
        Geodesic::GeodesicData Position(bool arcmode, double s12_a12, int outmask);

        void SetDistance(double s13);
        void SetArc(double a13);
        void GetSetDistance(bool arcmode, double s13_a13);
        bool Init();
        double Latitude();
        double Longitude();
        double Azimuth();
        GeoMath::Pair<double> AzimuthCosine();
        double EquatorialAzimuth();
        GeoMath::Pair<double> EquatorialAzimuthCosine();
        double EquatorialRadius();
        double Flattening();
        int Capabilities();
        bool Capabilities(int testcaps);
        double GenDistance(bool arcmode);
        double Distance();
        double Arc();
    };

}
