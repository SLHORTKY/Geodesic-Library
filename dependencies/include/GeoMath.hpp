#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

namespace Geodesic
{

    /// @brief simple unit conversion class
    class Units
    {
    public:
        static double toDeg(double radian = -0.0, double arcmin = -0.0, double arcsec = -0.0);
        static double toRad(double degree = -0.0, double arcmin = -0.0, double arcsec = -0.0);
        static double toArcmin(double degree, double radian, double arcsec);
        static double toArcsec(double degree, double radian, double arcmin);
        static double toKilometers(double meters, double miles, double feet, double nautical);
        static double toMeters(double kilometers, double miles, double feet, double nautical);
        static double toMiles(double kilometers, double meters, double feet, double nautical);
        static double toFeet(double kilometers, double meters, double miles, double nautical);
        static double toNautical(double kilometers, double meters, double miles, double feet);
        static double degToKilometers(double degrees, bool square = 0);
    };

    class GeoMath
    {
    public:
        template <typename T>
        struct Pair
        {
            Pair()
            {
                this->first = 0;
                this->second = 0;
            };
            Pair(T first, T second)
            {
                this->first = first;
                this->second = second;
            };
            T first = 0;
            T second = 0;
        };
        static constexpr int digits = 53;
        static double sq(double x) { return x * x; }
        /**
         * The inverse hyperbolic tangent function.  This is defined in terms of
         * Math.log1p(<i>x</i>) in order to maintain accuracy near <i>x</i> = 0.
         * In addition, the odd parity of the function is enforced.
         * <p>
         * @param x the argument.
         * @return atanh(<i>x</i>).
         **********************************************************************/
        static double IeeeRemainder(double a, double b)
        {
            return std::remainder(a, b);
        }
        static double atanh(double x)
        {
            double y = std::abs(x); // Enforce odd parity
            y = log1p(2 * y / (1 - y)) / 2;
            return x > 0 ? y : (x < 0 ? -y : x);
        }
        /**
         * Normalize a sine cosine pair.
         * <p>
         * @param p return parameter for normalized quantities with sinx<sup>2</sup>
         *   + cosx<sup>2</sup> = 1.
         * @param sinx the sine.
         * @param cosx the cosine.
         **********************************************************************/
        static void norm(Pair<double> &p, double sinx, double cosx)
        {
            double r = std::hypot(sinx, cosx);
            p.first = sinx / r;
            p.second = cosx / r;
        }
        /**
         * The error-free sum of two numbers.
         * <p>
         * @param u the first number in the sum.
         * @param v the second number in the sum.
         * @param p output Pair(<i>s</i>, <i>t</i>) with <i>s</i> = round(<i>u</i> +
         *   <i>v</i>) and <i>t</i> = <i>u</i> + <i>v</i> - <i>s</i>.
         * <p>
         * See D. E. Knuth, TAOCP, Vol 2, 4.2.2, Theorem B.
         **********************************************************************/
        static void sum(Pair<double> &p, double u, double v)
        {
            double s = u + v;
            double up = s - v;
            double vpp = s - up;
            up -= u;
            vpp -= v;
            double t = s != 0 ? 0.0 - (up + vpp) : s;
            // u + v =       s      + t
            //       = round(u + v) + t
            p.first = s;
            p.second = t;
        }
        /**
         * Evaluate a polynomial.
         * <p>
         * @param N the order of the polynomial.
         * @param p the coefficient array (of size <i>N</i> + <i>s</i> + 1 or more).
         * @param s starting index for the array.
         * @param x the variable.
         * @return the value of the polynomial.
         * <p>
         * Evaluate <i>y</i> = &sum;<sub><i>n</i>=0..<i>N</i></sub>
         * <i>p</i><sub><i>s</i>+<i>n</i></sub>
         * <i>x</i><sup><i>N</i>&minus;<i>n</i></sup>.  Return 0 if <i>N</i> &lt; 0.
         * Return <i>p</i><sub><i>s</i></sub>, if <i>N</i> = 0 (even if <i>x</i> is
         * infinite or a nan).  The evaluation uses Horner's method.
         **********************************************************************/
        static double polyval(int N, const std::vector<double> p, int s, double x)
        {
            double y = N < 0 ? 0 : p[s++];
            while (--N >= 0)
                y = y * x + p[s++];
            return y;
        }
        /**
         * Coarsen a value close to zero.
         * <p>
         * @param x the argument
         * @return the coarsened value.
         * <p>
         * This makes the smallest gap in <i>x</i> = 1/16 &minus; nextafter(1/16, 0)
         * = 1/2<sup>57</sup> for reals = 0.7 pm on the earth if <i>x</i> is an angle
         * in degrees.  (This is about 1000 times more resolution than we get with
         * angles around 90 degrees.)  We use this to avoid having to deal with near
         * singular cases when <i>x</i> is non-zero but tiny (e.g.,
         * 10<sup>&minus;200</sup>).  This converts &minus;0 to +0; however tiny
         * negative numbers get converted to &minus;0.
         **********************************************************************/
        static double AngRound(double x)
        {
            constexpr double z = 1 / 16.0;
            double y = std::abs(x);
            // compiler should not equate z- (z - y) to y
            y = y < z ? z - (z - y) : y;
            return std::copysign(y, x);
        }
        /**
         * Normalize an angle.
         * <p>
         * @param x the angle in degrees.
         * @return the angle reduced to the range [&minus;180&deg;, 180&deg;).
         * <p>
         * The range of <i>x</i> is unrestricted.
         **********************************************************************/
        static double AngNormalize(double x)
        {
            double y = IeeeRemainder(x, 360.0);
            return std::abs(y) == 180 ? std::copysign(180.0, x) : y;
        }
        /**
         * Normalize a latitude.
         * <p>
         * @param x the angle in degrees.
         * @return x if it is in the range [&minus;90&deg;, 90&deg;], otherwise
         *   return NaN.
         **********************************************************************/
        static double LatFix(double x)
        {
            return std::abs(x) > 90 ? std::numeric_limits<double>::quiet_NaN() : x;
        }
        /**
         * The exact difference of two angles reduced to [&minus;180&deg;, 180&deg;].
         * <p>
         * @param x the first angle in degrees.
         * @param y the second angle in degrees.
         * @param p output Pair(<i>d</i>, <i>e</i>) with <i>d</i> being the rounded
         *   difference and <i>e</i> being the error.
         * <p>
         * This computes <i>z</i> = <i>y</i> &minus; <i>x</i> exactly, reduced to
         * [&minus;180&deg;, 180&deg;]; and then sets <i>z</i> = <i>d</i> + <i>e</i>
         * where <i>d</i> is the nearest representable number to <i>z</i> and
         * <i>e</i> is the truncation error.  If <i>z</i> = &plusmn;0&deg; or
         * &plusmn;180&deg;, then the sign of <i>d</i> is given by the sign of
         * <i>y</i> &minus; <i>x</i>.  The maximum absolute value of <i>e</i> is
         * 2<sup>&minus;26</sup> (for doubles).
         **********************************************************************/
        static void AngDiff(Pair<double> &p, double x, double y)
        {
            sum(p, IeeeRemainder(-x, 360.0), IeeeRemainder(y, 360.0));
            sum(p, IeeeRemainder(p.first, 360.0), p.second);
            if (p.first == 0 || std::abs(p.first) == 180)
                p.first = std::copysign(p.first, p.second == 0 ? y - x : -p.second);
        }
        /**
         * Evaluate the sine and cosine function with the argument in degrees
         *
         * @param p return Pair(<i>s</i>, <i>t</i>) with <i>s</i> = sin(<i>x</i>) and
         *   <i>c</i> = cos(<i>x</i>).
         * @param x in degrees.
         * <p>
         * The results obey exactly the elementary properties of the trigonometric
         * functions, e.g., sin 9&deg; = cos 81&deg; = &minus; sin 123456789&deg;.
         **********************************************************************/
        static void sincosd(Pair<double> &p, double x)
        {
            double r;
            int q;
            r = IeeeRemainder(x, 360.0);
            q = (int)round(r / 90); // If r is NaN this returns 0
            r -= 90 * q;
            r = Geodesic::Units::toRad(r);
            // Possibly could call the gnu extension sincos
            double s = sin(r), c = cos(r);
            double sinx, cosx;
            switch (q & 3)
            {
            case 0:
                sinx = s;
                cosx = c;
                break;
            case 1:
                sinx = c;
                cosx = -s;
                break;
            case 2:
                sinx = -s;
                cosx = -c;
                break;
            default:
                sinx = -c;
                cosx = s;
                break; // case 3
            }
            if (sinx == 0)
                sinx = std::copysign(sinx, x);
            p.first = sinx;
            p.second = 0.0 + cosx;
        }
        /**
         * Evaluate the sine and cosine function with reduced argument plus correction
         *
         * @param p return Pair(<i>s</i>, <i>t</i>) with <i>s</i> =
         *   sin(<i>x</i> +  <i>t</i>) and <i>c</i> = cos(<i>x</i> + <i>t</i>).
         * @param x reduced angle in degrees.
         * @param t correction in degrees.
         * <p>
         * This is a variant of GeoMath.sincosd allowing a correction to the angle to
         * be supplied.  <i>x</i> x must be in [&minus;180&deg;, 180&deg;] and
         * <i>t</i> is assumed to be a <i>small</i> correction.  GeoMath.AngRound is
         * applied to the reduced angle to prevent problems with <i>x</i> + <i>t</i>
         * being extremely close but not exactly equal to one of the four cardinal
         * directions.
         **********************************************************************/
        static void sincosde(Pair<double> &p, double x, double t)
        {
            double r;
            int q;
            q = (int)round(x / 90); // If r is NaN this returns 0
            r = x - 90 * q;
            r = Geodesic::Units::toRad(GeoMath::AngRound(r + t));
            double s = sin(r), c = cos(r);
            double sinx, cosx;
            switch (q & 3)
            {
            case 0:
                sinx = s;
                cosx = c;
                break;
            case 1:
                sinx = c;
                cosx = -s;
                break;
            case 2:
                sinx = -s;
                cosx = -c;
                break;
            default:
                sinx = -c;
                cosx = s;
                break; // case 3
            }
            if (sinx == 0)
                sinx = std::copysign(sinx, x);
            p.first = sinx;
            p.second = 0.0 + cosx;
        }
        /**
         * Evaluate the atan2 function with the result in degrees
         *
         * @param y the sine of the angle
         * @param x the cosine of the angle
         * @return atan2(<i>y</i>, <i>x</i>) in degrees.
         * <p>
         * The result is in the range [&minus;180&deg; 180&deg;].  N.B.,
         * atan2d(&plusmn;0, &minus;1) = &plusmn;180&deg;.
         **********************************************************************/
        static double atan2d(double y, double x)
        {
            int q = 0;
            if (std::abs(y) > std::abs(x))
            {
                double t;
                t = x;
                x = y;
                y = t;
                q = 2;
            }
            if (x < 0)
            {
                x = -x;
                ++q;
            }

            double ang = Geodesic::Units::toDeg(std::atan2(y, x));
            switch (q)
            {
            case 1:
                ang = std::copysign(180.0, y) - ang;
                break;
            case 2:
                ang = 90 - ang;
                break;
            case 3:
                ang = -90 + ang;
                break;
            default:
                break;
            }
            return ang;
        }
        static double ulp(double x)
        {
            double next = std::nextafter(x, std::numeric_limits<double>::max());

            return next - x;
        }
        static double SinCosSeries(bool sinp, double sinx, double cosx, std::vector<double> c)
        {

            int k = c.size(), n = k - (sinp ? 1 : 0);
            double ar = 2 * (cosx - sinx) * (cosx + sinx),
                   y0 = (n & 1) != 0 ? c[--k] : 0, y1 = 0;
            n /= 2;

            while (n-- != 0)
            {

                y1 = ar * y0 - y1 + c[--k];
                y0 = ar * y1 - y0 + c[--k];
            }
            return sinp
                       ? 2 * sinx * cosx * y0 // sin(2 * x) * y0
                       : cosx * (y0 - y1);    // cos(x) * (y0 - y1)
        }

        static constexpr double EPSLN = 1.0e-10;
        static constexpr double HALF_PI = M_PI / 2.0;

        static constexpr double SPI = 3.14159265359;
        static constexpr double TWO_PI = M_PI * 2;

        static double sinH(double x)
        {
            double r = exp(x);
            r = (r - 1 / r) / 2;
            return r;
        }
        static double Hypot(double x, double y)
        {
            x = std::abs(x);
            y = std::abs(y);
            double a = std::max<double>(x, y);
            double b = std::min<double>(x, y) / (a > 0 ? a : 1);

            return a * sqrt(1 + pow(b, 2));
        }
        static double Log1py(double x)
        {
            double y = 1 + x;
            double z = y - 1;

            return z == 0 ? x : x * log(y) / z;
        }
        static double asinHy(double x)
        {
            double y = std::abs(x);
            y = Log1py(y * (1 + y / (Hypot(1, y) + 1)));

            return x < 0 ? -y : y;
        }
        static double gatg(double *pp, int length, double B)
        {
            double cos_2B = 2 * cos(2 * B);
            int i = length - 1;
            double h1 = pp[i];
            double h2 = 0;
            double h = 0;

            while (--i >= 0)
            {
                h = -h2 + cos_2B * h1 + pp[i];
                h2 = h1;
                h1 = h;
            }

            return (B + h * sin(2 * B));
        }
        static double clens(double *pp, int length, double arg_r)
        {
            double r = 2 * cos(arg_r);
            int i = length - 1;
            double hr1 = pp[i];
            double hr2 = 0;
            double hr = 1;

            while (--i >= 0)
            {
                hr = -hr2 + r * hr1 + pp[i];
                hr2 = hr1;
                hr1 = hr;
            }

            return sin(arg_r) * hr;
        }
        static double cosH(double x)
        {
            double r = exp(x);
            r = (r + 1 / r) / 2;
            return r;
        }
        static double *clens_cmplx(double *pp, int length, double arg_r, double arg_i)
        {
            double sin_arg_r = sin(arg_r);
            double cos_arg_r = cos(arg_r);
            double sinh_arg_i = sinH(arg_i);
            double cosh_arg_i = cosH(arg_i);
            double r = 2 * cos_arg_r * cosh_arg_i;
            double i = -2 * sin_arg_r * sinh_arg_i;
            int j = length - 1;
            double hr = pp[j];
            double hi1 = 0;
            double hr1 = 0;
            double hi = 0;
            double hr2;
            double hi2;

            while (--j >= 0)
            {
                hr2 = hr1;
                hi2 = hi1;
                hr1 = hr;
                hi1 = hi;
                hr = -hr2 + r * hr1 - i * hi1 + pp[j];
                hi = -hi2 + i * hr1 + r * hi1;
            }

            r = sin_arg_r * cosh_arg_i;
            i = cos_arg_r * sinh_arg_i;

            return new double[2]{r * hr - i * hi, r * hi + i * hr};
        }
        static double *pj_enfn(double es)
        {
            double C00 = 1.0;
            double C02 = 0.25;
            double C04 = 0.046875;
            double C06 = 0.01953125;
            double C08 = 0.01068115234375;
            double C22 = 0.75;
            double C44 = 0.46875;
            double C46 = 0.01302083333333333333;
            double C48 = 0.00712076822916666666;
            double C66 = 0.36458333333333333333;
            double C68 = 0.00569661458333333333;
            double C88 = 0.3076171875;

            double *en = new double[5];
            en[0] = C00 - es * (C02 + es * (C04 + es * (C06 + es * C08)));
            en[1] = es * (C22 - es * (C04 + es * (C06 + es * C08)));
            double t = es * es;
            en[2] = t * (C44 - es * (C46 + es * C48));
            t *= es;
            en[3] = t * (C66 - es * C68);
            en[4] = t * es * C88;
            return en;
        }
        static double pj_mlfn(double phi, double sphi, double cphi, double *en)
        {
            cphi *= sphi;
            sphi *= sphi;
            return (en[0] * phi - cphi * (en[1] + sphi * (en[2] + sphi * (en[3] + sphi * en[4]))));
        }

        static double pj_inv_mlfn(double arg, double es, double *en)
        {
            double MAX_ITER = 20;
            double k = 1 / (1 - es);
            double phi = arg;

            for (int i = MAX_ITER; i > 0; --i)
            {
                /* rarely goes over 2 iterations */
                double s = sin(phi);
                double t = 1 - es * s * s;

                t = (pj_mlfn(phi, s, cos(phi), en) - arg) * (t * sqrt(t)) * k;
                phi -= t;
                if (std::abs(t) < EPSLN)
                {
                    return phi;
                }
            }
            return phi;
        }
        static double sign(double value)
        {
            return value < 0 ? -1 : 1;
        }
        static double adjust_lon(double lon)
        {
            return (std::abs(lon) <= SPI) ? lon : (lon - (sign(lon) * TWO_PI));
        }
    };

}