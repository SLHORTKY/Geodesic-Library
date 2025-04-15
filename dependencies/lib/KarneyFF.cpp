#include "KarneyFF.hpp"
#include <limits>

#pragma region KarneyFF
double Geodesic::KarneyFF::tol0_ = GeoMath::ulp(1.0);
double Geodesic::KarneyFF::tol1_ = 200 * tol0_;
double Geodesic::KarneyFF::tol2_ = std::sqrt(tol0_);
double Geodesic::KarneyFF::tolb_ = tol0_;
double Geodesic::KarneyFF::xthresh_ = 1000 * tol2_;
double Geodesic::KarneyFF::tiny_ = std::sqrt(std::numeric_limits<double>::min());

Geodesic::KarneyFF Geodesic::kObj;

Geodesic::KarneyFF::KarneyFF()
{
    Geodesic::configureGeod(Geodesic::GEOD_TYPE::WGS_84);
    this->updateKarneyFFProperties();

    if (!(std::isfinite(_a) && _a > 0))
    {
        std::cout << "_a is +infinity";
    }
    if (!(std::isfinite(_b) && _b > 0))
    {
        std::cout << "_b is +infinity";
    }
}
void Geodesic::KarneyFF::updateKarneyFFProperties(KarneyFF &obj, GeodParams &params)
{
    obj._a = params.semi_major_axis;
    obj._f = params.flattening;
    obj._f1 = 1 - obj._f;
    obj._e2 = obj._f * (2 - obj._f);
    obj._ep2 = obj._e2 / GeoMath::sq(obj._f1); // e2 / (1 - e2)
    obj._n = obj._f / (2 - obj._f);
    obj._b = obj._a * obj._f1;
    obj._c2 = (GeoMath::sq(obj._a) + GeoMath::sq(obj._b) * (_e2 == 0 ? 1 : (obj._e2 > 0 ? GeoMath::atanh(std::sqrt(obj._e2)) : atan(std::sqrt(-obj._e2))) / std::sqrt(std::abs(obj._e2)))) / 2;
    obj._etol2 = 0.1 * tol2_ / std::sqrt(std::max<double>(0.001, std::abs(obj._f)) * std::min<double>(1.0, 1 - obj._f / 2) / 2);

    obj._A3x = std::vector<double>(nA3x_, 0);
    obj._C3x = std::vector<double>(nC3x_, 0);
    obj._C4x = std::vector<double>(nC4x_, 0);

    A3coeff();
    C3coeff();
    C4coeff();
}
void Geodesic::KarneyFF::Lambda12(Lambda12V &w, double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double salp1, double calp1, double slam120, double clam120, bool diffp, std::vector<double> &C1a, std::vector<double> &C2a, std::vector<double> &C3a, GeoMath::Pair<double> &p, LengthsV &v)
{
    if (sbet1 == 0 && calp1 == 0)
    {
        calp1 = -tiny_;
    }

    double // sin(alp1) * cos(bet1) = sin(alp0)
        salp0 = salp1 * cbet1,
        calp0 = std::hypot(calp1, salp1 * sbet1); // calp0 > 0

    double somg1, comg1, somg2, comg2, somg12, comg12;
    // tan(bet1) = tan(sig1) * cos(alp1)
    // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
    w.ssig1 = sbet1;
    somg1 = salp0 * sbet1;
    w.csig1 = comg1 = calp1 * cbet1;
    GeoMath::norm(p, w.ssig1, w.csig1);
    w.ssig1 = p.first;
    w.csig1 = p.second;
    // GeoMath.norm(somg1, comg1); -- don't need to normalize!

    // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
    // about this case, since this can yield singularities in the Newton
    // iteration.
    // sin(alp2) * cos(bet2) = sin(alp0)
    w.salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
    // calp2 = sqrt(1 - sq(salp2))
    //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
    // and subst for calp0 and rearrange to give (choose positive sqrt
    // to give alp2 in [0, pi/2]).
    w.calp2 = cbet2 != cbet1 || std::abs(sbet2) != -sbet1
                  ? sqrt(GeoMath::sq(calp1 * cbet1) + (cbet1 < -sbet1
                                                           ? (cbet2 - cbet1) * (cbet1 + cbet2)
                                                           : (sbet1 - sbet2) * (sbet1 + sbet2))) /
                        cbet2
                  : std::abs(calp1);
    // tan(bet2) = tan(sig2) * cos(alp2)
    // tan(omg2) = sin(alp0) * tan(sig2).
    w.ssig2 = sbet2;
    somg2 = salp0 * sbet2;
    w.csig2 = comg2 = w.calp2 * cbet2;
    GeoMath::norm(p, w.ssig2, w.csig2);
    w.ssig2 = p.first;
    w.csig2 = p.second;
    // GeoMath.norm(somg2, comg2); -- don't need to normalize!

    // sig12 = sig2 - sig1, limit to [0, pi]
    w.sig12 = atan2(std::max<double>(0.0, w.csig1 * w.ssig2 - w.ssig1 * w.csig2), w.csig1 * w.csig2 + w.ssig1 * w.ssig2);

    // omg12 = omg2 - omg1, limit to [0, pi]
    somg12 = std::max<double>(0.0, comg1 * somg2 - somg1 * comg2);
    comg12 = comg1 * comg2 + somg1 * somg2;
    // eta = omg12 - lam120
    double eta = std::atan2(somg12 * clam120 - comg12 * slam120, comg12 * clam120 + somg12 * slam120);

    double B312;
    double k2 = GeoMath::sq(calp0) * _ep2;
    w.eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
    C3f(w.eps, C3a);
    B312 = (GeoMath::SinCosSeries(true, w.ssig2, w.csig2, C3a) - GeoMath::SinCosSeries(true, w.ssig1, w.csig1, C3a));

    w.domg12 = -_f * A3f(w.eps) * salp0 * (w.sig12 + B312);
    w.lam12 = eta + w.domg12;

    if (diffp)
    {
        if (w.calp2 == 0)
        {
            w.dlam12 = -2 * _f1 * dn1 / sbet1;
        }
        else
        {
            Lengths(v, w.eps,
                    w.sig12, w.ssig1, w.csig1, dn1, w.ssig2, w.csig2, dn2,
                    cbet1, cbet2, GeodesicMask::REDUCEDLENGTH, C1a, C2a);
            w.dlam12 = v.m12b;
            w.dlam12 *= _f1 / (w.calp2 * cbet2);
        }
    }
}
void Geodesic::KarneyFF::Lengths(LengthsV &v, double eps, double sig12, double ssig1, double csig1, double dn1, double ssig2, double csig2, double dn2, double cbet1, double cbet2, int outmask, std::vector<double> &C1a, std::vector<double> &C2a)
{
    // Return m12b = (reduced length)/_b; also calculate s12b = distance/_b,
    // and m0 = coefficient of secular term in expression for reduced length.
    outmask &= GeodesicMask::OUT_MASK;
    // v to hold compute s12b, m12b, m0, M12, M21;

    double m0x = 0, J12 = 0, A1 = 0, A2 = 0;
    if ((outmask & (GeodesicMask::DISTANCE | GeodesicMask::REDUCEDLENGTH | GeodesicMask::GEODESICSCALE)) != 0)
    {
        A1 = A1m1f(eps);
        C1f(eps, C1a);
        if ((outmask & (GeodesicMask::REDUCEDLENGTH | GeodesicMask::GEODESICSCALE)) != 0)
        {
            A2 = A2m1f(eps);
            C2f(eps, C2a);
            m0x = A1 - A2;
            A2 = 1 + A2;
        }
        A1 = 1 + A1;
    }
    if ((outmask & GeodesicMask::DISTANCE) != 0)
    {
        double B1 = GeoMath::SinCosSeries(true, ssig2, csig2, C1a) - GeoMath::SinCosSeries(true, ssig1, csig1, C1a);
        // Missing a factor of _b
        v.s12b = A1 * (sig12 + B1);
        if ((outmask & (GeodesicMask::REDUCEDLENGTH | GeodesicMask::GEODESICSCALE)) != 0)
        {
            double B2 = GeoMath::SinCosSeries(true, ssig2, csig2, C2a) - GeoMath::SinCosSeries(true, ssig1, csig1, C2a);
            J12 = m0x * sig12 + (A1 * B1 - A2 * B2);
        }
    }
    else if ((outmask & (GeodesicMask::REDUCEDLENGTH | GeodesicMask::GEODESICSCALE)) != 0)
    {
        // Assume here that nC1_ >= nC2_
        for (int l = 1; l <= nC2_; ++l)
        {
            C2a[l] = A1 * C1a[l] - A2 * C2a[l];
        }
        J12 = m0x * sig12 + (GeoMath::SinCosSeries(true, ssig2, csig2, C2a) - GeoMath::SinCosSeries(true, ssig1, csig1, C2a));
    }
    if ((outmask & GeodesicMask::REDUCEDLENGTH) != 0)
    {
        v.m0 = m0x;
        v.m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;
    }
    if ((outmask & GeodesicMask::GEODESICSCALE) != 0)
    {
        double csig12 = csig1 * csig2 + ssig1 * ssig2;
        double t = _ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
        v.M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
        v.M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
    }
}
double Geodesic::KarneyFF::Astroid(double x, double y)
{
    double k;
    double p = GeoMath::sq(x),
           q = GeoMath::sq(y),
           r = (p + q - 1) / 6;
    if (!(q == 0 && r <= 0))
    {
        double             // Avoid possible division by zero when r = 0 by multiplying equations
                           // for s and t by r^3 and r, resp.
            S = p * q / 4, // S = r^3 * s
            r2 = GeoMath::sq(r),
            r3 = r * r2,
            // The discriminant of the quadratic equation for T3.  This is zero on
            // the evolute curve p^(1/3)+q^(1/3) = 1
            disc = S * (S + 2 * r3);
        double u = r;
        if (disc >= 0)
        {
            double T3 = S + r3;
            // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
            // of precision due to cancellation.  The result is unchanged because
            // of the way the T is used in definition of u.
            T3 += T3 < 0 ? -std::sqrt(disc) : std::sqrt(disc); // T3 = (r * t)^3
            // N.B. cbrt always returns the double root.  cbrt(-8) = -2.
            double T = std::cbrt(T3); // T = r * t
            // T can be zero; but then r2 / T -> 0.
            u += T + (T != 0 ? r2 / T : 0);
        }
        else
        {
            // T is complex, but the way u is defined the result is double.
            double ang = std::atan2(std::sqrt(-disc), -(S + r3));
            // There are three possible cube roots.  We choose the root which
            // avoids cancellation.  Note that disc < 0 implies that r < 0.
            u += 2 * r * std::cos(ang / 3);
        }
        double v = std::sqrt(GeoMath::sq(u) + q), // guaranteed positive
                                                  // Avoid loss of accuracy when u < 0.
            uv = u < 0 ? q / (v - u) : u + v,     // u+v, guaranteed positive
            w = (uv - q) / (2 * v);               // positive?
        // Rearrange expression for k to avoid loss of accuracy due to
        // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
        k = uv / (std::sqrt(uv + GeoMath::sq(w)) + w); // guaranteed positive
    }
    else
    { // q == 0 && r <= 0
        // y = 0 with |x| <= 1.  Handle this case directly.
        // for y small, positive root is k = abs(y)/sqrt(1-x^2)
        k = 0;
    }
    return k;
}
Geodesic::KarneyFF::InverseData Geodesic::KarneyFF::InverseInt(double lat1, double lon1, double lat2, double lon2, int outmask)
{
    InverseData result = InverseData();
    GeoMath::Pair<double> p = GeoMath::Pair<double>();
    GeodesicData r = GeodesicData();
    // Compute longitude difference (AngDiff does this carefully).  Result is
    // in [-180, 180] but -180 is only for west-going geodesics.  180 is for
    // east-going and meridional geodesics.
    r.lat1 = lat1 = GeoMath::LatFix(lat1);
    r.lat2 = lat2 = GeoMath::LatFix(lat2);
    // If really close to the equator, treat as on equator.
    lat1 = GeoMath::AngRound(lat1);
    lat2 = GeoMath::AngRound(lat2);
    double lon12, lon12s;
    GeoMath::AngDiff(p, lon1, lon2);
    lon12 = p.first;
    lon12s = p.second;
    if ((outmask & GeodesicMask::LONG_UNROLL) != 0)
    {
        r.lon1 = lon1;
        r.lon2 = (lon1 + lon12) + lon12s;
    }
    else
    {
        r.lon1 = GeoMath::AngNormalize(lon1);
        r.lon2 = GeoMath::AngNormalize(lon2);
    }
    // Make longitude difference positive.
    int lonsign = (int)copysign(1.0, lon12);
    lon12 *= lonsign;
    lon12s *= lonsign;
    double lam12 = Geodesic::Units::toRad(lon12), slam12, clam12;
    // Calculate sincos of lon12 + error (this applies AngRound internally).
    GeoMath::sincosde(p, lon12, lon12s);
    slam12 = p.first;
    clam12 = p.second;
    lon12s = (180 - lon12) - lon12s; // the supplementary longitude difference

    // Swap points so that point with higher (std::abs) latitude is point 1
    // If one latitude is a nan, then it becomes lat1.
    int swapp = std::abs(lat1) < std::abs(lat2) || lat2 != lat2 ? -1 : 1;
    if (swapp < 0)
    {
        lonsign *= -1;
        {
            double t = lat1;
            lat1 = lat2;
            lat2 = t;
        }
    }
    // Make lat1 <= 0
    int latsign = (int)std::copysign(1.0, -lat1);
    lat1 *= latsign;
    lat2 *= latsign;

    double sbet1, cbet1, sbet2, cbet2, s12x, m12x;
    s12x = m12x = std::numeric_limits<double>::quiet_NaN();

    GeoMath::sincosd(p, lat1);
    sbet1 = _f1 * p.first;
    cbet1 = p.second;

    GeoMath::norm(p, sbet1, cbet1);
    sbet1 = p.first;
    cbet1 = p.second;
    cbet1 = std::max<double>(tiny_, cbet1);

    GeoMath::sincosd(p, lat2);
    sbet2 = _f1 * p.first;
    cbet2 = p.second;
    // Ensure cbet2 = +epsilon at poles
    GeoMath::norm(p, sbet2, cbet2);
    sbet2 = p.first;
    cbet2 = p.second;
    cbet2 = std::max<double>(tiny_, cbet2);

    if (cbet1 < -sbet1)
    {
        if (cbet2 == cbet1)
        {
            sbet2 = std::copysign(sbet1, sbet2);
        }
    }
    else
    {
        if (std::abs(sbet2) == -sbet1)
        {
            cbet2 = cbet1;
        }
    }

    double dn1 = std::sqrt(1 + _ep2 * GeoMath::sq(sbet1)),
           dn2 = std::sqrt(1 + _ep2 * GeoMath::sq(sbet2));

    double a12, sig12, calp1, salp1, calp2, salp2;
    a12 = sig12 = calp1 = salp1 = calp2 = salp2 = std::numeric_limits<double>::quiet_NaN();
    // index zero elements of these arrays are unused
    std::vector<double> C1a(nC1_ + 1, 0);
    std::vector<double> C2a(nC2_ + 1, 0);
    std::vector<double> C3a(nC3_, 0);

    bool meridian = lat1 == -90 || slam12 == 0;
    LengthsV v = LengthsV();

    if (meridian)
    {

        // Endpoints are on a single full meridian, so the geodesic might lie on
        // a meridian.
        calp1 = clam12;
        salp1 = slam12; // Head to the target longitude
        calp2 = 1;
        salp2 = 0; // At the target we're heading north

        double // tan(bet) = tan(sig) * cos(alp)
            ssig1 = sbet1,
            csig1 = calp1 * cbet1,
            ssig2 = sbet2, csig2 = calp2 * cbet2;

        // sig12 = sig2 - sig1 (N.B. with Java max(+0.0, -0.0) -> +0.0)
        sig12 = std::atan2(std::max<double>(0.0, csig1 * ssig2 - ssig1 * csig2),
                           csig1 * csig2 + ssig1 * ssig2);
        Lengths(v, _n, sig12, ssig1, csig1, dn1,
                ssig2, csig2, dn2, cbet1, cbet2,
                outmask | GeodesicMask::DISTANCE | GeodesicMask::REDUCEDLENGTH,
                C1a, C2a);
        s12x = v.s12b;
        m12x = v.m12b;
        if ((outmask & GeodesicMask::GEODESICSCALE) != 0)
        {
            r.M12 = v.M12;
            r.M21 = v.M21;
        }

        if (sig12 < 1 || m12x >= 0)
        {
            // Need at least 2, to handle 90 0 90 180
            if (sig12 < 3 * tiny_ || // Prevent negative s12 or m12 for short lines
                (sig12 < tol0_ && (s12x < 0 || m12x < 0)))
            {
                sig12 = m12x = s12x = 0;
            }
            m12x *= _b;
            s12x *= _b;
            a12 = Geodesic::Units::toDeg(sig12);
        }
        else // m12 < 0, i.e., prolate and too close to anti-podal
        {
            meridian = false;
        }
    }

    double omg12 = std::numeric_limits<double>::quiet_NaN(), somg12 = 2, comg12 = std::numeric_limits<double>::quiet_NaN();
    if (!meridian && sbet1 == 0 && // and sbet2 == 0
                                   // Mimic the way Lambda12 works with calp1 = 0
        (_f <= 0 || lon12s >= _f * 180))
    {
        calp1 = calp2 = 0;
        salp1 = salp2 = 1;
        s12x = _a * lam12;
        sig12 = omg12 = lam12 / _f1;
        m12x = _b * sin(sig12);
        if ((outmask & GeodesicMask::GEODESICSCALE) != 0)
        {
            r.M12 = r.M21 = cos(sig12);
        }
        a12 = lon12 / _f1;
    }
    else if (!meridian)
    {
        double dnm;
        {
            InverseStartV s = InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                                           lam12, slam12, clam12,
                                           C1a, C2a, p, v);
            sig12 = s.sig12;
            salp1 = s.salp1;
            calp1 = s.calp1;
            salp2 = s.salp2;
            calp2 = s.calp2;
            dnm = s.dnm;
        }

        if (sig12 >= 0)
        {
            // Short lines (InverseStart sets salp2, calp2, dnm)
            s12x = sig12 * _b * dnm;
            m12x = GeoMath::sq(dnm) * _b * sin(sig12 / dnm);
            if ((outmask & GeodesicMask::GEODESICSCALE) != 0)
            {
                r.M12 = r.M21 = cos(sig12 / dnm);
            }
            a12 = Geodesic::Units::toDeg(sig12);
            omg12 = lam12 / (_f1 * dnm);
        }
        else
        {

            double ssig1, csig1, ssig2, csig2, eps, domg12;
            ssig1 = csig1 = ssig2 = csig2 = eps = domg12 = std::numeric_limits<double>::quiet_NaN();
            int numit = 0;
            // Bracketing range
            double salp1a = tiny_, calp1a = 1, salp1b = tiny_, calp1b = -1;
            Lambda12V w = Lambda12V();
            for (bool tripn = false, tripb = false;; ++numit)
            {

                double V, dV;

                Lambda12(w, sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                         slam12, clam12, numit < maxit1_, C1a, C2a, C3a, p, v);
                V = w.lam12;
                salp2 = w.salp2;
                calp2 = w.calp2;
                sig12 = w.sig12;
                ssig1 = w.ssig1;
                csig1 = w.csig1;
                ssig2 = w.ssig2;
                csig2 = w.csig2;
                eps = w.eps;
                domg12 = w.domg12;
                dV = w.dlam12;
                if (tripb ||                                     // Reversed test to allow escape with NaNs
                    !(std::abs(V) >= (tripn ? 8 : 1) * tol0_) || // Enough bisections to get accurate result
                    numit == maxit2_)
                {
                    break;
                }
                // Update bracketing values
                if (V > 0 && (numit > maxit1_ || calp1 / salp1 > calp1b / salp1b))
                {
                    salp1b = salp1;
                    calp1b = calp1;
                }
                else if (V < 0 && (numit > maxit1_ || calp1 / salp1 < calp1a / salp1a))
                {
                    salp1a = salp1;
                    calp1a = calp1;
                }
                if (numit < maxit1_ && dV > 0)
                {
                    double dalp1 = -V / dV;
                    if (std::abs(dalp1) < M_PI)
                    {
                        double sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
                               nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
                        if (nsalp1 > 0)
                        {
                            calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                            salp1 = nsalp1;
                            GeoMath::norm(p, salp1, calp1);
                            salp1 = p.first;
                            calp1 = p.second;
                            // In some regimes we don't get quadratic convergence because
                            // slope -> 0.  So use convergence conditions based on epsilon
                            // instead of sqrt(epsilon).
                            tripn = std::abs(V) <= 16 * tol0_;
                            continue;
                        }
                    }
                }

                salp1 = (salp1a + salp1b) / 2;
                calp1 = (calp1a + calp1b) / 2;
                GeoMath::norm(p, salp1, calp1);
                salp1 = p.first;
                calp1 = p.second;
                tripn = false;
                tripb = (std::abs(salp1a - salp1) + (calp1a - calp1) < tolb_ || std::abs(salp1 - salp1b) + (calp1 - calp1b) < tolb_);
            }
            {
                // Ensure that the reduced length and geodesic scale are computed in
                // a "canonical" way, with the I2 integral.
                int lengthmask = outmask | ((outmask & (GeodesicMask::REDUCEDLENGTH | GeodesicMask::GEODESICSCALE)) != 0
                                                ? GeodesicMask::DISTANCE
                                                : GeodesicMask::NONE);
                Lengths(v, eps, sig12,
                        ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
                        lengthmask, C1a, C2a);
                s12x = v.s12b;
                m12x = v.m12b;
                if ((outmask & GeodesicMask::GEODESICSCALE) != 0)
                {
                    r.M12 = v.M12;
                    r.M21 = v.M21;
                }
            }
            m12x *= _b;
            s12x *= _b;
            a12 = Geodesic::Units::toDeg(sig12);
            if ((outmask & GeodesicMask::AREA) != 0)
            {
                // omg12 = lam12 - domg12
                double sdomg12 = sin(domg12), cdomg12 = cos(domg12);
                somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                comg12 = clam12 * cdomg12 + slam12 * sdomg12;
            }
        }
    }

    if ((outmask & GeodesicMask::DISTANCE) != 0)
    {
        r.s12 = 0 + s12x; // Convert -0 to 0
    }
    if ((outmask & GeodesicMask::REDUCEDLENGTH) != 0)
    {
        r.m12 = 0 + m12x; // Convert -0 to 0
    }
    if ((outmask & GeodesicMask::AREA) != 0)
    {
        double // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
            salp0 = salp1 * cbet1,
            calp0 = std::hypot(calp1, salp1 * sbet1); // calp0 > 0
        double alp12;
        if (calp0 != 0 && salp0 != 0)
        {
            double // From Lambda12: tan(bet) = tan(sig) * cos(alp)
                ssig1 = sbet1,
                csig1 = calp1 * cbet1,
                ssig2 = sbet2, csig2 = calp2 * cbet2,
                k2 = GeoMath::sq(calp0) * _ep2,
                eps = k2 / (2 * (1 + std::sqrt(1 + k2)) + k2),
                // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
                A4 = GeoMath::sq(_a) * calp0 * salp0 * _e2;
            GeoMath::norm(p, ssig1, csig1);
            ssig1 = p.first;
            csig1 = p.second;
            GeoMath::norm(p, ssig2, csig2);
            ssig2 = p.first;
            csig2 = p.second;
            std::vector<double> C4a(nC4_);
            C4f(eps, C4a);
            double B41 = GeoMath::SinCosSeries(false, ssig1, csig1, C4a),
                   B42 = GeoMath::SinCosSeries(false, ssig2, csig2, C4a);
            r.S12 = A4 * (B42 - B41);
        }
        else // Avoid problems with indeterminate sig1, sig2 on equator
        {
            r.S12 = 0;
        }

        if (!meridian && somg12 == 2)
        {
            somg12 = sin(omg12);
            comg12 = cos(omg12);
        }

        if (!meridian && comg12 > -0.7071 && // Long difference not too big
            sbet2 - sbet1 < 1.75)
        { // Lat difference not too big
            // Use tan(Gamma/2) = tan(omg12/2)
            // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
            // with tan(x/2) = sin(x)/(1+cos(x))
            double domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
            alp12 = 2 * std::atan2(somg12 * (sbet1 * dbet2 + sbet2 * dbet1),
                                   domg12 * (sbet1 * sbet2 + dbet1 * dbet2));
        }
        else
        {
            // alp12 = alp2 - alp1, used in atan2 so no need to normalize
            double salp12 = salp2 * calp1 - calp2 * salp1,
                   calp12 = calp2 * calp1 + salp2 * salp1;
            // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
            // salp12 = -0 and alp12 = -180.  However this depends on the sign
            // being attached to 0 correctly.  The following ensures the correct
            // behavior.
            if (salp12 == 0 && calp12 < 0)
            {
                salp12 = tiny_ * calp1;
                calp12 = -1;
            }
            alp12 = std::atan2(salp12, calp12);
        }
        r.S12 += _c2 * alp12;
        r.S12 *= swapp * lonsign * latsign;
        // Convert -0 to 0
        r.S12 += 0;
    }

    // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
    if (swapp < 0)
    {
        {
            double t = salp1;
            salp1 = salp2;
            salp2 = t;
        }
        {
            double t = calp1;
            calp1 = calp2;
            calp2 = t;
        }
        if ((outmask & GeodesicMask::GEODESICSCALE) != 0)
        {
            double t = r.M12;
            r.M12 = r.M21;
            r.M21 = t;
        }
    }
    salp1 *= swapp * lonsign;
    calp1 *= swapp * latsign;
    salp2 *= swapp * lonsign;
    calp2 *= swapp * latsign;

    r.a12 = a12;
    result.g = r;
    result.salp1 = salp1;
    result.calp1 = calp1;
    result.salp2 = salp2;
    result.calp2 = calp2;
    return result;
}
Geodesic::KarneyFF::InverseStartV Geodesic::KarneyFF::InverseStart(double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double lam12, double slam12, double clam12, std::vector<double> &C1a, std::vector<double> &C2a, GeoMath::Pair<double> &p, LengthsV &v)
{
    InverseStartV w = InverseStartV();
    w.sig12 = -1; // Return value
    double        // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
        sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
        cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
    double sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
    bool shortline = cbet12 >= 0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
    double somg12, comg12;
    if (shortline)
    {
        double sbetm2 = GeoMath::sq(sbet1 + sbet2);
        // sin((bet1+bet2)/2)^2
        // =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
        sbetm2 /= sbetm2 + GeoMath::sq(cbet1 + cbet2);
        w.dnm = std::sqrt(1 + _ep2 * sbetm2);
        double omg12 = lam12 / (_f1 * w.dnm);
        somg12 = sin(omg12);
        comg12 = cos(omg12);
    }
    else
    {
        somg12 = slam12;
        comg12 = clam12;
    }

    w.salp1 = cbet2 * somg12;
    w.calp1 = comg12 >= 0
                  ? sbet12 + cbet2 * sbet1 * GeoMath::sq(somg12) / (1 + comg12)
                  : sbet12a - cbet2 * sbet1 * GeoMath::sq(somg12) / (1 - comg12);

    double ssig12 = std::hypot(w.salp1, w.calp1),
           csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

    if (shortline && ssig12 < _etol2)
    {
        // really short lines
        w.salp2 = cbet1 * somg12;
        w.calp2 = sbet12 - cbet1 * sbet2 * (comg12 >= 0 ? GeoMath::sq(somg12) / (1 + comg12) : 1 - comg12);
        GeoMath::norm(p, w.salp2, w.calp2);
        w.salp2 = p.first;
        w.calp2 = p.second;
        // Set return value
        w.sig12 = std::atan2(ssig12, csig12);
    }
    else if (std::abs(_n) > 0.1 || // Skip astroid calc if too eccentric
             csig12 >= 0 || ssig12 >= 6 * std::abs(_n) * M_PI * GeoMath::sq(cbet1))
    {
        // Nothing to do, zeroth order spherical approximation is OK
    }
    else
    {
        // Scale lam12 and bet2 to x, y coordinate system where antipodal point
        // is at origin and singular point is at y = 0, x = -1.
        double x, y, lamscale, betscale;
        double lam12x = std::atan2(-slam12, -clam12); // lam12 - pi
        if (_f >= 0)
        { // In fact f == 0 does not get here
            // x = dlong, y = dlat
            {
                double k2 = GeoMath::sq(sbet1) * _ep2,
                       eps = k2 / (2 * (1 + std::sqrt(1 + k2)) + k2);
                lamscale = _f * cbet1 * A3f(eps) * M_PI;
            }
            betscale = lamscale * cbet1;

            x = lam12x / lamscale;
            y = sbet12a / betscale;
        }
        else
        { // _f < 0
            // x = dlat, y = dlong
            double cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
                   bet12a = std::atan2(sbet12a, cbet12a);
            double m12b, m0;
            // In the case of lon12 = 180, this repeats a calculation made in
            // Inverse.
            Lengths(v, _n, M_PI + bet12a,
                    sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
                    cbet1, cbet2, GeodesicMask::REDUCEDLENGTH, C1a, C2a);
            m12b = v.m12b;
            m0 = v.m0;

            x = -1 + m12b / (cbet1 * cbet2 * m0 * M_PI);
            betscale = x < -0.01 ? sbet12a / x
                                 : -_f * GeoMath::sq(cbet1) * M_PI;
            lamscale = betscale / cbet1;
            y = lam12x / lamscale;
        }

        if (y > -tol1_ && x > -1 - xthresh_)
        {
            // strip near cut
            if (_f >= 0)
            {
                w.salp1 = std::min<double>(1.0, -x);
                w.calp1 = -std::sqrt(1 - GeoMath::sq(w.salp1));
            }
            else
            {
                w.calp1 = std::max<double>(x > -tol1_ ? 0.0 : -1.0, x);
                w.salp1 = std::sqrt(1 - GeoMath::sq(w.calp1));
            }
        }
        else
        {

            double k = Astroid(x, y);
            double omg12a = lamscale * (_f >= 0 ? -x * k / (1 + k) : -y * (1 + k) / k);
            somg12 = sin(omg12a);
            comg12 = -cos(omg12a);
            // Update spherical estimate of alp1 using omg12 instead of lam12
            w.salp1 = cbet2 * somg12;
            w.calp1 = sbet12a - cbet2 * sbet1 * GeoMath::sq(somg12) / (1 - comg12);
        }
    }
    // Sanity check on starting guess.  Backwards check allows NaN through.
    if (!(w.salp1 <= 0))
    {
        GeoMath::norm(p, w.salp1, w.calp1);
        w.salp1 = p.first;
        w.calp1 = p.second;
    }
    else
    {
        w.salp1 = 1;
        w.calp1 = 0;
    }
    return w;
}
Geodesic::GeodesicData Geodesic::KarneyFF::Inverse(double lat1, double lon1, double lat2, double lon2, int outmask)
{

    outmask &= GeodesicMask::OUT_MASK;
    InverseData result = InverseInt(lat1, lon1, lat2, lon2, outmask);
    GeodesicData r = result.g;
    if ((outmask & GeodesicMask::AZIMUTH) != 0)
    {
       r.azi1 = GeoMath::atan2d(result.salp1, result.calp1);
        r.azi2 = GeoMath::atan2d(result.salp2, result.calp2); 
    }
    return r;
}
Geodesic::GeodesicData Geodesic::KarneyFF::Inverse(double lat1, double lon1, double lat2, double lon2)
{
    return KarneyFF::Inverse(lat1, lon1, lat2, lon2, GeodesicMask::STANDARD);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::InverseLine(double lat1, double lon1, double lat2, double lon2)
{
    return InverseLine(lat1, lon1, lat2, lon2, GeodesicMask::ALL);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::InverseLine(double lat1, double lon1, double lat2, double lon2, int caps)
{
    InverseData result = InverseInt(lat1, lon1, lat2, lon2, 0);
    double salp1 = result.salp1, calp1 = result.calp1,
           azi1 = GeoMath::atan2d(salp1, calp1), a12 = result.g.a12;
    // Ensure that a12 can be converted to a distance
    if ((caps & (GeodesicMask::OUT_MASK & GeodesicMask::DISTANCE_IN)) != 0)
        caps |= GeodesicMask::DISTANCE;
    return Geodesic::GeodesicLine(*this, lat1, lon1, azi1, salp1, calp1, caps, true, a12);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::Line(double lat1, double lon1, double azi1)
{
    return Line(lat1, lon1, azi1, GeodesicMask::ALL);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::Line(double lat1, double lon1, double azi1, int caps)
{
    return GeodesicLine(*this, lat1, lon1, azi1, caps);
}
Geodesic::GeodesicData Geodesic::KarneyFF::Direct(double lat1, double lon1, double azi1, double s12)
{
    return Direct(lat1, lon1, azi1, false, s12, GeodesicMask::STANDARD);
}
Geodesic::GeodesicData Geodesic::KarneyFF::Direct(double lat1, double lon1, double azi1, double s12, int outmask)
{
    return Direct(lat1, lon1, azi1, false, s12, outmask);
}
Geodesic::GeodesicData Geodesic::KarneyFF::Direct(double lat1, double lon1, double azi1, bool arcmode, double s12_a12, int outmask)
{
    if (!arcmode)
        outmask |= GeodesicMask::DISTANCE_IN;
    return GeodesicLine(*this, lat1, lon1, azi1, outmask).Position(arcmode, s12_a12, outmask);
}
Geodesic::GeodesicData Geodesic::KarneyFF::ArcDirect(double lat1, double lon1, double azi1, double a12)
{
    return Direct(lat1, lon1, azi1, true, a12, GeodesicMask::STANDARD);
}
Geodesic::GeodesicData Geodesic::KarneyFF::ArcDirect(double lat1, double lon1, double azi1, double a12, int outmask)
{
    return Direct(lat1, lon1, azi1, true, a12, outmask);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::DirectLine(double lat1, double lon1, double azi1, double s12)
{
    return DirectLine(lat1, lon1, azi1, s12, GeodesicMask::ALL);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::DirectLine(double lat1, double lon1, double azi1, double s12, int caps)
{
    return GenDirectLine(lat1, lon1, azi1, false, s12, caps);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::ArcDirectLine(double lat1, double lon1, double azi1, double a12)
{
    return ArcDirectLine(lat1, lon1, azi1, a12, GeodesicMask::ALL);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::ArcDirectLine(double lat1, double lon1, double azi1, double a12, int caps)
{
    return GenDirectLine(lat1, lon1, azi1, true, a12, caps);
}
Geodesic::GeodesicLine Geodesic::KarneyFF::GenDirectLine(double lat1, double lon1, double azi1, bool arcmode, double s12_a12, int caps)
{
    azi1 = GeoMath::AngNormalize(azi1);
    double salp1, calp1;
    // Guard against underflow in salp0.  Also -0 is converted to +0.
    GeoMath::Pair<double> p = GeoMath::Pair<double>();
    GeoMath::sincosd(p, GeoMath::AngRound(azi1));
    salp1 = p.first;
    calp1 = p.second;
    // Automatically supply DISTANCE_IN if necessary
    if (!arcmode)
        caps |= GeodesicMask::DISTANCE_IN;
    return Geodesic::GeodesicLine(*this, lat1, lon1, azi1, salp1, calp1, caps, arcmode, s12_a12);
}

double Geodesic::KarneyFF::EllipsoidArea()
{
    return 4 * M_PI * _c2;
}
double Geodesic::KarneyFF::EquatorialRadius()
{
    return _a;
}
double Geodesic::KarneyFF::Flattening()
{
    return _f;
}
#pragma endregion

#pragma region GeodesicLine
const int Geodesic::GeodesicLine::nC1_ = Geodesic::KarneyFF::nC1_;
const int Geodesic::GeodesicLine::nC1p_ = Geodesic::KarneyFF::nC1p_;
const int Geodesic::GeodesicLine::nC2_ = Geodesic::KarneyFF::nC2_;
const int Geodesic::GeodesicLine::nC3_ = Geodesic::KarneyFF::nC3_;
const int Geodesic::GeodesicLine::nC4_ = Geodesic::KarneyFF::nC4_;

Geodesic::GeodesicLine::GeodesicLine()
{
    _caps = 0;
}
Geodesic::GeodesicLine::GeodesicLine(KarneyFF &g, double lat1, double lon1, double azi1)
{
    GeodesicLine(g, lat1, lon1, azi1, GeodesicMask::ALL);
}
Geodesic::GeodesicLine::GeodesicLine(KarneyFF &g, double lat1, double lon1, double azi1, int caps)
{
    azi1 = GeoMath::AngNormalize(azi1);
    double salp1, calp1;
    GeoMath::Pair<double> p = GeoMath::Pair<double>();

    GeoMath::sincosd(p, GeoMath::AngRound(azi1));
    salp1 = p.first;
    calp1 = p.second;
    LineInit(g, lat1, lon1, azi1, salp1, calp1, caps, p);
}
Geodesic::GeodesicLine::GeodesicLine(KarneyFF &g, double lat1, double lon1, double azi1, double salp1, double calp1, int caps, bool arcmode, double s13_a13)
{
    GeoMath::Pair<double> p = GeoMath::Pair<double>();
    LineInit(g, lat1, lon1, azi1, salp1, calp1, caps, p);
    GeodesicLine::GetSetDistance(arcmode, s13_a13);
}
void Geodesic::GeodesicLine::LineInit(KarneyFF &g, double lat1, double lon1, double azi1, double salp1, double calp1, int caps, GeoMath::Pair<double> &p)
{
    _a = g._a;
    _f = g._f;
    _b = g._b;
    _c2 = g._c2;
    _f1 = g._f1;
    // Always allow latitude and azimuth and unrolling the longitude
    _caps = caps | GeodesicMask::LATITUDE | GeodesicMask::AZIMUTH |
            GeodesicMask::LONG_UNROLL;

    _lat1 = GeoMath::LatFix(lat1);
    _lon1 = lon1;
    _azi1 = azi1;
    _salp1 = salp1;
    _calp1 = calp1;
    double cbet1, sbet1;
    GeoMath::sincosd(p, GeoMath::AngRound(_lat1));
    sbet1 = _f1 * p.first;
    cbet1 = p.second;
    // Ensure cbet1 = +epsilon at poles
    GeoMath::norm(p, sbet1, cbet1);
    sbet1 = p.first;
    cbet1 = std::max<double>(KarneyFF::tiny_, p.second);
    _dn1 = sqrt(1 + g._ep2 * GeoMath::sq(sbet1));

    // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    _salp0 = _salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = Math.hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    _calp0 = hypot(_calp1, _salp1 * sbet1);
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    _ssig1 = sbet1;
    _somg1 = _salp0 * sbet1;
    _csig1 = _comg1 = sbet1 != 0 || _calp1 != 0 ? cbet1 * _calp1 : 1;
    GeoMath::norm(p, _ssig1, _csig1);
    _ssig1 = p.first;
    _csig1 = p.second; // sig1 in (-pi, pi]
    // GeoMath.norm(_somg1, _comg1); -- don't need to normalize!

    _k2 = GeoMath::sq(_calp0) * g._ep2;
    double eps = _k2 / (2 * (1 + sqrt(1 + _k2)) + _k2);

    if ((_caps & GeodesicMask::CAP_C1) != 0)
    {
        _A1m1 = KarneyFF::A1m1f(eps);
        _C1a = std::vector<double>(nC1_ + 1);
        KarneyFF::C1f(eps, _C1a);
        _B11 = GeoMath::SinCosSeries(true, _ssig1, _csig1, _C1a);
        double s = sin(_B11), c = cos(_B11);
        // tau1 = sig1 + B11
        _stau1 = _ssig1 * c + _csig1 * s;
        _ctau1 = _csig1 * c - _ssig1 * s;
        // Not necessary because C1pa reverts C1a
        //    _B11 = -SinCosSeries(true, _stau1, _ctau1, _C1pa, nC1p_);
    }

    if ((_caps & GeodesicMask::CAP_C1p) != 0)
    {
        _C1pa = std::vector<double>(nC1p_ + 1);
        KarneyFF::C1pf(eps, _C1pa);
    }

    if ((_caps & GeodesicMask::CAP_C2) != 0)
    {
        _C2a = std::vector<double>(nC2_ + 1);
        _A2m1 = KarneyFF::A2m1f(eps);
        KarneyFF::C2f(eps, _C2a);
        _B21 = GeoMath::SinCosSeries(true, _ssig1, _csig1, _C2a);
    }

    if ((_caps & GeodesicMask::CAP_C3) != 0)
    {
        _C3a = std::vector<double>(nC3_);
        g.C3f(eps, _C3a);
        _A3c = -_f * _salp0 * g.A3f(eps);
        _B31 = GeoMath::SinCosSeries(true, _ssig1, _csig1, _C3a);
    }

    if ((_caps & GeodesicMask::CAP_C4) != 0)
    {
        _C4a = std::vector<double>(nC4_);
        g.C4f(eps, _C4a);
        // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
        _A4 = GeoMath::sq(_a) * _calp0 * _salp0 * g._e2;
        _B41 = GeoMath::SinCosSeries(false, _ssig1, _csig1, _C4a);
    }
}
Geodesic::GeodesicData Geodesic::GeodesicLine::Position(double s12)
{
    return Position(false, s12, GeodesicMask::STANDARD);
}
Geodesic::GeodesicData Geodesic::GeodesicLine::Position(double s12, int outmask)
{
    return Position(false, s12, outmask);
}
Geodesic::GeodesicData Geodesic::GeodesicLine::ArcPosition(double a12)
{
    return Position(true, a12, GeodesicMask::STANDARD);
}
Geodesic::GeodesicData Geodesic::GeodesicLine::ArcPosition(double a12, int outmask)
{
    return Position(true, a12, outmask);
}
Geodesic::GeodesicData Geodesic::GeodesicLine::Position(bool arcmode, double s12_a12, int outmask)
{
    outmask &= _caps & GeodesicMask::OUT_MASK;
    Geodesic::GeodesicData r = Geodesic::GeodesicData();
    if (!(Init() &&
          (arcmode ||
           (_caps & (GeodesicMask::OUT_MASK & GeodesicMask::DISTANCE_IN)) != 0)))
        // Uninitialized or impossible distance calculation requested
        return r;
    r.lat1 = _lat1;
    r.azi1 = _azi1;
    r.lon1 = ((outmask & GeodesicMask::LONG_UNROLL) != 0) ? _lon1 : GeoMath::AngNormalize(_lon1);

    // Avoid warning about uninitialized B12.
    double sig12, ssig12, csig12, B12 = 0, AB1 = 0;
    if (arcmode)
    {
        // Interpret s12_a12 as spherical arc length
        r.a12 = s12_a12;
        sig12 = Geodesic::Units::toRad(s12_a12);
        GeoMath::Pair<double> p = GeoMath::Pair<double>();
        GeoMath::sincosd(p, s12_a12);
        ssig12 = p.first;
        csig12 = p.second;
    }
    else
    {
        // Interpret s12_a12 as distance
        r.s12 = s12_a12;
        double
            tau12 = s12_a12 / (_b * (1 + _A1m1)),
            s = sin(tau12),
            c = cos(tau12);
        // tau2 = tau1 + tau12
        B12 = -GeoMath::SinCosSeries(true,
                                     _stau1 * c + _ctau1 * s,
                                     _ctau1 * c - _stau1 * s,
                                     _C1pa);
        sig12 = tau12 - (B12 - _B11);
        ssig12 = sin(sig12);
        csig12 = cos(sig12);
        if (std::abs(_f) > 0.01)
        {

            double
                ssig2 = _ssig1 * csig12 + _csig1 * ssig12,
                csig2 = _csig1 * csig12 - _ssig1 * ssig12;
            B12 = GeoMath::SinCosSeries(true, ssig2, csig2, _C1a);
            double serr = (1 + _A1m1) * (sig12 + (B12 - _B11)) - s12_a12 / _b;
            sig12 = sig12 - serr / sqrt(1 + _k2 * GeoMath::sq(ssig2));
            ssig12 = sin(sig12);
            csig12 = cos(sig12);
            // Update B12 below
        }
        r.a12 = Geodesic::Units::toDeg(sig12);
    }

    double ssig2, csig2, sbet2, cbet2, salp2, calp2;
    // sig2 = sig1 + sig12
    ssig2 = _ssig1 * csig12 + _csig1 * ssig12;
    csig2 = _csig1 * csig12 - _ssig1 * ssig12;
    double dn2 = sqrt(1 + _k2 * GeoMath::sq(ssig2));
    if ((outmask & (GeodesicMask::DISTANCE | GeodesicMask::REDUCEDLENGTH |
                    GeodesicMask::GEODESICSCALE)) != 0)
    {
        if (arcmode || std::abs(_f) > 0.01)
            B12 = GeoMath::SinCosSeries(true, ssig2, csig2, _C1a);
        AB1 = (1 + _A1m1) * (B12 - _B11);
    }
    // sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = _calp0 * ssig2;
    // Alt: cbet2 = Math.hypot(csig2, salp0 * ssig2);
    cbet2 = hypot(_salp0, _calp0 * csig2);
    if (cbet2 == 0)
        // I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
        cbet2 = csig2 = KarneyFF::tiny_;
    // tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = _salp0;
    calp2 = _calp0 * csig2; // No need to normalize

    if ((outmask & GeodesicMask::DISTANCE) != 0 && arcmode)
        r.s12 = _b * ((1 + _A1m1) * sig12 + AB1);

    if ((outmask & GeodesicMask::LONGITUDE) != 0)
    {
        // tan(omg2) = sin(alp0) * tan(sig2)
        double somg2 = _salp0 * ssig2, comg2 = csig2, // No need to normalize
            E = std::copysign(1, _salp0);             // east or west going?
        // omg12 = omg2 - omg1
        double omg12 = ((outmask & GeodesicMask::LONG_UNROLL) != 0)
                           ? E * (sig12 - (atan2(ssig2, csig2) - atan2(_ssig1, _csig1)) + (atan2(E * somg2, comg2) - atan2(E * _somg1, _comg1)))
                           : atan2(somg2 * _comg1 - comg2 * _somg1,
                                   comg2 * _comg1 + somg2 * _somg1);
        double lam12 = omg12 + _A3c *
                                   (sig12 + (GeoMath::SinCosSeries(true, ssig2, csig2, _C3a) - _B31));
        double lon12 = Geodesic::Units::toDeg(lam12);
        r.lon2 = ((outmask & GeodesicMask::LONG_UNROLL) != 0) ? _lon1 + lon12 : GeoMath::AngNormalize(r.lon1 + GeoMath::AngNormalize(lon12));
    }

    if ((outmask & GeodesicMask::LATITUDE) != 0)
        r.lat2 = GeoMath::atan2d(sbet2, _f1 * cbet2);

    if ((outmask & GeodesicMask::AZIMUTH) != 0)
        r.azi2 = GeoMath::atan2d(salp2, calp2);

    if ((outmask &
         (GeodesicMask::REDUCEDLENGTH | GeodesicMask::GEODESICSCALE)) != 0)
    {
        double
            B22 = GeoMath::SinCosSeries(true, ssig2, csig2, _C2a),
            AB2 = (1 + _A2m1) * (B22 - _B21),
            J12 = (_A1m1 - _A2m1) * sig12 + (AB1 - AB2);
        if ((outmask & GeodesicMask::REDUCEDLENGTH) != 0)
            // Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
            // accurate cancellation in the case of coincident points.
            r.m12 = _b * ((dn2 * (_csig1 * ssig2) - _dn1 * (_ssig1 * csig2)) - _csig1 * csig2 * J12);
        if ((outmask & GeodesicMask::GEODESICSCALE) != 0)
        {
            double t = _k2 * (ssig2 - _ssig1) * (ssig2 + _ssig1) / (_dn1 + dn2);
            r.M12 = csig12 + (t * ssig2 - csig2 * J12) * _ssig1 / _dn1;
            r.M21 = csig12 - (t * _ssig1 - _csig1 * J12) * ssig2 / dn2;
        }
    }

    if ((outmask & GeodesicMask::AREA) != 0)
    {
        double
            B42 = GeoMath::SinCosSeries(false, ssig2, csig2, _C4a);
        double salp12, calp12;
        if (_calp0 == 0 || _salp0 == 0)
        {
            // alp12 = alp2 - alp1, used in atan2 so no need to normalize
            salp12 = salp2 * _calp1 - calp2 * _salp1;
            calp12 = calp2 * _calp1 + salp2 * _salp1;
        }
        else
        {
            salp12 = _calp0 * _salp0 *
                     (csig12 <= 0 ? _csig1 * (1 - csig12) + ssig12 * _ssig1 : ssig12 * (_csig1 * ssig12 / (1 + csig12) + _ssig1));
            calp12 = GeoMath::sq(_salp0) + GeoMath::sq(_calp0) * _csig1 * csig2;
        }
        r.S12 = _c2 * atan2(salp12, calp12) + _A4 * (B42 - _B41);
    }

    return r;
}

void Geodesic::GeodesicLine::SetDistance(double s13)
{
    _s13 = s13;
    Geodesic::GeodesicData g = Position(false, _s13, 0);
    _a13 = g.a12;
}
void Geodesic::GeodesicLine::SetArc(double a13)
{
    _a13 = a13;
    Geodesic::GeodesicData g = Position(true, _a13, GeodesicMask::DISTANCE);
    _s13 = g.s12;
}
void Geodesic::GeodesicLine::GetSetDistance(bool arcmode, double s13_a13)
{
    if (arcmode)
        SetArc(s13_a13);
    else
        SetDistance(s13_a13);
}
bool Geodesic::GeodesicLine::Init()
{
    return _caps != 0;
}
double Geodesic::GeodesicLine::Latitude()
{
    return Init() ? _lat1 : std::numeric_limits<double>::quiet_NaN();
}
double Geodesic::GeodesicLine::Longitude()
{
    return Init() ? _lon1 : std::numeric_limits<double>::quiet_NaN();
    0;
}
double Geodesic::GeodesicLine::Azimuth()
{
    return Init() ? _azi1 : std::numeric_limits<double>::quiet_NaN();
}
Geodesic::GeoMath::Pair<double> Geodesic::GeodesicLine::AzimuthCosine()
{
    return GeoMath::Pair<double>{
        Init() ? _salp1 : std::numeric_limits<double>::quiet_NaN(),
        Init() ? _calp1 : std::numeric_limits<double>::quiet_NaN(),
    };
}
double Geodesic::GeodesicLine::EquatorialAzimuth()
{
    return Init() ? GeoMath::atan2d(_salp0, _calp0) : std::numeric_limits<double>::quiet_NaN();
}
Geodesic::GeoMath::Pair<double> Geodesic::GeodesicLine::EquatorialAzimuthCosine()
{
    return GeoMath::Pair<double>{
        Init() ? _salp0 : std::numeric_limits<double>::quiet_NaN(),
        Init() ? _calp0 : std::numeric_limits<double>::quiet_NaN(),
    };
}
double Geodesic::GeodesicLine::EquatorialRadius()
{
    return Init() ? _a : std::numeric_limits<double>::quiet_NaN();
}
double Geodesic::GeodesicLine::Flattening()
{
    return Init() ? _f : std::numeric_limits<double>::quiet_NaN();
}
int Geodesic::GeodesicLine::Capabilities()
{
    return _caps;
}
bool Geodesic::GeodesicLine::Capabilities(int testcaps)
{
    testcaps &= GeodesicMask::OUT_ALL;
    return (_caps & testcaps) == testcaps;
}
double Geodesic::GeodesicLine::GenDistance(bool arcmode)
{
    return Init() ? (arcmode ? _a13 : _s13) : std::numeric_limits<double>::quiet_NaN();
}
double Geodesic::GeodesicLine::Distance()
{
    return GenDistance(false);
}
double Geodesic::GeodesicLine::Arc()
{
    return GenDistance(true);
}
#pragma endregion

Geodesic::GeodParams Geodesic::geodParameters;
#pragma region Geodesic namespace
void Geodesic::configureGeod(GEOD_TYPE geod_type, GeodParams &geodParams)
{
    switch (geod_type)
    {
    case GEOD_TYPE::WGS_84:
        geodParams.semi_major_axis = 6378137.0;
        geodParams.flattening = 1 / 298.257223563;
        break;
    case GEOD_TYPE::GRS_80:
        geodParams.semi_major_axis = 6378137.0;
        geodParams.flattening = 1 / 298.257222101;
        break;
    case GEOD_TYPE::AIRY_1830:
        geodParams.semi_major_axis = 6377563.396;
        geodParams.flattening = 1 / 299.3249646;
        break;
    case GEOD_TYPE::INTL_1924:
        geodParams.semi_major_axis = 6378388.0;
        geodParams.flattening = 1 / 297.0;
        break;
    case GEOD_TYPE::CLARKE_1880:
        geodParams.semi_major_axis = 6378249.145;
        geodParams.flattening = 1 / 293.465;
        break;
    case GEOD_TYPE::GRS_67:
        geodParams.semi_major_axis = 6378160.0;
        geodParams.flattening = 1 / 298.25;
        break;

    default: // WGS84
        geodParams.semi_major_axis = 6378137.0;
        geodParams.flattening = 1 / 298.257223563;
        break;
    }
    geodParams.geod_type = geod_type;
    geodParams.semi_minor_axis = geodParams.semi_major_axis * (1 - geodParams.flattening);
    geodParams.eccentricity = sqrt(1 - pow(geodParams.semi_minor_axis / geodParams.semi_major_axis, 2));

    kObj.updateKarneyFFProperties(kObj, geodParams);
}
void Geodesic::configureGeod(double semi_major_axis, double flattening, GeodParams &geodParams)
{
    geodParams.geod_type = GEOD_TYPE::CUSTOM;
    geodParams.semi_major_axis = semi_major_axis;
    geodParams.flattening = flattening;
    geodParams.semi_minor_axis = geodParams.semi_major_axis * (1 - geodParams.flattening);
    geodParams.eccentricity = sqrt(1 - pow(geodParams.semi_minor_axis / geodParams.semi_major_axis, 2));
    kObj.updateKarneyFFProperties(kObj, geodParams);
}
#pragma endregion

#pragma region GEODESIC UNIT CONVERSIONS
double Geodesic::Units::toDeg(double radian, double arcmin, double arcsec)
{
    double degree = radian * 180.0 / M_PI;
    degree += arcmin / 60.0;
    degree += arcsec / 3600.0;
    return degree;
}
double Geodesic::Units::toRad(double degree, double arcmin, double arcsec)
{
    degree += arcmin / 60.0;
    degree += arcsec / 3600.0;
    return degree / 180.0 * M_PI;
}
double Geodesic::Units::toArcmin(double degree, double radian, double arcsec)
{
    degree += radian * 180 / M_PI;
    degree += arcsec / 3600;
    return degree * 60;
}
double Geodesic::Units::toArcsec(double degree, double radian, double arcmin)
{
    degree += radian * 180 / M_PI;
    degree += arcmin / 60;
    return degree * 3600;
}
double Geodesic::Units::toKilometers(double meters, double miles, double feet, double nautical)
{
    return meters / 1000 + feet / 3280.8399 + nautical / 0.54 + miles * 1.609344;
}
double Geodesic::Units::toMeters(double kilometers, double miles, double feet, double nautical)
{
    return (kilometers + toKilometers(0, miles, feet, nautical)) * 1000;
}
double Geodesic::Units::toMiles(double kilometers, double meters, double feet, double nautical)
{
    return (kilometers + toKilometers(meters, 0, feet, nautical)) / 1.609344;
}
double Geodesic::Units::toFeet(double kilometers, double meters, double miles, double nautical)
{
    return (miles + toMiles(kilometers, meters, 0, nautical)) * 5280;
}
double Geodesic::Units::toNautical(double kilometers, double meters, double miles, double feet)
{
    return (kilometers + toKilometers(meters, miles, feet, 0)) / 1.852;
}
double Geodesic::Units::degToKilometers(double degrees, bool square)
{
    const double rate = 110.567;
    if (square = 1)
        return degrees * rate * rate;
    return degrees * rate;
}
#pragma endregion

