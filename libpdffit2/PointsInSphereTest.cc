/***********************************************************************
* Short Title: unit tests for PointsInSphere module
*
* Comments:
*
* $Id: PointsInSphereTest.cc,v 1.1 2006/03/21 23:03:36 juhas Exp $
*
* <license text>
***********************************************************************/

#include <algorithm>
#include <stdexcept>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "PointsInSphere.h"
using namespace std;
using namespace CppUnit;
using namespace NS_POINTSINSPHERE;

const double eps = 1.0e-12;

struct vidxgroup
{
    double vijk[4];
    vidxgroup(double v, int* ijk)
    {
	vijk[0] = v;
	for (size_t i = 0; i != 3; ++i) { vijk[i+1] = ijk[i]; }
    }
    vidxgroup(double v, int i, int j, int k)
    {
	vijk[0] = v; vijk[1] = i; vijk[2] = j; vijk[3] = k;
    }
};

bool operator<(const vidxgroup &x, const vidxgroup &y)
{
    return (x.vijk[0] < y.vijk[0] - eps) ||
	    lexicographical_compare(x.vijk+1, x.vijk+4, y.vijk+1, y.vijk+4);
}

bool operator==(const vidxgroup &x, const vidxgroup &y)
{
    bool eq = (fabs(x.vijk[0] - y.vijk[0]) < eps) &&
	    equal(x.vijk+1, x.vijk+4, y.vijk+1);
    return eq;
}

ostream& operator<<(ostream &s, const vidxgroup &x)
{
    return s << "<" << x.vijk[0] << ";" << int(x.vijk[1])
	<< ',' << int(x.vijk[2]) << ',' << int(x.vijk[3]) << '>';
}


////////////////////////////////////////////////////////////////////////
// PointsInSphereTest
////////////////////////////////////////////////////////////////////////

class PointsInSphereTest : public TestFixture
{

    CPPUNIT_TEST_SUITE(PointsInSphereTest);
    CPPUNIT_TEST(testCubic);
    CPPUNIT_TEST(testOrthorombic);
    CPPUNIT_TEST(testHexagonal);
    CPPUNIT_TEST(testFCC);
    CPPUNIT_TEST_SUITE_END();
private:
    LatticeParameters* latpar;
public:
    void setUp()
    {
	latpar = new LatticeParameters(1, 1, 1, 90, 90, 90);
    }
    void tearDown()
    {
	delete latpar;
    }
    int count(double Rmin, double Rmax)
    {
	int c = 0;
	for (   PointsInSphere sph(Rmin, Rmax, *latpar);
		not sph.finished(); sph.next(), ++c  )
	{ }
	return c;
    }
    vector<vidxgroup> sortedPoints(double Rmin, double Rmax)
    {
	vector<vidxgroup> ridx;
	for (   PointsInSphere sph(Rmin, Rmax, *latpar);
		not sph.finished(); sph.next()  )
	{
	    ridx.push_back(vidxgroup(sph.r(), sph.mno));
	}
	sort(ridx.begin(), ridx.end());
	return ridx;
    }
public:
    void testCubic()
    {
	latpar->a = latpar->b = latpar->c = 1.0;
	latpar->alpha = latpar->beta = latpar->gamma = 90.0;
	latpar->update();
	CPPUNIT_ASSERT_EQUAL(0, count(0.0, 0.0));
	CPPUNIT_ASSERT_EQUAL(0, count(eps, 0.5));
	CPPUNIT_ASSERT_EQUAL(0, count(1.0 + eps, 1.1));
	CPPUNIT_ASSERT_EQUAL(1, count(0.0, eps));
	CPPUNIT_ASSERT_EQUAL(7, count(0.0, 1 + eps));
	CPPUNIT_ASSERT_EQUAL(19, count(0.0, sqrt(2.0) + eps));
	CPPUNIT_ASSERT_EQUAL(12, count(1.0 + eps, sqrt(2.0) + eps));
    }
    void testOrthorombic()
    {
	latpar->a = 1.0; latpar->b = 2.0; latpar->c = 3.0;
	latpar->alpha = latpar->beta = latpar->gamma = 90.0;
	latpar->update();
	CPPUNIT_ASSERT_EQUAL(3, count(0.0, 1.1));
	CPPUNIT_ASSERT_EQUAL(4, count(1.9, 2.1));
	vidxgroup ep[] = {
	    vidxgroup(0, 0, 0, 0),
	    vidxgroup(1, -1, 0, 0),
	    vidxgroup(1, 1, 0, 0),
	    vidxgroup(2, -2, 0, 0),
	    vidxgroup(2, 0, -1, 0),
	    vidxgroup(2, 0, 1, 0),
	    vidxgroup(2, 2, 0, 0),
	    vidxgroup(sqrt(5.0), -1, -1, 0),
	    vidxgroup(sqrt(5.0), -1, 1, 0),
	    vidxgroup(sqrt(5.0), 1, -1, 0),
	    vidxgroup(sqrt(5.0), 1, 1, 0),
	    vidxgroup(sqrt(8.0), -2, -1, 0),
	    vidxgroup(sqrt(8.0), -2, 1, 0),
	    vidxgroup(sqrt(8.0), 2, -1, 0),
	    vidxgroup(sqrt(8.0), 2, 1, 0),
	    vidxgroup(3, -3, 0, 0),
	    vidxgroup(3, 0, 0, -1),
	    vidxgroup(3, 0, 0, 1),
	    vidxgroup(3, 3, 0, 0),
	};
	vector<vidxgroup> exp_pts(ep, ep + sizeof(ep)/sizeof(vidxgroup));
	vector<vidxgroup> act_pts = sortedPoints(0.0, 3.0+eps);
	CPPUNIT_ASSERT_EQUAL(exp_pts.size(), act_pts.size());
	for (size_t i = 0; i != exp_pts.size(); ++i)
	{
	    CPPUNIT_ASSERT_EQUAL(exp_pts[i], act_pts[i]);
	}
    }
    void testHexagonal()
    {
	latpar->a = 1.0; latpar->b = 1.0; latpar->c = 2.0;
	latpar->alpha = latpar->beta = 90.0; latpar->gamma = 120.0;
	latpar->update();
	CPPUNIT_ASSERT_EQUAL(7, count(0.0, 1+eps));
	vidxgroup ep[] = {
	    vidxgroup(0, 0, 0, 0), 
	    vidxgroup(1, -1, -1, 0), 
	    vidxgroup(1, -1, 0, 0), 
	    vidxgroup(1, 0, -1, 0), 
	    vidxgroup(1, 0, 1, 0), 
	    vidxgroup(1, 1, 0, 0), 
	    vidxgroup(1, 1, 1, 0), 
	    vidxgroup(sqrt(3.0), -2, -1, 0), 
	    vidxgroup(sqrt(3.0), -1, -2, 0), 
	    vidxgroup(sqrt(3.0), -1, 1, 0), 
	    vidxgroup(sqrt(3.0), 1, -1, 0), 
	    vidxgroup(sqrt(3.0), 1, 2, 0), 
	    vidxgroup(sqrt(3.0), 2, 1, 0), 
	    vidxgroup(2, -2, -2, 0), 
	    vidxgroup(2, -2, 0, 0), 
	    vidxgroup(2, 0, -2, 0), 
	    vidxgroup(2, 0, 0, -1), 
	    vidxgroup(2, 0, 0, 1), 
	    vidxgroup(2, 0, 2, 0), 
	    vidxgroup(2, 2, 0, 0), 
	    vidxgroup(2, 2, 2, 0), 
	};
	vector<vidxgroup> exp_pts(ep, ep + sizeof(ep)/sizeof(vidxgroup));
	vector<vidxgroup> act_pts = sortedPoints(0.0, 2.0+eps);
	CPPUNIT_ASSERT_EQUAL(exp_pts.size(), act_pts.size());
	for (size_t i = 0; i != exp_pts.size(); ++i)
	{
	    CPPUNIT_ASSERT_EQUAL(exp_pts[i], act_pts[i]);
	}
    }
    void testFCC()
    {
	latpar->a = latpar->b = latpar->c = sqrt(0.5);
	latpar->alpha = latpar->beta = latpar->gamma = 60.0;
	latpar->update();
	CPPUNIT_ASSERT_EQUAL(13, count(0.0, sqrt(0.5)+eps));
	CPPUNIT_ASSERT_EQUAL(19, count(0.0, 1.0+eps));
    }
};


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(PointsInSphereTest);

/***********************************************************************
* Too see what people have been up to just run:
*   cvs log
***********************************************************************/
