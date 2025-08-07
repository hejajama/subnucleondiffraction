#include "unit_test_framework.hpp"
#include "interpolation.hpp"

TEST(interpolate_square)
{
    double x[5] = {0, 1, 2, 3, 4};
    double y[5] = {0, 1, 4, 9, 16};
    Interpolator inter(x,y,5);
    inter.SetMethod(INTERPOLATE_SPLINE);
    inter.Initialize();
    
    ASSERT_ALMOST_EQUAL(inter.Evaluate(0),0,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(1),1,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(2),4,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(3),9,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(4),16,1e-6);

    // Same test with std::vector
    std::vector<double> xvec(x, x + 5);
    std::vector<double> yvec(y, y + 5);
    Interpolator inter2(xvec, yvec);
    inter2.SetMethod(INTERPOLATE_SPLINE);
    inter2.Initialize();
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(0),0,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(1),1,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(2),4,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(3),9,1e-6);
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(4),16,1e-6);
}

TEST(interpolate_out_of_range)
{
    double x[5] = {0, 1, 2, 3, 4};
    double y[5] = {0, 1, 4, 9, 16};
    Interpolator inter(x,y,5);
    inter.SetMethod(INTERPOLATE_SPLINE);
    inter.SetFreeze(true); // Freeze below and above the given data
    inter.Initialize();
    
    ASSERT_ALMOST_EQUAL(inter.Evaluate(-1),0,1e-6);
    ASSERT_ALMOST_EQUAL(inter.Evaluate(5),16,1e-6);

    Interpolator inter2(x,y,5);
    inter2.Initialize();
    try {
        inter2.Evaluate(-1);
        ASSERT_TRUE(false); // Should not reach here
    } catch (const std::out_of_range& e) {
        ASSERT_TRUE(true); // Exception caught as expected
    }
    
}

TEST(interpolate_copy_constructor)
{
    std::vector<double> x = {0, 1, 2, 3, 4};
    std::vector<double> y = {0, 1, 4, 9, 16};
    Interpolator inter(x, y);
    inter.Initialize();

    Interpolator inter2 (inter);

    inter.Clear();
    inter.Clear();

    ASSERT_ALMOST_EQUAL(inter2.Evaluate(3),9,1e-6);

    Interpolator* inter3 = new Interpolator(x,y);

    Interpolator inter4(*inter3);

    inter3->Clear();
    delete inter3;

    ASSERT_ALMOST_EQUAL(inter4.Evaluate(3),9,1e-6);

    // Same test but construct interpolator using x and y arrays, as 
    // in this case Interpolator does not reserve memory for x and y data
    double x_arr[5] = {0, 1, 2, 3, 4};
    double y_arr[5] = {0, 1, 4, 9, 16};
    Interpolator inter5(x_arr, y_arr, 5);
    inter5.Initialize();

    Interpolator inter6(inter5);

    inter5.Clear();

    ASSERT_ALMOST_EQUAL(inter6.Evaluate(3), 9, 1e-6);


}

TEST(interpolate_assignment_operator)
{
    std::vector<double> x = {0, 1, 2, 3, 4};
    std::vector<double> y = {0, 1, 4, 9, 16};
    Interpolator inter(x, y);
    inter.Initialize();

    Interpolator inter2;

    inter2 = inter;


    inter.Clear();
    ASSERT_ALMOST_EQUAL(inter2.Evaluate(3), 9, 1e-6);

    Interpolator* inter3 = new Interpolator(x, y);
    inter3->Initialize();

    Interpolator inter4;
    inter4 = *inter3;

    inter3->Clear();
    delete inter3;

    ASSERT_ALMOST_EQUAL(inter4.Evaluate(3), 9, 1e-6);

    // Same test but construct interpolator using x and y arrays
    double x_arr[5] = {0, 1, 2, 3, 4};
    double y_arr[5] = {0, 1, 4, 9, 16};
    Interpolator inter5(x_arr, y_arr, 5);
    inter5.Initialize();

    Interpolator inter6;
    inter6 = inter5;

    inter5.Clear();

    ASSERT_ALMOST_EQUAL(inter6.Evaluate(3), 9, 1e-6);
}