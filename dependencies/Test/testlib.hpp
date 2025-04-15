#pragma once
#include <iostream>
#include <cmath>


static int operation_counter = 0;
template <typename T>
bool assertEquals(T first, T second, double tol = 1e-12)
{
    auto a_dif = std::abs(first - second);
    bool value = a_dif <= tol;
    return value;
}
template <>
bool assertEquals(double first, double second, double tol)
{
    auto a_dif = std::abs(first - second);
    bool value = a_dif <= tol;
    return value;
}
bool assertTrue(bool expression)
{
    return expression ? expression : false;
}
bool assertFalse(bool expression)
{
    return expression ? !expression : true;
}
static bool equiv(double x, double y)
{
    // Test for equivalence
    return (std::isnan(x) && std::isnan(y)) ||
           (x == y && std::copysign(1, x) == std::copysign(1, y));
}
template <typename T>
void reprVector(std::vector<T> values)
{
    for (size_t i = 0; i < values.size(); i++)
    {
        std::cout << values[i] << " ";
    }
    std::cout << "" << std::endl;
}
void logResult(std::string function_name, bool success_result, int error_count, std::vector<int> errorInds , std::string Test_tag = "NOT SPECIFIED TEST")
{
    std::cout << "TEST TAG : "<< Test_tag << std::endl;
    std::cout << "FUNCTION UNDER TESTING : " << function_name << std::endl;
    std::cout << "SUCCESS_RESULT : " << (success_result ? "PASSED" : "FAILED") << std::endl;
    std::cout << "FAILED_CASE_COUNT : " << error_count << std::endl;
    std::cout << (error_count > 0 ? "FAILED_CASE_IND : " : "");
    reprVector(errorInds);

    std::cout << "" << std::endl;
}

